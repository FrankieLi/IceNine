"""Shared pieces of the run-time profiling (Task 3): paths, the isolation record, the stage
instrumentation of the reconstructor and of the network stage (all by patching from scripts via
scripts/nn_hybrid/stage_timer.py; nothing in icenine/ changes), the interleaving order and small
statistics helpers.

Import this module first: it pins OMP / MKL / torch to one thread (the single-worker timing
condition) before numpy / torch are imported by the study scripts."""

import json
import os
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Sequence, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
CACHE = HERE / "cache"  # gitignored
OUT = ICENINE_PY / "benchmarks" / "profiling"
for _p in (
    HERE,
    ICENINE_PY / "scripts" / "nn_hybrid",
    ICENINE_PY / "scripts" / "coarse_proxy",
    ICENINE_PY / "scripts" / "findoptimal_robustness",
    ICENINE_PY / "scripts",
    ICENINE_PY / "benchmarks",
):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import stage_timer as ST  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import stats as shared_stats  # noqa: E402

torch.set_num_threads(1)

Q8_LOCAL = "local"


# ---------------------------------------------------------------------------
# Isolation record
# ---------------------------------------------------------------------------


def _run(cmd: Sequence[str]) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=20).stdout.strip()
    except Exception as e:  # pragma: no cover - environment dependent
        return f"unavailable: {e!r}"


def isolation_record(note: str = "") -> Dict[str, Any]:
    """Load average, power source, thermal state, thread settings and the busiest other
    processes, as the plan's isolation conditions ask. Called at the start and end of a run."""
    ps = _run(["ps", "-Ao", "pcpu,pid,comm", "-r"]).splitlines()[1:8]
    busy = []
    for line in ps:
        parts = line.split(None, 2)
        if len(parts) == 3 and float(parts[0]) >= 5.0 and int(parts[1]) != os.getpid():
            busy.append(f"{parts[0]}% {Path(parts[2]).name}")
    batt = _run(["pmset", "-g", "batt"]).splitlines()
    return dict(
        note=note,
        time=time.strftime("%Y-%m-%d %H:%M:%S"),
        loadavg=list(os.getloadavg()),
        power=batt[0] if batt else "",
        battery=batt[1].strip() if len(batt) > 1 else "",
        thermal=_run(["pmset", "-g", "therm"]).replace("\n", " | "),
        processes_over_5pct_cpu=busy,
        OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS"),
        MKL_NUM_THREADS=os.environ.get("MKL_NUM_THREADS"),
        torch_threads=torch.get_num_threads(),
        cpu=_run(["sysctl", "-n", "machdep.cpu.brand_string"]),
        cores=_run(["sysctl", "-n", "hw.perflevel0.logicalcpu"])
        + " P + "
        + _run(["sysctl", "-n", "hw.perflevel1.logicalcpu"])
        + " E",
    )


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=float))


# ---------------------------------------------------------------------------
# Interleaving and statistics
# ---------------------------------------------------------------------------


def rotated(pipes: Sequence[str], index: int) -> List[str]:
    """The order of the pipelines for one case: a rotation by the case index (so every pipeline
    runs first, second, ... equally often: the ABAB interleaving of the plan, generalised)."""
    k = index % len(pipes)
    return list(pipes[k:]) + list(pipes[:k])


def wilson(k: int, n: int, z: float = 1.96) -> Tuple[float, float, float]:
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    lo, hi = shared_stats.wilson(k, n, z)
    return k / n, lo, hi


def pct(x: Sequence[float], q: float) -> float:
    return float(np.percentile(np.asarray(x, dtype=float), q)) if len(x) else float("nan")


def summarize_times(x: Sequence[float]) -> Dict[str, float]:
    """median and the 10-90% range of per-run times."""
    return dict(n=len(x), median=pct(x, 50), p10=pct(x, 10), p90=pct(x, 90), mean=float(np.mean(x)))


# ---------------------------------------------------------------------------
# Stage instrumentation (patching only; results of the patched functions are untouched)
# ---------------------------------------------------------------------------


class Instrument:
    """Installs the plan's stage wrappers on one StageTimer.

    reconstructor: AdaptiveVoxelReconstructor.reconstruct_voxel ("reconstruct"),
        icenine.reconstructor.run_discrete_search_spaced (per level; patched where the reconstructor
        looks the name up), MCOptimizer.optimize (quick MC per level / FindOptimal / other MC by
        max_mc_steps), MCOptimizer.variance_minimizing_optimize ("variance"), and
        VoxelCostFunction.evaluate (calls and time by role: local / global / mc).
    network (net=True): perturbation_sweep.prepare_nominal ("prepare"), the harness renderers
        render_windows / render_distractor_windows / make_realistic_dataset, decode_windows
        ("decode"), Gauss-Newton extract_measurements and CentroidGaussNewton.solve.
    Use `with Instrument(T, rec_params).installed():` or call install()/uninstall()."""

    def __init__(self, timer: ST.StageTimer, max_mc_steps: int, net: bool = True):
        self.T = timer
        self.max_steps = int(max_mc_steps)
        self.net = net
        self.tag = "none"  # L0..L3 inside reconstruct_voxel, "f1" inside the F1 post hoc
        self.level = -1
        self.roles: Dict[int, str] = {}  # id(cost fn instance) -> role
        self._installed = False

    def mark(self, obj: Any, role: str) -> None:
        self.roles[id(obj)] = role

    def _eval_name(self, a: tuple, k: dict) -> str:
        obj = a[0]
        role = self.roles.get(id(obj))
        if role is None:
            role = "local" if getattr(obj, "pixel_radius", 1) == 0 else "global"
        return "evaluate_" + role

    def install(self) -> None:
        import icenine.reconstructor as RC
        from icenine.cost_functions import VoxelCostFunction
        from icenine.orientation_search import MCOptimizer

        T = self.T

        def total_name(a: tuple, k: dict) -> str:
            self.level, self.tag = -1, "none"
            return "reconstruct"

        def discrete_name(a: tuple, k: dict) -> str:
            self.level += 1
            self.tag = f"L{self.level}"
            return f"discrete_{self.tag}"

        def mc_name(a: tuple, k: dict) -> str:
            steps = k.get("max_mc_steps")
            if steps == 10:
                return f"quick_mc_{self.tag}"
            return "find_optimal" if steps == self.max_steps else "mc_other"

        T.patch(RC.AdaptiveVoxelReconstructor, "reconstruct_voxel", total_name)
        T.patch(RC, "run_discrete_search_spaced", discrete_name)
        T.patch(MCOptimizer, "optimize", mc_name)
        T.patch(MCOptimizer, "variance_minimizing_optimize", "variance")
        T.patch(VoxelCostFunction, "evaluate", self._eval_name)
        if self.net:
            import icenine.orientation_baselines as OB
            import icenine.orientation_eval as OE
            import perturbation_sweep as ps

            T.patch(ps, "prepare_nominal", "prepare")
            T.patch(OE, "render_windows", "render_physics")
            T.patch(OE, "render_distractor_windows", "render_distractors")
            T.patch(OE, "make_realistic_dataset", "render_realism")
            T.patch(OE, "decode_windows", "decode")
            T.patch(OB, "extract_measurements", "gn_extract")
            T.patch(OB.CentroidGaussNewton, "solve", "gn_solve")
        self._installed = True

    def uninstall(self) -> None:
        self.T.uninstall()
        self._installed = False

    def installed(self) -> Any:
        import contextlib

        @contextlib.contextmanager
        def cm() -> Any:
            self.install()
            try:
                yield self
            finally:
                self.uninstall()

        return cm()

    # per-run bookkeeping -------------------------------------------------------------------

    def eval_counts(self, before: Dict[str, Any]) -> Dict[str, int]:
        c0 = before["calls"]
        return {
            r: int(self.T.calls.get("evaluate_" + r, 0) - c0.get("evaluate_" + r, 0))
            for r in ("global", "local", "mc")
        }


def stage_delta(T: ST.StageTimer, before: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
    """Inclusive and exclusive seconds and call counts per stage since `before` (a snapshot)."""
    out: Dict[str, Dict[str, float]] = {}
    for kind in ("inclusive", "exclusive"):
        out[kind] = {k: round(v, 6) for k, v in T.delta_since(before, kind).items()}
    calls = {}
    for k, v in T.calls.items():
        d = v - before["calls"].get(k, 0)
        if d:
            calls[k] = int(d)
    out["calls"] = calls  # type: ignore[assignment]
    return out


def reduced_err_deg(R: np.ndarray, R_true: np.ndarray) -> float:
    """Cubic-reduced misorientation (deg), the error of the E0 / sweep studies."""
    from csl import reduced_misorientation_deg

    return float(reduced_misorientation_deg(np.asarray(R, dtype=np.float64), R_true))


def clone_ctx_with_passes(ctx: SimpleNamespace, passes: int) -> SimpleNamespace:
    """A copy of the worker context whose args.passes is `passes` (H1 = one net pass)."""
    c = SimpleNamespace(**vars(ctx))
    c.args = SimpleNamespace(**{**vars(ctx.args), "passes": passes})
    return c


def optional(x: Optional[Any], default: Any) -> Any:
    return default if x is None else x


def run_pool(
    func: Any, init: Any, wargs: Dict[str, Any], todo: List[Any], workers: int, label: str
) -> None:
    """Spawn pool with this study's own worker initialiser (instrumentation installed)."""
    import multiprocessing as mp

    t0 = time.time()
    with mp.get_context("spawn").Pool(workers, initializer=init, initargs=(wargs,)) as pool:
        for k, res in enumerate(pool.imap(func, todo), 1):  # ordered: case order = run order
            print(f"  [{label} {k}/{len(todo)}] {res} ({time.time() - t0:.0f}s total)", flush=True)
