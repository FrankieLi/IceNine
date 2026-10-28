"""S2 Step 0a: a no-start search (``reconstruct_voxel``) split by caller, and the same search with
the fused prototype in place of ``VoxelCostFunction.evaluate``.

Callers (wall time of the section, and the time and count of the ``evaluate`` calls made inside
it): ``screen`` (the level screen: ``evaluate`` of the pixel-radius cost, or ``BatchedScreen``'s
``ScreenEvaluator``), ``recost`` (the local cost of the screen passers inside
``run_discrete_search_spaced``), ``quick_mc`` (the per-level MC, sequential or ``BatchedQuickMC``),
``find_optimal`` (the optimizer runs of ``refine_from_candidates``), ``variance``
(``variance_minimizing_optimize``), ``final`` (the last ``evaluate``), ``rest`` (wall minus the
sections). The batched quick-MC cost (``BatchedCostEvaluator.evaluate_many``) is timed on its own
(orientations scored, seconds). Wrappers as in ``s4/step0_timing.py`` and S1's ``search_split.py``:
the library is patched in this process only.

Arms: ``default`` (the config's keys), ``bqmc`` (``BatchedQuickMC`` on, CPU), ``fused`` and
``fused_bqmc`` (the same with every covered ``evaluate`` call going through
``s2fused.FusedEval``; uncovered cost functions keep torch). A flip in the fused cost can change
the search path, so evaluation counts and results are reported per arm, side by side.

Sets: ``three`` = ThreeVoxels voxel 0 with ``Example2.Simulation.config`` (S1's search split);
``e0`` = the first ``--n-voxels`` E0 clean voxels (``ReconstructQ8`` physics, seed 0 of each voxel).
Single worker, preflight-gated (``--allow-busy`` labels the run "contended"); arms interleaved.

S2 Step 5 (``--lib``): the ``fused`` and ``fused_bqmc`` arms turn the library key ``FusedCost``
on instead of patching in the prototype (``evaluate`` and the batched quick-MC cost then run the
library's fused entries); ``evaluate_fused_n`` counts the calls of covered key-on cost functions
and ``batched_quick_mc.orientations_fused`` the orientations of the batched fused entry.

  uv run python scripts/speed_parity/s2/s2_search.py --set three --out benchmarks/s2_step0
  uv run python scripts/speed_parity/s2/s2_search.py --set e0 --n-voxels 10 --out ...
"""

import argparse
import contextlib
import io
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
torch.set_num_threads(1)

import s2fused as SF  # noqa: E402
from s2_timing import gate  # noqa: E402

from icenine import batched_cost as BC  # noqa: E402
from icenine import batched_screen as BS  # noqa: E402
from icenine import cost_functions as CF  # noqa: E402
from icenine import orientation_search as OS  # noqa: E402
from icenine import reconstructor as REC  # noqa: E402

CALLERS = ("screen", "recost", "quick_mc", "find_optimal", "variance", "final", "other")
ARMS = {
    "default": (False, False),
    "bqmc": (True, False),
    "fused": (False, True),
    "fused_bqmc": (True, True),
}


class Split:
    """Section stack, section walls and evaluate time per caller of one run."""

    def __init__(self) -> None:
        self.stack: List[str] = []
        self.wall: Dict[str, float] = defaultdict(float)
        self.ev_s: Dict[str, float] = defaultdict(float)
        self.ev_n: Dict[str, int] = defaultdict(int)
        self.fused_n = 0
        self.batched_s = 0.0
        self.batched_n = 0
        self.batched_fused_n = 0  # orientations the batched cost scored with the fused entry

    def caller(self, cost_fn: Any) -> str:
        st = self.stack
        if "variance" in st:
            return "variance"
        if "refine" in st:
            return "find_optimal" if "opt" in st else "final"
        if "discrete" in st:
            return "screen" if cost_fn.pixel_radius > 0 else "recost"
        if "opt" in st or "bqmc" in st:
            return "quick_mc"
        return "other"

    def section_label(self, name: str) -> str:
        if name == "opt":
            return "find_optimal" if "refine" in self.stack else "quick_mc"
        return {"bqmc": "quick_mc", "discrete": "discrete", "variance": "variance"}.get(name, name)


STATE: Dict[str, Any] = dict(split=None, fused=False)


def _section(name: str, fn: Callable[..., Any]) -> Callable[..., Any]:
    def wrapped(*a: Any, **k: Any) -> Any:
        sp: Optional[Split] = STATE["split"]
        if sp is None:
            return fn(*a, **k)
        label = sp.section_label(name)
        nested = name in sp.stack or (name == "opt" and "variance" in sp.stack)
        sp.stack.append(name)
        t0 = time.perf_counter()
        try:
            return fn(*a, **k)
        finally:
            dt = time.perf_counter() - t0
            sp.stack.pop()
            if not nested:
                sp.wall[label] += dt

    return wrapped


_orig_evaluate = CF.VoxelCostFunction.evaluate


def _evaluate(self: Any, orientation: Any, voxel_vertices: Any, phase_index: int = 0) -> Any:
    sp: Optional[Split] = STATE["split"]
    fe = None
    if STATE["fused"]:
        fe = self.__dict__.get("_s2_fused", False)
        if fe is False:
            try:
                fe = SF.FusedEval(self)
            except SF.Unsupported:
                fe = None
            self.__dict__["_s2_fused"] = fe
    t0 = time.perf_counter()
    if fe is not None:
        r = fe.evaluate(orientation, voxel_vertices, phase_index)
    else:
        r = _orig_evaluate(self, orientation, voxel_vertices, phase_index)
    if sp is not None:
        c = sp.caller(self)
        sp.ev_s[c] += time.perf_counter() - t0
        sp.ev_n[c] += 1
        # the library path (``--lib``): the key on, and the cost function covered by the fused path
        cov = self.__dict__.get("_s2_lib_covered")
        if cov is None:
            cov = bool(getattr(self, "fused_cost", False)) and (
                CF.fused_unsupported_reason(self) is None
            )
            self.__dict__["_s2_lib_covered"] = cov
        sp.fused_n += fe is not None or cov
    return r


def _screen(fn: Callable[..., Any]) -> Callable[..., Any]:
    def wrapped(*a: Any, **k: Any) -> Any:
        sp: Optional[Split] = STATE["split"]
        t0 = time.perf_counter()
        try:
            return fn(*a, **k)
        finally:
            if sp is not None:
                sp.ev_s["screen"] += time.perf_counter() - t0

    return wrapped


_orig_many = BC.BatchedCostEvaluator.evaluate_many


def _many(self: Any, orientations: Any, *a: Any, **k: Any) -> Any:
    sp: Optional[Split] = STATE["split"]
    t0 = time.perf_counter()
    n_fused = self.counts.get("orient_fused", 0)
    r = _orig_many(self, orientations, *a, **k)
    if sp is not None:
        sp.batched_s += time.perf_counter() - t0
        sp.batched_n += len(orientations)
        sp.batched_fused_n += self.counts.get("orient_fused", 0) - n_fused
    return r


def install() -> None:
    CF.VoxelCostFunction.evaluate = _evaluate  # type: ignore[method-assign]
    BC.BatchedCostEvaluator.evaluate_many = _many  # type: ignore[method-assign]
    REC.run_discrete_search_spaced = _section(  # type: ignore[assignment]
        "discrete", REC.run_discrete_search_spaced
    )
    A = REC.AdaptiveVoxelReconstructor
    A._batched_quick_mc = _section("bqmc", A._batched_quick_mc)  # type: ignore[method-assign]
    A.refine_from_candidates = _section(  # type: ignore[method-assign]
        "refine", A.refine_from_candidates
    )
    for cls in (OS.MCOptimizer, OS.CMAOptimizer, OS.RiemannianAdamOptimizer):
        cls.optimize = _section("opt", cls.optimize)  # type: ignore[method-assign]
    OS.MCOptimizer.variance_minimizing_optimize = _section(  # type: ignore[method-assign]
        "variance", OS.MCOptimizer.variance_minimizing_optimize
    )
    for name in ("overlaps_exact", "overlaps_interp"):
        setattr(BS.ScreenEvaluator, name, _screen(getattr(BS.ScreenEvaluator, name)))


def timed_run(rec: Any, verts: Any, phase: int, rng: Any, truth: Optional[np.ndarray]) -> Dict:
    sp = Split()
    STATE["split"] = sp
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            res = rec.reconstruct_voxel(verts, phase, rng=rng)
    finally:
        STATE["split"] = None
    wall = time.perf_counter() - t0
    err = None
    if truth is not None:
        import stats

        err = float(stats.misorientation_deg_cubic(res.orientation[None], truth[None])[0])
    disc = sp.wall.get("discrete", 0.0)
    callers = {}
    for c in CALLERS:
        callers[c] = dict(eval_s=sp.ev_s.get(c, 0.0), eval_n=sp.ev_n.get(c, 0))
    # section walls: the discrete stage is screen + recost + its own machinery
    sections = dict(
        discrete_s=disc,
        quick_mc_s=sp.wall.get("quick_mc", 0.0),
        find_optimal_s=sp.wall.get("find_optimal", 0.0),
        variance_s=sp.wall.get("variance", 0.0),
        refine_s=sp.wall.get("refine", 0.0),
    )
    sections["final_and_refine_machinery_s"] = (
        sections["refine_s"] - sections["find_optimal_s"] - sections["variance_s"]
    )
    sections["rest_s"] = wall - disc - sections["quick_mc_s"] - sections["refine_s"]
    g, loc, tot = rec.last_eval_counts
    return dict(
        wall_s=wall,
        callers=callers,
        sections=sections,
        evaluate_s=sum(sp.ev_s.get(c, 0.0) for c in CALLERS if c != "screen"),
        evaluate_n=sum(sp.ev_n.values()),
        evaluate_fused_n=sp.fused_n,
        screen_s=sp.ev_s.get("screen", 0.0),
        batched_quick_mc=dict(
            s=sp.batched_s, orientations=sp.batched_n, orientations_fused=sp.batched_fused_n
        ),
        evals=dict(global_=g, local=loc, total=tot),
        cost=float(res.cost),
        err_deg=err,
    )


def cases(a: argparse.Namespace) -> Any:
    from s3_step0 import E0Set, ThreeVoxels

    if a.set == "three":
        S3 = ThreeVoxels()
        verts, truth, phase = S3.case()
        rec = S3.reconstructor(a.config)
        yield dict(
            name="three0", rec=rec, verts=verts, phase=phase, truth=truth,
            rng=lambda: np.random.default_rng(1),
        )  # fmt: skip
        return
    SE = E0Set()
    for vpos in range(a.n_voxels):
        verts, truth, phase, rec = SE.load(vpos)
        yield dict(
            name=f"e0_{vpos}", rec=rec, verts=verts, phase=phase, truth=truth,
            rng=lambda vpos=vpos: SE.run_seed(vpos, 0),
        )  # fmt: skip


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--set", choices=["three", "e0"], required=True)
    ap.add_argument("--config", default="Example2.Simulation.config")
    ap.add_argument("--n-voxels", type=int, default=10)
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--allow-busy", action="store_true")
    ap.add_argument(
        "--lib",
        action="store_true",
        help="S2 Step 5: the fused arms turn the library key FusedCost on (no prototype)",
    )
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    info = gate(a.allow_busy)
    (out / f"preflight_search_{a.set}{a.tag}.json").write_text(json.dumps(info, indent=1) + "\n")
    install()
    arms = [x for x in a.arms.split(",") if x]
    runs: List[Dict[str, Any]] = []
    keys: Dict[str, Any] = {}
    for case in cases(a):
        rec = case["rec"]
        sp = rec.setup.search_params
        keys = dict(
            batched_screen=sp.batched_screen,
            cpp_overlap_rule=sp.cpp_overlap_rule,
            cpp_edge_rules=sp.cpp_edge_rules,
            local_optimizer=sp.local_optimizer,
            use_hybrid_optimizer=sp.use_hybrid_optimizer,
            min_sin_eta=getattr(rec, "min_sin_eta", None),
        )
        sp.batched_quick_mc, sp.search_device = False, "cpu"
        sp.fused_cost = False
        STATE["fused"] = False
        timed_run(rec, case["verts"], case["phase"], case["rng"](), None)  # warm-up, untimed
        for rep in range(a.repeats):
            order = arms if rep % 2 == 0 else arms[::-1]
            for arm in order:
                sp.batched_quick_mc, fused = ARMS[arm]
                # --lib: the library key (the reconstructor builds its cost functions with it);
                # otherwise the Step 0 prototype replaces evaluate
                sp.fused_cost = fused and a.lib
                STATE["fused"] = fused and not a.lib
                r = timed_run(rec, case["verts"], case["phase"], case["rng"](), case["truth"])
                r.update(case=case["name"], arm=arm, repeat=rep)
                runs.append(r)
                ev = r["callers"]
                print(
                    f"{case['name']} rep {rep} {arm:10s} wall {r['wall_s']:6.2f}s"
                    f"  evaluate {r['evaluate_s']:6.2f}s n={r['evaluate_n']}"
                    f" (fused {r['evaluate_fused_n']})  screen {r['screen_s']:.2f}s"
                    f"  bqmc {r['batched_quick_mc']['s']:.2f}s"
                    f"  "
                    + " ".join(f"{c}={ev[c]['eval_n']}" for c in CALLERS)
                    + f"  err {r['err_deg']:.3f} cost {r['cost']:.4f}",
                    flush=True,
                )
        sp.batched_quick_mc, STATE["fused"] = False, False
        sp.fused_cost = False
    res = dict(
        label=info["label"],
        lib=a.lib,
        preflight=info,
        torch=torch.__version__,
        torch_threads=torch.get_num_threads(),
        set=a.set,
        config=a.config if a.set == "three" else "E0 (ReconstructQ8 physics)",
        keys=keys,
        arms=arms,
        repeats=a.repeats,
        runs=runs,
    )
    (out / f"search_{a.set}{a.tag}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
