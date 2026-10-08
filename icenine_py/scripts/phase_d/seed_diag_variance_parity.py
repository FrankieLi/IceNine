"""C++ vs Python VarianceMinimizing stage on one voxel of the Phase D full-sample images.

Subcommands (run from icenine_py/):
  costs --cpp-log LOG --voxel V     cost Python assigns to the start orientations that a patched C++
                                    build printed (lines "VM <9 matrix entries> cost <c>") vs C++
  trace --cpp-log LOG --voxel V     Python variance stage started at the same orientation, N seeds:
                                    steps, runs, fraction of runs below the variance threshold, the
                                    restart start-cost distribution (C++ log: same quantities)
The C++ log comes from a scratch copy of Src/ with debug prints (not committed); see
MIGRATION_HISTORY "C++ vs Python variance stage".
"""

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import seed_diag as SD  # noqa: E402


def parse_cpp(log: Path) -> Tuple[List[Tuple[np.ndarray, float]], List[Dict[str, float]]]:
    """VM lines (start matrix, cost) and VT lines (per-run radius, steps, cost, var, ext) of the
    first voxel in the log (the log may hold several copies of the voxel, in order)."""
    vm, vt = [], []
    for line in log.read_text().splitlines():
        if line.startswith("VM "):
            t = line.split()
            vm.append((np.array([float(x) for x in t[1:10]]).reshape(3, 3), float(t[11])))
        elif line.startswith("VT "):
            t = line.split()
            vt.append(
                dict(
                    radius_deg=float(t[3]),
                    steps=int(t[5]),
                    cost=float(t[7]),
                    var=float(t[11]),
                    ext=int(t[13]),
                    total=int(t[15]),
                    max=int(t[17]),
                )
            )
    return vm, vt


def cpp_stages(vt: List[Dict[str, float]]) -> List[List[Dict[str, float]]]:
    """Split the VT lines into stages (a new stage starts when `total` decreases)."""
    out: List[List[Dict[str, float]]] = []
    prev = 1 << 60
    for r in vt:
        if r["total"] <= prev:
            out.append([])
        out[-1].append(r)
        prev = r["total"]
    return out


def setup(voxel: int):
    os.chdir(SD.EX)
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile

    cfg = ConfigFile.from_file(str(SD.HERE / "configs" / "ReconstructPhaseD_mc_clean.config"))
    stack = SD.EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
    s = RC.setup_reconstruction(cfg, exp_data=ExperimentalData.from_binary_memmap(stack, 180, 2))
    rec = RC.AdaptiveVoxelReconstructor(s)
    v = s.sample.get_mic().voxels[voxel]
    verts = RC._get_voxel_vertices(v)
    local = rec._make_local_cost_fn()
    truth = MicFile.read(SD.MIC).voxels[voxel].orientation
    return rec, local, verts, v.phase, np.asarray(truth, dtype=np.float64)


def legacy_variance_fn(rev: str = "HEAD"):
    """The variance_minimizing_optimize of git revision `rev` (before the restart fix)."""
    import subprocess

    src = subprocess.run(
        ["git", "show", f"{rev}:icenine_py/icenine/orientation_search.py"],
        capture_output=True,
        text=True,
        check=True,
        cwd=str(HERE.parents[1]),
    ).stdout
    ns: Dict[str, Any] = {
        "__name__": "icenine._legacy_orientation_search",
        "__package__": "icenine",
    }
    exec(compile(src, "legacy_orientation_search.py", "exec"), ns)
    return ns["MCOptimizer"].variance_minimizing_optimize


def final_box_width(rec: Any) -> float:
    """Box width of the variance stage: as `_finish` computes it from the final diameter."""
    p = rec.params
    diameter = p.local_grid_radius / 1.5 ** (p.max_local_resolution + 1)
    return max(diameter / 3.0, math.radians(0.2)) / (2**p.min_local_resolution)


def cmd_costs(a: argparse.Namespace) -> None:
    rec, local, verts, phase, truth = setup(a.voxel)
    vm, _ = parse_cpp(Path(a.cpp_log))
    rows = []
    for R, c_cpp in vm[: a.n]:
        c_py = local.evaluate(R.astype(np.float32), verts, phase).cost
        rows.append((c_cpp, float(c_py)))
    arr = np.array(rows)
    d = arr[:, 1] - arr[:, 0]
    res = dict(
        voxel=a.voxel,
        n=len(rows),
        max_abs_diff=float(np.abs(d).max()),
        mean_abs_diff=float(np.abs(d).mean()),
        n_exact=int((np.abs(d) < 1e-6).sum()),
        cpp_cost_range=[float(arr[:, 0].min()), float(arr[:, 0].max())],
        first5=[[round(x, 6), round(y, 6)] for x, y in rows[:5]],
    )
    print(json.dumps(res, indent=1))
    Path(a.out).write_text(json.dumps(res, indent=1) + "\n")


def cmd_trace(a: argparse.Namespace) -> None:
    rec, local, verts, phase, truth = setup(a.voxel)
    vm, vt = parse_cpp(Path(a.cpp_log))
    R0 = vm[0][0]
    box = final_box_width(rec)
    out: Dict[str, Any] = dict(voxel=a.voxel, box_deg=math.degrees(box), n_seeds=a.n_seeds)
    stages = cpp_stages(vt)
    out["cpp"] = [summarise([r for r in st]) for st in stages]
    runs = []
    for seed in range(a.n_seeds):
        mc, _ = rec._make_optimizers(local, verts, phase, np.random.default_rng(seed))
        log: List[Dict[str, float]] = []
        total = [0]
        orig = mc._zero_temp_with_variance

        def traced(init, radius, n_steps, _o=orig, _l=log, _t=total):  # type: ignore
            r = _o(init, radius, n_steps)
            _t[0] += n_steps
            _l.append(
                dict(
                    radius_deg=math.degrees(radius),
                    steps=n_steps,
                    cost=float(r[1]),
                    var=float(r[2]),
                    ext=int(r[2] > 0.02**2),
                    total=_t[0],
                    start_cost=(
                        float(local.evaluate(init, verts, phase).cost) if a.start_costs else -1
                    ),
                )
            )
            if _t[0] >= a.max_total:
                raise StopIteration
            return r

        mc._zero_temp_with_variance = traced  # type: ignore[method-assign]
        if a.legacy:
            import types

            mc.variance_minimizing_optimize = types.MethodType(  # type: ignore[method-assign]
                legacy_variance_fn(a.legacy), mc
            )
        try:
            res = mc.variance_minimizing_optimize(
                R0.astype(np.float32), box, rec.params.max_mc_steps, 2, 0.0, 0.02**2
            )
            capped = False
        except StopIteration:
            capped = True
        s = summarise(log)
        s["capped"] = capped
        runs.append(s)
        print(seed, s, flush=True)
    out["python"] = runs
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")


def summarise(log: List[Dict[str, float]]) -> Dict[str, Any]:
    v = np.array([r["var"] for r in log])
    c = np.array([r["cost"] for r in log])
    return dict(
        n_runs=len(log),
        total_steps=int(log[-1]["total"]) if log else 0,
        frac_ext=float((v > 0.02**2).mean()),
        mean_run_cost=float(c.mean()),
        median_run_cost=float(np.median(c)),
        q25_q75_run_cost=[float(x) for x in np.quantile(c, [0.25, 0.75])],
        var_median=float(np.median(v)),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("costs")
    c.add_argument("--cpp-log", required=True)
    c.add_argument("--voxel", type=int, required=True)
    c.add_argument("--n", type=int, default=160)
    c.add_argument("--out", required=True)
    t = sub.add_parser("trace")
    t.add_argument("--cpp-log", required=True)
    t.add_argument("--voxel", type=int, required=True)
    t.add_argument("--n-seeds", type=int, default=8)
    t.add_argument("--max-total", type=int, default=30000)
    t.add_argument("--start-costs", action="store_true")
    t.add_argument("--legacy", default="", help="git revision whose variance stage to run")
    t.add_argument("--out", required=True)
    a = ap.parse_args()
    {"costs": cmd_costs, "trace": cmd_trace}[a.cmd](a)


if __name__ == "__main__":
    main()
