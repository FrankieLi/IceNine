"""S3 Step 4 timings: the lockstep quick MC (BatchedQuickMC) against the sequential one.

``search``     whole no-start searches (``reconstruct_voxel``) of ThreeVoxels voxel 0 or of the first
               E0 clean voxels, in arms ``off`` (sequential quick MC), ``cpu`` (BatchedQuickMC on,
               SearchDevice cpu) and ``mps`` (BatchedQuickMC on, SearchDevice mps). Records wall time,
               the quick-MC time of every level (recorder events ``discrete`` to ``quick_mc``), the
               candidates per level, the evaluations, the batched cost's own split (torch / host
               prep / stage D) and the result. Repeats interleave the arms, so a drifting machine
               load hits all arms alike; an untimed warm-up run of the ``off`` arm touches the images
               first.
``crossover``  the batched cost alone, per orientation, on the CPU compact variant and on the
               fixed-shape variant on MPS, over a range of batch sizes (checks ``GPU_MIN_BATCH``).

Both call ``preflight.require_quiet`` first; ``--allow-busy`` records the preflight and continues
(timings then labelled "contended"). Single torch thread. The key-on arms draw a different random
stream than the key-off arm (draw option B), so their searches and evaluation counts differ; the
quick-MC time is compared per evaluation as well as per level.

  uv run python scripts/speed_parity/s3_timing.py search --set three --out DIR --allow-busy
  uv run python scripts/speed_parity/s3_timing.py search --set e0 --n-voxels 10 --out DIR ...
  uv run python scripts/speed_parity/s3_timing.py crossover --out DIR --allow-busy
"""

import argparse
import contextlib
import io
import json
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(ICENINE_PY))
sys.path.insert(0, str(ICENINE_PY / "scripts" / "common"))
sys.path.insert(0, str(HERE))
torch.set_num_threads(1)  # single-worker timings

import preflight  # noqa: E402

from s3_step0 import E0Set, ThreeVoxels, orientation_pool  # noqa: E402

from icenine import batched_cost as BC  # noqa: E402

ARMS = {"off": (False, "cpu"), "cpu": (True, "cpu"), "mps": (True, "mps")}


def gate(a: argparse.Namespace) -> Dict[str, Any]:
    info = preflight.preflight()
    try:
        preflight.require_quiet(info=info)
        info["quiet"] = True
    except preflight.MachineBusyError as exc:
        if not a.allow_busy:
            raise
        info["quiet"] = False
        info["not_quiet_reason"] = str(exc)
        print("WARNING: machine not quiet, timings are 'contended':", exc, flush=True)
    info["label"] = "single-worker" if info["quiet"] else "contended"
    return info


def timed_run(
    rec: Any, verts: Any, phase: int, rng: np.random.Generator, R_true: Optional[np.ndarray]
) -> Dict[str, Any]:
    """One reconstruct_voxel with the quick-MC time and candidate count of each level."""
    stamp: Dict[str, Any] = dict(t=None, n=None)
    levels: List[Dict[str, Any]] = []

    def recorder(name: str, d: Dict[str, Any]) -> None:
        if name == "discrete":
            stamp["t"] = time.perf_counter()
            stamp["n"] = len(d["R"])
        elif name == "quick_mc":
            levels.append(dict(n=stamp["n"], quick_mc_s=time.perf_counter() - stamp["t"]))

    rec.recorder = recorder
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            res = rec.reconstruct_voxel(verts, phase, rng=rng)
    finally:
        rec.recorder = None
    wall = time.perf_counter() - t0
    g, loc, tot = rec.last_eval_counts
    err = None
    if R_true is not None:
        from scipy.spatial.transform import Rotation

        err = float(np.degrees(Rotation.from_matrix(res.orientation @ R_true.T).magnitude()))
    evs = list(rec._quick_mc_evaluators.values())
    split = None
    if evs and evs[0] is not None:
        e = evs[0]
        split = dict(timing=dict(e.timing), counts=dict(e.counts))
    n_cand = sum(lv["n"] for lv in levels)
    q_s = sum(lv["quick_mc_s"] for lv in levels)
    return dict(
        wall_s=wall,
        levels=levels,
        quick_mc_s=q_s,
        quick_mc_evals=12 * n_cand,
        quick_mc_us_per_eval=1e6 * q_s / (12 * n_cand) if n_cand else None,
        evals=dict(global_=g, local=loc, total=tot),
        cost=float(res.cost),
        err_deg=err,
        batched_split=split,
    )


def cmd_search(a: argparse.Namespace) -> None:
    out = Path(a.out).resolve()  # before the example set chdirs
    out.mkdir(parents=True, exist_ok=True)
    info = gate(a)
    arms = [x for x in a.arms.split(",") if x]
    if "mps" in arms and not BC.device_available("mps"):
        arms.remove("mps")
        print("MPS not available: skipping the mps arm", flush=True)

    def make_cases() -> Any:
        """Yield one case at a time: E0Set.load re-attaches the shared reconstructor's images."""
        if a.set == "three":
            S3 = ThreeVoxels()
            verts, truth, phase = S3.case()
            yield dict(
                name=f"ThreeVoxels voxel 0 ({a.config})",
                rec=S3.reconstructor(a.config),
                verts=verts,
                phase=phase,
                truth=truth,
                seeds=list(range(a.seeds)),
                rng=lambda seed: np.random.default_rng(seed),
            )
            return
        SE = E0Set()
        for vpos in range(a.n_voxels):
            verts, truth, phase, rec = SE.load(vpos)
            yield dict(
                name=f"E0 clean voxel {vpos}",
                rec=rec,
                verts=verts,
                phase=phase,
                truth=truth,
                seeds=list(range(a.seeds)),
                rng=lambda seed, vpos=vpos: SE.run_seed(vpos, seed),
            )

    runs: List[Dict[str, Any]] = []
    for case in make_cases():
        rec = case["rec"]
        sp = rec.setup.search_params
        sp.batched_quick_mc, sp.search_device = False, "cpu"
        print(f"[{case['name']}] warm-up (untimed, key off)", flush=True)
        timed_run(rec, case["verts"], case["phase"], case["rng"](case["seeds"][0]), None)
        if "mps" in arms:  # first MPS use compiles kernels: one untimed key-on mps run
            sp.batched_quick_mc, sp.search_device = True, "mps"
            timed_run(rec, case["verts"], case["phase"], case["rng"](case["seeds"][0]), None)
        for rep in range(a.repeats):
            for seed in case["seeds"]:
                for arm in arms:
                    sp.batched_quick_mc, sp.search_device = ARMS[arm]
                    r = timed_run(
                        rec, case["verts"], case["phase"], case["rng"](seed), case["truth"]
                    )
                    r.update(case=case["name"], arm=arm, repeat=rep, seed=seed)
                    runs.append(r)
                    print(
                        f"  {case['name']} rep {rep} seed {seed} {arm:>3}: wall {r['wall_s']:.2f}s,"
                        f" quick MC {r['quick_mc_s']:.2f}s ({r['quick_mc_us_per_eval']:.0f} us/eval,"
                        f" cand/level {[lv['n'] for lv in r['levels']]}), evals {r['evals']['total']},"
                        f" err {r['err_deg']:.3f} deg",
                        flush=True,
                    )
        sp.batched_quick_mc, sp.search_device = False, "cpu"
    res = dict(
        preflight=info,
        torch=torch.__version__,
        torch_threads=torch.get_num_threads(),
        label=info["label"],
        set=a.set,
        arms=arms,
        repeats=a.repeats,
        runs=runs,
    )
    (out / f"timing_search_{a.set}{a.tag}.json").write_text(
        json.dumps(res, indent=1, default=float)
    )
    (out / f"preflight_search_{a.set}{a.tag}.json").write_text(json.dumps(info, indent=1))


def cmd_crossover(a: argparse.Namespace) -> None:
    out = Path(a.out).resolve()  # before the example set chdirs
    out.mkdir(parents=True, exist_ok=True)
    info = gate(a)
    S = E0Set()
    verts, truth, phase, rec = S.load(a.vpos)
    cf = rec._make_local_cost_fn()
    pool = orientation_pool(truth, 1024)
    sizes = [int(x) for x in a.sizes.split(",")]
    devs = ["cpu"] + (["mps"] if BC.device_available("mps") else [])
    evs = {d: BC.BatchedCostEvaluator(cf, d) for d in devs}
    rows: List[Dict[str, Any]] = []
    for B in sizes:
        row: Dict[str, Any] = dict(B=B)
        for label, dev, variant in (
            ("cpu_compact", "cpu", "compact"),
            ("mps_fixed", "mps", "fixed"),
            ("cpu_fixed", "cpu", "fixed"),
        ):
            if dev not in evs:
                continue
            ev = evs[dev]
            ev.evaluate_many(pool[:B], verts, phase, variant=variant)  # warm-up
            ts = []
            for _ in range(a.reps):
                t0 = time.perf_counter()
                ev.evaluate_many(pool[:B], verts, phase, variant=variant)
                ts.append(time.perf_counter() - t0)
            row[label + "_us_per_orient"] = 1e6 * statistics.median(ts) / B
        rows.append(row)
        print(row, flush=True)
    cross = [
        r["B"]
        for r in rows
        if "mps_fixed_us_per_orient" in r
        and r["mps_fixed_us_per_orient"] < r["cpu_compact_us_per_orient"]
    ]
    res = dict(
        preflight=info,
        label=info["label"],
        rows=rows,
        gpu_min_batch=BC.GPU_MIN_BATCH,
        mps_faster_from_B=min(cross) if cross else None,
        voxel=a.vpos,
    )
    (out / "timing_crossover.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "preflight_crossover.json").write_text(json.dumps(info, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("search")
    s.add_argument("--set", choices=["three", "e0"], required=True)
    s.add_argument("--config", default="ReconstructQ8.config")
    s.add_argument("--n-voxels", type=int, default=10)
    s.add_argument("--seeds", type=int, default=1)
    s.add_argument("--repeats", type=int, default=2)
    s.add_argument("--arms", default="off,cpu,mps")
    s.add_argument("--tag", default="")
    s.add_argument("--out", required=True)
    s.add_argument("--allow-busy", action="store_true")
    c = sub.add_parser("crossover")
    c.add_argument("--vpos", type=int, default=0)
    c.add_argument("--sizes", default="1,2,4,8,16,32,48,64,96,128,192,256,512")
    c.add_argument("--reps", type=int, default=9)
    c.add_argument("--out", required=True)
    c.add_argument("--allow-busy", action="store_true")
    a = ap.parse_args()
    {"search": cmd_search, "crossover": cmd_crossover}[a.cmd](a)


if __name__ == "__main__":
    main()
