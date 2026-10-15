"""Task E measurement 4: milliseconds per evaluation of each smooth arm against the strict cost.

One image set at a time, one process, preflight-gated (``preflight.require_quiet``; without a
quiet machine the script stops unless ``--contended`` is given, and then labels every number
"contended"). Orientations: 500 within 5 deg of the truth and 500 uniformly random, fixed seed.

* strict: ``VoxelCostFunction.evaluate`` (complete cost, pixel radius 0, all reflections, lenient
  rule), one call per orientation, as the quick MC calls it.
* each arm: ``SmoothCost.evaluate_many`` at the level-0 widths, in batches of 256 (how the ranker
  calls it: the whole list of a level at once) and one orientation per call.
* the one-off build of a ``SmoothCost`` (lit-pixel trees; observed rays for E4), per image set.
Each loop is timed three times; the median is reported.

  uv run python scripts/eval_reduction/smooth_timing.py --set full
  uv run python scripts/eval_reduction/smooth_timing.py --set e0
Writes benchmarks/eval_reduction/E/timing/{timing_<set>.json, preflight_<set>.json}.
"""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import preflight  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "E" / "timing"
ARMS: Dict[str, str] = {
    "S": "S",
    "E1": "E1",
    "E2": "E2",
    "E3": "E3",
    "E4": "E4",
    "E1+E2": "E1+E2",
    "E1+E2+E3": "E1+E2+E3",
    "E2+E3": "E2+E3",
    "E3+E4": "E3+E4",
    "E1+E2+E3+E4": "E1+E2+E3+E4",
}
L0_DIAM = math.radians(5.0)


def med3(fn: Any) -> float:
    ts = []
    for _ in range(3):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def context(which: str) -> Dict[str, Any]:
    import cost_edges as CE
    import cost_edges_run as RUN

    if which == "full":
        RUN.init_full()
        v = 15901
    else:
        RUN.init_e0()
        v = RUN.e0_voxels()[0]
    base, Rt, verts, phase = RUN.case_data(v)
    return dict(
        base=base, Rt=Rt, verts=verts, phase=phase, voxel=v, cf=CE.make_cf(base, "full", "old")
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="full", choices=["full", "e0"])
    ap.add_argument("--contended", action="store_true")
    a = ap.parse_args()
    import torch

    torch.set_num_threads(1)
    from icenine.smooth_cost import SmoothCost, SmoothSpec

    OUT.mkdir(parents=True, exist_ok=True)
    pf = preflight.preflight()
    try:
        preflight.require_quiet(info=pf)
        label = "single-worker"
    except preflight.MachineBusyError as e:
        if not a.contended:
            print(f"not quiet: {e}")
            sys.exit(2)
        label = "contended"
        pf = dict(pf, not_quiet=str(e)[:400])
    (OUT / f"preflight_{a.set}.json").write_text(json.dumps(pf, indent=1) + "\n")
    ctx = context(a.set)
    b = ctx["base"]
    rng = np.random.default_rng(12345)
    near = (
        Rotation.from_rotvec(
            rng.normal(size=(500, 3)) / np.sqrt(3) * np.radians(rng.uniform(0, 5, size=(500, 1)))
        ).as_matrix()
        @ ctx["Rt"]
    )
    far = Rotation.random(500, random_state=12345).as_matrix()
    Rs = np.concatenate([near, far])
    R32 = Rs.astype(np.float32)
    cf, verts, ph = ctx["cf"], ctx["verts"], ctx["phase"]
    res: Dict[str, Any] = dict(
        set=a.set, voxel=ctx["voxel"], label=label, n_orientations=len(Rs), level0_diameter_deg=5.0
    )
    for R in R32[:20]:
        cf.evaluate(R, verts, ph)  # warm-up (frame table)
    t = med3(lambda: [cf.evaluate(R, verts, ph) for R in R32])
    res["strict_ms_per_eval"] = 1000.0 * t / len(Rs)
    t0 = time.perf_counter()
    sc = SmoothCost(
        b.simulator, b.detector_list, b.range_map, b.exp_data, b.sample, b.structure_list,
        eta_limit=b.eta_limit,
    )  # fmt: skip
    sc._ensure_trees()
    res["build_trees_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    sc.observed_rays()
    sc.attached_rays(verts)
    res["build_rays_s"] = time.perf_counter() - t0
    arms: Dict[str, Any] = {}
    for name, s in ARMS.items():
        sp = SmoothSpec.parse(s)
        sc.evaluate_many(Rs[:64], verts, ph, sp, diameter_rad=L0_DIAM)  # warm-up
        tb = med3(lambda: sc.evaluate_many(Rs, verts, ph, sp, diameter_rad=L0_DIAM, chunk=256))
        t1 = med3(
            lambda: [
                sc.evaluate_many(R[None], verts, ph, sp, diameter_rad=L0_DIAM) for R in Rs[:200]
            ]
        )
        arms[name] = dict(
            ms_per_eval_batch256=1000.0 * tb / len(Rs),
            ms_per_eval_single=1000.0 * t1 / 200,
            ratio_batch256_to_strict=(1000.0 * tb / len(Rs)) / res["strict_ms_per_eval"],
        )
        print(name, arms[name], flush=True)
    res["arms"] = arms
    res["preflight_after"] = preflight.preflight()
    (OUT / f"timing_{a.set}.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "preflight_after"}, indent=1))


if __name__ == "__main__":
    main()
