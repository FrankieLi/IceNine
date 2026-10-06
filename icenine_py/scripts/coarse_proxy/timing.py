#!/usr/bin/env python3
"""Cost accounting of the proxy on ONE worker: OMP / MKL / torch threads = 1, nothing else running.

Median over 1000 harvested candidates (25 random cases x 40 candidates, seeded) of
  * one Q_max-8 local evaluation (pixel radius 0): the unit of "evaluation equivalents";
  * one Q_max-5 pixel-radius-3 global evaluation (the coarse level-0 cost);
  * the F-lowQ passes for Q = 4 and 5 (geometry pass only, and including their two low-Q costs);
  * the full 62-feature E2 pass (geometry + two Q_max-8 costs);
  * one model prediction (GBT, batch of 200, per candidate).
Output benchmarks/coarse_proxy/timing.json.

  uv run python scripts/coarse_proxy/timing.py --set lowq5+c8 --target reg
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

C, M = B.C, B.M
import features as F  # noqa: E402
import models as MD  # noqa: E402

N_CASES, PER_CASE = 25, 40


def time_per_call(fn: Callable[[], Any], reps: int = 1) -> float:
    """Seconds per call, averaged over `reps` calls (the medians are taken by the caller)."""
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t0) / reps


def main() -> None:
    import torch
    import joblib
    from icenine.cost_functions import VoxelCostFunction
    from optimizer_sweep import voxel_context

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--set", default="lowq5+c8")
    ap.add_argument("--target", default="reg")
    a = ap.parse_args()
    torch.set_num_threads(1)
    C.init_worker(C.worker_args())
    W = C.get_worker()
    D = B.load_e2_aligned()
    D.update(B.load_per_task(B.CACHE / "lowq", ["X4", "X5"]))
    rng = np.random.default_rng(0)
    ts = B.tasks()
    cases = rng.choice(len(ts), N_CASES, replace=False)
    t: Dict[str, list] = {
        k: []
        for k in ("local_q8", "global_q5_r3", "lowq4_geom", "lowq4_full", "lowq5_geom",
                  "lowq5_full", "e2_full")
    }  # fmt: skip
    sizes = [len(np.load(C.CACHE_DIR / "e2" / f"v{v}_{var}.npz")["err"]) for v, _, var in ts]
    off = np.concatenate([[0], np.cumsum(sizes)])
    for ci in cases:
        v, vpos, var = ts[ci]
        keys = np.load(C.CACHE_DIR / "images" / f"v{v}_{var}.npz")["keys"]
        C.attach(keys)
        vctx = voxel_context(W.ctx, v)
        vert, ph = vctx.vertices, vctx.voxel.phase
        fe = {q: F.FeatureExtractor(W.local_fn, W.ctx.geo, ph, q_max=q) for q in (4.0, 5.0)}
        fe_full = F.FeatureExtractor(W.local_fn, W.ctx.geo, ph)
        for f in [*fe.values(), fe_full]:
            f.set_image(keys)
        g5 = VoxelCostFunction(
            simulator=W.local_fn.simulator, detector_list=W.local_fn.detector_list,
            range_map=W.local_fn.range_map, exp_data=W.local_fn.exp_data, sample=W.local_fn.sample,
            structure_list=W.local_fn.structure_list, mode="hard", eta_limit=W.local_fn.eta_limit,
            pixel_radius=3, max_q=5.0, min_sin_eta=W.local_fn.min_sin_eta,
        )  # fmt: skip
        idx = np.nonzero(D["source"][off[ci] : off[ci + 1]] == 0)[0]
        for i in rng.choice(idx, PER_CASE, replace=False):
            R = D["R"][off[ci] + i]
            Rf = R.astype(np.float32)
            for f in [*fe.values(), fe_full]:
                f.features(R, vert, ph)  # warm-up
            t["local_q8"].append(time_per_call(lambda: W.local_fn.evaluate(Rf, vert, ph)))
            t["global_q5_r3"].append(time_per_call(lambda: g5.evaluate(Rf, vert, ph)))
            for q in (4, 5):
                f = fe[float(q)]
                t[f"lowq{q}_geom"].append(
                    time_per_call(lambda: f.features(R, vert, ph, with_cost=False))
                )
                t[f"lowq{q}_full"].append(time_per_call(lambda: f.features(R, vert, ph)))
            t["e2_full"].append(time_per_call(lambda: fe_full.features(R, vert, ph)))
    med = {k: float(np.median(v)) for k, v in t.items()}
    # model prediction, batch 200, per candidate
    X = MD.feature_matrix(a.set, D)[:200]
    m = joblib.load(MD.model_path(a.set, a.target, 0))
    reps = [time_per_call(lambda: MD.predict_score(a.target, m, X)) for _ in range(20)]
    med["model_predict_per_candidate"] = float(np.median(reps)) / 200
    unit = med["local_q8"]
    out = dict(
        n_candidates=len(t["local_q8"]), seconds_median=med,
        equivalents_of_local_q8={k: v / unit for k, v in med.items()},
        set=a.set, target=a.target, threads=1,
    )  # fmt: skip
    B.OUT.mkdir(parents=True, exist_ok=True)
    C.save_json(B.OUT / "timing.json", out)
    for k, v in med.items():
        print(f"{k:28s} {v * 1e3:8.3f} ms   {v / unit:7.2f} local-Q8 evaluations")


if __name__ == "__main__":
    main()
