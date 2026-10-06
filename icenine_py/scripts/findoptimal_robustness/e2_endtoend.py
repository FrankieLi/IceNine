#!/usr/bin/env python3
"""E2 end to end: the classifier used inside the reconstruction.

rerank  (a) rank_key hook: at the end of every level the candidates are ordered by the
        classifier's score (probability of "within 3 deg of the truth") instead of the post-MC
        local cost before the top 1/4 is kept. Full reconstruct_voxel runs, seed 0, same images and
        rng as E0. The model of the fold NOT containing the voxel is used (voxel-disjoint).
final   (b) post hoc on the F1 data: among the original answer and the refined CSL relatives
        (Sigma <= 29, the F1 candidate set) the classifier's top score is returned instead of the
        lowest local cost. No new search.

uv run python scripts/findoptimal_robustness/e2_endtoend.py rerank --model GBT --workers 10
uv run python scripts/findoptimal_robustness/e2_endtoend.py final --model GBT --workers 10
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import e2_models as M  # noqa: E402
import features as F  # noqa: E402

_STATE = {}


def _prep(vidx, variant, kind):
    from optimizer_sweep import voxel_context

    W = C.get_worker()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    phase = vctx.voxel.phase
    if "fe" not in _STATE:
        _STATE["fe"] = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase)
    fe = _STATE["fe"]
    fe.set_image(keys)
    return W, vctx, fe


def _model(fold, kind):
    import joblib

    key = (fold, kind)
    if key not in _STATE:
        _STATE[key] = joblib.load(C.CACHE_DIR / "models" / f"fold{fold}_{kind}.joblib")
    return _STATE[key]


def task_rerank(item):
    vidx, vpos, variant, kind, path = item
    t_start = time.time()
    W, vctx, fe = _prep(vidx, variant, kind)
    fold = int(M.fold_of(np.array([vpos]))[0])
    model = _model(fold, kind)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    n_scored = [0]

    def rank_key(level, cands):
        X = np.stack([fe.features(c.orientation, vertices, phase) for c in cands])
        n_scored[0] += len(cands)
        return -M.score(model, X, True)

    W.rec.rank_key = rank_key
    t0 = time.perf_counter()
    with C.quiet():
        res = W.rec.reconstruct_voxel(vertices, phase, rng=C.run_seed(vpos, 0))
    dt = time.perf_counter() - t0
    W.rec.rank_key = None
    g, loc, _ = W.rec.last_eval_counts
    np.savez_compressed(
        path, R_final=np.asarray(res.orientation, float), cost_final=float(res.cost), runtime=dt,
        evals_global=g, evals_local=loc, n_scored=n_scored[0], R_true=vctx.R_true,
    )  # fmt: skip
    return f"rerank {kind} voxel {vidx} {variant} {time.time() - t_start:.0f}s"


def task_final(item):
    vidx, vpos, variant, kind, path = item
    t_start = time.time()
    W, vctx, fe = _prep(vidx, variant, kind)
    fold = int(M.fold_of(np.array([vpos]))[0])
    model = _model(fold, kind)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    d = np.load(C.CACHE_DIR / "e0" / f"v{vidx}_{variant}.npz")
    f1 = np.load(C.CACHE_DIR / "f1" / f"v{vidx}_{variant}.npz")
    out = {}
    for s in range(C.N_SEEDS):
        cands = [d[f"s{s}_R_final"]] + list(f1[f"s{s}_ref_R"])
        costs = [float(d[f"s{s}_cost_final"])] + list(f1[f"s{s}_ref_cost"])
        X = np.stack([fe.features(R, vertices, phase) for R in cands])
        sc = M.score(model, X, True)
        out[f"s{s}_R"] = np.stack(cands)
        out[f"s{s}_cost"] = np.array(costs)
        out[f"s{s}_score"] = sc
    np.savez_compressed(path, R_true=vctx.R_true, **out)
    return f"final {kind} voxel {vidx} {variant} {time.time() - t_start:.0f}s"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("mode", choices=["rerank", "final"])
    ap.add_argument("--model", default="GBT")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    cache = C.CACHE_DIR / f"e2_{a.mode}_{a.model}"
    cache.mkdir(parents=True, exist_ok=True)
    its = []
    for vpos, v in enumerate(vox):
        for var in C.VARIANTS:
            ok = (C.CACHE_DIR / "e0" / f"v{v}_{var}.npz").exists()
            if a.mode == "final":
                ok = ok and (C.CACHE_DIR / "f1" / f"v{v}_{var}.npz").exists()
            out = cache / f"v{v}_{var}.npz"
            if ok and not out.exists():
                its.append((v, vpos, var, a.model, str(out)))
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    C.run_pool(task_rerank if a.mode == "rerank" else task_final, its, a.workers, a.mode)


if __name__ == "__main__":
    main()
