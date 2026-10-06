#!/usr/bin/env python3
"""E2 dataset: candidates with one-pass features and labels.

Per (voxel, variant) task, features (features.py) and labels (error to the truth, Brandon CSL
class) of
  harvest   every candidate of E0 seeds 0 and 1: the post-quick-MC candidates of all 4 levels
            (source 0; level, rank, kept flag recorded) and FindOptimal's per-candidate results
            (source 1);
  synthetic the truth perturbed by residuals (source 2) and every distinct Sigma <= 29 relative of
            the truth perturbed by residuals (source 3); the residual angle is drawn from the
            pooled E0 error of the nearest candidate to the truth at levels 0-2 (levels with a
            candidate within 6 deg), axis isotropic.

  uv run python scripts/findoptimal_robustness/e2_dataset.py run --workers 10
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import csl  # noqa: E402
import features as F  # noqa: E402

N_BASIN_SYN = 24
HARVEST_SEEDS = (0, 1)
_FE = {}


def residual_pool(cache: Path, files, max_files=400) -> np.ndarray:
    out = []
    for f in files[:max_files]:
        d = np.load(f)
        if "unbuildable" in d.files:
            continue
        for s in range(C.N_SEEDS):
            for L in range(3):
                k = f"s{s}_L{L}_qmc_R"
                if k in d.files:
                    e = float(C.err_deg(d[k], d["R_true"]).min())
                    if e < 6.0:
                        out.append(e)
    return np.array(out)


def perturb(R, angles, rng):
    ax = rng.normal(size=(len(angles), 3))
    ax /= np.linalg.norm(ax, axis=1, keepdims=True)
    dR = Rotation.from_rotvec(np.radians(angles)[:, None] * ax).as_matrix()
    return dR @ R  # sample-frame perturbation of the orientation, as the search's local grid


def label(R_true, R):
    err = float(C.err_deg(R, R_true))
    if err < 1.0:
        return err, "basin", 0.0
    c = csl.csl_classify(R_true, R, max_sigma=29)
    return err, c["label"], float(c["deviation"])


def task(item):
    vidx, vpos, variant, path, pool = item
    W = C.get_worker()
    t_start = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    from optimizer_sweep import voxel_context

    vctx = voxel_context(W.ctx, vidx)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    if "fe" not in _FE:
        _FE["fe"] = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase)
    fe = _FE["fe"]
    fe.set_image(keys)
    d = np.load(C.CACHE_DIR / "e0" / f"v{vidx}_{variant}.npz")
    R_true = d["R_true"]
    rng = np.random.default_rng([30_000 + vpos, C.VARIANTS.index(variant)])
    Rs, meta = [], []

    def add(R, source, **kw):
        Rs.append(np.asarray(R, dtype=np.float64))
        meta.append(dict(source=source, **kw))

    for s in HARVEST_SEEDS:
        for L in range(4):
            k = f"s{s}_L{L}_qmc_R"
            if k not in d.files:
                continue
            R, cost, nk = d[k], d[f"s{s}_L{L}_qmc_cost"], int(d[f"s{s}_L{L}_qmc_n_keep"])
            for i in range(len(R)):
                add(R[i], 0, seed=s, level=L, rank=i, kept=int(i < nk), cost=cost[i])
        for i, (R, c) in enumerate(zip(d[f"s{s}_find_R_out"], d[f"s{s}_find_cost"])):
            add(R, 1, seed=s, level=4, rank=i, kept=1, cost=c)
    basin = perturb(R_true, rng.choice(pool, N_BASIN_SYN), rng)
    for R in basin:
        add(R, 2, seed=-1, level=-1, rank=-1, kept=-1, cost=np.nan)
    rel, _ = csl.csl_relatives(R_true, max_sigma=29)
    relp = perturb_each(rel, rng.choice(pool, len(rel)), rng)
    for R in relp:
        add(R, 3, seed=-1, level=-1, rank=-1, kept=-1, cost=np.nan)
    X = np.stack([fe.features(R, vertices, phase) for R in Rs])
    err, lab, dev = [], [], []
    for R in Rs:
        e, l, dv = label(R_true, R)
        err.append(e)
        lab.append(l)
        dev.append(dv)
    out = dict(
        X=X,
        err=np.array(err),
        csl_label=np.array(lab),
        csl_dev=np.array(dev),
        vidx=np.full(len(Rs), vidx),
        vpos=np.full(len(Rs), vpos),
        variant=np.full(len(Rs), C.VARIANTS.index(variant)),
        R=np.stack(Rs),
    )
    for k in ("source", "seed", "level", "rank", "kept", "cost"):
        out[k] = np.array([m[k] for m in meta])
    # reflections the truth shares with its relatives (for the shared-reflection dependence)
    out["truth_features"] = fe.features(R_true, vertices, phase)
    np.savez_compressed(path, **out)
    return f"voxel {vidx} {variant}: {len(Rs)} candidates in {time.time() - t_start:.0f}s"


def perturb_each(Rs, angles, rng):
    ax = rng.normal(size=(len(angles), 3))
    ax /= np.linalg.norm(ax, axis=1, keepdims=True)
    dR = Rotation.from_rotvec(np.radians(angles)[:, None] * ax).as_matrix()
    return dR @ Rs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    e0 = C.CACHE_DIR / "e0"
    files = [e0 / f"v{v}_{var}.npz" for v in vox for var in C.VARIANTS]
    pool = residual_pool(e0, [f for f in files if f.exists()])
    print(
        f"residual pool: n={len(pool)}, quantiles {np.round(np.quantile(pool, [0.1, 0.5, 0.9]), 2)}"
    )
    cache = C.CACHE_DIR / "e2"
    cache.mkdir(parents=True, exist_ok=True)
    its = [
        (v, vpos, var, str(cache / f"v{v}_{var}.npz"), pool)
        for vpos, v in enumerate(vox)
        for var in C.VARIANTS
        if (e0 / f"v{v}_{var}.npz").exists() and not (cache / f"v{v}_{var}.npz").exists()
    ]
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    if its:  # only when work runs: a no-op invocation must not overwrite the committed pool
        np.save(C.OUT_DIR / "residual_pool.npy", pool)
    C.run_pool(task, its, a.workers, "e2")


if __name__ == "__main__":
    main()
