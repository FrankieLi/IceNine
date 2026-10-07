#!/usr/bin/env python3
"""
Phase A1 (finisher/MC study): the cost landscape around the truth, clean vs realistic.

Cases are exactly the T5 cases (scripts/finisher_diagnosis/diagnose.build_items: 200 H3 + 100 H0).
For each case and BOTH variants of the same voxel/radius/direction (variant_batches yields clean and
realistic), the VoxelCostFunction is evaluated at the truth rotated by r about N_DIR random axes for
a radius grid r (including the finisher's own final MC step), plus at the T5 finisher result. Raw
costs are cached per task; scripts/cost_sensitivity/summary.py analyses them.

Usage (from icenine_py/):
  uv run python scripts/cost_sensitivity/landscape.py pilot
  uv run python scripts/cost_sensitivity/landscape.py run --workers 10
"""

import argparse
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "finisher_diagnosis"))
sys.path.insert(0, str(HERE.parent / "nn_hybrid"))
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent))

import diagnose as D  # noqa: E402
import findoptimal_sweep as fs  # noqa: E402
import run as nnrun  # noqa: E402

CACHE_DIR = HERE / "cache"
N_DIR = 400
RADII_DEG = np.array([0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1])
VARIANTS = ["clean", "all"]  # "all" = realistic


def unit_directions(n: int, seed: int) -> np.ndarray:
    """n random unit vectors (n, 3), deterministic in seed."""
    ax = np.random.default_rng(seed).normal(size=(n, 3))
    return ax / np.linalg.norm(ax, axis=1, keepdims=True)


def rotate_about(R: np.ndarray, radii_deg: np.ndarray, dirs: np.ndarray) -> np.ndarray:
    """R rotated (left-multiplied) by radius x direction: (len(radii), len(dirs), 3, 3)."""
    rv = np.radians(radii_deg)[:, None, None] * dirs[None]
    rots = Rotation.from_rotvec(rv.reshape(-1, 3)).as_matrix()
    return (rots @ R).reshape(len(radii_deg), len(dirs), 3, 3)


def task(item: Tuple[Any, ...]) -> Tuple[int, int, float]:
    vidx, vpos, ri, dirs, pipe, path, ref_nroi, ref_fail, raw_R, _raw_start = item
    assert fs._W is not None
    ctx = fs._W.ctx
    a = ctx.args
    t0 = time.time()
    box, step = D.final_box(fs._W.rec)
    radii = np.concatenate([RADII_DEG, [math.degrees(step)]])
    res: Dict[Tuple[int, int], Dict[str, np.ndarray]] = {}
    for vb in nnrun.variant_batches(ctx, vidx, vpos, ri, ref_nroi, ref_fail, VARIANTS):
        for jj, j in enumerate(dirs):
            fs.attach_images(nnrun.case_keys(ctx, vb, j))
            R_true = np.asarray(vb.vctx.R_true, dtype=np.float64)
            seed = nnrun.b_seed(a, vpos, ri, j, vb.vi)
            u = unit_directions(N_DIR, seed + 31)
            pts = rotate_about(R_true, radii, u)
            cost = np.array(
                [[D._cost(pts[k, m], vb.vctx) for m in range(N_DIR)] for k in range(len(radii))]
            )
            res[(vb.vi, jj)] = dict(
                cost=cost,
                dir=u,
                cost_true=np.array(D._cost(R_true, vb.vctx)),
                cost_res=np.array(D._cost(raw_R[j], vb.vctx)),
                ang_res=np.array(D.angle_deg(raw_R[j], R_true)),
            )
    out: Dict[str, np.ndarray] = dict(
        dirs=np.array(dirs), radii=radii, step_deg=np.array(math.degrees(step)),
        vidx=np.array(vidx), ri=np.array(ri), pipe=np.array(pipe),
    )  # fmt: skip
    for k in res[(0, 0)]:
        out[k] = np.stack(
            [np.stack([res[(vi, jj)][k] for jj in range(len(dirs))]) for vi in range(len(VARIANTS))]
        )  # (variant, case, ...)
    np.savez(path, **out)
    return vidx, ri, time.time() - t0


def build_items(cache: Path, only_first: int) -> Tuple[List[Tuple[Any, ...]], Dict[str, Any]]:
    items, wargs = D.build_items(cache, only_first=only_first)
    return items, wargs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run"])
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()
    import multiprocessing as mp

    cache = CACHE_DIR / args.cmd
    cache.mkdir(parents=True, exist_ok=True)
    items, wargs = build_items(cache, 2 if args.cmd == "pilot" else 0)
    todo = [it for it in items if not Path(it[5]).exists()]
    print(
        f"{len(items)} tasks ({len(items) - len(todo)} cached), {args.workers} workers", flush=True
    )
    t0 = time.time()
    with mp.get_context("spawn").Pool(
        min(args.workers, max(len(todo), 1)), initializer=fs.init_worker, initargs=(wargs,)
    ) as pool:
        for k, (vidx, ri, secs) in enumerate(pool.imap_unordered(task, todo), 1):
            print(
                f"  [{k}/{len(todo)}] voxel {vidx} r#{ri} {secs:.0f}s ({time.time()-t0:.0f}s)",
                flush=True,
            )
    print(f"finished; wall {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
