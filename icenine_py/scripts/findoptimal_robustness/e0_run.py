#!/usr/bin/env python3
"""E0: instrumented full reconstructions (reconstruct_voxel from scratch) of the study's voxels.

One task = one (voxel, variant): the case images are built once, saved (compressed pixel keys) and
N_SEEDS reconstruct runs are recorded with the recorder hook. Resumable: a finished task leaves
cache/e0/v{voxel}_{variant}.npz. Seed 0 uses experiment A's rng, so for the sweep's first 50
voxels it reproduces findoptimal_sweep experiment A exactly.

  uv run python scripts/findoptimal_robustness/e0_run.py select
  uv run python scripts/findoptimal_robustness/e0_run.py pilot
  uv run python scripts/findoptimal_robustness/e0_run.py run --workers 10 [--variants clean all]
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

N_SPARE = 8  # extra voxels queued in case some cannot be built


def task(item):
    vidx, vpos, variant, n_seeds, path = item
    W = C.get_worker()
    t_start = time.time()
    case = C.build_case(vidx, vpos, variant)
    if case is None:
        np.savez(path, unbuildable=True)
        return f"voxel {vidx} {variant}: unbuildable"
    keys, R_true, vctx = case
    img_dir = C.CACHE_DIR / "images"
    img_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(img_dir / f"v{vidx}_{variant}.npz", keys=keys)
    C.attach(keys)
    rec = W.rec
    vertices, phase = vctx.vertices, vctx.voxel.phase
    cost_true = W.local_fn.evaluate(R_true.astype(np.float32), vertices, phase).cost
    out = dict(R_true=R_true, cost_true=cost_true, n_pixels=len(keys), vidx=vidx, vpos=vpos)
    for s in range(n_seeds):
        recorder = C.Recorder()
        rec.recorder = recorder
        t0 = time.perf_counter()
        with C.quiet():
            res = rec.reconstruct_voxel(vertices, phase, rng=C.run_seed(vpos, s))
        dt = time.perf_counter() - t0
        rec.recorder = None
        g, loc, _ = rec.last_eval_counts
        out[f"s{s}_R_final"] = np.asarray(res.orientation, dtype=np.float64)
        out[f"s{s}_cost_final"] = float(res.cost)
        out[f"s{s}_runtime"] = dt
        out[f"s{s}_evals_global"] = g
        out[f"s{s}_evals_local"] = loc
        for k, v in recorder.to_arrays().items():
            out[f"s{s}_{k}"] = v
    np.savez_compressed(path, **out)
    return f"voxel {vidx} {variant} done in {time.time() - t_start:.0f}s"


def items(info, variants, n_seeds, n_vox):
    cache = C.CACHE_DIR / "e0"
    cache.mkdir(parents=True, exist_ok=True)
    vox = [int(v) for v in info["voxel_indices"]][:n_vox]
    return [
        (v, vpos, var, n_seeds, str(cache / f"v{v}_{var}.npz"))
        for vpos, v in enumerate(vox)
        for var in variants
    ]


def load_info():
    p = C.OUT_DIR / "voxels.npz"
    if not p.exists():
        raise SystemExit("run `select` first")
    return dict(np.load(p))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["select", "pilot", "run"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--variants", nargs="*", default=C.VARIANTS)
    ap.add_argument("--n-seeds", type=int, default=C.N_SEEDS)
    ap.add_argument("--n-voxels", type=int, default=C.N_VOXELS + N_SPARE)
    a = ap.parse_args()
    C.OUT_DIR.mkdir(parents=True, exist_ok=True)
    if a.cmd == "select":
        info = C.select_voxels(C.N_VOXELS + N_SPARE)
        np.savez(C.OUT_DIR / "voxels.npz", **info)
        print("saved voxels.npz;", len(info["voxel_indices"]), "voxels")
        return
    info = load_info()
    if a.cmd == "pilot":
        its = items(info, a.variants, a.n_seeds, 4)
        for it in its:
            Path(it[4]).unlink(missing_ok=True)
        t0 = time.time()
        C.run_pool(task, its, a.workers, "pilot")
        print(f"pilot wall {time.time() - t0:.0f}s for {len(its)} tasks x {a.n_seeds} seeds")
        return
    its = items(info, a.variants, a.n_seeds, a.n_voxels)
    todo = [it for it in its if not Path(it[4]).exists()]
    print(f"{len(its)} tasks, {len(todo)} to run, {a.workers} workers", flush=True)
    C.run_pool(task, todo, a.workers, "e0")


if __name__ == "__main__":
    main()
