#!/usr/bin/env python3
"""E0: instrumented full reconstructions (reconstruct_voxel from scratch) of the study's voxels.

One task = one (voxel, variant): the case images are built once, saved (compressed pixel keys) and
N_SEEDS reconstruct runs are recorded with the recorder hook. Resumable: a finished task leaves
cache/e0{tag}/v{voxel}_{variant}.npz. Seed 0 uses experiment A's rng, so for the sweep's first 50
voxels it reproduces findoptimal_sweep experiment A exactly.

  uv run python scripts/findoptimal_robustness/e0_run.py select
  uv run python scripts/findoptimal_robustness/e0_run.py pilot
  uv run python scripts/findoptimal_robustness/e0_run.py run --workers 10 [--variants clean all]

S3 A/B (the lockstep quick MC, config key BatchedQuickMC, draw option B): ``--batched-quick-mc``
turns the key on (``--search-device`` picks the device of its batched cost) and ``--tag _s3on``
keeps its results in cache/e0_s3on, apart from the key-off arm's cache/e0. A case's detector images
are built once and saved in cache/images; a later run (the other arm) reuses them.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

N_SPARE = 8  # extra voxels queued in case some cannot be built


def get_case(vidx, vpos, variant):
    """C.build_case, but the pixel keys of a case already built (cache/images) are reused: only the
    voxel context is rebuilt. Returns (keys, R_true, vctx) or None (unbuildable)."""
    img = C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz"
    if img.exists():
        vctx = C.osw.voxel_context(C.get_worker().ctx, vidx)
        return np.load(img)["keys"], vctx.R_true, vctx
    case = C.build_case(vidx, vpos, variant)
    if case is not None:
        img.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(img, keys=case[0])
    return case


def task(item):
    vidx, vpos, variant, n_seeds, path, opts = item
    W = C.get_worker()
    t_start = time.time()
    case = get_case(vidx, vpos, variant)
    if case is None:
        np.savez(path, unbuildable=True)
        return f"voxel {vidx} {variant}: unbuildable"
    keys, R_true, vctx = case
    C.attach(keys)
    rec = W.rec
    sp = rec.setup.search_params
    sp.batched_quick_mc = bool(opts["batched_quick_mc"])
    sp.search_device = str(opts["search_device"])
    vertices, phase = vctx.vertices, vctx.voxel.phase
    cost_true = W.local_fn.evaluate(R_true.astype(np.float32), vertices, phase).cost
    out = dict(R_true=R_true, cost_true=cost_true, n_pixels=len(keys), vidx=vidx, vpos=vpos)
    out["batched_quick_mc"] = sp.batched_quick_mc
    out["search_device"] = sp.search_device
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


def items(info, variants, n_seeds, n_vox, tag="", opts=None):
    cache = C.CACHE_DIR / f"e0{tag}"
    cache.mkdir(parents=True, exist_ok=True)
    opts = opts or dict(batched_quick_mc=False, search_device="cpu")
    vox = [int(v) for v in info["voxel_indices"]][:n_vox]
    return [
        (v, vpos, var, n_seeds, str(cache / f"v{v}_{var}.npz"), opts)
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
    ap.add_argument("--batched-quick-mc", action="store_true", help="S3: BatchedQuickMC on")
    ap.add_argument("--search-device", default="cpu", choices=["cpu", "mps", "cuda"])
    ap.add_argument("--tag", default="", help="results go to cache/e0<tag> (e.g. _s3on)")
    a = ap.parse_args()
    opts = dict(batched_quick_mc=a.batched_quick_mc, search_device=a.search_device)
    C.OUT_DIR.mkdir(parents=True, exist_ok=True)
    if a.cmd == "select":
        info = C.select_voxels(C.N_VOXELS + N_SPARE)
        np.savez(C.OUT_DIR / "voxels.npz", **info)
        print("saved voxels.npz;", len(info["voxel_indices"]), "voxels")
        return
    info = load_info()
    if a.cmd == "pilot":
        its = items(info, a.variants, a.n_seeds, 4, a.tag, opts)
        for it in its:
            Path(it[4]).unlink(missing_ok=True)
        t0 = time.time()
        C.run_pool(task, its, a.workers, "pilot")
        print(f"pilot wall {time.time() - t0:.0f}s for {len(its)} tasks x {a.n_seeds} seeds")
        return
    its = items(info, a.variants, a.n_seeds, a.n_voxels, a.tag, opts)
    todo = [it for it in its if not Path(it[4]).exists()]
    print(f"{len(its)} tasks, {len(todo)} to run, {a.workers} workers", flush=True)
    C.run_pool(task, todo, a.workers, "e0")


if __name__ == "__main__":
    main()
