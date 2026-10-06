#!/usr/bin/env python3
"""E1 fixes that change the search itself (F1b, F2, F3): full reconstruct_voxel re-runs with a
knob changed, on the same images and with the same rng seeds as E0 (so a run differs from its E0
twin only through the knob).

  F2a    keep the best 1/2 (not 1/4) per level, FindOptimal cap 60 candidates
  F2b    keep the union of the top 1/4 by post-MC local cost and the top 1/4 by discrete score
  F3a    coarse levels use Q_max 8 from level 0 (n_q_max starts at 8, not 5)
  F3b    coarse cost pixel_radius 1 (not 3)
  F3c    both F3a and F3b
  F1b    before each level's quick MC add the CSL relatives (Sigma <= 11) of the 3 best candidates
         by discrete score

  uv run python scripts/findoptimal_robustness/fixes_run.py run --fix F2a --seeds 0 --workers 10
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import csl  # noqa: E402

FIXES = {
    "F2a": dict(keep_fraction=0.5, max_discrete_candidates=60),
    "F2b": dict(keep_union_discrete=True),
    "F3a": dict(n_q_start_offset=3.0),
    "F3b": dict(global_pixel_radius=1),
    "F3c": dict(n_q_start_offset=3.0, global_pixel_radius=1),
    "F1b": dict(csl_expand=3),
}


def make_expander(top_k: int, max_sigma: int = 11):
    from icenine.orientation_search import SearchCandidate

    def expand(level, candidates):
        order = sorted(range(len(candidates)), key=lambda i: candidates[i].cost)[:top_k]
        extra = []
        for i in order:
            rel, _ = csl.csl_relatives(
                np.asarray(candidates[i].orientation, dtype=np.float64), max_sigma=max_sigma
            )
            extra += [SearchCandidate(orientation=r.astype(np.float32), cost=1.0) for r in rel]
        return extra

    return expand


def apply_knobs(rec, knobs):
    rec.keep_fraction, rec.keep_union_discrete = 0.25, False
    rec.n_q_start_offset, rec.global_pixel_radius = 0.0, 3
    rec.extra_candidates = None
    rec.params.max_discrete_candidates = 30
    for k, v in knobs.items():
        if k == "csl_expand":
            rec.extra_candidates = make_expander(v)
        elif k == "max_discrete_candidates":
            rec.params.max_discrete_candidates = v
        else:
            setattr(rec, k, v)


def task(item):
    fix, vidx, vpos, variant, seeds, path = item
    W = C.get_worker()
    t_start = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    from optimizer_sweep import voxel_context

    vctx = voxel_context(W.ctx, vidx)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    apply_knobs(W.rec, FIXES[fix])
    out = {}
    for s in seeds:
        rec = C.Recorder()
        W.rec.recorder = rec
        t0 = time.perf_counter()
        with C.quiet():
            res = W.rec.reconstruct_voxel(vertices, phase, rng=C.run_seed(vpos, s))
        dt = time.perf_counter() - t0
        W.rec.recorder = None
        g, loc, _ = W.rec.last_eval_counts
        out[f"s{s}_R_final"] = np.asarray(res.orientation, dtype=np.float64)
        out[f"s{s}_cost_final"] = float(res.cost)
        out[f"s{s}_runtime"] = dt
        out[f"s{s}_evals_global"] = g
        out[f"s{s}_evals_local"] = loc
        arr = rec.to_arrays()
        for k, v in arr.items():
            if k.startswith("L") and k.endswith("_qmc_R") or k.endswith("_qmc_cost"):
                out[f"s{s}_{k}"] = v
            if k.startswith("L") and k.endswith("_qmc_n_keep"):
                out[f"s{s}_{k}"] = v
    apply_knobs(W.rec, {})
    np.savez_compressed(path, **out)
    return f"{fix} voxel {vidx} {variant} done in {time.time() - t_start:.0f}s"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--fix", required=True, choices=sorted(FIXES))
    ap.add_argument("--seeds", type=int, nargs="+", default=[0])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--n-voxels", type=int, default=C.N_VOXELS)
    a = ap.parse_args()
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    cache = C.CACHE_DIR / "fix" / a.fix
    cache.mkdir(parents=True, exist_ok=True)
    vox = [int(v) for v in info["voxel_indices"]]
    its = []
    for vpos, v in enumerate(vox[: a.n_voxels + 8]):
        for var in C.VARIANTS:
            e0 = C.CACHE_DIR / "e0" / f"v{v}_{var}.npz"
            if not e0.exists() or "unbuildable" in np.load(e0).files:
                continue
            tag = "".join(str(s) for s in a.seeds)
            out = cache / f"v{v}_{var}_s{tag}.npz"
            if not out.exists():
                its.append((a.fix, v, vpos, var, a.seeds, str(out)))
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    C.run_pool(task, its, a.workers, a.fix)


if __name__ == "__main__":
    main()
