#!/usr/bin/env python3
"""F1: CSL-relative check after the search (post hoc on the E0 final answers).

For every E0 run's final answer g: all distinct Sigma relatives (Sigma <= 29) are generated, each
gets the coarse quick MC (the settings of the last level: 10 steps x 5 restarts, local cost), the
best few by local cost go through FindOptimal + VarianceMinimizing (refine_from_candidates) and
the lowest final local cost among {g, refined relatives} is returned. The selection variants
(Sigma <= 11 / <= 29, top-k) are derived in summarize_fixes.py from the stored per-relative results.

  uv run python scripts/findoptimal_robustness/f1_run.py run --workers 10
"""

import argparse
import math
import sys
import time
from pathlib import Path
from typing import Any, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import csl  # noqa: E402

TOP_K = 3
SIG_SMALL = 11


def refine_one(
    W: Any, rng: np.random.Generator, vertices: np.ndarray, phase: int, R: np.ndarray
) -> Tuple[np.ndarray, float, int]:
    from icenine.orientation_search import MCOptimizer, SearchCandidate

    lf = W.local_fn
    mc = MCOptimizer(cost_fn=lf, voxel_vertices=vertices, phase_index=phase, rng=rng)
    n0 = lf.eval_count
    with C.quiet():
        res = W.rec.refine_from_candidates(
            [SearchCandidate(orientation=R.astype(np.float32), cost=1.0)],
            vertices,
            phase,
            local_cost_fn=lf,
            mc_optimizer=mc,
        )
    return np.asarray(res.orientation, dtype=np.float64), float(res.cost), lf.eval_count - n0


def quick_mc(
    W: Any,
    rng: np.random.Generator,
    vertices: np.ndarray,
    phase: int,
    R: np.ndarray,
    diameter: float,
) -> Tuple[np.ndarray, float]:
    from icenine.orientation_search import MCOptimizer

    p = W.rec.params
    lf = W.local_fn
    mc = MCOptimizer(cost_fn=lf, voxel_vertices=vertices, phase_index=phase, rng=rng)
    radius = max(diameter / 3.0, math.radians(0.2))
    box = radius / (2**p.min_local_resolution)
    res = mc.optimize(
        initial_orientation=R.astype(np.float32),
        angular_box_side=box,
        angular_step=box * p.mc_radius_scale_factor,
        max_mc_steps=10,
        max_restarts=5,
        max_convergence_cost=p.max_convergence_cost,
    )
    return np.asarray(res.orientation, dtype=np.float64), float(res.cost)


def task(item: Tuple[Any, ...]) -> str:
    vidx, vpos, variant, seeds, e0_path, path, src_path = item
    W = C.get_worker()
    t_start = time.time()
    d = np.load(e0_path)
    if "unbuildable" in d.files:
        return f"voxel {vidx} {variant}: skipped"
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    case_ctx = C.get_worker().ctx
    from optimizer_sweep import voxel_context

    vctx = voxel_context(case_ctx, vidx)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    lf = W.local_fn
    out = {}
    for s in seeds:
        rng = np.random.default_rng([20_000 + vpos, s])
        if src_path:  # seed 0: v{v}_{var}.npz, seeds > 0: v{v}_{var}_s{seed}.npz
            sp = src_path if s == 0 else src_path.replace(".npz", f"_s{s}.npz")
            g = np.load(sp)["R_final"]
        else:
            g = d[f"s{s}_R_final"]
        diameter = float(d[f"s{s}_L3_disc_diameter"])
        rel, labels = csl.csl_relatives(g, max_sigma=29)
        n0 = lf.eval_count
        post, cost, rel_ev = [], [], []
        for Rk in rel:
            e_before = lf.eval_count
            Rp, cp = quick_mc(W, rng, vertices, phase, Rk, diameter)
            post.append(Rp)
            cost.append(cp)
            rel_ev.append(lf.eval_count - e_before)
        post, cost = np.array(post), np.array(cost)
        n_quick = lf.eval_count - n0
        sig = np.array([int("".join(ch for ch in lab if ch.isdigit())) for lab in labels])
        pick = set()
        for sel in (sig <= SIG_SMALL, sig <= 29):
            idx = np.nonzero(sel)[0]
            pick.update(idx[np.argsort(cost[idx])[:TOP_K]].tolist())
        ref_idx = sorted(pick)
        ref_R, ref_cost, ref_ev = [], [], []
        for i in ref_idx:
            R_, c_, e_ = refine_one(W, rng, vertices, phase, post[i])
            ref_R.append(R_)
            ref_cost.append(c_)
            ref_ev.append(e_)
        out[f"s{s}_rel_label"] = np.array(labels)
        out[f"s{s}_rel_sigma"] = sig
        out[f"s{s}_rel_R0"] = rel
        out[f"s{s}_rel_R_post"] = post
        out[f"s{s}_rel_cost_post"] = cost
        out[f"s{s}_rel_evals"] = np.array(rel_ev)
        out[f"s{s}_ref_idx"] = np.array(ref_idx)
        out[f"s{s}_ref_R"] = np.array(ref_R)
        out[f"s{s}_ref_cost"] = np.array(ref_cost)
        out[f"s{s}_ref_evals"] = np.array(ref_ev)
        out[f"s{s}_quick_evals"] = n_quick
    np.savez_compressed(path, **out)
    return f"voxel {vidx} {variant} F1 done in {time.time() - t_start:.0f}s"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--n-seeds", type=int, default=C.N_SEEDS)
    ap.add_argument(
        "--source",
        default="",
        help="cache dir of runs with R_final (seed 0 only), e.g. e2_rerank_GBT; default E0 answers",
    )
    ap.add_argument(
        "--source-seeds",
        type=int,
        nargs="+",
        default=[0],
        help="seeds of the --source runs (file v*_s<seed>.npz for seed > 0); output of a seed set "
        "other than [0] is v*_s<seeds>.npz with s<seed>_ keys (default [0]: unchanged)",
    )
    ap.add_argument(
        "--cache-root",
        default="",
        help="directory holding the --source dir and receiving f1_<source> "
        "(default: this study's cache, so existing behaviour is unchanged)",
    )
    ap.add_argument("--limit", type=int, default=0, help="only the first N tasks (testing)")
    a = ap.parse_args()
    root = Path(a.cache_root).resolve() if a.cache_root else C.CACHE_DIR
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    cache = root / ("f1" if not a.source else f"f1_{a.source}")
    cache.mkdir(parents=True, exist_ok=True)
    vox = [int(v) for v in info["voxel_indices"]]
    its = []
    for vpos, v in enumerate(vox):
        for var in C.VARIANTS:
            e0 = C.CACHE_DIR / "e0" / f"v{v}_{var}.npz"
            tag = (
                ""
                if (not a.source or a.source_seeds == [0])
                else "_s" + "".join(str(x) for x in a.source_seeds)
            )
            out = cache / f"v{v}_{var}{tag}.npz"
            src = root / a.source / f"v{v}_{var}.npz" if a.source else None
            src_ok = src is None or all(
                Path(str(src).replace(".npz", "" if x == 0 else f"_s{x}.npz")).exists()
                for x in a.source_seeds
            )
            if e0.exists() and not out.exists() and src_ok:
                its.append(
                    (
                        v,
                        vpos,
                        var,
                        tuple(a.source_seeds) if src else tuple(range(a.n_seeds)),
                        str(e0),
                        str(out),
                        str(src) if src else "",
                    )
                )
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    C.run_pool(task, its, a.workers, "f1")


if __name__ == "__main__":
    main()
