#!/usr/bin/env python3
"""Labels of the Q_max-8 cost proxy and their validation.

REGRESSION TARGET  y_bcost ("basin cost"): the Q_max-8 local cost a candidate would reach if the
search refined it. For every (voxel, variant) case, with its detector images attached:
  * candidate within 3 deg of the truth (reduced misorientation): the cost at the truth after
    `refine_from_candidates([truth])` (FindOptimal + VarianceMinimizing), one number per case;
  * else within 3 deg of an exact Sigma <= 29 relative of the truth: the cost after a SHORT
    refinement of that relative (quick MC of the last level, then one MC of up to 200 steps with the
    final box, 0.329 deg), one number per relative (about 270 per case);
  * else (censored): the candidate's own cost (harvested candidates: the post-quick-MC cost;
    FindOptimal results: their refined cost; synthetic candidates: no label, NaN).
  The candidate's OWN cost is never the target of a candidate in the truth's basin or at a
  relative: that is the ranking that fails (pruning recall 0.900). Rationale: in all 347 wrong runs
  of the FindOptimal study the cost at the truth was lower than the trap's, so a perfect y_bcost
  ranks the truth's basin first.
  Validation of the censored branch: 5 randomly drawn censored harvested candidates per case
  (2,000 in total) are refined with the same short procedure; `validate` reports the Spearman
  correlation of their y_bcost (own post-quick-MC cost) with the refined cost.
CLASSIFICATION LABELS  y3 = error < 3 deg, y1 = error < 1 deg, and the level-matched label
  y_lm = y3 for candidates of levels 0-2 (and synthetic ones), y1 for level 3 and FindOptimal
  results.

  uv run python scripts/coarse_proxy/labels.py refine --workers 10
  uv run python scripts/coarse_proxy/labels.py build
  uv run python scripts/coarse_proxy/labels.py validate
"""

import argparse
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

C = B.C
import csl  # noqa: E402
import f1_run  # noqa: E402

NEAR_DEG = 3.0
N_VAL_PER_TASK = 5


def reduced_angle_matrix(Ra: np.ndarray, Rb: np.ndarray) -> np.ndarray:
    """(n, m) cubic-symmetry-reduced angles (deg) between orientations Ra (n,3,3) and Rb (m,3,3);
    equals csl.reduced_misorientation_deg(Ra[i], Rb[j]) (via max over S of tr(Ra S Rb^T))."""
    n, m = len(Ra), len(Rb)
    P = np.einsum("nik,skl->nsil", Ra, csl.OPS).reshape(n * 24, 9)
    T = (P @ Rb.reshape(m, 9).T).reshape(n, 24, m).max(axis=1)
    return np.degrees(np.arccos(np.clip((T - 1.0) / 2.0, -1.0, 1.0)))


def short_refine(
    W: Any, rng: np.random.Generator, vertices: Any, phase: int, R: np.ndarray
) -> Tuple[np.ndarray, float, int]:
    """Quick MC of the last level (10 steps x 5 restarts, diameter of level 3), then one MC of up
    to max_mc_steps with the final box. Returns (R, cost, local evaluations)."""
    from icenine.orientation_search import MCOptimizer

    p = W.rec.params
    lf = W.local_fn
    n0 = lf.eval_count
    d3 = p.local_grid_radius / 1.5**p.max_local_resolution  # diameter at level 3
    Rq, _ = f1_run.quick_mc(W, rng, vertices, phase, R, d3)
    d_end = p.local_grid_radius / 1.5 ** (p.max_local_resolution + 1)
    box = max(d_end / 3.0, math.radians(0.2)) / (2**p.min_local_resolution)
    mc = MCOptimizer(cost_fn=lf, voxel_vertices=vertices, phase_index=phase, rng=rng)
    res = mc.optimize(
        initial_orientation=Rq.astype(np.float32),
        angular_box_side=box,
        angular_step=box * p.mc_radius_scale_factor,
        max_mc_steps=p.max_mc_steps,
        max_restarts=p.successive_restarts,
        max_convergence_cost=p.max_convergence_cost,
    )
    return np.asarray(res.orientation, dtype=np.float64), float(res.cost), lf.eval_count - n0


def task(item: Tuple[Any, ...]) -> str:
    vidx, vpos, variant, path = item
    from optimizer_sweep import voxel_context

    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    e2 = np.load(C.CACHE_DIR / "e2" / f"v{vidx}_{variant}.npz")
    R_true = np.load(C.CACHE_DIR / "e0" / f"v{vidx}_{variant}.npz")["R_true"]
    rng = np.random.default_rng([40_000 + vpos, C.VARIANTS.index(variant)])
    Rt_full, c_full, ev_full = f1_run.refine_one(W, rng, vertices, phase, R_true)
    Rt_short, c_short, ev_short = short_refine(W, rng, vertices, phase, R_true)
    rel, labels = csl.csl_relatives(R_true, max_sigma=29)
    post, cost, evs = [], [], []
    for Rk in rel:
        Rp, cp, ev = short_refine(W, rng, vertices, phase, Rk)
        post.append(Rp)
        cost.append(cp)
        evs.append(ev)
    # validation draw: censored harvested candidates (source 0, > 3 deg from truth and relatives)
    R, err, src = e2["R"], e2["err"], e2["source"]
    near_rel = reduced_angle_matrix(R, rel).min(axis=1)
    censored = np.nonzero((src == 0) & (err >= NEAR_DEG) & (near_rel >= NEAR_DEG))[0]
    pick = np.sort(rng.choice(censored, size=min(N_VAL_PER_TASK, len(censored)), replace=False))
    val_cost, val_ev = [], []
    for i in pick:
        _, cv, ev = short_refine(W, rng, vertices, phase, R[i])
        val_cost.append(cv)
        val_ev.append(ev)
    sig = np.array([int("".join(ch for ch in lab if ch.isdigit())) for lab in labels])
    np.savez_compressed(
        path,
        truth_full_R=Rt_full, truth_full_cost=c_full, truth_full_evals=ev_full,
        truth_short_cost=c_short, truth_short_evals=ev_short,
        rel_R0=rel, rel_sigma=sig, rel_post_R=np.array(post), rel_post_cost=np.array(cost),
        rel_evals=np.array(evs), val_idx=pick, val_cost=np.array(val_cost),
        val_evals=np.array(val_ev), seconds=time.time() - t0,
    )  # fmt: skip
    return f"voxel {vidx} {variant}: {len(rel)} relatives, {time.time() - t0:.0f}s"


def build_labels(
    e2: Dict[str, np.ndarray], lab: Dict[str, Any]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(y_bcost, category, distance to nearest relative) for one case. category: 0 = within 3 deg
    of the truth, 1 = of a relative, 2 = censored (own cost), 3 = no label (synthetic, far)."""
    err, src, own = e2["err"], e2["source"], e2["cost"]
    d_rel = reduced_angle_matrix(e2["R"], lab["rel_R0"])
    k = d_rel.argmin(axis=1)
    near = d_rel[np.arange(len(k)), k]
    y = np.full(len(err), np.nan)
    cat = np.full(len(err), 3, dtype=np.int8)
    harvest = src <= 1
    y[harvest] = own[harvest]
    cat[harvest] = 2
    is_rel = near < NEAR_DEG
    y[is_rel] = lab["rel_post_cost"][k[is_rel]]
    cat[is_rel] = 1
    is_truth = err < NEAR_DEG
    y[is_truth] = float(lab["truth_full_cost"])
    cat[is_truth] = 0
    return y, cat, near


def cmd_build() -> None:
    ys, cats, nears = [], [], []
    for v, _, var in B.tasks():
        e2 = dict(np.load(C.CACHE_DIR / "e2" / f"v{v}_{var}.npz"))
        lab = dict(np.load(B.CACHE / "labels" / f"v{v}_{var}.npz"))
        y, cat, near = build_labels(e2, lab)
        ys.append(y)
        cats.append(cat)
        nears.append(near)
    D = B.M.load_dataset()
    y, cat = np.concatenate(ys), np.concatenate(cats)
    err, lvl = D["err"], D["level"]
    y3, y1 = (err < 3.0).astype(np.int8), (err < 1.0).astype(np.int8)
    y_lm = np.where((lvl == 3) | (lvl == 4), y1, y3).astype(np.int8)
    np.savez_compressed(
        B.CACHE / "labels" / "y.npz", y_bcost=y, category=cat, near_rel=np.concatenate(nears),
        y3=y3, y1=y1, y_lm=y_lm,
    )  # fmt: skip
    names = ["truth basin", "relative basin", "censored (own cost)", "unlabelled"]
    for c in range(4):
        print(f"category {c} {names[c]:22s}: {int((cat == c).sum())}")


def cmd_validate() -> None:
    from scipy.stats import spearmanr

    own, ref, lvl, per_case = [], [], [], []
    D = B.M.load_dataset()
    off = 0
    for v, _, var in B.tasks():
        lab = np.load(B.CACHE / "labels" / f"v{v}_{var}.npz")
        n = len(np.load(C.CACHE_DIR / "e2" / f"v{v}_{var}.npz")["err"])
        idx = lab["val_idx"]
        c_own = D["cost"][off + idx]
        own.append(c_own)
        ref.append(lab["val_cost"])
        lvl.append(D["level"][off + idx])
        if len(idx) >= 3:
            per_case.append(spearmanr(c_own, lab["val_cost"])[0])
        off += n
    own, ref, lvl = np.concatenate(own), np.concatenate(ref), np.concatenate(lvl)
    rho = float(spearmanr(own, ref)[0])
    res = dict(
        n=int(len(own)), spearman_pooled=rho,
        spearman_within_case_mean=float(np.nanmean(per_case)),
        median_cost_own=float(np.median(own)), median_cost_refined=float(np.median(ref)),
        median_drop=float(np.median(own - ref)), frac_refined_lower=float(np.mean(ref < own)),
        by_level={
            int(lv): dict(
                n=int((lvl == lv).sum()),
                spearman=float(spearmanr(own[lvl == lv], ref[lvl == lv])[0]),
            )
            for lv in np.unique(lvl)
            if (lvl == lv).sum() > 30
        },
        verdict="regression label kept" if rho >= 0.8 else "FALL BACK to the y3/y1 classifiers",
    )  # fmt: skip
    B.OUT.mkdir(parents=True, exist_ok=True)
    C.save_json(B.OUT / "label_validation.json", res)
    print(res)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["refine", "build", "validate"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    cache = B.CACHE / "labels"
    cache.mkdir(parents=True, exist_ok=True)
    if a.cmd == "build":
        return cmd_build()
    if a.cmd == "validate":
        return cmd_validate()
    its = [
        (v, vpos, var, str(cache / f"v{v}_{var}.npz"))
        for v, vpos, var in B.tasks()
        if not (cache / f"v{v}_{var}.npz").exists()
    ]
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    C.run_pool(task, its, a.workers, "labels")


if __name__ == "__main__":
    main()
