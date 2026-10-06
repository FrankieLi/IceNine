#!/usr/bin/env python3
"""E2: can a one-pass classifier tell the truth basin from traps better than the local cost?

Dataset: cache/e2/*.npz (e2_dataset.py). Models: logistic regression and histogram gradient-boosted
trees, with and without the two cost features, trained to predict "within 3 deg of the truth".
Everything is split by GRAIN: 4 folds of about 50 voxels (both variants of a voxel and all voxels of
a grain together; the 200 voxels come from 158 grains); a model is evaluated only on grains it
never saw. Predictions of the 4 held-out folds are pooled.

Evaluation sets (held-out voxels):
  A contested   harvested candidates that matter for the search: level-0..2 candidates with a
                post-MC cost rank below 2 x n_keep, all level-3 candidates and FindOptimal's
                results; basin (< 1 deg) against trap (> 3 deg); ROC-AUC
  B synthetic   synthetic perturbed truths against perturbed Sigma<=29 relatives of the truth
  C csl-near    basin against harvested candidates within 3 deg of an exact Sigma<=29 relative
  D pruning     per (voxel, variant, seed, level 0-2) with a basin(3) candidate present: is the
                best such candidate among the n_keep best by the score? (the search's actual
                pruning decision; the baseline is the local cost)
  E final       per run, among FindOptimal's returned candidates (>= 1 basin and >= 1 trap): is the
                top-scored one a basin candidate?

  uv run python scripts/findoptimal_robustness/e2_models.py run [--save-models]
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

FOLDS = 4
SIZES = (25, 50, 100, 150)
COST_COLS = 2  # the last two features are cost_local, cost_global3


def load_dataset() -> Dict[str, np.ndarray]:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    parts = []
    for v in vox:
        for var in C.VARIANTS:
            f = C.CACHE_DIR / "e2" / f"v{v}_{var}.npz"
            if f.exists():
                d = np.load(f)
                parts.append({k: d[k] for k in d.files if k not in ("R", "truth_features")})
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def fold_of(vpos: np.ndarray, seed: int = 0) -> np.ndarray:
    """Fold of each voxel position, GRAIN-disjoint: all voxels of one grain (same orientation, hence
    the same CSL relatives) are in one fold; grains are assigned greedily, largest first, to the
    fold with the fewest voxels (ties broken by a seeded shuffle)."""
    grain = np.load(C.OUT_DIR / "voxels.npz")["voxel_grain_id"]
    n = len(grain)
    ids, counts = np.unique(grain, return_counts=True)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(ids))
    order = order[np.argsort(-counts[order], kind="stable")]
    load = np.zeros(FOLDS, dtype=int)
    fold_of_grain = {}
    for k in order:
        f = int(np.argmin(load))
        fold_of_grain[int(ids[k])] = f
        load[f] += counts[k]
    per_voxel = np.array([fold_of_grain[int(g)] for g in grain])
    return per_voxel[np.asarray(vpos)]


def make_model(kind: str):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    if kind == "LR":
        return make_pipeline(
            StandardScaler(), LogisticRegression(C=1.0, max_iter=3000, class_weight="balanced")
        )
    return HistGradientBoostingClassifier(
        max_depth=4, learning_rate=0.1, max_iter=200, l2_regularization=1.0, random_state=0,
        class_weight="balanced",
    )  # fmt: skip


def cols(with_cost: bool, n_feat: int) -> slice:
    return slice(0, n_feat) if with_cost else slice(0, n_feat - COST_COLS)


def fit(kind: str, X: np.ndarray, err: np.ndarray, with_cost: bool):
    m = make_model(kind)
    m.fit(X[:, cols(with_cost, X.shape[1])], (err < 3.0).astype(int))
    return m


def score(m, X: np.ndarray, with_cost: bool) -> np.ndarray:
    return m.predict_proba(X[:, cols(with_cost, X.shape[1])])[:, 1]


def eval_sets(D: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Boolean masks (n,) of the evaluation populations."""
    src, lvl, rank = D["source"], D["level"], D["rank"]
    err = D["err"]
    n_keep = {0: 32, 1: 8, 2: 2}  # typical (E0 median) keep counts; real ones differ per run
    contested = (
        (src == 1)
        | ((src == 0) & (lvl == 3))
        | ((src == 0) & (lvl <= 2) & (rank < 2 * np.array([n_keep.get(int(l), 2) for l in lvl])))
    )
    harvest = (src == 0) | (src == 1)
    basin = err < 1.0
    trap = err > 3.0
    return dict(
        A=(contested & (basin | trap)),
        B=((src == 2) | (src == 3)),
        C=(harvest & (basin | (trap & (D["csl_dev"] <= 3.0) & (D["csl_label"] != "none")))),
    )


def groups(D: Dict[str, np.ndarray], kind: str) -> List[np.ndarray]:
    """Index arrays of the groups of evaluation D (pruning levels 0-2) or E (FindOptimal results)."""
    idx = np.arange(len(D["err"]))
    out = []
    key_cols = ("vpos", "variant", "seed", "level")
    if kind == "D":
        sel = (D["source"] == 0) & (D["level"] <= 2)
    else:
        sel = D["source"] == 1
    keys = np.stack([D[c][sel] for c in key_cols], axis=1)
    ii = idx[sel]
    _, inv = np.unique(keys, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    order = np.argsort(inv, kind="stable")
    bounds = np.flatnonzero(np.diff(inv[order])) + 1
    return [ii[g] for g in np.split(order, bounds)]


def pruning_recall(D, scores, grp, n_keep_fn) -> Tuple[float, int]:
    """Fraction of level groups with a basin(3) candidate whose best-ranked basin candidate (by
    `scores`, higher better) is within the n_keep best of the group."""
    ok, n = 0, 0
    for g in grp:
        e = D["err"][g]
        if not (e < 3.0).any():
            continue
        nk = max(1, len(g) // 4)
        order = np.argsort(-scores[g], kind="stable")
        ranks = np.empty(len(g), dtype=int)
        ranks[order] = np.arange(len(g))
        ok += int(ranks[e < 3.0].min() < nk)
        n += 1
    return ok / max(n, 1), n


def final_precision(D, scores, grp) -> Tuple[float, int]:
    ok, n = 0, 0
    for g in grp:
        e = D["err"][g]
        if (e < 1.0).any() and (e > 3.0).any():
            ok += int(e[np.argmax(scores[g])] < 1.0)
            n += 1
    return ok / max(n, 1), n


def auc(y: np.ndarray, s: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    if y.min() == y.max():
        return float("nan")
    return float(roc_auc_score(y, s))


def metrics(
    D: Dict[str, np.ndarray], scores: Dict[str, np.ndarray], mask_override=None
) -> Dict[str, Any]:
    sets = eval_sets(D)
    out: Dict[str, Any] = {}
    gD, gE = groups(D, "D"), groups(D, "E")
    for name, s in scores.items():
        r: Dict[str, Any] = {}
        for k in ("A", "B", "C"):
            m = sets[k]
            if k == "B":
                y = (D["source"][m] == 2).astype(int)
            else:
                y = (D["err"][m] < 1.0).astype(int)
            r[f"auc_{k}"] = auc(y, s[m])
            r[f"n_{k}"] = int(m.sum())
        r["prune_recall"], r["n_prune_groups"] = pruning_recall(D, s, gD, None)
        r["final_precision"], r["n_final_groups"] = final_precision(D, s, gE)
        out[name] = r
    return out


def cv_scores(
    D, with_variants_train=("clean", "all"), test_variant=None, kinds=("LR", "GBT"), save=False
):
    """Pooled held-out scores of all models over the voxel folds. Returns dict name -> scores."""
    fold = fold_of(D["vpos"])
    X = D["X"]
    var = D["variant"]
    train_var = np.isin(var, [C.VARIANTS.index(v) for v in with_variants_train])
    S = {"cost": -X[:, -COST_COLS]}
    for kind in kinds:
        for wc in (True, False):
            S[f"{kind}{'' if wc else '-nocost'}"] = np.full(len(X), np.nan)
    for f in range(FOLDS):
        tr = (fold != f) & train_var
        te = fold == f
        for kind in kinds:
            for wc in (True, False):
                m = fit(kind, X[tr], D["err"][tr], wc)
                S[f"{kind}{'' if wc else '-nocost'}"][te] = score(m, X[te], wc)
                if save:
                    import joblib

                    (C.CACHE_DIR / "models").mkdir(parents=True, exist_ok=True)
                    joblib.dump(
                        m,
                        C.CACHE_DIR / "models" / f"fold{f}_{kind}{'' if wc else '-nocost'}.joblib",
                    )
    return S


def subset(D, mask):
    return {k: v[mask] for k, v in D.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--save-models", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    D = load_dataset()
    n = len(D["err"])
    print(
        f"dataset: {n} candidates, {D['X'].shape[1]} features, {len(np.unique(D['vpos']))} voxels"
    )
    results: Dict[str, Any] = {"n_candidates": n}
    lines: List[str] = []

    def fmt(res):
        o = []
        for name, r in res.items():
            o.append(
                f"  {name:12s} AUC A {r['auc_A']:.3f} (n={r['n_A']})  B {r['auc_B']:.3f}  C {r['auc_C']:.3f} (n={r['n_C']})   "
                f"prune recall {r['prune_recall']:.3f} (n={r['n_prune_groups']})  final precision {r['final_precision']:.3f} (n={r['n_final_groups']})"
            )
        return "\n".join(o)

    # (1) trained on both variants, tested on held-out voxels, by variant of the test data
    S = cv_scores(D, save=a.save_models)
    for tv in ("both", "clean", "all"):
        m = np.ones(n, bool) if tv == "both" else D["variant"] == C.VARIANTS.index(tv)
        res = metrics(subset(D, m), {k: v[m] for k, v in S.items()})
        results[f"cv_train_both_test_{tv}"] = res
        lines.append(f"[train both variants, held-out voxels, test on {tv}]\n" + fmt(res))
    np.savez_compressed(C.CACHE_DIR / "e2_cv_scores.npz", **{f"score_{k}": v for k, v in S.items()})
    # ROC curves on set A (pooled held-out)
    from sklearn.metrics import roc_curve

    A = eval_sets(D)["A"]
    roc = {}
    for k, v in S.items():
        fpr, tpr, _ = roc_curve((D["err"][A] < 1.0).astype(int), v[A])
        roc[f"{k}_fpr"], roc[f"{k}_tpr"] = fpr, tpr
    np.savez_compressed(C.OUT_DIR / "e2_roc.npz", **roc)
    # (2) cross-variant: train on one variant's voxels, test on the other variant (held-out voxels)
    for trv, tev in (("clean", "all"), ("all", "clean")):
        S2 = cv_scores(D, with_variants_train=(trv,))
        m = D["variant"] == C.VARIANTS.index(tev)
        res = metrics(subset(D, m), {k: v[m] for k, v in S2.items()})
        results[f"train_{trv}_test_{tev}"] = res
        lines.append(f"[train {trv} voxels only, test {tev} on held-out voxels]\n" + fmt(res))
    # (3) AUC by negative class (basin vs each CSL class), cost vs GBT, set A+harvest
    harvest = D["source"] <= 1
    lab = D["csl_label"]
    by_class = {}
    for cl in sorted(set(lab[harvest & (D["err"] > 3.0)].tolist())):
        neg = harvest & (D["err"] > 3.0) & (lab == cl) & (D["csl_dev"] <= 3.0)
        pos = harvest & (D["err"] < 1.0)
        if neg.sum() < 20:
            continue
        m = neg | pos
        y = (D["err"][m] < 1.0).astype(int)
        by_class[cl] = dict(
            n_neg=int(neg.sum()),
            auc_cost=auc(y, S["cost"][m]),
            auc_gbt=auc(y, S["GBT"][m]),
            auc_lr=auc(y, S["LR"][m]),
        )
    results["auc_by_csl_class"] = by_class
    lines.append(
        "[AUC basin vs candidates within 3 deg of an exact CSL relative, by class]\n"
        + "\n".join(
            f"  Sigma{k:5s} n_neg {v['n_neg']:6d}  cost {v['auc_cost']:.3f}  LR {v['auc_lr']:.3f}  GBT {v['auc_gbt']:.3f}"
            for k, v in by_class.items()
        )
    )
    # (4) learning curve
    fold = fold_of(D["vpos"])
    rng = np.random.default_rng(1)
    lc: Dict[str, Any] = {}
    for size in SIZES:
        acc = {
            k: [] for k in ("GBT_auc", "GBT_prune", "LR_auc", "LR_prune", "GBT_final", "LR_final")
        }
        for f in range(FOLDS):
            te = fold == f
            tr_vox = np.unique(D["vpos"][fold != f])
            for rep in range(3 if size < 150 else 1):
                pick = rng.choice(tr_vox, size=min(size, len(tr_vox)), replace=False)
                tr = np.isin(D["vpos"], pick)
                for kind in ("GBT", "LR"):
                    m = fit(kind, D["X"][tr], D["err"][tr], True)
                    sc = score(m, D["X"][te], True)
                    r = metrics(subset(D, te), {kind: sc})[kind]
                    acc[f"{kind}_auc"].append(r["auc_A"])
                    acc[f"{kind}_prune"].append(r["prune_recall"])
                    acc[f"{kind}_final"].append(r["final_precision"])
        lc[size] = {k: [float(np.mean(v)), float(np.std(v)), len(v)] for k, v in acc.items()}
        lines.append(
            f"[learning curve, {size} train voxels] "
            + "  ".join(f"{k} {v[0]:.3f}+-{v[1]:.3f}" for k, v in lc[size].items())
        )
    # baseline (cost) on the same test folds for reference
    base = metrics(D, {"cost": S["cost"]})["cost"]
    lc["cost_baseline"] = base
    results["learning_curve"] = lc
    lines.append(
        f"[cost baseline, all held-out] AUC A {base['auc_A']:.3f}  prune recall {base['prune_recall']:.3f}  final precision {base['final_precision']:.3f}"
    )
    # feature importance (permutation-free): LR coefficients of the last fold are unstable; report GBT impurity-free proxy via drop in AUC? skipped
    C.save_json(C.OUT_DIR / "e2_results.json", results)
    (C.OUT_DIR / "e2_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
