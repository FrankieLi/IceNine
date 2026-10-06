#!/usr/bin/env python3
"""Offline models of the Q_max-8 cost proxy (grain-disjoint 4-fold evaluation).

Feature sets (X columns of the E2 table, low-Q tables of lowq_dataset.py):
  cost8      the local Q_max-8 cost of the candidate itself (the current ranking; not trained)
  cost_q4/5  the raw local cost at max_q = 4 / 5 (not trained)
  hand       the untrained hand-made key: mean of overall hit0, hit_any, hit1, hit3 (E2 table)
  cache4/5   F-cache: log_n_pairs and the families with |q| <= Q of the cached full E2 table
             (no new pass; a diagnostic, since it needs the full pass to exist)
  lowq4/5    F-lowQ: FeatureExtractor(q_max=Q): a pass over |q| <= Q reflections only + the two cost
             columns at max_q = Q
  lowq4+c8 / lowq5+c8   F-lowQ plus the free post-quick-MC Q_max-8 cost (the deployable variant)
  e2full     the 62 E2 columns (reference; one full pass + two Q_max-8 costs)
Targets: clf3 (GBT classifier, error < 3 deg), clf_lm (level-matched: < 3 deg at levels 0-2 and for
synthetic candidates, < 1 deg at level 3 and for FindOptimal results), reg (GBT regressor of the
basin cost y_bcost, labels.py). Trained on both variants mixed; evaluated on grains never seen.

Metrics (the E2 sets; harvested vs synthetic reported separately because the synthetic candidates
are built from the truth):
  AUC  A contested harvested candidates, B synthetic perturbed truths vs perturbed relatives,
       C basin vs harvested candidates near a CSL relative;
  D    pruning recall per level at keep 1/4 and 1/8 (is the best basin candidate - error < 3 deg,
       < 1 deg at level 3 - among the n_keep = max(1, int(n * frac)) best of its level group?);
  E    final precision among FindOptimal's returned candidates; Spearman of the score with -y_bcost.

  uv run python scripts/coarse_proxy/models.py run
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))  # before common sets 1

import argparse  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Dict, List, Optional, Tuple  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

C, M = B.C, B.M
import features as F  # noqa: E402

FRACS = (0.25, 0.125)
N_FEAT_FULL = 62
# (kind, Q, add Q_max-8 cost)
SETS: Dict[str, Tuple[str, Optional[int], bool]] = {
    "cache4": ("cache", 4, False),
    "cache5": ("cache", 5, False),
    "lowq4": ("lowq", 4, False),
    "lowq5": ("lowq", 5, False),
    "lowq4+c8": ("lowq", 4, True),
    "lowq5+c8": ("lowq", 5, True),
    "e2full": ("full", None, False),
}
TARGETS = ("clf3", "clf_lm", "reg")
UNTRAINED = ("cost8", "cost_q4", "cost_q5", "hand")
N_FAM = {4: 2, 5: 3}


def load_all() -> Dict[str, np.ndarray]:
    D = M.load_dataset()
    low = B.load_per_task(B.CACHE / "lowq", ["X4", "X5"])
    D["X4"], D["X5"] = low["X4"], low["X5"]
    y = np.load(B.CACHE / "labels" / "y.npz")
    for k in y.files:
        D[k] = y[k]
    assert len(D["X4"]) == len(D["err"]) == len(D["y_bcost"])
    return D


def feature_matrix(name: str, D: Dict[str, np.ndarray]) -> np.ndarray:
    kind, q, c8 = SETS[name]
    X = D["X"]
    if kind == "full":
        return X
    if kind == "cache":
        cols = [0] + list(range(5, 5 + 3 * N_FAM[q]))
        out = X[:, cols]
    else:
        out = D[f"X{q}"]
    return np.hstack([out, X[:, -2:-1]]) if c8 else out


def untrained_score(name: str, D: Dict[str, np.ndarray]) -> np.ndarray:
    """Higher = better."""
    if name == "cost8":
        return -D["X"][:, -2]
    if name == "cost_q4":
        return -D["X4"][:, -2]
    if name == "cost_q5":
        return -D["X5"][:, -2]
    return D["X"][:, 1:5].mean(axis=1)  # hand: hit0, hit_any, hit1, hit3


def make_target_model(target: str) -> Any:
    from sklearn.ensemble import HistGradientBoostingRegressor

    if target == "reg":
        return HistGradientBoostingRegressor(
            max_depth=4, learning_rate=0.1, max_iter=200, l2_regularization=1.0, random_state=0
        )
    return M.make_model("GBT")


def fit_target(target: str, X: np.ndarray, D: Dict[str, np.ndarray], tr: np.ndarray) -> Any:
    if target == "reg":
        ok = tr & np.isfinite(D["y_bcost"])
        return make_target_model(target).fit(X[ok], D["y_bcost"][ok])
    y = D["y3"] if target == "clf3" else D["y_lm"]
    return make_target_model(target).fit(X[tr], y[tr])


def predict_score(target: str, model: Any, X: np.ndarray) -> np.ndarray:
    """Higher = better."""
    if target == "reg":
        return -model.predict(X)
    return model.predict_proba(X)[:, 1]


def model_path(name: str, target: str, fold: int) -> Path:
    return B.CACHE / "models" / f"{name}__{target}__fold{fold}.joblib"


def cv_scores(D: Dict[str, np.ndarray], save: bool = True) -> Dict[str, np.ndarray]:
    import joblib

    fold = M.fold_of(D["vpos"])
    S: Dict[str, np.ndarray] = {n: untrained_score(n, D) for n in UNTRAINED}
    (B.CACHE / "models").mkdir(parents=True, exist_ok=True)
    for name in SETS:
        X = feature_matrix(name, D)
        for target in TARGETS:
            key = f"{name}|{target}"
            S[key] = np.full(len(X), np.nan)
            for f in range(M.FOLDS):
                t0 = time.time()
                tr = fold != f
                m = fit_target(target, X, D, tr)
                S[key][fold == f] = predict_score(target, m, X[fold == f])
                if save:
                    joblib.dump(m, model_path(name, target, f))
                print(f"  {key} fold {f}: {time.time() - t0:.0f}s", flush=True)
    return S


# -- metrics ------------------------------------------------------------------------------------
def level_groups(D: Dict[str, np.ndarray]) -> Tuple[List[np.ndarray], np.ndarray]:
    """Index arrays of the harvested quick-MC groups (voxel, variant, seed, level 0-3) and the
    level of each."""
    sel = (D["source"] == 0) & (D["level"] <= 3)
    keys = np.stack([D[c][sel] for c in ("vpos", "variant", "seed", "level")], axis=1)
    ii = np.arange(len(D["err"]))[sel]
    _, inv = np.unique(keys, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    order = np.argsort(inv, kind="stable")
    parts = np.split(order, np.flatnonzero(np.diff(inv[order])) + 1)
    grp = [ii[p] for p in parts]
    return grp, np.array([int(D["level"][g[0]]) for g in grp])


def recall(
    D: Dict[str, np.ndarray], s: np.ndarray, grp: List[np.ndarray], frac: float, thr: np.ndarray
) -> Tuple[float, int]:
    """Fraction of groups holding a candidate below its threshold (error, deg) whose best-scored
    such candidate is among the max(1, int(n * frac)) best of the group."""
    ok = n = 0
    for g, t in zip(grp, thr):
        hit = D["err"][g] < t
        if not hit.any():
            continue
        nk = max(1, int(len(g) * frac))
        ranks = np.empty(len(g), dtype=int)
        ranks[np.argsort(-s[g], kind="stable")] = np.arange(len(g))
        ok += int(ranks[hit].min() < nk)
        n += 1
    return ok / max(n, 1), n


def evaluate(
    D: Dict[str, np.ndarray], S: Dict[str, np.ndarray], mask: np.ndarray
) -> Dict[str, Any]:
    """Metrics of every score set on the candidates selected by `mask` (a variant subset)."""
    from scipy.stats import spearmanr

    Ds = {k: v[mask] for k, v in D.items()}
    sets = M.eval_sets(Ds)
    grp, glev = level_groups(Ds)
    thr = np.where(glev == 3, 1.0, 3.0)
    gE = M.groups(Ds, "E")
    lab = (Ds["category"] != 3) & (Ds["source"] <= 1)
    out: Dict[str, Any] = {}
    for name, s_full in S.items():
        s = s_full[mask]
        r: Dict[str, Any] = {}
        for k in ("A", "B", "C"):
            m = sets[k]
            y = (Ds["source"][m] == 2) if k == "B" else (Ds["err"][m] < 1.0)
            r[f"auc_{k}"] = M.auc(y.astype(int), s[m])
            r[f"n_{k}"] = int(m.sum())
        for frac in FRACS:
            tag = {0.25: "4", 0.125: "8"}[frac]
            for lv in range(4):
                sel = [g for g, lg in zip(grp, glev) if lg == lv]
                t = thr[glev == lv]
                r[f"recall_L{lv}_1/{tag}"], r[f"n_L{lv}"] = recall(Ds, s, sel, frac, t)
            sel = [g for g, lg in zip(grp, glev) if lg <= 2]
            r[f"recall_L012_1/{tag}"], r["n_L012"] = recall(Ds, s, sel, frac, thr[glev <= 2])
        r["final_precision"], r["n_final"] = M.final_precision(Ds, s, gE)
        r["spearman_ybcost"] = float(spearmanr(s[lab], -Ds["y_bcost"][lab])[0])
        out[name] = r
    return out


def decisions(res: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """D1-D3 from the pooled (both variants) held-out metrics."""
    rec = "recall_L012_1/4"
    rec8 = "recall_L012_1/8"
    out: Dict[str, Any] = {}
    # D1: F-cache Q5 GBT (classifier, E2 protocol) pruning recall < 0.95?
    c5 = res["cache5|clf3"][rec]
    c5r = res["cache5|reg"][rec]
    out["D1"] = dict(
        cache5_clf3_recall=c5, cache5_reg_recall=c5r, e2full_clf3_recall=res["e2full|clf3"][rec],
        low_q_signal_mostly_lost=bool(c5 < 0.95),
    )  # fmt: skip
    # deployable candidates: F-lowQ + cost8 (and F-lowQ alone as the ablation)
    dep = [k for k in res if k.startswith("lowq") and k.endswith("+c8|" + k.split("|")[1])]
    best = max(dep, key=lambda k: res[k][rec])
    e2 = res["e2full|clf3"][rec]
    out["D1"].update(
        best_deployable=best, best_deployable_recall=res[best][rec], e2_gap=e2 - res[best][rec],
        adds_nothing_beyond_e2=bool(e2 - res[best][rec] > 0.01),
    )  # fmt: skip
    # D2: recall at keep 1/8 >= 0.98 for the best deployable model
    best8 = max(dep, key=lambda k: res[k][rec8])
    out["D2"] = dict(model=best8, recall_1_8=res[best8][rec8], holds=bool(res[best8][rec8] >= 0.98))
    # D3: regression vs classifier on the best feature set by recall D (1/4)
    sets = sorted({k.split("|")[0] for k in dep})
    d3 = {}
    for fs in sets:
        reg, clf = res[f"{fs}|reg"][rec], max(res[f"{fs}|clf3"][rec], res[f"{fs}|clf_lm"][rec])
        d3[fs] = dict(reg=reg, best_clf=clf, carry="reg" if reg >= clf - 0.005 else "clf")
    out["D3"] = d3
    return out


def table(res: Dict[str, Dict[str, Any]]) -> str:
    head = (
        f"{'model':18s} {'AUC A':>6s} {'B':>6s} {'C':>6s} | recall 1/4: L0    L1    L2    L3@1  "
        f"L0-2 | recall 1/8: L0    L1    L2    L3@1  L0-2 |  E     rho"
    )
    rows = [head]
    for k, r in res.items():
        rows.append(
            f"{k:18s} {r['auc_A']:6.3f} {r['auc_B']:6.3f} {r['auc_C']:6.3f} |"
            + "".join(f" {r[f'recall_L{i}_1/4']:5.3f}" for i in range(4))
            + f" {r['recall_L012_1/4']:5.3f} |"
            + "".join(f" {r[f'recall_L{i}_1/8']:5.3f}" for i in range(4))
            + f" {r['recall_L012_1/8']:5.3f} | {r['final_precision']:5.3f}"
            + f" {r['spearman_ybcost']:6.3f}"
        )
    return "\n".join(rows)


def retrain(tag: str, name: str, target: str) -> None:
    """Retrain once with the candidates harvested from the proxy run `tag` (endtoend.py harvest)
    added to the training rows of the folds that do not contain their voxel; saved in
    models_retrain/ (same file names as models/)."""
    import joblib

    D = load_all()
    fold = M.fold_of(D["vpos"])
    X = feature_matrix(name, D)
    files = sorted((B.CACHE / "e2e" / f"{tag}_harvest").glob("*.npz"))
    H = [np.load(f) for f in files]
    Xh = np.vstack([h["X"] for h in H])
    err_h = np.concatenate([h["err"] for h in H])
    lvl_h = np.concatenate([h["level"] for h in H])
    fold_h = np.concatenate([np.full(len(h["err"]), int(h["fold"])) for h in H])
    y3h = (err_h < 3.0).astype(np.int8)
    ylmh = np.where(lvl_h == 3, err_h < 1.0, err_h < 3.0).astype(np.int8)
    ybh = np.concatenate([h["y_bcost"] for h in H])
    out = B.CACHE / "models_retrain"
    out.mkdir(parents=True, exist_ok=True)
    for f in range(M.FOLDS):
        tr, trh = fold != f, fold_h != f
        Dc = dict(
            y_bcost=np.concatenate([D["y_bcost"][tr], ybh[trh]]),
            y3=np.concatenate([D["y3"][tr], y3h[trh]]),
            y_lm=np.concatenate([D["y_lm"][tr], ylmh[trh]]),
        )
        Xc = np.vstack([X[tr], Xh[trh]])
        m = fit_target(target, Xc, Dc, np.ones(len(Xc), bool))
        joblib.dump(m, out / f"{name}__{target}__fold{f}.joblib")
        print(f"retrained fold {f}: {len(Xc)} rows ({int(trh.sum())} harvested)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "retrain"])
    ap.add_argument("--reuse-scores", action="store_true", help="skip training, reuse scores.npz")
    ap.add_argument("--tag", default="p_i")
    ap.add_argument("--set", default="lowq5+c8")
    ap.add_argument("--target", default="reg")
    a = ap.parse_args()
    if a.cmd == "retrain":
        return retrain(a.tag, a.set, a.target)
    t0 = time.time()
    D = load_all()
    print(f"dataset {len(D['err'])} candidates; labelled {int(np.isfinite(D['y_bcost']).sum())}")
    sp = B.CACHE / "scores.npz"
    if a.reuse_scores and sp.exists():
        z = np.load(sp)
        S = {k.replace("__", "|"): z[k] for k in z.files}
    else:
        S = cv_scores(D)
        np.savez_compressed(sp, **{k.replace("|", "__"): v for k, v in S.items()})
    out: Dict[str, Any] = {}
    lines: List[str] = []
    for tv in ("both", "clean", "all"):
        mask = np.ones(len(D["err"]), bool)
        if tv != "both":
            mask = D["variant"] == C.VARIANTS.index(tv)
        res = evaluate(D, S, mask)
        out[tv] = res
        lines.append(f"[held-out grains, test on {tv}]\n" + table(res))
    out["decisions"] = decisions(out["both"])
    B.OUT.mkdir(parents=True, exist_ok=True)
    C.save_json(B.OUT / "offline_metrics.json", out)
    (B.OUT / "offline_metrics.txt").write_text("\n\n".join(lines) + "\n")
    print("\n\n".join(lines))
    print(out["decisions"])
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
