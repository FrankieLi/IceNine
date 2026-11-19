"""AR5 step 25, phase 0 (gates s0 and s1): stall point D on stored J5D rows. Read only.

Reads, per disclosed sample, ``J5D_<s>/{e1_log,mc_traces}.json.gz``, ``J5D_<s>/run.json`` and
``stored_J5SD_<s>/repair.json`` under the cache directory given as an argument. Writes the three
plan tables, the per-sample AUC/saving numbers and the s1 verdict. No truth is read, no run is
made. Stall point D, rescued and never step follow the plan text ("AR5 step 25", Terms).
The logistic model and the depth-3 tree are fitted here with numpy/scipy (scikit-learn is not in
the locked environment); the boosted tree is NOT reproduced.

Usage: uv run python scripts/autoresearch/ar5/es_sizing.py <cache_dir> <out_dir>
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from doc_tables import write_tables  # noqa: E402
from stats import wilson  # noqa: E402

TAU = 0.95
GUARD = 0.80
TRAIN = ("held3", "held4", "held5", "held6")
SAMPLES = TRAIN + ("held", "mid", "grad", "big")
COLS = [
    "h1",
    "h2",
    "gap",
    "cost1",
    "cma_max",
    "cma_mean",
    "lift_max",
    "lift_best",
    "ng2_best",
    "lastg_best",
    "gains_best",
    "n_above80",
    "reseed",
    "evD",
]
Row = Dict[str, float]


def _guard(name: str) -> None:
    if "held2" in name or "fine7" in name or "held7" in name:
        raise SystemExit("refused: sample not allowed in this step")


def j5sd_total(cache: Path, s: str) -> Tuple[int, Dict[str, Any]]:
    st = json.load(open(cache / f"J5D_{s}" / "run.json"))["stats"]
    rep = json.load(open(cache / f"stored_J5SD_{s}" / "repair.json"))
    seed = sum(
        int(st[k])
        for k in ("n_evals_seed", "n_evals_reseed", "n_evals_neighbor", "n_evals_revisit")
    )
    return seed + int(rep["pass_evals"]) + int(rep["n_evals"]), st


def extract(cache: Path, s: str) -> Tuple[List[Row], Dict[str, Any]]:
    """Rows of the seed steps that reach D, plus per-sample counts."""
    _guard(s)
    e1 = json.load(gzip.open(cache / f"J5D_{s}" / "e1_log.json.gz", "rt"))
    mc = json.load(gzip.open(cache / f"J5D_{s}" / "mc_traces.json.gz", "rt"))
    total, st = j5sd_total(cache, s)
    bym = {(m["step"] - 1, m["j"]): m for m in mc}
    rows: List[Row] = []
    bad_map = 0
    for i, e in enumerate(e1):
        evD = e["costed"]
        tried = {f[0] for f in e["fits"]}
        stateD: List[Tuple[float, float]] = []
        cma_hit: Dict[int, float] = {}
        mc_info: Dict[int, Tuple[float, float, int, float, int]] = {}
        found = False
        for f in e["fits"]:
            pos, j, kind, hit, cost, ev = f[:6]
            if kind == "C":
                evD += ev
                stateD.append((cost, hit))
                cma_hit[pos] = hit
                found |= hit >= TAU
                continue
            m = bym.get((i, j))
            if m is None or m["voxel"] != e["voxel"] or m["evals"] != ev:
                bad_map += 1
                continue
            tr = m["trace"]
            if m["checks"]:  # first PH(500) check
                c_ev, c_hit = m["checks"][0][0], m["checks"][0][1]
                tr = [t for t in tr if t[0] <= c_ev]
                cD = tr[-1][1] if tr else cost
                evD += min(ev, c_ev)
                stateD.append((cD, c_hit))
                h = c_hit
            else:  # ended before any check
                evD += ev
                stateD.append((cost, hit))
                h = hit
            mc_info[pos] = (
                h,
                0.0,
                sum(1 for t in tr if t[0] > 208),
                tr[-1][0] if tr else 0,
                len(tr),
            )
            found |= h >= TAU
        if found or len(tried) < 5:
            continue
        hs = sorted((h for _, h in stateD), reverse=True)
        ch = [cma_hit[p] for p in sorted(cma_hit)]
        takenD = min(stateD)
        takenJ = min((f[4], f[3]) for f in e["fits"])
        lifts = [mc_info[p][0] - cma_hit[p] for p in mc_info if p in cma_hit]
        bp = max(mc_info, key=lambda p: mc_info[p][0]) if mc_info else 0
        bi = mc_info.get(bp, (0, 0, 0, 0, 0))
        rows.append(
            dict(
                y=int(e["stop"] >= 0),
                reseed=int(e["reseed"]),
                h1=hs[0],
                h2=hs[1],
                gap=hs[0] - hs[1],
                cost1=takenD[0],
                cma_max=max(ch),
                cma_mean=float(np.mean(ch)),
                lift_max=max(lifts) if lifts else 0,
                lift_best=bi[0] - cma_hit.get(bp, 0),
                ng2_best=bi[2],
                lastg_best=bi[3],
                gains_best=bi[4],
                n_above80=sum(h >= GUARD for h in hs),
                evD=evD,
                after=e["evals"] - evD,
                hitD=takenD[1],
                hitJ=takenJ[1],
            )
        )
    info = dict(
        steps=len(e1),
        j5sd=total,
        bad_map=bad_map,
        j5d_seed_fail=int(st.get("index_fallback_seed_test", 0)),
    )
    return rows, info


# ---------------------------------------------------------------- small models (numpy/scipy)
def auc(y: np.ndarray, p: np.ndarray) -> float:
    """P(random positive scores above random negative), ties count half (rank sum)."""
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(p, kind="mergesort")
    ps = p[order]
    ranks = np.empty(len(p))
    i = 0
    while i < len(p):
        j = i
        while j + 1 < len(p) and ps[j + 1] == ps[i]:
            j += 1
        ranks[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def fit_logistic(X: np.ndarray, y: np.ndarray, C: float = 1.0) -> np.ndarray:
    """L2 logistic regression (penalty 0.5|w|^2 + C*sum loss, free intercept), L-BFGS."""
    A = np.hstack([X, np.ones((len(X), 1))])

    def f(w: np.ndarray) -> Tuple[float, np.ndarray]:
        z = A @ w
        loss = C * np.sum(np.logaddexp(0, z) - y * z) + 0.5 * np.sum(w[:-1] ** 2)
        g = C * A.T @ (1 / (1 + np.exp(-z)) - y)
        g[:-1] += w[:-1]
        return loss, g

    return minimize(
        f, np.zeros(A.shape[1]), jac=True, method="L-BFGS-B", options=dict(maxiter=2000)
    ).x


def predict_logistic(w: np.ndarray, X: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-(np.hstack([X, np.ones((len(X), 1))]) @ w)))


def fit_tree(X: np.ndarray, y: np.ndarray, depth: int = 3, min_leaf: int = 40) -> Any:
    """CART (gini), depth <= 3, at least 40 per leaf. A node is (feature, threshold, left,
    right); a leaf is the fraction of rescued steps."""
    n = len(y)
    if depth == 0 or n < 2 * min_leaf or y.min() == y.max():
        return float(y.mean())
    best: Optional[Tuple[float, int, float]] = None
    for k in range(X.shape[1]):
        o = np.argsort(X[:, k], kind="mergesort")
        xs, cs = X[o, k], np.cumsum(y[o])
        for i in range(min_leaf, n - min_leaf + 1):
            if xs[i - 1] == xs[i]:
                continue
            pl, pr = cs[i - 1] / i, (cs[-1] - cs[i - 1]) / (n - i)
            imp = (i * 2 * pl * (1 - pl) + (n - i) * 2 * pr * (1 - pr)) / n
            if best is None or imp < best[0] - 1e-15:
                best = (imp, k, (xs[i - 1] + xs[i]) / 2)
    p = y.mean()
    if best is None or best[0] >= 2 * p * (1 - p) - 1e-12:
        return float(p)
    _, k, t = best
    m = X[:, k] <= t
    return (
        k,
        t,
        fit_tree(X[m], y[m], depth - 1, min_leaf),
        fit_tree(X[~m], y[~m], depth - 1, min_leaf),
    )


def predict_tree(node: Any, X: np.ndarray) -> np.ndarray:
    if not isinstance(node, tuple):
        return np.full(len(X), node)
    k, t, left, right = node
    out = np.empty(len(X))
    m = X[:, k] <= t
    out[m] = predict_tree(left, X[m])
    out[~m] = predict_tree(right, X[~m])
    return out


# ---------------------------------------------------------------- evaluation
def mat(rows: Sequence[Row]) -> Tuple[np.ndarray, np.ndarray]:
    return (
        np.array([[r[c] for c in COLS] for r in rows], float),
        np.array([r["y"] for r in rows], int),
    )


def logistic_fold(train: Sequence[Row], test: Sequence[Row], total: int) -> Dict[str, float]:
    Xa, ya = mat(train)
    Xb, yb = mat(test)
    mu, sd = Xa.mean(0), Xa.std(0) + 1e-9
    w = fit_logistic((Xa - mu) / sd, ya)
    pa, pb = predict_logistic(w, (Xa - mu) / sd), predict_logistic(w, (Xb - mu) / sd)
    t = float(np.quantile(pa[ya == 1], 0.10))  # stops at most ~10% of training rescued steps
    stop = (pb < t) & (np.array([r["hitD"] for r in test]) >= GUARD)
    after = np.array([r["after"] for r in test])
    return dict(
        auc=auc(yb, pb),
        save=100 * after[stop].sum() / total,
        resc_stopped=int((stop & (yb == 1)).sum()),
        rescued=int(yb.sum()),
        t=t,
    )


def tree_fold(train: Sequence[Row], test: Sequence[Row]) -> float:
    Xa, ya = mat(train)
    Xb, yb = mat(test)
    return auc(yb, predict_tree(fit_tree(Xa, ya), Xb))


def main(cache_dir: str, out_dir: str) -> None:
    cache, out = Path(cache_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    data: Dict[str, List[Row]] = {}
    info: Dict[str, Dict[str, Any]] = {}
    for s in SAMPLES:
        data[s], info[s] = extract(cache, s)
    t1 = [
        "| sample | steps | reach D | rescued | never | after D, never | after D, rescued "
        "| seed test if all stop |",
        "|---|---|---|---|---|---|---|---|",
    ]
    t3 = ["| sample | R0 saving | lost 0.95 ends | changed ends |", "|---|---|---|---|"]
    for s in SAMPLES:
        r, inf = data[s], info[s]
        y = np.array([x["y"] for x in r])
        aft = np.array([x["after"] for x in r])
        hD = np.array([x["hitD"] for x in r])
        hJ = np.array([x["hitJ"] for x in r])
        tot = inf["j5sd"]
        new_fail = int(((hD < GUARD) & (hJ >= GUARD)).sum())
        t1.append(
            f"| {s} | {inf['steps']:,} | {len(r)} | {int(y.sum())} | {int((1 - y).sum())} "
            f"| {100 * aft[y == 0].sum() / tot:.2f}% | {100 * aft[y == 1].sum() / tot:.2f}% "
            f"| {new_fail} ({inf['j5d_seed_fail']}) |"
        )
        if s in TRAIN + ("held", "mid"):
            stop = hD >= GUARD
            lost = int((stop & (y == 1)).sum())
            lo, hi = wilson(lost, inf["steps"])
            t3.append(
                f"| {s} | {100 * aft[stop].sum() / tot:.2f}% | {lost} "
                f"({100 * lo:.2f}% to {100 * hi:.2f}%) | {int((stop & (hJ != hD)).sum())} |"
            )
    t2 = [
        "| sample | AUC logistic | AUC tree | AUC boosted | saving (logistic) | rescued stopped |",
        "|---|---|---|---|---|---|",
    ]
    verdict: Dict[str, Dict[str, float]] = {}
    for s in SAMPLES[:-1]:
        train = sum((data[x] for x in TRAIN if x != s), [])
        lg = logistic_fold(train, data[s], info[s]["j5sd"])
        ta = tree_fold(train, data[s])
        verdict[s] = dict(auc_lr=lg["auc"], auc_tree=ta, save=lg["save"])
        t2.append(
            f"| {s} | {lg['auc']:.3f} | {ta:.3f} | not reproduced | {lg['save']:.2f}% "
            f"| {lg['resc_stopped']} of {lg['rescued']} |"
        )
    need = TRAIN + ("held", "mid")
    ok = all(verdict[s]["auc_lr"] >= 0.70 and verdict[s]["save"] >= 4.0 for s in need)
    s1 = "PASS" if ok else "FAIL"
    write_tables(
        out / "tables.md",
        {"es_sizing": "\n".join(t1), "es_predictor": "\n".join(t2), "es_r0": "\n".join(t3)},
    )
    (out / "s1_verdict.json").write_text(
        json.dumps(dict(s1=s1, per_sample=verdict, info=info), indent=1, default=float)
    )
    print("\n".join(t1 + [""] + t2 + [""] + t3 + ["", f"s1: {s1}"]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
