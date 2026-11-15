"""AR5 step 18 (task ar5r): the accounting, threshold procedure, gate and pass rules of the learned
skip of BFS local optimizations (pure; numpy and the standard library only).

Terms are those of ``learnskip_features``. **Wrong skip**: a skipped fit that ``J5`` accepted.
**Changed final answer**: a skipped fit that ``J5`` rejected and whose end was the voxel's final
answer before the pass (step 17's review). **Net saving**: evaluations of skipped fits that ``J5``
rejected, minus every check (one counted evaluation per fit). A model call counts zero.

Threshold procedure (frozen): skip a fit if its score (the predicted chance that ``J5`` accepts it)
is under the threshold. Each training sample is scored by a model trained on the other training
samples (out of fold). Candidates: the quantiles at 0, 0.1%, ..., 80% of all out-of-fold scores.
Take the one with the largest sum of net shares over the training samples among those where every
training sample has wrong skips at most 0.1% of its fits and changed final answers at most 2; none:
no skip. A tie in the sum takes the smaller threshold (the first in ascending order).
"""

from fractions import Fraction
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

import fitskip_rules as FR

SAMPLES = FR.SAMPLES  # held3, held, mid, grad, big
EXPECT = FR.EXPECT  # stored counts of arm J5 (totals with the pass, fit counts)
FITS = dict(held3=72_458, held=76_911, mid=52_027, grad=53_398, big=32_543)
WRONG_MAX_SHARE = Fraction(1, 1000)
FC_MAX = 2
QUANTILES = np.linspace(0.0, 0.8, 801)
FEATURE_COUNT = 16
MODEL_SPEC = dict(
    library="scikit-learn 1.6.1 (training only; the trees are exported to JSON)",
    model="HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, random_state=0), "
    "all else at the 1.6.1 defaults (31 leaves, 20 rows a leaf, no L2, 255 bins, early stopping "
    "on over 10,000 rows: a 10% validation split drawn with seed 0, stop after 10 rounds with no "
    "loss gain)",
    max_iter=300,
    learning_rate=0.1,
    random_state=0,
    label="J5 accepted the fit; one row per fit of J5r, unweighted",
    score="the predicted chance that the fit is accepted",
)

# the plan's sizing table: per sample the learned skip (net share %, wrong skips, allowed wrong
# skips, changed final answers), S(0.3, both) (net share %, wrong, changed) and the ceiling %
PLAN_TABLE: Dict[str, Dict[str, Any]] = dict(
    held3=dict(
        share=5.46, wrong=30, allowed=72, fc=1, s_share=3.98, s_wrong=21, s_fc=1, ceil=13.31
    ),
    held=dict(share=6.16, wrong=70, allowed=76, fc=0, s_share=4.21, s_wrong=34, s_fc=0, ceil=14.90),
    mid=dict(share=5.79, wrong=46, allowed=52, fc=0, s_share=3.62, s_wrong=26, s_fc=0, ceil=16.75),
    grad=dict(
        share=11.83, wrong=43, allowed=53, fc=0, s_share=8.74, s_wrong=24, s_fc=0, ceil=20.31
    ),
    big=dict(share=5.38, wrong=17, allowed=32, fc=0, s_share=3.20, s_wrong=6, s_fc=0, ceil=14.19),
)
SIZING_GATE_MODEL_MID_THRESHOLD = 0.00957  # the plan's number for the gate model of ``mid``
SIZING_FINAL_THRESHOLD = 0.01112
SIZING_FINAL_ROUNDS = 246

GATE: Dict[str, Any] = dict(
    g0="prediction 1 (the tracked script reproduces the plan's table to the evaluation) and J5r "
    "has 72,458 / 76,911 / 52,027 / 53,398 / 32,543 fits",
    g1="on each held-out sample: wrong skips at most 0.1% of its fits, changed final answers at "
    "most 2, net saving at least 4.5% (9% on grad)",
    g1_net_min=dict(held3=0.045, held=0.045, mid=0.045, grad=0.09, big=0.045),
    g2="the live J5L run on mid with the gate model, up to its first wrong skip in run order, "
    "equals the replay fit by fit (features to 1e-9, scores, skips, ends, evaluations), and over "
    "the whole run if it has none",
    g3="that run's whole mid layer with its pass passes rules (1) to (5)",
    fails="stop, report, no tuning, no other model or threshold, no whole layer",
)
JUDGE: Dict[str, Any] = dict(
    r1_total_fewer=0.03,
    r2_broken_minus_fixed_max=5,
    r3_median_error_higher_max_deg=0.001,
    r4_wrong_seeds_more_max=1,
    r4_grid_within=1,
    r4_reseeds_more_max_share=0.02,
    r5_vs_s8=dict(fewer=0.20, broken_minus_fixed_max=5, median_error_higher_max_deg=0.005),
    r5_note="if J5 fails (5) on a sample, J5L is judged there by (1) to (4)",
    pass_samples=("held4", "held3", "held", "mid", "grad", "big"),
    step_passes="J5L passes on held4 and on the five disclosed samples; fails on held4: no key",
)
PREDICTIONS = (
    "1. The tracked script reproduces the plan's table to the evaluation. 2. Live mid (gate "
    "model): matches its replay to its first wrong skip; 4.5% to 7% under J5. 3. J5L (final "
    "model) against J5, with the pass: held4 4% to 7% fewer; held3, held, mid, big 5% to 7.5%; "
    "grad 10% to 13.5%; broken minus fixed -2 to 4; median error up at most 0.0005 degree; "
    "reseeds at most 2% more. 4. mid passes rule 5 (about 1 in 2)."
)


class SampleData:
    """One sample's features and labels (rows in run order)."""

    def __init__(
        self,
        name: str,
        X: np.ndarray,
        accepted: np.ndarray,
        evals: np.ndarray,
        final_before: np.ndarray,
    ) -> None:
        self.name = name
        self.X = X
        self.accepted = np.asarray(accepted, bool)
        self.evals = np.asarray(evals, np.int64)
        self.final_before = np.asarray(final_before, bool)
        self.total = EXPECT[name]["total"]


def account(skip: np.ndarray, d: SampleData) -> Dict[str, Any]:
    """The counts of a skip vector (bool, run order) on a sample, as the replay of step 17: a
    skipped fit that J5 rejected saves its evaluations; a wrong skip earns nothing; every fit is
    checked once."""
    skip = np.asarray(skip, bool)
    rej = skip & ~d.accepted
    saved = int(d.evals[rej].sum())
    n = len(d.accepted)
    net = saved - n
    return dict(
        n=n,
        skipped=int(skip.sum()),
        wrong=int((skip & d.accepted).sum()),
        changed=int((rej & d.final_before).sum()),
        saved=saved,
        checks=n,
        net=net,
        share=Fraction(net, d.total),
    )


def within_limits(a: Dict[str, Any]) -> bool:
    """Wrong skips at most 0.1% of the sample's fits; changed final answers at most 2."""
    return bool(Fraction(a["wrong"]) <= WRONG_MAX_SHARE * a["n"] and a["changed"] <= FC_MAX)


def choose_threshold(
    oof: Dict[str, np.ndarray], data: Dict[str, SampleData]
) -> Dict[str, Optional[float]]:
    """The frozen procedure on out-of-fold scores ``oof[sample]``. Returns ``threshold`` (None: no
    skip), ``share_sum`` (float) and ``n_candidates``."""
    allsc = np.concatenate([oof[s] for s in oof])
    grid = np.unique(np.quantile(allsc, QUANTILES))
    best_t: Optional[float] = None
    best_sum: Optional[Fraction] = None
    for t in grid:
        tot = Fraction(0)
        ok = True
        for s in oof:
            a = account(oof[s] < t, data[s])
            if not within_limits(a):
                ok = False
                break
            tot += a["share"]
        if ok and (best_sum is None or tot > best_sum):
            best_t, best_sum = float(t), tot
    return dict(
        threshold=best_t,
        share_sum=float(best_sum) if best_sum is not None else None,
        n_candidates=int(len(grid)),
    )


def skip_vector(scores: np.ndarray, threshold: Optional[float]) -> np.ndarray:
    """Skip where the score is under the threshold; none when there is no threshold."""
    if threshold is None:
        return np.zeros(len(scores), bool)
    return np.asarray(scores) < threshold


def gate_g1(per_sample: Dict[str, Dict[str, Any]]) -> Dict[str, bool]:
    """Per held-out sample: within the limits and net saving at least the plan's share."""
    return {
        s: bool(within_limits(a) and a["share"] >= Fraction(str(GATE["g1_net_min"][s])))
        for s, a in per_sample.items()
    }


def matches_plan_table(
    s: str, a: Dict[str, Any], s_acct: Dict[str, Any], ceil: Fraction
) -> Dict[str, bool]:
    """The held-out result and the S(0.3, both) result against the plan's table (shares at two
    decimals in percent, the counts as integers)."""
    t = PLAN_TABLE[s]

    def pct(x: Fraction) -> float:
        return round(float(x) * 100, 2)

    return dict(
        share=pct(a["share"]) == t["share"],
        wrong=a["wrong"] == t["wrong"],
        allowed=int(WRONG_MAX_SHARE * a["n"]) == t["allowed"],
        changed=a["changed"] == t["fc"],
        s_share=pct(s_acct["share"]) == t["s_share"],
        s_wrong=s_acct["wrong"] == t["s_wrong"],
        s_changed=s_acct["changed"] == t["s_fc"],
        ceiling=pct(ceil) == t["ceil"],
    )


def ceiling(d: SampleData) -> Fraction:
    """Skip every J5-rejected fit free: its evaluations as a share of the J5 total (no checks)."""
    return Fraction(int(d.evals[~d.accepted].sum()), d.total)


def s_rule_skip(d: SampleData, t: float = 0.3) -> np.ndarray:
    """Step 17's S(t, both) on the features: start hit ratio under t times best (columns 4, 3)."""
    return d.X[:, 3] < t * d.X[:, 2]


def first_wrong_skip(skip: np.ndarray, accepted: np.ndarray) -> Optional[int]:
    """Index (run order) of the first skipped fit that J5 accepted, else None."""
    w = np.where(np.asarray(skip, bool) & np.asarray(accepted, bool))[0]
    return int(w[0]) if len(w) else None


def judge(
    j5l: Dict[str, Any], j5: Dict[str, Any], s8: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """The pass rules (1) to (5) on summaries made by the step's summary script. Each summary has
    ``total`` (evaluations with the pass), ``broken`` and ``fixed`` against the other arm are given
    in ``j5l['vs_j5']`` and ``j5l['vs_s8']`` (dicts with ``broken``, ``fixed``), ``median``
    (degree),
    ``wrong_seeds``, ``grid``, ``reseeds``. Pure; not used before the runs exist."""
    r: Dict[str, Any] = {}
    fewer = 1 - j5l["total"] / j5["total"]
    r["r1"] = bool(fewer >= JUDGE["r1_total_fewer"])
    r["r1_fewer"] = fewer
    b, f = j5l["vs_j5"]["broken"], j5l["vs_j5"]["fixed"]
    r["r2"] = bool(b - f <= JUDGE["r2_broken_minus_fixed_max"])
    r["r3"] = bool(j5l["median"] - j5["median"] <= JUDGE["r3_median_error_higher_max_deg"])
    r["r4"] = bool(
        j5l["wrong_seeds"] - j5["wrong_seeds"] <= JUDGE["r4_wrong_seeds_more_max"]
        and j5l["grid"] <= j5["grid"] + JUDGE["r4_grid_within"]
        and j5l["reseeds"] <= j5["reseeds"] * (1 + JUDGE["r4_reseeds_more_max_share"])
    )
    if s8 is not None:
        c = JUDGE["r5_vs_s8"]
        b5, f5 = j5l["vs_s8"]["broken"], j5l["vs_s8"]["fixed"]
        r["r5"] = bool(
            1 - j5l["total"] / s8["total"] >= c["fewer"]
            and b5 - f5 <= c["broken_minus_fixed_max"]
            and j5l["median"] - s8["median"] <= c["median_error_higher_max_deg"]
        )
        j5_b, j5_f = j5["vs_s8"]["broken"], j5["vs_s8"]["fixed"]
        j5_r5 = bool(
            1 - j5["total"] / s8["total"] >= c["fewer"]
            and j5_b - j5_f <= c["broken_minus_fixed_max"]
            and j5["median"] - s8["median"] <= c["median_error_higher_max_deg"]
        )
        r["j5_passes_r5"] = j5_r5
        r["passes"] = bool(r["r1"] and r["r2"] and r["r3"] and r["r4"] and (r["r5"] or not j5_r5))
    else:
        r["passes"] = bool(r["r1"] and r["r2"] and r["r3"] and r["r4"])
    return r
