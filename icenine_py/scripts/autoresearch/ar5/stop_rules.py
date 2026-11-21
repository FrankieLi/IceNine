"""AR5 step 1: the pure rule functions (numpy only; no truth, no images).

The try test of the boundary pass (P1, P2, P3), the choice of a try's end (rule ``hit`` of the
pass), the pick, the hit-ratio flag, the ladder accounting and the cut of the vote count. Every
function takes arrays and returns arrays or numbers; the scripts ``stop_*.py`` read the stored
data and call them. Tie rules are stated in each docstring. Tests: ``tests/test_ar5_stop_rules.py``.
"""

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts" / "common"))
import stats  # noqa: E402

TOL = 1e-9  # a hit-ratio tie within this is a tie
SEED_TEST = 0.80  # the seed test: a result with hit ratio >= this is accepted
MAX_LOSS = 0.05  # a rule may lose at most this share of the net gain of the pass
D_VALUES = (0.0, 0.02, 0.05, 0.10, 0.20)  # P1 settings (fixed in the plan)
M_VALUES = (1, 2)  # P2 settings
THRESHOLDS = tuple(round(0.5 + 0.05 * i, 2) for i in range(10))  # 0.50 .. 0.95


# ---- Part 1: the try test -----------------------------------------------------------------
def p1_keep(start_hit: np.ndarray, old_hit: np.ndarray, d: float) -> np.ndarray:
    """P1: a try is skipped when its start hit ratio is more than ``d`` under the voxel's hit
    ratio (``old_hit``, one value per try). A start exactly ``d`` under is kept."""
    return np.asarray(start_hit) >= np.asarray(old_hit) - d - TOL


def p2_keep(voxel: np.ndarray, start_hit: np.ndarray, start_cost: np.ndarray, m: int) -> np.ndarray:
    """P2: per voxel keep the ``m`` tries with the highest start hit ratio. Ties: the lower start
    cost, then the earlier try."""
    voxel = np.asarray(voxel)
    keep = np.zeros(len(voxel), bool)
    order = np.lexsort((np.arange(len(voxel)), start_cost, -np.asarray(start_hit), voxel))
    v_sorted = voxel[order]
    first = np.r_[True, v_sorted[1:] != v_sorted[:-1]]
    start_of = np.maximum.accumulate(np.where(first, np.arange(len(order)), 0))
    rank = np.arange(len(order)) - start_of
    keep[order[rank < m]] = True
    return keep


def rule_keep(
    kind: str,
    voxel: np.ndarray,
    start_hit: np.ndarray,
    start_cost: np.ndarray,
    old_hit: np.ndarray,
    d: Optional[float] = None,
    m: Optional[int] = None,
) -> np.ndarray:
    """The mask of kept tries for ``kind`` in P1, P2, P3 (P3 = both)."""
    keep = np.ones(len(voxel), bool)
    if kind in ("P1", "P3"):
        assert d is not None
        keep &= p1_keep(start_hit, old_hit, d)
    if kind in ("P2", "P3"):
        assert m is not None
        keep &= p2_keep(voxel, start_hit, start_cost, m)
    if kind not in ("P1", "P2", "P3"):
        raise KeyError(kind)
    return keep


def choose_try(
    voxel: np.ndarray,
    old_hit: np.ndarray,
    old_cost: np.ndarray,
    t_hit: np.ndarray,
    t_cost: np.ndarray,
    keep: np.ndarray,
) -> np.ndarray:
    """Per voxel (length ``len(old_hit)``) the index of the try whose end replaces the voxel's
    end, or -1 (the old end stays). Rule ``hit`` of the pass: highest hit ratio, then lowest cost,
    then the old end first, then the earlier try. Only kept tries compete."""
    voxel = np.asarray(voxel)
    keep = np.asarray(keep, bool)
    out = np.full(len(old_hit), -1, dtype=np.int64)
    idx = np.where(keep)[0]
    if len(idx) == 0:
        return out
    order = idx[np.lexsort((idx, t_cost[idx], -t_hit[idx], voxel[idx]))]  # best first per voxel
    v_sorted = voxel[order]
    first = np.r_[True, v_sorted[1:] != v_sorted[:-1]]
    best = order[first]
    bv = voxel[best]
    # the old end wins ties (hit and cost both equal)
    better = (t_hit[best] > old_hit[bv] + 0.0) | (
        (t_hit[best] == old_hit[bv]) & (t_cost[best] < old_cost[bv])
    )
    out[bv[better]] = best[better]
    return out


def pass_evals(t_evals: np.ndarray, keep: np.ndarray) -> int:
    """Evaluations of the pass under a try test: one start evaluation per try (kept or not), plus
    the evaluations of the kept tries."""
    return int(len(t_evals) + np.asarray(t_evals)[np.asarray(keep, bool)].sum())


def net_gain(wrong_before: np.ndarray, wrong_after: np.ndarray) -> int:
    """Wrong voxels before the pass minus wrong voxels after it."""
    return int(np.sum(wrong_before) - np.sum(wrong_after))


def paired_vs_full(
    wrong_before: np.ndarray, wrong_full: np.ndarray, wrong_rule: np.ndarray
) -> Dict[str, int]:
    """Voxels fixed (wrong after the full pass, right under the rule), broken (right after the
    full pass, wrong under the rule) and wrong under both. ``wrong_before`` is not used in the
    split; it is accepted so that callers pass the three maps together."""
    full, rule = np.asarray(wrong_full, bool), np.asarray(wrong_rule, bool)
    return dict(
        fixed=int(np.sum(full & ~rule)),
        broken=int(np.sum(~full & rule)),
        both_wrong=int(np.sum(full & rule)),
    )


def pick_setting(
    rows: Sequence[Dict[str, Any]], full_net_gain: float, max_loss: float = MAX_LOSS
) -> Optional[Dict[str, Any]]:
    """The row with the fewest evaluations among those with ``net_gain >= (1 - max_loss) *
    full_net_gain``; ties go to the earlier row. None when no row qualifies."""
    ok = [r for r in rows if r["net_gain"] >= (1.0 - max_loss) * full_net_gain - TOL]
    if not ok:
        return None
    return min(ok, key=lambda r: r["evals"])  # min keeps the first of equal keys


def pick_d(
    rows: Sequence[Dict[str, Any]], full_net_gain: float, max_loss: float = MAX_LOSS
) -> Dict[str, Any]:
    """The P1 row P3 takes its ``d`` from: the pick of ``pick_setting``; if no P1 setting loses
    at most ``max_loss``, the row with the highest net gain (the earlier on a tie)."""
    r = pick_setting(rows, full_net_gain, max_loss)
    if r is not None:
        return r
    return max(rows, key=lambda x: x["net_gain"])


def try_test_passes(
    evals_rule: float, evals_full: float, gain_rule: float, gain_full: float
) -> bool:
    """The pass rule of the try test on ``mid``: at least half of the pass's evaluations cut and
    at most 5% of its net gain lost."""
    return bool(
        evals_rule <= 0.5 * evals_full + TOL and gain_rule >= (1.0 - MAX_LOSS) * gain_full - TOL
    )


# ---- Part 2: the ladder -------------------------------------------------------------------
def ladder_case(
    r1: Dict[str, float], r3: Dict[str, float], grid_evals: float, test: float = SEED_TEST
) -> Dict[str, Any]:
    """One case through the ladder R1, R3, grid. A rung passes when its end hit ratio is at least
    ``test``. Evaluations add up over the rungs run; if both fail, the full grid search is paid
    on top (its own level 0 is paid again; the library cannot continue after level 0)."""
    ev = float(r1["evals"])
    if r1["hit"] >= test - TOL:
        return dict(rung="R1", evals=ev)
    ev += float(r3["evals"])
    if r3["hit"] >= test - TOL:
        return dict(rung="R3", evals=ev)
    return dict(rung="grid", evals=ev + float(grid_evals))


def ladder_verdict(
    ladder_evals: float, grid_evals: float, right_ladder: int, right_grid: int
) -> Dict[str, bool]:
    """The pass rule of Part 2: at most two thirds of the grid fallbacks' evaluations, and at
    most one case fewer right than the grid search."""
    ev_ok = ladder_evals <= (2.0 / 3.0) * grid_evals + TOL
    right_ok = right_ladder >= right_grid - 1
    return dict(evals_ok=bool(ev_ok), right_ok=bool(right_ok), passes=bool(ev_ok and right_ok))


# ---- Part 3: is a voxel finished? -----------------------------------------------------------
def flag_counts(hit: np.ndarray, wrong: np.ndarray, thr: float) -> Dict[str, int]:
    """A voxel is flagged when its hit ratio is under ``thr`` (strictly, beyond ``TOL``)."""
    flagged = np.asarray(hit) < thr - TOL
    wrong = np.asarray(wrong, bool)
    return dict(
        flagged=int(flagged.sum()),
        wrong_flagged=int((flagged & wrong).sum()),
        wrong_unflagged=int((~flagged & wrong).sum()),
        n=int(len(wrong)),
        n_wrong=int(wrong.sum()),
    )


def pick_threshold(
    hits: Sequence[np.ndarray],
    wrongs: Sequence[np.ndarray],
    thresholds: Sequence[float] = THRESHOLDS,
    recall: float = 0.8,
) -> Optional[float]:
    """The lowest threshold whose flag catches at least ``recall`` of the wrong voxels, pooled
    over the samples given. None when no threshold does."""
    n_wrong = sum(int(np.sum(w)) for w in wrongs)
    for t in sorted(thresholds):
        caught = sum(flag_counts(h, w, t)["wrong_flagged"] for h, w in zip(hits, wrongs))
        if n_wrong > 0 and caught >= recall * n_wrong - TOL:
            return float(t)
    return None


def hit_after(
    old_hit: np.ndarray, voxel: np.ndarray, t_hit: np.ndarray, choice: np.ndarray
) -> np.ndarray:
    """Hit ratio of every voxel after the pass: the taken try's end, else the old one."""
    out = np.asarray(old_hit, dtype=float).copy()
    took = np.where(np.asarray(choice) >= 0)[0]
    out[took] = np.asarray(t_hit)[np.asarray(choice)[took]]
    return out


def sample_unflagged(flagged: np.ndarray, boundary: np.ndarray, every: int = 20) -> np.ndarray:
    """Every ``every``-th unflagged boundary voxel in index order (the ``every``-th, the
    ``2*every``-th ...)."""
    un = np.where(~np.asarray(flagged, bool) & np.asarray(boundary, bool))[0]
    return un[every - 1 :: every]


# ---- the vote count (H3) --------------------------------------------------------------------
def best_cut(wrong_counts: np.ndarray, right_counts: np.ndarray) -> Dict[str, Any]:
    """The cut on a count (``count >= cut`` predicts wrong) that maximizes the share of wrong
    voxels at or above it minus the share of right voxels at or above it (Youden's J). Ties: the
    lowest cut. Chosen after the result and labelled so on the page."""
    w, r = np.asarray(wrong_counts), np.asarray(right_counts)
    best: Optional[Dict[str, Any]] = None
    for c in np.unique(np.concatenate([w, r])):
        wg, rg = int((w >= c).sum()), int((r >= c).sum())
        j = wg / max(len(w), 1) - rg / max(len(r), 1)
        if best is None or j > best["j"] + 1e-12:
            best = dict(cut=float(c), wrong_ge=wg, right_ge=rg, j=float(j))
    assert best is not None
    best.update(n_wrong=int(len(w)), n_right=int(len(r)))
    return best


def separates(wrong_ge: int, n_wrong: int, right_ge: int, n_right: int) -> bool:
    """True when the Wilson 95% intervals of the share at or above the cut do not overlap for the
    wrong and the right voxels."""
    lo_w, hi_w = stats.wilson(wrong_ge, n_wrong)
    lo_r, hi_r = stats.wilson(right_ge, n_right)
    return bool(lo_w > hi_r or lo_r > hi_w)


def wilson_row(k: int, n: int) -> List[float]:
    lo, hi = stats.wilson(k, n)
    return [k, n, lo, hi]
