"""AR5 step 7: the pure rules of the skip gate (numpy only; no truth, no images).

The key ``BoundaryPassSkipHitRatio h``: a boundary voxel whose own hit ratio (the hit ratio of its
end before the pass) is at least ``h`` gets no try, and keeps its end. The *replay* removes the
skipped voxels' tries from a stored pass; it is exact because the tries of a voxel do not depend
on one another (each starts from a neighbour's orientation in the map as it was before the pass)
and the winner rule ``hit`` reads only the voxel's own old end and its own tries.

The decision rule is the plan's "AR5 step 7" rule, frozen in ``ar5/gate/frozen.json``. Tests:
``tests/test_ar5_gate.py``.
"""

from typing import Any, Dict, Tuple

import numpy as np

H_S = 0.95  # frozen skip cut
MIN_CUT = 0.5  # the pass must be cut by at least 50%
LOSS_SHARE_PCT = 5  # broken minus fixed at most 5% of the pass's net gain (rounded up), >= 1
MEDIAN_TOL_FULL = 0.001  # median error at most 0.001 degree above the full pass
HOUSE_FEWER_PCT = 20  # house rules against S8: 20% fewer evaluations
HOUSE_BROKEN_MINUS_FIXED = 5
HOUSE_MEDIAN_TOL = 0.005
EPS = 1e-12


def skip_mask(own_hit: np.ndarray, h: float) -> np.ndarray:
    """True where the own hit ratio is at least ``h`` (the voxel gets no try)."""
    return np.asarray(own_hit, dtype=float) >= h


def winner_of(
    old_hit: float,
    old_cost: float,
    t_hit: np.ndarray,
    t_cost: np.ndarray,
) -> int:
    """Rule ``hit``: index of the winning try, or -1 for the old end. Highest hit ratio, then lowest
    cost, then the earliest (the old end first)."""
    cost = np.concatenate([[old_cost], t_cost]).astype(float)
    hit = np.concatenate([[old_hit], t_hit]).astype(float)
    return int(np.lexsort((np.arange(len(cost)), cost, -hit))[0]) - 1


def replay(
    old_hit: np.ndarray,
    old_cost: np.ndarray,
    t_voxel: np.ndarray,
    t_hit: np.ndarray,
    t_cost: np.ndarray,
    t_evals: np.ndarray,
    h: float,
) -> Dict[str, Any]:
    """The stored pass with the tries of every voxel with own hit ratio >= ``h`` removed.

    ``old_hit`` / ``old_cost`` (n,) are the scores of the old ends; ``t_*`` (m,) the stored tries.
    Returns ``kept`` (m,) bool, ``evals`` and ``evals_full`` (tries only; the old-end scores are
    not in them), ``cut``, try and voxel counts, ``skipped`` (n,) bool (boundary voxels skipped)
    and ``winner`` (n,): index into the tries of the end kept at each voxel (-1: the old end; also
    for voxels with no try or a skipped voxel)."""
    n = len(old_hit)
    tv = np.asarray(t_voxel, dtype=np.int64)
    skip_v = skip_mask(old_hit, h)
    kept = ~skip_v[tv]
    boundary = np.zeros(n, dtype=bool)
    boundary[tv] = True
    winner = np.full(n, -1, dtype=np.int64)
    order = np.argsort(tv, kind="stable")
    sv = tv[order]
    starts = np.flatnonzero(np.r_[True, sv[1:] != sv[:-1]])
    ends = np.r_[starts[1:], len(sv)]
    for a, b in zip(starts, ends):
        v = int(sv[a])
        if skip_v[v]:
            continue
        idx = order[a:b]  # try indices of voxel v, in the order made
        w = winner_of(float(old_hit[v]), float(old_cost[v]), t_hit[idx], t_cost[idx])
        winner[v] = -1 if w < 0 else int(idx[w])
    evals_full = int(np.asarray(t_evals, dtype=np.int64).sum())
    evals = int(np.asarray(t_evals, dtype=np.int64)[kept].sum())
    return dict(
        kept=kept,
        evals=evals,
        evals_full=evals_full,
        cut=1.0 - evals / evals_full if evals_full else 0.0,
        tries_kept=int(kept.sum()),
        tries_skipped=int((~kept).sum()),
        voxels_skipped=int((skip_v & boundary).sum()),
        boundary_voxels=int(boundary.sum()),
        skipped=skip_v & boundary,
        winner=winner,
    )


def fixed_broken(wrong_ref: np.ndarray, wrong_new: np.ndarray) -> Tuple[int, int]:
    """(fixed, broken) of ``new`` against ``ref``, paired by voxel: fixed = wrong under ref only."""
    a = np.asarray(wrong_ref, dtype=bool)
    b = np.asarray(wrong_new, dtype=bool)
    return int((a & ~b).sum()), int((~a & b).sum())


def net_gain(wrong_before: int, wrong_after: int) -> int:
    """Wrong voxels before the pass minus after it."""
    return int(wrong_before) - int(wrong_after)


def allowed_loss(net: int) -> int:
    """Largest allowed broken-minus-fixed: 5% of the net gain, rounded up (integer arithmetic),
    at least 1 voxel."""
    return max(1, (LOSS_SHARE_PCT * int(net) + 99) // 100)


def gate_rule(
    cut: float,
    wrong_before: int,
    wrong_full: int,
    fixed: int,
    broken: int,
    med_full: float,
    med_skip: float,
    total_skip: int,
    total_s8: int,
    fixed_s8: int,
    broken_s8: int,
    med_s8: float,
) -> Dict[str, Any]:
    """The plan's decision rule on one sample: the skip passes if the pass is cut by at least 50%;
    broken minus fixed against the full pass is within ``allowed_loss``; median error is at most
    0.001 degree higher than the full pass; and the run with the skipped pass meets the house rules
    against S8 (at least 20% fewer evaluations, broken minus fixed at most 5, median at most 0.005
    degree higher)."""
    allowed = allowed_loss(net_gain(wrong_before, wrong_full))
    out: Dict[str, Any] = dict(
        allowed_loss=allowed,
        pass_cut=bool(cut >= MIN_CUT - EPS),
        pass_loss=bool(broken - fixed <= allowed),
        pass_median=bool(med_skip - med_full <= MEDIAN_TOL_FULL + EPS),
        pass_s8_evals=bool(100 * int(total_skip) <= (100 - HOUSE_FEWER_PCT) * int(total_s8)),
        pass_s8_wrong=bool(broken_s8 - fixed_s8 <= HOUSE_BROKEN_MINUS_FIXED),
        pass_s8_median=bool(med_skip - med_s8 <= HOUSE_MEDIAN_TOL + EPS),
    )
    out["passes"] = all(v for k, v in out.items() if k.startswith("pass_"))
    return out
