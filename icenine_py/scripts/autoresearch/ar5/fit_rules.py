"""AR5 step 2: the pure rule functions (numpy only; no truth, no images).

The early stop of the seed step, the paired split, the pick rule and mid's pass rule of Part 1
(a seed refinement whose step matches the error of its start), and the pass rule of run S10.
The scripts ``fit_*.py`` read the stored data and call them. Tests: ``tests/test_ar5_fit_rules.py``.
All comparisons that decide a rule use integers or an explicit 1e-12 slack, so the edge of a rule
does not depend on float rounding.
"""

from typing import Dict, Optional, Sequence, Tuple

import numpy as np

SEED_TEST = 0.80  # the seed test: a result with hit ratio >= this is accepted
TAU = 0.95  # the early stop of the seed step (run S8)
LOSS_FRAC = 0.01  # Part 1: an arm may lose at most this share of all steps against M
MID_FALL = 0.80  # Part 1 pass rule on mid: seed-step evaluations fall by at least this
S10_FALL = 0.20  # Part 2: at least this share fewer evaluations than S8 with the pass
S10_NET = 5  # Part 2: broken minus fixed voxels, at most
S10_MED = 0.005  # Part 2: median error at most this many degrees above S8's
EPS = 1e-12


def stop_position(hits: Sequence[float], tau: float = TAU) -> Optional[int]:
    """Position of the first end with hit ratio >= ``tau`` (the seed step stops there), or None."""
    for i, h in enumerate(hits):
        if h >= tau:
            return i
    return None


def early_stop(hits: Sequence[float], costs: Sequence[float], tau: float = TAU) -> Tuple[int, bool]:
    """(position taken, stopped early) among the ends of the refinements, in refinement order.
    The first end with hit ratio >= ``tau`` is kept (later ends are not refined). If none, the
    lowest final cost of all ends; a cost tie takes the earlier end."""
    p = stop_position(hits, tau)
    if p is not None:
        return p, True
    return int(np.argmin(np.asarray(costs, dtype=float))), False  # argmin: first minimum


def seed_test(hit: float) -> bool:
    """The seed test (hit ratio of the seed result at least 0.80)."""
    return bool(hit >= SEED_TEST)


def paired_split(a: Sequence[bool], b: Sequence[bool]) -> Dict[str, int]:
    """Paired counts of two boolean flags per step: both, only a, only b, neither."""
    a_ = np.asarray(a, dtype=bool)
    b_ = np.asarray(b, dtype=bool)
    return dict(
        both=int((a_ & b_).sum()),
        only_a=int((a_ & ~b_).sum()),
        only_b=int((~a_ & b_).sum()),
        neither=int((~a_ & ~b_).sum()),
    )


def fall(x: float, ref: float) -> float:
    """Share by which ``x`` is below ``ref`` (0.0 when ``ref`` is 0)."""
    return 0.0 if ref == 0 else float(ref - x) / float(ref)


def _no_worse(arm: int, ref: int, n_steps: int) -> bool:
    """The arm's count is at most 1% of all steps below the reference's (integer test)."""
    return (ref - arm) * 100 <= n_steps  # (ref - arm) / n_steps <= 0.01


def pick_arm(arms: Dict[str, Dict[str, int]], n_steps: int, ref: str = "M") -> str:
    """The cheapest arm (total evaluations; a tie keeps the arm listed first) for which (a) the
    steps that pass the seed test are at most 1% of ``n_steps`` fewer than with ``ref`` and (b)
    the steps that end within 1 degree of the truth are at most 1% of ``n_steps`` fewer. Counts
    are net (arm total against reference total). ``ref`` itself is returned when no other arm
    qualifies. ``arms[name]`` holds ``n_pass``, ``n_right`` and ``evals``."""
    r = arms[ref]
    ok = [
        name
        for name, c in arms.items()
        if name != ref
        and _no_worse(c["n_pass"], r["n_pass"], n_steps)
        and _no_worse(c["n_right"], r["n_right"], n_steps)
    ]
    if not ok:
        return ref
    return min(ok, key=lambda k: (arms[k]["evals"], list(arms).index(k)))


def mid_pass(arm: Dict[str, int], ref: Dict[str, int], n_steps: int) -> Dict[str, object]:
    """Pass rule of Part 1 on mid: (a) and (b) as in ``pick_arm``, and (c) the seed-step
    evaluations fall by at least 80% against ``ref``."""
    a = _no_worse(arm["n_pass"], ref["n_pass"], n_steps)
    b = _no_worse(arm["n_right"], ref["n_right"], n_steps)
    f = fall(arm["evals"], ref["evals"])
    c = (ref["evals"] - arm["evals"]) * 5 >= ref["evals"] * 4 and ref["evals"] > 0
    return dict(a=a, b=b, c=bool(c), evals_fall=f, passes=bool(a and b and c))


def s10_pass(
    evals: int, evals_s8: int, broken: int, fixed: int, med: float, med_s8: float
) -> Dict[str, object]:
    """Pass rule of S10 on one sample: (a) at least 20% fewer evaluations than S8 (both with the
    boundary pass); (b) broken minus fixed at most 5 voxels; (c) median error at most 0.005 degree
    above S8's."""
    a = (evals_s8 - evals) * 5 >= evals_s8 and evals_s8 > 0
    b = (broken - fixed) <= S10_NET
    c = med <= med_s8 + S10_MED + EPS
    return dict(
        fewer=fall(evals, evals_s8),
        a=bool(a),
        b=bool(b),
        c=bool(c),
        net_broken=int(broken - fixed),
        med_diff=float(med - med_s8),
        passes=bool(a and b and c),
    )
