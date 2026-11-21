"""AR5 step 14: the pure numbers and rules of the verdict (no study imports).

Used by ``polish_stats.py`` and ``polish_summary.py``, tested in ``tests/test_ar5_polish.py``. The
rules are the plan's ("Decision rules", fixed before any run). Arm ``J`` (``H`` with rule ``G2``)
against arm ``H`` on a sample (aligned draws; counts after the boundary pass):
(1) total with the pass: ``held2`` at least 6% fewer, ``held`` at least 5% fewer, ``mid`` and
    ``grad`` at most 0.5% more, ``big`` identical.
(2) broken minus fixed at most 5. (3) median error at most 0.001 degree above H's.
(4) at most 1 more wrong seed (a taken seed end wrong where H's was right; matched or unmatched).
(5) the house rule against S8 with the pass: at least 20% fewer evaluations, broken minus fixed at
    most 5, median error at most 0.005 degree above S8's. If H fails (5) on the sample, J is
    judged there by (1) to (4). The step passes if J passes on all five samples.
"""

from typing import Any, Dict

EPS = 1e-9
NET = 5
MED = 0.001
HOUSE_CUT = 0.20
HOUSE_NET = 5
HOUSE_MED = 0.005
SAMPLES = ("held2", "held", "mid", "grad", "big")
MIN_FEWER = {"held2": 6.0, "held": 5.0}  # percent fewer, total with the pass
MAX_MORE = {"mid": 0.5, "grad": 0.5}  # percent more allowed


def pct_fewer(new: float, old: float) -> float:
    """Percent fewer: 100 * (old - new) / old (0 when old is 0)."""
    return 0.0 if old == 0 else 100.0 * (old - new) / old


def house(evals: int, evals_s8: int, broken: int, fixed: int, med: float, med_s8: float) -> Dict:
    a = (evals_s8 - evals) * 5 >= evals_s8 and evals_s8 > 0
    b = (broken - fixed) <= HOUSE_NET
    c = med <= med_s8 + HOUSE_MED + EPS
    return dict(a=bool(a), b=bool(b), c=bool(c), passes=bool(a and b and c))


def rules(
    sample: str,
    total_j: int,
    total_h: int,
    net_vs_h: int,
    med_j: float,
    med_h: float,
    extra_wrong_seeds: int,
    house_j: Dict[str, bool],
    house_h: Dict[str, bool],
) -> Dict[str, Any]:
    cut = pct_fewer(total_j, total_h)
    if sample in MIN_FEWER:
        r1 = cut >= MIN_FEWER[sample] - EPS
    elif sample in MAX_MORE:
        r1 = cut >= -MAX_MORE[sample] - EPS
    else:
        r1 = total_j == total_h
    r2 = net_vs_h <= NET
    r3 = med_j <= med_h + MED + EPS
    r4 = extra_wrong_seeds <= 1
    r5 = bool(house_j["passes"])
    h_fails_5 = not house_h["passes"]
    ok = bool(r1 and r2 and r3 and r4 and (r5 or h_fails_5))
    return dict(
        r1=bool(r1),
        r2=bool(r2),
        r3=bool(r3),
        r4=bool(r4),
        r5=r5,
        h_fails_5=bool(h_fails_5),
        total_cut_pct=cut,
        passes=ok,
    )


def step_verdict(per: Dict[str, bool]) -> bool:
    """The step passes if J passes on every sample."""
    return all(per[s] for s in SAMPLES)
