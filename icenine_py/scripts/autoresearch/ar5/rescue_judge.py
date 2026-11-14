"""AR5 step 9, phase A step 6: the pure numbers and rules of the verdict (no study imports).

Used by ``rescue_stats.py`` and ``rescue_summary.py`` and tested in ``tests/test_ar5_rescue.py``.
The decision rules are the plan's, frozen before any run (plan, "Decision rules"). Where a rule is
unclear the stricter reading is applied and said so.

Rules for sample ``s`` against arm ``A`` (all counts after the boundary pass unless said):
(1) ``held`` and ``mid``: fewer grid fallbacks than A and a total (with the pass) at least 10%
    fewer; ``grad``: a total fewer; ``big``: identical (evaluations and map).
(2) broken minus fixed against A at most 5 (after the pass).
(3) median error at most 0.001 degree above A's (after the pass).
(4) at most 1 wrong rescued seed.
(5) the house rule against S8 with the pass: at least 20% fewer evaluations, broken minus fixed at
    most 5, median error at most 0.005 degree above S8's. If A fails (5) on the sample, G is
    judged there by (1) to (4). Stricter reading: "A fails (5)" is taken from the same three tests.
The step passes if G passes on ``held`` and on all disclosed samples (``mid``, ``grad``, ``big``).
"""

from typing import Any, Dict

EPS = 1e-9
NET_A = 5
MED_A = 0.001
MIN_CUT = 0.10
HOUSE_CUT = 0.20
HOUSE_NET = 5
HOUSE_MED = 0.005


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
    total_g: int,
    total_a: int,
    fallbacks_g: int,
    fallbacks_a: int,
    net_vs_a: int,
    med_g: float,
    med_a: float,
    wrong_rescued: int,
    identical_to_a: bool,
    house_g: Dict[str, bool],
    house_a: Dict[str, bool],
) -> Dict[str, Any]:
    cut = pct_fewer(total_g, total_a)
    if sample in ("held", "mid"):
        r1 = fallbacks_g < fallbacks_a and cut >= 100 * MIN_CUT - EPS
    elif sample == "grad":
        r1 = total_g < total_a
    else:
        r1 = bool(identical_to_a)
    r2 = net_vs_a <= NET_A
    r3 = med_g <= med_a + MED_A + EPS
    r4 = wrong_rescued <= 1
    r5 = bool(house_g["passes"])
    a_fails_5 = not house_a["passes"]
    ok = bool(r1 and r2 and r3 and r4 and (r5 or a_fails_5))
    return dict(
        r1=bool(r1),
        r2=bool(r2),
        r3=bool(r3),
        r4=bool(r4),
        r5=r5,
        a_fails_5=bool(a_fails_5),
        total_cut_pct=cut,
        passes=ok,
    )


def step_verdict(per: Dict[str, bool]) -> bool:
    """The step passes if G passes on held and on every disclosed sample."""
    return all(per[s] for s in ("held", "mid", "grad", "big"))
