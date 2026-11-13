"""AR5 step 12: the pure numbers and rules of the verdict (no study imports).

Used by ``stallhit_stats.py`` and ``stallhit_summary.py``, tested in ``tests/test_ar5_stallhit.py``.
The decision rules are the plan's, frozen before any run ("Decision rules"). Rules for sample ``s``
of arm ``H`` against arm ``G`` (aligned draws; counts after the boundary pass unless said):
(1) ``held2``: total at least 2% fewer and seed-step evaluations at least 5% fewer; ``mid``: total
    at least 5% fewer; ``big``: at least 2.5%; ``held``: fewer (more than 0%); ``grad``: at most 1%
    more.
(2) broken minus fixed at most 5. (3) median error at most 0.001 degree above G's.
(4) at most 1 more wrong seed (a taken seed end wrong where G's was right; matched or unmatched
    steps) and grid searches within 1 of G's.
(5) the house rule against S8 with the pass: at least 20% fewer evaluations, broken minus fixed at
    most 5, median error at most 0.005 degree above S8's. If G fails (5) on the sample, H is judged
    there by (1) to (4). The step passes if H passes on ``held2``, ``held`` and all disclosed
    samples; if it fails on ``held2``: no key.
"""

from typing import Any, Dict

EPS = 1e-9
NET = 5
MED = 0.001
HOUSE_CUT = 0.20
HOUSE_NET = 5
HOUSE_MED = 0.005
CUT_MIN = {"held2": 2.0, "mid": 5.0, "big": 2.5}  # percent fewer, total with the pass
SEED_MIN_HELD2 = 5.0
GRAD_MORE = 1.0  # percent more allowed
SAMPLES = ("held2", "held", "mid", "grad", "big")


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
    total_h: int,
    total_g: int,
    seed_h: int,
    seed_g: int,
    net_vs_g: int,
    med_h: float,
    med_g: float,
    extra_wrong_seeds: int,
    grid_h: int,
    grid_g: int,
    house_h: Dict[str, bool],
    house_g: Dict[str, bool],
) -> Dict[str, Any]:
    cut = pct_fewer(total_h, total_g)
    seed_cut = pct_fewer(seed_h, seed_g)
    if sample == "held2":
        r1 = cut >= CUT_MIN["held2"] - EPS and seed_cut >= SEED_MIN_HELD2 - EPS
    elif sample == "held":
        r1 = total_h < total_g
    elif sample == "grad":
        r1 = cut >= -GRAD_MORE - EPS
    else:
        r1 = cut >= CUT_MIN[sample] - EPS
    r2 = net_vs_g <= NET
    r3 = med_h <= med_g + MED + EPS
    r4 = extra_wrong_seeds <= 1 and abs(grid_h - grid_g) <= 1
    r5 = bool(house_h["passes"])
    g_fails_5 = not house_g["passes"]
    ok = bool(r1 and r2 and r3 and r4 and (r5 or g_fails_5))
    return dict(
        r1=bool(r1),
        r2=bool(r2),
        r3=bool(r3),
        r4=bool(r4),
        r5=r5,
        g_fails_5=bool(g_fails_5),
        total_cut_pct=cut,
        seed_cut_pct=seed_cut,
        passes=ok,
    )


def step_verdict(per: Dict[str, bool]) -> bool:
    """The step passes if H passes on held2, held and every disclosed sample."""
    return all(per[s] for s in SAMPLES)
