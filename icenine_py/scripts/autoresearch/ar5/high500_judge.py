"""AR5 step 16: the pure numbers and rules of the verdict (no study imports).

The rules are ``high500_rules.JUDGE`` (frozen), repeated here as code. Arm ``J5`` against arm ``J``
on a sample (aligned draws; counts after the boundary pass):
(1) total with the pass: at least 4% fewer on ``held3``, 5% on ``held`` and ``mid``, 1% on ``big``;
    on ``grad`` at most 0.5% more.
(2) broken minus fixed at most 5.  (3) median error at most 0.001 degree above J's.
(4) at most 1 more wrong seed (a taken seed end wrong where J's was right; matched or unmatched),
    and grid searches within 1 of J's.
(5) the house rule against S8 with the pass: at least 20% fewer evaluations, broken minus fixed at
    most 5, median error at most 0.005 degree above S8's. If J fails (5) on the sample, J5 is judged
    there by (1) to (4). The step passes if J5 passes on all five samples.
"""

from typing import Any, Dict

EPS = 1e-9
FEWER_MIN = {"held3": 4.0, "held": 5.0, "mid": 5.0, "big": 1.0}  # percent fewer
MORE_MAX = {"grad": 0.5}  # percent more allowed
NET = 5
MED = 0.001
HOUSE_NET = 5
HOUSE_MED = 0.005
SAMPLES = ("held3", "held", "mid", "grad", "big")


def pct_fewer(new: float, old: float) -> float:
    return 0.0 if old == 0 else 100.0 * (old - new) / old


def house(evals: int, evals_s8: int, broken: int, fixed: int, med: float, med_s8: float) -> Dict:
    a = (evals_s8 - evals) * 5 >= evals_s8 and evals_s8 > 0
    b = (broken - fixed) <= HOUSE_NET
    c = med <= med_s8 + HOUSE_MED + EPS
    return dict(a=bool(a), b=bool(b), c=bool(c), passes=bool(a and b and c))


def rules(
    sample: str,
    total_j5: int,
    total_j: int,
    net_vs_j: int,
    med_j5: float,
    med_j: float,
    extra_wrong_seeds: int,
    grid_j5: int,
    grid_j: int,
    house_j5: Dict[str, bool],
    house_j: Dict[str, bool],
) -> Dict[str, Any]:
    cut = pct_fewer(total_j5, total_j)
    if sample in FEWER_MIN:
        r1 = cut >= FEWER_MIN[sample] - EPS
    else:
        r1 = cut >= -MORE_MAX[sample] - EPS
    r2 = net_vs_j <= NET
    r3 = med_j5 <= med_j + MED + EPS
    r4 = extra_wrong_seeds <= 1 and abs(grid_j5 - grid_j) <= 1
    r5 = bool(house_j5["passes"])
    j_fails_5 = not house_j["passes"]
    ok = bool(r1 and r2 and r3 and r4 and (r5 or j_fails_5))
    return dict(
        r1=bool(r1),
        r2=bool(r2),
        r3=bool(r3),
        r4=bool(r4),
        r5=r5,
        j_fails_5=bool(j_fails_5),
        total_cut_pct=cut,
        passes=ok,
    )


def step_verdict(per: Dict[str, bool]) -> bool:
    return all(per[s] for s in SAMPLES)
