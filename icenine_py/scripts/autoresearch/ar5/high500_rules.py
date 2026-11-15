"""AR5 step 16 (task ar5p): rule ``PH(500)``, the high stop alone at a 500-evaluation stall in the
seed step's MC refinement (pure; numpy, ``band_rules`` and ``stallhit_rules`` only).

``PH(k)`` is step 12's rule (frozen there), here with ``k = 500`` in place of ``PH(3000)``. On one
MC refinement of the seed step (the MC stage then the variance stage, one run) a **stall** is ``k``
strict-cost evaluations with no strict gain of the best cost. At a stall there is one **stall
check**: the hit ratio of the running best (the first lowest-cost orientation) is evaluated once.
At 0.95 or more the refinement stops and keeps that orientation (the check is the final evaluation,
counted). Under 0.95 it goes on (the check counted as one evaluation) and checks again ``k``
evaluations after the next strict gain. There is no low side of any kind and no other ``k``. The
refinement has its own generator, so a stop moves no other random draw.

The meter, the decision and the walk are ``band_rules`` with ``low=None`` (step 15's code with no
low side); this file fixes the constants, the gate numbers of the plan and the decision rules.
"""

from typing import Any, Dict, List

import band_rules as BR
import stallhit_rules as SR

TAU = BR.TAU  # 0.95: the seed step's stop value
K = 500  # the stall size of this step: the only form
PH_K = SR.K  # 3,000: arm J's own stall size, replaced by K in arm J5 (never present beside it)
SEED_TEST = SR.SEED_TEST  # 0.80: unchanged
RECORD_K = 500  # the stall size of step 15's recorded events that this step walks

RULE = (
    "Rule PH(500): in the seed step, the MC refinement of a candidate (MC then variance stage) has "
    "a stall check after 500 strict-cost evaluations with no strict gain of its best cost: the hit "
    "ratio of the running best (the first lowest-cost orientation) is evaluated once; at 0.95 or "
    "more the refinement stops and the check is its final evaluation (counted); else it goes on "
    "(the check counted as one evaluation) and checks again 500 evaluations after the next gain. "
    "No low side. It replaces PH(3000). CMA fits, the five candidates, the 0.95 stop, the seed "
    "test 0.80, rules G1 and G2b and the grid search are unchanged. Arm J5 = arm J (step 14) with "
    "PH(500)."
)

# The stored counts of arm J (plan stop rule): total with the pass, MC evaluations of the seed step
EXPECT: Dict[str, Dict[str, int]] = {
    "held": dict(total=24_457_773, mc=11_381_609),
    "mid": dict(total=11_295_158, mc=3_795_130),
    "grad": dict(total=9_175_867, mc=1_078_188),
    "big": dict(total=3_635_256, mc=424_450),
}

# Prediction 1 (plan table): moved steps, path changes, net saving of the arithmetic walk of every
# step (every check counted; growth after a changed voxel is not replayed)
PRED1: Dict[str, Dict[str, int]] = {
    "held": dict(moved=8, path_change=0, saving=1_584_152),
    "mid": dict(moved=6, path_change=1, saving=740_633),
    "grad": dict(moved=1, path_change=0, saving=134_974),
    "big": dict(moved=1, path_change=0, saving=88_144),
}


def decide(hit: float) -> str:
    """``"high"`` (stop: the best is at 0.95 or more) or ``"go"``; never a low stop."""
    return BR.decide(hit, TAU, None)


def refine_cuts(refs: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
    """``PH(500)`` on every recorded MC refinement of one seed step (``refs``: the
    ``refinements`` of a recording row of step 15: ``pos``, ``evals`` as J ran it, ``cont_j``
    J's own checks that went on, ``final_cost``, ``events`` by stall size). Returns per position a
    ``band_rules.cut_entry`` plus ``base`` (J's evaluations) and ``w`` (the walk)."""
    out: Dict[int, Dict[str, Any]] = {}
    for r in refs:
        w = BR.walk_events(
            r["events"].get(str(RECORD_K), []),
            int(r["evals"]) - int(r["cont_j"]),
            TAU,
            None,
        )
        c = BR.cut_entry(w, float(r["final_cost"]))
        c["base"] = int(r["evals"])
        c["w"] = w
        out[int(r["pos"])] = c
    return out


def path_change_refs(refs: List[Dict[str, Any]]) -> List[int]:
    """Positions of the refinements whose stop moves the end and whose own end (``end_hit``) was
    under 0.95: there the stopped end (at 0.95 or more) and J's end are at different sides of the
    stop value, so the seed step's path can change. A moved end where J's end is also at 0.95 or
    more keeps the step's taken candidate (only the orientation differs a little)."""
    cuts = refine_cuts(refs)
    return [
        int(r["pos"])
        for r in refs
        if cuts[int(r["pos"])]["kind"] == "moved" and float(r["end_hit"]) < TAU
    ]


# ---- the decision rules of the plan, fixed before any run
GATE: Dict[str, Any] = dict(
    samples=("held", "mid", "grad", "big"),
    g0="prediction 1: the stored counts, every recorded step equal to the stored one, and the "
    "table's moved steps, path changes and net savings (to the evaluation)",
    g1_moved_max=dict(held=10, mid=10, grad=2, big=2),
    g1_path_change_max=2,
    g2="on the live replayed steps: no taken end wrong where J's is right, no extra failed seed "
    "test, no extra grid search",
    g3_saving_min=0.055,  # net saving, every check counted, share of J's run with the pass
    g3_samples=("held", "mid"),
    fails="stop and report: no tuning, no other k, no whole layer",
)
JUDGE: Dict[str, Any] = dict(
    r1_total_fewer=dict(held3=0.04, held=0.05, mid=0.05, big=0.01),
    r1_grad_more_max=0.005,
    r2_broken_minus_fixed_max=5,
    r3_median_error_higher_max_deg=0.001,
    r4_wrong_seeds_more_max=1,
    r4_grid_within=1,
    r5_vs_s8=dict(fewer=0.20, broken_minus_fixed_max=5, median_error_higher_max_deg=0.005),
    pass_samples=("held3", "held", "mid", "grad", "big"),
)
