"""AR5 step 24 (task ar5x): rule ``BF`` (ball first). Pure: numpy only; no image, no reconstructor,
no truth. The constants, the decision rules and predictions of the plan (frozen), and the pure
arithmetic of the two-pass seed step.

Arm ``J5DB`` is arm ``J5D`` (step 22: ``J5`` with rule ``DEF``) with rule ``BF``; ``J5SDB`` is
``J5DB`` with step 21's screen pass and wide repair (imported unchanged).

**Rule ``BF``.** The seed step runs in two passes.

1. Ball pass. For the best 5 candidates in turn: the CMA fit as now; if its end is under 0.95, the
   MC refinement's first stage as now, its kept state is stored, and the ball stage (step 23's rule
   ``CB``: ``ball_rules.ball_stage``) runs on a copy of the refinement's generator. The pass stops
   at the first end at 0.95 or more, and the step takes it.
2. Fallback, only if no end of the ball pass reaches 0.95. ``J5D``'s seed step runs as it always
   does (``PH(500)``, ``DEF``, candidate order, pauses and resumes), with two changes: each CMA fit
   is taken from the ball pass, not run again, and each MC refinement starts its stock second stage
   from its kept state, not from the candidate. Its result is ``J5D``'s result for that step,
   evaluation by evaluation.

Script-layer note (phase A): a kept state is reproduced by running the refinement again from its
candidate (the refinement has its own generator, ``default_rng(0)``, so the first stage is the
same); its first stage is CHECKED equal to the kept one and is not counted twice.
"""

from typing import Any, Dict, List, Optional

import ball_rules as BR

TAU = BR.TAU
CMA_EVALS = BR.CMA_EVALS
FAIL_HALVE = BR.FAIL_HALVE
FLOOR_DEG = BR.FLOOR_DEG
CAP_PROPOSALS = BR.CAP_PROPOSALS
ARM_DB = "J5DB"
ARM_SDB = "J5SDB"
NEW_SAMPLE = "held7"
GENERATOR_SEED = 7  # held7 = held with generator seed 7: a new sample, not an extra random seed
DISCLOSED = BR.DISCLOSED
LIVE_SAMPLES = ("mid", "held", "grad", "big")  # the disclosed samples run live
HELD_TYPES = ("held5", "held4", "held3", "held", "held6", "held7")
GATE_LIVE = ("mid", "held")
LIVE_REST = ("grad", "big")

RULE = (
    "Rule BF (ball first): the seed step runs in two passes. (1) Ball pass: for the best 5 "
    "candidates in turn, the CMA fit as now; if its end is under 0.95, the MC refinement's first "
    "stage as now, then its kept state is stored (best orientation, cost, evaluation count, "
    "generator state), then the ball stage (step 23's rule CB: start at the first stage's best, "
    "radius the box, halve after 10 proposals in a row with no gain, end under 0.01 degree or at "
    "300 proposals) runs on a copy of the refinement's generator. The pass stops at the first end "
    "at 0.95 or more, and the step takes it. (2) Fallback, only if no end of the ball pass reaches "
    "0.95: J5D's seed step runs as it always does (PH(500), DEF, candidate order, pauses and "
    "resumes), with two changes: each CMA fit is taken from the ball pass, not run again, and each "
    "MC refinement starts its stock second stage from its kept state, not from the candidate. Its "
    "result is J5D's result for that step, evaluation by evaluation. If no end reaches 0.95, the "
    "lowest-cost end of the fallback is taken, as in J5D. The ball stage's constants (10, 0.01 "
    "degree, 300, the box as start radius) are step 23's, not changed and not tuned. No new "
    "threshold, no hit ratio inside the stages, no neighbour, no truth, no crystal symmetry. "
    "Unchanged: costing, CMA fits, 0.80, 0.95, G1, G2b, the grid search with its own final "
    "refinement, BFS, the screen pass and the wide repair. Every evaluation of both passes counts, "
    "including the ball stage in steps that fall back. Arm J5DB = J5D with BF; arm J5SDB = J5DB "
    "with step 21's screen pass and wide repair."
)

PREDICTIONS: Dict[str, str] = {
    "1": "Replay on mid: the 1,453 steps that end in the ball pass equal step 23's stored rows "
    "(stop, taken end, evaluations); the 173 fallback steps equal J5D's stored steps (fits, stop, "
    "taken end); the index stage is 2,099,764 evaluations exactly; 0 lost ends; 0 checks and 0 "
    "pauses in the ball pass.",
    "2": "The tracked sizing script reproduces the first table exactly on all eight samples.",
    "3": "Live J5SDB against J5SD: 4.5% to 10% fewer evaluations on the held types, 5% to 9.5% on "
    "mid, 2% to 5% on grad, 1.5% to 4% on big; broken minus fixed -3 to +3; grid searches within "
    "2 of J5SD's; seed-test failures within 10 of J5SD's; median error within 0.0005 degree.",
    "4": "held7: 4.5% to 10% fewer evaluations than J5SD; broken minus fixed -3 to +3; median "
    "error within 0.0005 degree; the house rule against S8 passes.",
}

GATE: Dict[str, Any] = dict(
    g0="prediction 1, exactly (mid)",
    g1="prediction 2, exactly, with lost 0.95 ends at most 3% of the steps on every sample (by "
    "construction 0; any lost end is a bug and stops the step) and the index stage at least 10% "
    "fewer evaluations than J5D's on the held types and mid",
    g1_fewer_min=0.10,
    g1_fewer_samples=("held5", "held4", "held3", "held", "held6", "mid"),
    g1_lost_max=0.03,
    g2="live J5SDB on mid and held: at least 4% fewer evaluations than J5SD, broken minus fixed at "
    "most 3, median error at most 0.001 degree higher, grid searches at most 2 more than J5SD's, "
    "seed-test failures at most 10 more",
    g2_samples=GATE_LIVE,
    g2_fewer_min=0.04,
    g2_broken_minus_fixed_max=3,
    g2_median_error_higher_max_deg=0.001,
    g2_grid_more_max=2,
    g2_seed_test_more_max=10,
    fails="stop and report: no tuning (no change of the ball's constants or of the fallback), no "
    "held7",
)

JUDGE: Dict[str, Any] = dict(
    r1_fewer_min=dict(
        held5=0.04,
        held4=0.04,
        held3=0.04,
        held=0.04,
        held6=0.04,
        held7=0.04,
        mid=0.04,
        grad=0.015,
        big=0.01,
    ),
    r2_broken_minus_fixed_max=3,
    r3_median_error_higher_max_deg=0.001,
    r4_vs_s8=dict(fewer=0.20, broken_minus_fixed_max=5, median_error_higher_max_deg=0.005),
    passes="J5SDB passes on a sample if (1) to (4) hold against J5SD (and S8 for 4)",
    step_passes="J5SDB passes on held7 and on the four disclosed samples run live (mid, held, "
    "grad, big)",
    fails_on_held7="no key value; study merge",
)

CAPS_MIN: Dict[str, int] = dict(S8=480, rest=300)


# ---- pure arithmetic of the two-pass step --------------------------------------------------------
def fewer(new: float, old: float) -> float:
    """The share fewer (positive: ``new`` is smaller than ``old``)."""
    return BR.fewer(new, old)


def judge_vs_j5sd(
    sample: str, fewer_: float, broken_minus_fixed: int, median_higher: float
) -> Dict[str, Any]:
    """Rules (1) to (3) of the plan against ``J5SD``."""
    j = JUDGE
    r1 = bool(fewer_ >= j["r1_fewer_min"][sample])
    r2 = bool(broken_minus_fixed <= j["r2_broken_minus_fixed_max"])
    r3 = bool(median_higher <= j["r3_median_error_higher_max_deg"])
    return dict(r1=r1, r2=r2, r3=r3, ok=bool(r1 and r2 and r3))


def judge_vs_s8(fewer_: float, broken_minus_fixed: int, median_higher: float) -> Dict[str, Any]:
    """Rule (4), the house rule against ``S8``."""
    return BR.judge_vs_s8(fewer_, broken_minus_fixed, median_higher)


def fallback_extra(fallback_fits: List[List[Any]], stage1_by_j: Dict[int, int]) -> Dict[str, int]:
    """The evaluations the fallback pays beyond the ball pass: every MC refinement it ran costs its
    own count minus its first stage (the first stage is already counted in the ball pass; the CMA
    fits are taken from the ball pass and are not in the fallback's count at all).

    ``fallback_fits``: the fallback step's ``fits`` rows ``[pos, j, kind, hit, cost, evals, R]``;
    ``stage1_by_j``: the first-stage evaluations of the ball pass's refinements by candidate."""
    extra = 0
    n_ref = 0
    for f in fallback_fits:
        if f[2] != "M":
            continue
        j = int(f[1])
        assert j in stage1_by_j, f"the fallback ran a refinement the ball pass never ran (j={j})"
        extra += int(f[5]) - int(stage1_by_j[j])
        n_ref += 1
    return dict(extra=int(extra), refinements=n_ref)


def step_evals(ball_evals: int, ball_stop: int, extra: Optional[int]) -> int:
    """The step's evaluations: the ball pass's (costing, CMA fits, first stages, ball stages), plus
    the fallback's extra when the ball pass found no end at 0.95."""
    if ball_stop >= 0:
        assert not extra, "a step that ends in the ball pass has no fallback"
        return int(ball_evals)
    return int(ball_evals) + int(extra or 0)


def lost_end(j5d_stop: int, bf_stop: int) -> bool:
    """A seed step that loses its 0.95 end: J5D had one (stop >= 0), BF has none."""
    return bool(j5d_stop >= 0 and bf_stop < 0)


def feature_row(
    index: int,
    ball_refs: List[Dict[str, Any]],
    cma_fits: List[Dict[str, Any]],
    found: bool,
    miso_deg: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """One feature row of a fallback step, no extra evaluation: the five ball ends' and CMA ends'
    own hit ratios and costs, their spread in angle, and whether the fallback found 0.95."""
    return dict(
        index=int(index),
        ball_hit=[float(r["hit"]) for r in ball_refs],
        ball_cost=[float(r["cost"]) for r in ball_refs],
        cma_hit=[float(c["hit"]) for c in cma_fits],
        cma_cost=[float(c["cost"]) for c in cma_fits],
        spread_deg=None if miso_deg is None else [float(x) for x in miso_deg],
        fallback_found=bool(found),
    )
