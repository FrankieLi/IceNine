"""AR5 step 23 (task ar5w): rule ``CB`` (converging ball) for the second MC stage of the seed step's
MC refinement. Pure: numpy only; no image, no reconstructor, no truth.

Arm ``J5DC`` is arm ``J5D`` (step 22: ``J5`` with rule ``DEF``) with rule ``CB``; ``J5SDC`` is
``J5DC`` with step 21's screen pass and wide repair (imported unchanged).

**Terms.** *Stage 1*: the library's ``FindOptimal`` stage of ``refine_from_candidates`` (blocks of
31 samples, at most 208 evaluations). *Stage 2*: the library's ``variance_minimizing_optimize``
(the variance stage). *Box*: ``max(diameter / 3, 0.2 deg) / 2^MinLocalResolution`` (0.329 degree
in these runs), the ``search_box_side`` the library passes to stage 2. *Ball of radius r*: a
rotation about a uniformly random axis by the angle r times the cube root of a uniform number, so
that the points fill a ball of true angle r evenly. *Gain*: a strict drop of the running best
cost.

**Rule ``CB``.** Stage 1 runs as now. Stage 2 is replaced. It evaluates the stage-1 best (one
counted evaluation, as the stock stage does), then makes proposals: each is the current best turned
(left product, as the stock stage composes) by a rotation drawn in a ball of radius r, r starting
at the box. A proposal with a strictly lower cost becomes the best and the new centre (r is not
changed by a gain). After 10 proposals in a row with no gain r halves and the count restarts. The
stage ends when r is under 0.01 degree (``floor``; tested first), or after 300 proposals
(``cap``). Then the library's final evaluation, as now. The draws come from the refinement's own
generator, continuing from stage 1: per proposal exactly four numbers, three normals (the axis,
normalised) then one uniform (the radius). No hit ratio, no neighbour, no truth, no crystal
symmetry in the stage.

The constants 10, 0.01 degree and 300 are set before any run and are not tuned.
"""

import math
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

# ---- constants (frozen) -------------------------------------------------------------------------
FAIL_HALVE = 10  # proposals in a row with no gain before r halves
FLOOR_DEG = 0.01  # the stage ends when r falls under this
CAP_PROPOSALS = 300  # ... or after this many proposals
ARM_DC = "J5DC"
ARM_SDC = "J5SDC"
NEW_SAMPLE = "held7"
GENERATOR_SEED = 7  # held7 = held with generator seed 7: a new sample, not an extra random seed
TAU = 0.95
CMA_EVALS = 201  # one CMA fit of the seed step
DISCLOSED = ("held5", "held4", "held3", "held", "mid", "grad", "big", "held6")
LIVE_SAMPLES = ("mid", "held", "held6", "grad", "big")  # the disclosed samples run live
HELD_TYPES = ("held5", "held4", "held3", "held", "held6", "held7")
GATE_LIVE = ("mid", "held")
LIVE_REST = ("held6", "grad", "big")  # J5SDC live after the gate

RULE = (
    "Rule CB (converging ball): in the seed step's MC refinement (library refine_from_candidates "
    "on one candidate) stage 1 (FindOptimal) runs as now. Stage 2 (the variance stage) is "
    "replaced: "
    "it evaluates the stage-1 best (one counted evaluation), then proposes the current best turned "
    "by a rotation drawn in a ball of radius r (axis uniform on the sphere, angle r times the cube "
    "root of a uniform number), r starting at the box (0.329 degree); a strictly lower cost "
    "becomes the best and the new centre; after 10 proposals in a row with no gain r halves; the "
    "stage ends when r is under 0.01 degree (floor, tested first) or after 300 proposals (cap); "
    "then the final "
    "evaluation, as now. The refinement's own generator continues from stage 1 (three normals, "
    "one uniform per proposal). No hit ratio, no neighbour, no truth, no crystal symmetry in the "
    "stage. Unchanged: 20 costed and 5 refined candidates, CMA fits, 0.80, 0.95, PH(500) and DEF "
    "(inert under CB: zero checks, zero pauses or stop), G1, G2b, the grid search with its own "
    "final refinement, BFS, the screen pass and the wide repair. Arm J5DC = J5D with CB; arm "
    "J5SDC = J5DC with step 21's screen pass and wide repair."
)

PREDICTIONS: Dict[str, str] = {
    "1": "With CB the first stage of every refinement equals J5D's, evaluation by evaluation "
    "(mid); with CB off the replay reproduces J5D's seed steps on mid (ends, evaluations).",
    "2": "Replay, eight disclosed samples: the seed step's index stage 40% to 60% fewer "
    "evaluations than J5D's on the held types and mid, 20% to 45% on grad and big; 250 to 450 "
    "evaluations per refinement (J5D: 793 to 1,341); the second stage ends at the floor in at "
    "least 90% of refinements; steps with an end at 0.95 within 3% of J5D's count.",
    "3": "Live J5SDC against J5SD: 15% to 28% fewer evaluations on the held types, 10% to 20% on "
    "mid, 3% to 8% on grad, 1% to 4% on big; broken minus fixed -3 to +4; grid searches within 2 "
    "of J5D's; median error within 0.0005 degree.",
    "4": "held7: 15% to 28% fewer evaluations than J5SD; broken minus fixed -3 to +4; median error "
    "within 0.0005 degree; the house rule against S8 passes.",
}

GATE: Dict[str, Any] = dict(
    g0="prediction 1, exactly (mid)",
    g1="replay on all eight: the index stage at least 30% fewer evaluations than J5D's on the "
    "held types and mid; steps that lose their 0.95 end (J5D had one, CB has none) at most 3% of "
    "the replayed steps",
    g1_fewer_min=0.30,
    g1_fewer_samples=("held5", "held4", "held3", "held", "held6", "mid"),
    g1_lost_max=0.03,
    g2="live J5SDC on mid and held: at least 10% fewer evaluations than J5SD, broken minus fixed "
    "at most 4, median error at most 0.001 degree higher, grid searches at most 2 more than "
    "J5D's",
    g2_samples=GATE_LIVE,
    g2_fewer_min=0.10,
    g2_broken_minus_fixed_max=4,
    g2_median_error_higher_max_deg=0.001,
    g2_grid_more_max=2,
    fails="stop and report: no tuning (no change of 10, 0.01 degree, 300 or the start radius), "
    "no held7",
)

JUDGE: Dict[str, Any] = dict(
    r1_fewer_min=dict(
        held5=0.10,
        held4=0.10,
        held3=0.10,
        held=0.10,
        held6=0.10,
        held7=0.10,
        mid=0.10,
        grad=0.02,
        big=0.01,
    ),
    r2_broken_minus_fixed_max=4,
    r3_median_error_higher_max_deg=0.001,
    r4_vs_s8=dict(fewer=0.20, broken_minus_fixed_max=5, median_error_higher_max_deg=0.005),
    passes="J5SDC passes on a sample if (1) to (4) hold against J5SD (and S8 for 4)",
    step_passes="J5SDC passes on held7 and on mid, held, held6, grad and big run live",
    fails_on_held7="no key value; study merge",
)

CAPS_MIN: Dict[str, int] = dict(S8=480, rest=300)


# ---- the ball and the stage (pure) -------------------------------------------------------------
def rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    """The rotation matrix about the unit ``axis`` by ``angle`` (radians)."""
    x, y, z = float(axis[0]), float(axis[1]), float(axis[2])
    k = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + math.sin(angle) * k + (1.0 - math.cos(angle)) * (k @ k)


def draw_ball(rng: np.random.Generator, r: float) -> Tuple[np.ndarray, float]:
    """One rotation in a ball of radius ``r`` (radians): axis uniform on the sphere (three
    normals, normalised), angle ``r`` times the cube root of one uniform. Returns the matrix and
    the angle. Exactly four numbers are drawn."""
    v = rng.normal(size=3)
    u = float(rng.uniform())
    n = float(np.linalg.norm(v))
    angle = r * u ** (1.0 / 3.0)
    axis = v / n if n > 0.0 else np.array([1.0, 0.0, 0.0])
    return rodrigues(axis, angle), angle


def ball_stage(
    evaluate: Callable[[np.ndarray], Tuple[float, Any]],
    R0: np.ndarray,
    r0: float,
    rng: np.random.Generator,
    fail_halve: int = FAIL_HALVE,
    floor: float = math.radians(FLOOR_DEG),
    cap: int = CAP_PROPOSALS,
) -> Dict[str, Any]:
    """The converging ball (stage 2 of rule CB).

    ``evaluate(R)`` returns ``(cost, info)`` and is called once on ``R0`` (the start) and once per
    proposal: ``1 + proposals`` calls in all. Returns ``R`` and ``cost`` and ``info`` of the best,
    ``proposals``, ``gains``, ``halvings``, ``end`` (``"floor"`` or ``"cap"``), ``r_final``
    (radians), ``evals`` (the calls made) and ``start_cost``."""
    best_R = np.asarray(R0, dtype=np.float64).copy()
    best_cost, best_info = evaluate(best_R)
    start_cost = float(best_cost)
    r = float(r0)
    fails = gains = halvings = proposals = 0
    end: Optional[str] = None
    if r < floor:
        end = "floor"
    while end is None:
        dR, _ = draw_ball(rng, r)
        trial = dR @ best_R
        c, info = evaluate(trial)
        proposals += 1
        if c < best_cost:
            best_R, best_cost, best_info = trial, c, info
            gains += 1
            fails = 0
        else:
            fails += 1
            if fails >= fail_halve:
                r *= 0.5
                halvings += 1
                fails = 0
        if r < floor:
            end = "floor"
        elif proposals >= cap:
            end = "cap"
    return dict(
        R=best_R,
        cost=float(best_cost),
        info=best_info,
        proposals=proposals,
        gains=gains,
        halvings=halvings,
        end=end,
        r_final=r,
        evals=1 + proposals,
        start_cost=start_cost,
    )


def floor_proposals_min(r0: float, floor: float = math.radians(FLOOR_DEG)) -> int:
    """The fewest proposals that can end a stage at the floor (every proposal fails)."""
    n, r = 0, float(r0)
    while r >= floor:
        r *= 0.5
        n += FAIL_HALVE
    return n


def fewer(new: float, old: float) -> float:
    """The share fewer (positive: ``new`` is smaller than ``old``)."""
    return 1.0 - float(new) / float(old)


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
    h = JUDGE["r4_vs_s8"]
    a = bool(fewer_ >= h["fewer"])
    b = bool(broken_minus_fixed <= h["broken_minus_fixed_max"])
    c = bool(median_higher <= h["median_error_higher_max_deg"])
    return dict(fewer=a, broken=b, median=c, ok=bool(a and b and c))


def lost_end(j5d_stop: int, cb_stop: int) -> bool:
    """A seed step that loses its 0.95 end: J5D had one (stop >= 0), CB has none."""
    return bool(j5d_stop >= 0 and cb_stop < 0)


def split_refinement(ev: int, stage1: int, proposals: int) -> Dict[str, int]:
    """One refinement's evaluations: stage 1, stage 2 (the start evaluation and the proposals) and
    the rest (the final evaluation, 1)."""
    s2 = 1 + int(proposals)
    return dict(stage1=int(stage1), stage2=s2, final=int(ev) - int(stage1) - s2)


def first_stage_equal(
    stage1: int, cb_trace: List[List[float]], j5d_trace: List[List[float]]
) -> bool:
    """The gains of stage 1 (evaluation numbers up to ``stage1``) equal, evaluation by evaluation
    and cost by cost: rule CB changes nothing before stage 2."""
    a = [[int(n), float(c)] for n, c in cb_trace if int(n) <= stage1]
    b = [[int(n), float(c)] for n, c in j5d_trace if int(n) <= stage1]
    return bool(a == b)
