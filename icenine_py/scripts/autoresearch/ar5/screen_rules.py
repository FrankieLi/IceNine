"""AR5 step 21 (task ar5u): score each pass candidate once, and let step 20's repair do the fitting
(pure: numpy and the pure functions of step 20's and step 19's frozen modules; no image, no
reconstructor of its own, no truth).

Arm ``J5S`` = arm ``J5`` with screen tries in the boundary pass and step 20's repair on the wide
flag. Everything that decides a voxel is here or in the two modules it imports and does not copy:

* ``cmarepair_rules`` (step 20): the repair (``repair``: CMA 120, first step 0.05 degree, own end
  and every candidate of the snapshot, per-fit generators ``fit_rng``, keep rule ``hit``) and the
  step 19 functions it imports (``contested``, ``after_pass_state``, ``other_end_hit``,
  ``repair_candidates``).
* **Screen try**: one boundary-pass try (the CMA fit of the pass, first step 0.05 degree) with a
  budget of ONE evaluation: the candidate's start is scored and is the try's end. The library's
  final evaluation is counted too: 2 evaluations a try. ``Reconstructor.local_optimization``
  refuses a budget under 2 in its override (``cma_max_evals``), so the script layer sets the
  budget on the optimizer that the library makes (``CMAOptimizer`` itself takes 1): the library
  code that runs is the library's. Phase B's key must take the budget 1 in the library itself.
* **Wide flag**: a voxel is flagged if it is contested (step 19's flag, with the screen ends in
  its margin), or if its ratio is under 0.95 and one of its pass candidates (try ends) lies more
  than 1 degree from its kept end.

The pass's draws are step 8's aligned generators (``make_rng(v, "pass", k)``, k the number of tries
made at the voxel so far), as ``J5``'s pass. A one-evaluation try draws one number it does not use.
"""

from typing import Any, Dict, List, Sequence, Tuple

import cmarepair_rules as CR
import numpy as np

# ---- this step's constants -----------------------------------------------------------------------
ARM = "J5S"
SCREEN_MAX_EVALS = 1  # the optimizer's budget of a screen try (the start only)
EVALS_PER_SCREEN = SCREEN_MAX_EVALS + 1  # ... plus the library's final evaluation
FULL_MAX_EVALS = 60  # J5's own pass tries (the gate's re-run)
PASS_SIGMA0_DEG = 0.05  # first step of a try (J5's pass)
RATIO_MAX = CR.RATIO_MAX  # 0.95: the wide flag's ratio limit (the contested flag's, too)
WIDE_DEG = CR.FAR_DEG  # 1 degree: a candidate this far from the kept end makes a voxel wide
WRONG_DEG = CR.WRONG_DEG
SAMPLES = ("held4", "held3", "held", "mid", "grad", "big")  # the plan's order

# ---- step 19's and step 20's functions, imported (not copied) -----------------------------------
repair = CR.repair
fit_rng = CR.fit_rng
fit_kwargs = CR.fit_kwargs
contested = CR.contested
after_pass_state = CR.after_pass_state
other_end_hit = CR.other_end_hit
repair_candidates = CR.repair_candidates
wrong_mask = CR.wrong_mask
fixed_broken = CR.fixed_broken
misorientation_deg = CR.misorientation_deg
REPAIR_MAX_EVALS = CR.CMA_MAX_EVALS
EVALS_PER_FIT = CR.EVALS_PER_FIT

RULE = (
    "Arm J5S = arm J5 with screen tries in the boundary pass and step 20's repair on the wide "
    "flag. BFS and seeds as J5 (J5's stored map before the pass). The pass as the library's, every "
    "try a screen try: the candidate's start is scored (one evaluation, the library's final "
    "evaluation counted: 2 a try) and is the try's end; the keep rule 'hit' unchanged; each try's "
    "generator as J5's pass (make_rng(v, 'pass', k)). Then step 20's repair, imported unchanged "
    "(CMA 120 evaluations, first step 0.05 degree, one fit from the voxel's own end and one from "
    "every candidate of the snapshot, per-fit generators default_rng([0, v, c + 1]), keep rule "
    "'hit'), on the wide flag: a voxel is flagged if it is contested (hit ratio under 0.95 of its "
    "neighbours' median AND margin under 0.15; step 19's flag, the screen ends counting as the "
    "other ends) or if its hit ratio is under 0.95 and one of its pass candidates lies more than "
    "1 degree from its kept end. No voxel is skipped by its own hit ratio: every boundary voxel "
    "scores all its candidates. No truth is used."
)

# ---- stored counts (plan, "What the record says") -----------------------------------------------
# step 20's tables: J5's run with the pass, J5C's total, S8's total (evaluations)
TOTALS: Dict[str, Tuple[int, int, int]] = dict(
    held4=(27_350_255, 27_487_469, 68_033_600),
    held3=(23_503_385, 23_650_279, 67_408_448),
    held=(22_939_812, 23_063_958, 59_923_523),
    mid=(10_589_493, 10_666_328, 26_443_074),
    grad=(9_041_686, 9_158_451, 18_000_230),
    big=(3_567_098, 3_586_700, 8_340_428),
)
# g0: step 20's repair evaluations (J5C), which the full-tries replay reproduces
G0_J5C_EVALS: Dict[str, int] = dict(
    held4=137_214, held3=146_894, held=124_146, mid=76_835, grad=116_765, big=19_602
)
# g1 = prediction 1: the sizing, reproduced exactly
P1_PASS_EVALS: Dict[str, int] = dict(
    held4=181_534, held3=161_864, held=170_478, mid=100_058, grad=120_710, big=37_242
)
P1_FLAGGED: Dict[str, int] = dict(
    held4=1_869, held3=1_787, held=1_832, mid=1_187, grad=1_445, big=366
)
P1_REPAIR_EVALS: Dict[str, int] = dict(
    held4=1_617_286, held3=1_512_258, held=1_547_106, mid=759_275, grad=776_336, big=162_987
)
P1_TOTAL: Dict[str, int] = dict(
    held4=23_612_288,
    held3=20_240_655,
    held=19_457_817,
    mid=8_397_057,
    grad=6_257_077,
    big=2_631_446,
)
P1_FIXED: Dict[str, int] = dict(held4=3, held3=4, held=9, mid=3, grad=2, big=3)
P1_BROKEN: Dict[str, int] = dict(held4=2, held3=2, held=3, mid=3, grad=2, big=1)

GATE: Dict[str, Any] = dict(
    g0="the full-budget re-run of J5's pass tries (60 evaluations, the gain trace kept) "
    "reproduces J5's stored tries (count, evaluations, cost, hit ratio, end) and the after-pass "
    "map, and the full-tries replay (step 20's repair on its flag) reproduces J5C: repair "
    "evaluations 137,214 / 146,894 / 124,146 / 76,835 / 116,765 / 19,602 and J5C's wrong voxels "
    "(held4, held3, held, mid, grad, big)",
    g1="prediction 1, exactly: pass evaluations 181,534 / 161,864 / 170,478 / 100,058 / 120,710 / "
    "37,242; flagged 1,869 / 1,787 / 1,832 / 1,187 / 1,445 / 366; repair evaluations 1,617,286 / "
    "1,512,258 / 1,547,106 / 759,275 / 776,336 / 162,987; totals 23,612,288 / 20,240,655 / "
    "19,457,817 / 8,397,057 / 6,257,077 / 2,631,446; fixed / broken against J5C 3/2, 4/2, 9/3, "
    "3/3, 2/2, 3/1",
    g2="on each disclosed sample: at least 10% fewer evaluations than J5C, broken minus fixed "
    "against J5C at most 0, median error not higher, and on mid broken minus fixed against S8 at "
    "most +1",
    g2_min_saving_pct=10.0,
    g2_max_net_vs_j5c=0,
    mid_net_vs_s8_max=1,
)
JUDGE: Dict[str, Any] = dict(
    r1="at least 8% fewer evaluations than J5C",
    r2="broken minus fixed against J5C at most 2",
    r3="median error at most 0.001 degree higher than J5C's",
    r4="house rule against S8 with the pass: at least 20% fewer evaluations, broken minus fixed "
    "at most 5, median error at most 0.005 degree higher",
    step="J5S passes on held5 and on the six disclosed samples",
    min_saving_vs_j5c_pct=8.0,
    max_net_vs_j5c=2,
    max_median_rise_vs_j5c_deg=0.001,
    min_cut_vs_s8_pct=20.0,
    max_net_vs_s8=5,
    max_median_rise_vs_s8_deg=0.005,
    stop_line=dict(min_saving_vs_j5c_pct=8.0, max_net_vs_j5c=2, max_repair_share_of_j5_pct=10.0),
)
PREDICTIONS = {
    "1": "the tracked scripts reproduce the sizing exactly on (held4, held3, held, mid, grad, "
    "big): pass evaluations 181,534 / 161,864 / 170,478 / 100,058 / 120,710 / 37,242; "
    "flagged 1,869 / "
    "1,787 / 1,832 / 1,187 / 1,445 / 366; repair evaluations 1,617,286 / 1,512,258 / 1,547,106 / "
    "759,275 / 776,336 / 162,987; totals 23,612,288 / 20,240,655 / 19,457,817 / 8,397,057 / "
    "6,257,077 / 2,631,446; fixed / broken against J5C 3/2, 4/2, 9/3, 3/3, 2/2, 3/1",
    "2": "held5: 11% to 18% fewer evaluations than J5C; flagged 1,500 to 2,200; repair 5% to 8% of "
    "J5; against J5C fixed 2 to 12, broken 0 to 4; median error 0.001 to 0.004 degree lower",
    "3": "held5, against S8: J5C broken minus fixed -35 to -5; J5S within 3 of J5C",
}


# ---- the pass ----------------------------------------------------------------------------------


def pass_try(
    hooks: Any, bfs: Any, v: int, R0: np.ndarray, rng: np.random.Generator, budget: int
) -> Dict[str, Any]:
    """One try of the boundary pass at voxel ``v`` from ``R0``: the CMA fit of the pass (first step
    0.05 degree) with ``budget`` evaluations (the start's own included), plus the library's final
    evaluation. ``budget == FULL_MAX_EVALS`` is ``boundary_common.try_fit`` (J5's tries);
    ``budget == SCREEN_MAX_EVALS`` is a screen try."""
    rec = hooks.rec
    vox = hooks.mic.voxels[v]
    old = rec.params.local_optimizer
    rec.params.local_optimizer = "cma"
    kw: Dict[str, Any] = {}
    patched = False
    if budget < 2:
        maker = rec._make_cma_optimizer

        def make(*a: Any, **k: Any) -> Any:
            opt = maker(*a, **k)
            opt.max_evals = int(budget)  # CMAOptimizer takes 1; the reconstructor's override not
            return opt

        rec._make_cma_optimizer = make
        patched = True
    else:
        kw["cma_max_evals"] = int(budget)
    try:
        r = rec.local_optimization(
            voxel_vertices=hooks.vertices(v),
            phase_index=vox.phase,
            initial_orientation=np.asarray(R0, np.float64),
            rng=rng,
            cma_sigma0_deg=PASS_SIGMA0_DEG,
            **kw,
        )
    finally:
        rec.params.local_optimizer = old
        if patched:
            del rec._make_cma_optimizer
    hit, conf = bfs._ratios(r.overlap_info)
    return dict(
        end=np.asarray(r.orientation, np.float64),
        cost=float(r.cost),
        hit=float(hit),
        conf=float(conf),
        evals=int(rec.last_local_optimization_evals),
    )


# ---- the wide flag -----------------------------------------------------------------------------


def wide_flag(
    flag: np.ndarray,
    ratio: np.ndarray,
    t_voxel: np.ndarray,
    t_end: np.ndarray,
    R_after: np.ndarray,
    sym: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """(wide (n,) bool, extra (n,) bool). ``flag``/``ratio``: step 19's ``contested`` result on the
    state after the pass. A voxel is wide if contested, or if its ratio is under ``RATIO_MAX`` and
    one of its pass candidates (the try ends, rows ``t_voxel``/``t_end``) is more than ``WIDE_DEG``
    from its kept end ``R_after``. ``extra``: the wide voxels that are not contested. The
    contested part is the input flag, unchanged."""
    flag = np.asarray(flag, dtype=bool)
    n = len(flag)
    has_far = np.zeros(n, dtype=bool)
    if len(t_voxel):
        tv = np.asarray(t_voxel, dtype=np.int64)
        far = np.asarray(misorientation_deg(t_end, R_after[tv], sym), dtype=float) > WIDE_DEG
        has_far[tv[far]] = True
    wide = flag | ((np.asarray(ratio, dtype=float) < RATIO_MAX) & has_far)
    return wide, wide & ~flag


# ---- accounting --------------------------------------------------------------------------------


def total_evals(j5_total: int, stored_pass_evals: int, pass_evals: int, repair_evals: int) -> int:
    """The arm's evaluations: ``J5``'s run with its pass, minus the pass's tries as stored, plus the
    new pass's tries and the repair's (fits and the scores of voxels with no stored score)."""
    return int(j5_total) - int(stored_pass_evals) + int(pass_evals) + int(repair_evals)


def saving_pct(base: int, total: int) -> float:
    return 100.0 * (int(base) - int(total)) / int(base)


def judge(
    total: int,
    j5c_total: int,
    s8_total: int,
    net_vs_j5c: int,
    net_vs_s8: int,
    median: float,
    median_j5c: float,
    median_s8: float,
) -> Dict[str, Any]:
    """The plan's four rules of arm ``J5S`` on one sample. ``net_*``: broken minus fixed."""
    r1 = saving_pct(j5c_total, total)
    r4 = saving_pct(s8_total, total)
    ok1 = r1 >= JUDGE["min_saving_vs_j5c_pct"]
    ok2 = net_vs_j5c <= JUDGE["max_net_vs_j5c"]
    ok3 = median <= median_j5c + JUDGE["max_median_rise_vs_j5c_deg"]
    ok4 = (
        r4 >= JUDGE["min_cut_vs_s8_pct"]
        and net_vs_s8 <= JUDGE["max_net_vs_s8"]
        and median <= median_s8 + JUDGE["max_median_rise_vs_s8_deg"]
    )
    return dict(
        saving_vs_j5c_pct=r1,
        cut_vs_s8_pct=r4,
        ok1=bool(ok1),
        ok2=bool(ok2),
        ok3=bool(ok3),
        ok4=bool(ok4),
        ok=bool(ok1 and ok2 and ok3 and ok4),
    )


def gate_g2(
    sample: str,
    total: int,
    j5c_total: int,
    fixed: int,
    broken: int,
    median: float,
    median_j5c: float,
    net_vs_s8: int,
) -> Dict[str, Any]:
    """The gate's test g2 for one disclosed sample (fixed / broken against J5C)."""
    save = saving_pct(j5c_total, total)
    ok_save = save >= GATE["g2_min_saving_pct"]
    ok_net = (broken - fixed) <= int(GATE["g2_max_net_vs_j5c"])
    ok_med = median <= median_j5c
    out: Dict[str, Any] = dict(
        rule="g2",
        sample=sample,
        saving_vs_j5c_pct=save,
        fixed=fixed,
        broken=broken,
        median=median,
        median_j5c=median_j5c,
        ok_saving=bool(ok_save),
        ok_net=bool(ok_net),
        ok_median=bool(ok_med),
    )
    ok = ok_save and ok_net and ok_med
    if sample == "mid":
        out["net_vs_s8"] = net_vs_s8
        out["ok_s8"] = bool(net_vs_s8 <= int(GATE["mid_net_vs_s8_max"]))
        ok = ok and out["ok_s8"]
    out["ok"] = bool(ok)
    return out


def pass_evals_of(tries_evals: Sequence[int]) -> int:
    return int(np.sum(np.asarray(tries_evals, dtype=np.int64))) if len(tries_evals) else 0


def screen_flags(
    R: np.ndarray,
    cost: np.ndarray,
    hit: np.ndarray,
    old_hit: np.ndarray,
    Rb: np.ndarray,
    t_voxel: np.ndarray,
    t_end: np.ndarray,
    t_hit: np.ndarray,
    nbrs: Sequence[Sequence[int]],
    sym: np.ndarray,
) -> Dict[str, Any]:
    """The flags of arm ``J5S`` on the state after the pass (``R``, ``cost``, ``hit`` from
    ``after_pass_state`` on the screen tries): the contested voxels, the wide set, the extra
    voxels (listed apart) and the candidates of the wide set on the snapshot ``R``."""
    other = other_end_hit(R, Rb, old_hit, t_voxel, t_end, t_hit, sym)
    flag, ratio, margin = contested(hit, other, nbrs)
    wide, extra = wide_flag(flag, ratio, t_voxel, t_end, R, sym)
    flagged: List[int] = [int(v) for v in np.where(wide)[0]]
    cands = repair_candidates(R, nbrs, sym, flagged)
    return dict(
        other=other,
        contested=flag,
        ratio=ratio,
        margin=margin,
        wide=wide,
        extra=extra,
        flagged=flagged,
        cands=cands,
        fits=sum(len(cands[v]) + 1 for v in flagged),
    )
