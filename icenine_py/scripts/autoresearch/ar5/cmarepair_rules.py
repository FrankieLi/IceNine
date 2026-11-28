"""AR5 step 20 (task ar5t): the repair of contested voxels with a 120-evaluation CMA fit in place of
step 19's MC fit (pure: numpy and the pure functions of step 19's frozen ``repair_rules``; no image,
no reconstructor, no truth).

The flag, the state after the pass, the candidates, the keep rule and the scoring helpers are step
19's frozen functions, imported from ``repair_rules`` and not copied (a change there changes the
hash this step freezes). What this step changes is the fit and its generators:

* **CMA fit**: one local CMA run of the library (``optimizer="cma"``), first step 0.05 degree (the
  pass's try step), at most 120 evaluations with the start's own, plus the library's final
  evaluation: 121 in all.
* **Per-fit generator**: the fit at voxel ``v`` from start voxel ``c`` (``c = -1`` for the voxel's
  own end) draws from ``default_rng([0, v, c + 1])``. No stop, budget or fit depends on another
  fit, so no draw moves when a fit is added, dropped or cut.
* **repair**: after the pass, on a snapshot of the map: for each flagged voxel in index order, one
  CMA fit from its own end and one from each candidate of the snapshot; keep, among the current
  end and these ends, the highest hit ratio, then the lowest cost, then the current end
  (``boundary_pass.keep_hit``). A voxel's hit ratio is the stored one; a missing one costs one
  counted evaluation. Every evaluation is counted. No truth is used.

The fit itself is injected (``fit(v, R0, rng) -> (end, cost, hit, evals)``), so the rule is tested
without images; ``cmarepair_run.py`` supplies the reconstructor's CMA fit.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import repair_rules as R19  # step 19's frozen rule: the flag and the keep rule are imported

from icenine.boundary_pass import keep_hit

# ---- step 19's frozen functions and constants, imported (not copied) ---------------------------
RATIO_MAX = R19.RATIO_MAX
MARGIN_MAX = R19.MARGIN_MAX
FAR_DEG = R19.FAR_DEG
CAND_DEG = R19.CAND_DEG
WRONG_DEG = R19.WRONG_DEG
after_pass_state = R19.after_pass_state
other_end_hit = R19.other_end_hit
contested = R19.contested
repair_candidates = R19.repair_candidates
wrong_mask = R19.wrong_mask
fixed_broken = R19.fixed_broken
RepairResult = R19.RepairResult
misorientation_deg = R19.misorientation_deg

# ---- this step's fit ---------------------------------------------------------------------------
REPAIR_SEED = 0  # first entry of every per-fit generator
REPAIR_OPTIMIZER = "cma"
CMA_SIGMA0_DEG = 0.05  # first step: the pass's try step
CMA_MAX_EVALS = 120  # at most this many evaluations with the start's own ...
EVALS_PER_FIT = CMA_MAX_EVALS + 1  # ... plus the library's final evaluation
NO_OTHER = R19.NO_OTHER

SAMPLES = R19.SAMPLES  # held3, held, mid, grad, big (disclosed)
ARM = "J5C"

RULE = (
    "Arm J5C = arm J5 plus a repair after its boundary pass. Contested voxel: step 19's frozen "
    "flag (after-pass hit ratio under 0.95 of the median hit ratio of its neighbours, radius 2 "
    "side lengths, AND margin under 0.15; a voxel with no other pass end is never flagged), "
    "imported from repair_rules. For each, in index order: one CMA fit (the library's local CMA "
    "run, optimizer='cma', first step 0.05 degree, at most 120 evaluations with the start's own, "
    "plus the final evaluation: 121 in all) from its own end and one from each candidate of the "
    "after-pass snapshot (neighbour orientations more than 0.5 degree from its own, merged "
    "within 0.5 degree); each fit with its own generator default_rng([0, v, c + 1]) (c = -1 for "
    "the voxel's own end); keep, among the current end and these ends, the highest hit ratio, then "
    "the lowest cost, then the current end. Every evaluation counted. No truth is used."
)

# stored counts (plan stop rule), step 19's
J5_TOTAL: Dict[str, int] = dict(R19.J5_TOTAL)
J5_WRONG: Dict[str, int] = dict(R19.J5_WRONG)
FLAGGED: Dict[str, int] = dict(R19.FLAGGED)
FITS: Dict[str, int] = dict(R19.FITS)
J5L = dict(R19.J5L)

# prediction 1 (plan): the sizing, reproduced exactly by the tracked scripts
P1_EVALS: Dict[str, int] = dict(held3=146_894, held=124_146, mid=76_835, grad=116_765, big=19_602)
P1_FIXED: Dict[str, int] = dict(held3=14, held=18, mid=15, grad=31, big=14)
P1_BROKEN: Dict[str, int] = dict(held3=0, held=0, mid=0, grad=0, big=0)
P1_J5L = dict(flagged=145, fits=650, evals=78_650, fixed=18, broken=0)

# gate (plan "Decision rules", fixed before any run)
GATE: Dict[str, Any] = dict(
    g0="step 19's flags and J5 totals reproduce: flagged 194 / 183 / 142 / 238 / 56 and fits "
    "1,214 / 1,026 / 635 / 965 / 162 (held3, held, mid, grad, big); J5L on mid 145 / 650; J5 wrong "
    "54 / 67 / 35 / 56 / 21; J5 totals 23,503,385 / 22,939,812 / 10,589,493 / 9,041,686 / "
    "3,567,098",
    g1="prediction 1, exactly: repair evaluations 146,894 / 124,146 / 76,835 / 116,765 / 19,602; "
    "fixed / broken 14/0, 18/0, 15/0, 31/0, 14/0; J5L 145 / 650, 78,650, 18/0",
    g2="on each disclosed sample: repair evaluations at most 1.5% of J5 (2.0% on grad), fixed "
    "minus broken against J5 at least 3, and on mid broken minus fixed against S8 at most +1",
    g2_min_net_fixed=3,
    max_share_pct=1.5,
    max_share_pct_grad=2.0,
    mid_net_vs_s8_max=1,
)
# per-sample pass of arm J5C (disclosed samples and held4)
JUDGE: Dict[str, Any] = dict(
    r1="repair evaluations at most 1.5% of J5's run with the pass (2.0% on grad)",
    r2="broken minus fixed against J5 at most 0",
    r3="median error at most 0.001 degree higher than J5's",
    r4="house rule against S8 with the pass: at least 20% fewer evaluations, broken minus fixed at "
    "most 5, median error at most 0.005 degree higher",
    step="J5C passes on held4 and on the five disclosed samples, and on mid broken minus fixed "
    "against S8 is at most +1",
    max_share_pct=1.5,
    max_share_pct_grad=2.0,
    max_net_broken_vs_j5=0,
    max_median_rise_vs_j5_deg=0.001,
    min_cut_vs_s8_pct=20.0,
    max_net_broken_vs_s8=5,
    max_median_rise_vs_s8_deg=0.005,
    mid_net_vs_s8_max=1,
)

PREDICTIONS = {
    "1": "the tracked scripts reproduce the sizing exactly: flagged / fits 194 / 1,214, 183 / "
    "1,026, 142 / 635, 238 / 965, 56 / 162; repair evaluations 146,894 / 124,146 / 76,835 / "
    "116,765 / 19,602; fixed / broken 14/0, 18/0, 15/0, 31/0, 14/0; J5L 145 / 650, 78,650, 18/0",
    "2": "held4: flagged 150 to 250; repair cost 0.4% to 1.0% of J5; fixed against J5 8 to 25, "
    "broken 0 to 2; median error not higher by more than 0.0001 degree",
    "3": "held4, against S8: J5 broken minus fixed -10 to +5; J5C -30 to -5",
}


def max_share(sample: str) -> float:
    """The cap (percent of J5's run with the pass) on the repair's evaluations."""
    return float(GATE["max_share_pct_grad"] if sample == "grad" else GATE["max_share_pct"])


def fit_kwargs() -> Dict[str, Any]:
    """The keyword arguments of ``Reconstructor.local_optimization`` for one repair fit."""
    return dict(
        optimizer=REPAIR_OPTIMIZER, cma_sigma0_deg=CMA_SIGMA0_DEG, cma_max_evals=CMA_MAX_EVALS
    )


def fit_rng(v: int, c: Optional[int]) -> np.random.Generator:
    """The generator of the fit at voxel ``v`` from start ``c`` (None: the voxel's own end)."""
    cc = -1 if c is None else int(c)
    return np.random.default_rng([REPAIR_SEED, int(v), cc + 1])


# ---- the repair ----------------------------------------------------------------------------

FitFn = Any  # fit(v, R0, rng) -> (end, cost, hit, evals)
ScoreFn = Any  # score(v, R) -> (cost, hit)


def repair(
    R: np.ndarray,
    cost: np.ndarray,
    hit: np.ndarray,
    flagged: Sequence[int],
    cands: Dict[int, List[int]],
    fit: FitFn,
    score: Optional[ScoreFn] = None,
) -> RepairResult:
    """The repair of the plan on a map after the pass. ``R``, ``cost``, ``hit``: the state after
    the pass (not changed; copies are returned). ``flagged``: the contested voxels (taken in index
    order); ``cands``: their candidate voxels read from the snapshot ``R``. ``fit(v, R0, rng)``:
    one CMA fit of voxel ``v`` from ``R0`` -> (end, cost, hit, evaluations); ``rng`` is the fit's
    own ``fit_rng(v, c)``. ``score(v, R)`` is called, and one evaluation counted, for a flagged
    voxel whose stored hit or cost is NaN."""
    snap = np.asarray(R, dtype=np.float64).copy()
    out = RepairResult(
        R=snap.copy(),
        cost=np.asarray(cost, dtype=float).copy(),
        hit=np.asarray(hit, dtype=float).copy(),
        flagged=sorted(int(v) for v in flagged),
    )
    for v in out.flagged:
        if np.isnan(out.hit[v]) or np.isnan(out.cost[v]):
            if score is None:
                raise ValueError(f"voxel {v} has no stored hit ratio and no score function")
            out.cost[v], out.hit[v] = score(v, snap[v])
            out.n_evals_score += 1
        starts: List[Tuple[Optional[int], np.ndarray]] = [(None, snap[v])]
        starts += [(int(c), snap[int(c)]) for c in cands.get(v, [])]
        rows: List[Dict[str, Any]] = []
        for c, R0 in starts:
            end, ec, eh, ev = fit(v, np.asarray(R0, dtype=np.float64), fit_rng(v, c))
            rows.append(
                dict(
                    voxel=v,
                    start=c,
                    end=np.asarray(end, dtype=np.float64),
                    cost=float(ec),
                    hit=float(eh),
                    evals=int(ev),
                )
            )
        out.fits += rows
        out.n_evals_fits += sum(r["evals"] for r in rows)
        k = keep_hit((out.cost[v], out.hit[v]), [(r["cost"], r["hit"]) for r in rows])
        if k > 0:
            new = rows[k - 1]
            out.R[v] = new["end"]
            out.cost[v] = new["cost"]
            out.hit[v] = new["hit"]
            out.changed.append(v)
    return out


def gate_g2(sample: str, fixed: int, broken: int, share_pct: float) -> Dict[str, Any]:
    """The gate's test g2 for one disclosed sample (cost and net fixed against J5)."""
    cap = max_share(sample)
    ok_net = (fixed - broken) >= int(GATE["g2_min_net_fixed"])
    ok_cost = share_pct <= cap
    return dict(
        rule="g2",
        sample=sample,
        fixed=fixed,
        broken=broken,
        share_pct=share_pct,
        cap_pct=cap,
        ok_net=bool(ok_net),
        ok_cost=bool(ok_cost),
        ok=bool(ok_net and ok_cost),
    )
