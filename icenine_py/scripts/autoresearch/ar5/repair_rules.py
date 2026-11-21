"""AR5 step 19 (task ar5s): the repair of contested voxels after the boundary pass (pure: numpy and
two pure functions of the library; no image, no reconstructor, no truth).

After the pass (arm ``J5``), flag the **contested** voxels and refit them with the run's MC local
optimization. All terms are the plan's (``PLAN.md``, "AR5 step 19"):

* **hit ratio** of a voxel: the one its end has after the pass (the old end's score or the kept
  try's); **ratio**: that over the median hit ratio of the voxel's neighbours (radius 2 side
  lengths); **margin**: that hit ratio minus the best hit among the voxel's *other* pass ends (the
  old end and the try ends more than 1 degree from the kept end); a voxel with none is never
  flagged.
* **contested**: ratio under ``RATIO_MAX`` (0.95) and margin under ``MARGIN_MAX`` (0.15), both
  strict. No truth is used.
* **candidates**: as the pass (``icenine.boundary_pass.boundary_candidates``): the neighbours'
  orientations more than 0.5 degree from the voxel's own, merged within 0.5 degree, read from a
  snapshot of the map after the pass.
* **repair**: for each contested voxel, in index order: one MC local optimization (``optimizer=
  "mc"``) from its own end and one from each candidate of the snapshot; keep, among the current end
  and these ends, the highest hit ratio, then the lowest cost, then the current end
  (``boundary_pass.keep_hit``). Its own generator ``default_rng(0)`` (one for the whole repair, so
  no other draw moves). Every evaluation is counted. A voxel's hit ratio is the stored one; a
  missing one (NaN) costs one counted evaluation (``score``).

The fit itself is injected (``fit(v, R0, rng) -> (end, cost, hit, evals)``), so the rule is tested
without images; ``repair_run.py`` supplies the reconstructor's MC fit.
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from icenine.boundary_pass import boundary_candidates, keep_hit
from icenine.ray_indexing import misorientation_deg

RATIO_MAX = 0.95  # flag: ratio under this
MARGIN_MAX = 0.15  # flag: margin under this
FAR_DEG = 1.0  # a pass end counts as "other" when more than this from the kept end
CAND_DEG = 0.5  # the pass's boundary / merge distance
WRONG_DEG = 1.0  # right: at most this from the truth (symmetry-reduced)
REPAIR_SEED = 0  # the repair's own generator
REPAIR_OPTIMIZER = "mc"  # the run's classic variance-minimizing MC fit
NO_OTHER = 1.0  # margin of a voxel with no other pass end (never flagged: 1.0 >= MARGIN_MAX)

SAMPLES = ("held3", "held", "mid", "grad", "big")  # the plan's order
ARM = "J5R"

RULE = (
    "Arm J5R = arm J5 plus a repair after its boundary pass. Contested voxel: after-pass hit ratio "
    "under 0.95 of the median hit ratio of its neighbours (radius 2 side lengths) AND margin under "
    "0.15 (its hit ratio minus the best hit among its other pass ends: the old end and the try "
    "ends more than 1 degree from the kept end; a voxel with none is never flagged). For each, in "
    "index order: one MC local optimization (the run's classic variance-minimizing fit, "
    "optimizer='mc') from its own end and one from each candidate of the after-pass snapshot "
    "(neighbour orientations more than 0.5 degree from its own, merged within 0.5 degree); keep, "
    "among the current end and these ends, the highest hit ratio, then the lowest cost, then the "
    "current end. Its own generator default_rng(0). Every evaluation counted. No truth is used."
)

# stored counts (plan stop rule): J5 total evaluations with the pass; J5 wrong voxels after it
J5_TOTAL: Dict[str, int] = dict(
    held3=23_503_385, held=22_939_812, mid=10_589_493, grad=9_041_686, big=3_567_098
)
J5_WRONG: Dict[str, int] = dict(held3=54, held=67, mid=35, grad=56, big=21)
# the plan's table (g0): flagged voxels and repair fits (a fit per start: own end plus candidates)
FLAGGED: Dict[str, int] = dict(held3=194, held=183, mid=142, grad=238, big=56)
FITS: Dict[str, int] = dict(held3=1_214, held=1_026, mid=635, grad=965, big=162)
J5L = dict(total=9_990_663, wrong=39, flagged=145, fits=650)  # step 18's live run of mid
# D-20261115-N7: neighbours are the pass's own (mic.get_neighbors, radius 2 side lengths), the
# angle is the library's misorientation_deg; the plan's first counts (193 / 183 / 143 / 238 / 56
# and 1,217 / 1,035 / 643 / 967 / 161) used a cKDTree ball and an arccos angle and are superseded.

# gate (plan "Decision rules", fixed before any run)
GATE = dict(
    g0="prediction 1 (amended, D-20261115-N7): flagged 194 / 183 / 142 / 238 / 56 and fits "
    "1,214 / 1,026 / 635 / 965 / 162 (held3, held, mid, grad, big); J5L on mid 145 / 650; "
    "J5 wrong 54 / 67 / 35 / 56 / 21; J5 totals as stored",
    g1="J5R on mid: fixed minus broken against J5 at least 3; repair evaluations at most 2.5% "
    "of J5",
    g2="J5R on held3, held, grad, big: broken minus fixed against J5 at most 0; repair evaluations "
    "at most 2.5% of J5 (3.5% on grad)",
    g1_min_net_fixed=3,
    g2_max_net_broken=0,
    max_share_pct=2.5,
    max_share_pct_grad=3.5,
)
# per-sample pass of arm J5R (disclosed samples and held4)
JUDGE = dict(
    r1="repair evaluations at most 2.5% of J5's run with the pass (3.5% on grad)",
    r2="broken minus fixed against J5 at most 0",
    r3="median error at most 0.001 degree higher than J5's",
    r4="house rule against S8 with the pass: at least 20% fewer evaluations, broken minus fixed at "
    "most 5, median error at most 0.005 degree higher",
    step="J5R passes on held4 and on the five disclosed samples, and on mid broken minus fixed "
    "against S8 is at most +1",
    max_share_pct=2.5,
    max_share_pct_grad=3.5,
    max_net_broken_vs_j5=0,
    max_median_rise_vs_j5_deg=0.001,
    min_cut_vs_s8_pct=20.0,
    max_net_broken_vs_s8=5,
    max_median_rise_vs_s8_deg=0.005,
    mid_net_vs_s8_max=1,
)

PREDICTIONS = {
    "1": "the scripts reproduce the table's flagged voxels and fits, and the stored counts",
    "2": "MC fits at contested voxels take 100 to 250 evaluations (median); repair cost 0.5% to 2% "
    "of J5 on held3, held, mid, big, 1% to 3.5% on grad",
    "3": "fixed against J5: held3 6 to 14, held 7 to 16, mid 4 to 10, grad 10 to 23, big 2 to 4; "
    "broken 0 to 2 per sample; median error moves at most 0.0005 degree",
    "4": "mid after the repair: broken minus fixed against S8 from -6 to +1",
    "5": "held4: repair cost under 2% of J5; fixed minus broken against J5 at least 3",
    "6": "J5L repaired, on mid: -7 to +2 against S8",
}


def max_share(sample: str) -> float:
    """The cap (percent of J5's run with the pass) on the repair's evaluations."""
    return float(GATE["max_share_pct_grad"] if sample == "grad" else GATE["max_share_pct"])


# ---- the state after the pass --------------------------------------------------------------


def after_pass_state(
    R_before: np.ndarray,
    old_score: np.ndarray,
    t_voxel: np.ndarray,
    t_end: np.ndarray,
    t_cost: np.ndarray,
    t_hit: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(R after (n, 3, 3), cost (n,), hit (n,)) of every voxel after the pass, from the stored
    tries: the old end's score (``old_score`` rows: cost, hit, confidence) and, per voxel, its tries
    (rows grouped by voxel in order made). The kept end is rule ``hit`` (``keep_hit``)."""
    R = R_before.copy()
    cost = np.asarray(old_score[:, 0], dtype=float).copy()
    hit = np.asarray(old_score[:, 1], dtype=float).copy()
    order = np.argsort(t_voxel, kind="stable")
    tv = np.asarray(t_voxel)[order]
    if len(tv):
        edges = np.flatnonzero(np.diff(tv)) + 1
        lo = np.concatenate([[0], edges])
        hi = np.concatenate([edges, [len(tv)]])
        for a, b in zip(lo, hi):
            v = int(tv[a])
            rows = order[a:b]
            k = keep_hit((cost[v], hit[v]), [(float(t_cost[r]), float(t_hit[r])) for r in rows])
            if k > 0:
                r = rows[k - 1]
                R[v] = t_end[r]
                cost[v] = float(t_cost[r])
                hit[v] = float(t_hit[r])
    return R, cost, hit


def other_end_hit(
    R_after: np.ndarray,
    R_before: np.ndarray,
    old_hit: np.ndarray,
    t_voxel: np.ndarray,
    t_end: np.ndarray,
    t_hit: np.ndarray,
    sym: np.ndarray,
) -> np.ndarray:
    """(n,) best hit among each voxel's other pass ends: the try ends more than ``FAR_DEG`` from the
    kept end, and the old end when it is more than ``FAR_DEG`` from the kept end; -1 for none."""
    n = len(R_after)
    other = np.full(n, -1.0)
    if len(t_voxel):
        tv = np.asarray(t_voxel)
        far = np.asarray(misorientation_deg(t_end, R_after[tv], sym), dtype=float) > FAR_DEG
        np.maximum.at(other, tv[far], np.asarray(t_hit, dtype=float)[far])
    oldfar = np.asarray(misorientation_deg(R_before, R_after, sym), dtype=float) > FAR_DEG
    return np.where(oldfar, np.maximum(other, np.asarray(old_hit, dtype=float)), other)


def margin_of(hit: np.ndarray, other: np.ndarray) -> np.ndarray:
    """Margin: hit minus the best other end's hit; ``NO_OTHER`` where there is none."""
    return np.where(other < 0, NO_OTHER, hit - other)


def ratio_of(hit: np.ndarray, nbrs: Sequence[Sequence[int]]) -> np.ndarray:
    """Hit over the median hit of the neighbours (inf when a voxel has no neighbour)."""
    med = np.array(
        [np.median(hit[np.asarray(x, dtype=np.int64)]) if len(x) else np.nan for x in nbrs]
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        r = hit / med
    return np.where(np.isnan(med), np.inf, r)


def contested(
    hit: np.ndarray, other: np.ndarray, nbrs: Sequence[Sequence[int]]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(flag (n,) bool, ratio, margin): ratio under ``RATIO_MAX``, margin under ``MARGIN_MAX``."""
    ratio = ratio_of(hit, nbrs)
    margin = margin_of(hit, other)
    return (ratio < RATIO_MAX) & (margin < MARGIN_MAX), ratio, margin


def repair_candidates(
    R: np.ndarray,
    nbrs: Sequence[Sequence[int]],
    sym: np.ndarray,
    voxels: Sequence[int],
) -> Dict[int, List[int]]:
    """{voxel: candidate voxels} of ``voxels`` on the snapshot ``R``, as the pass does it
    (``boundary_candidates``); a voxel that is no boundary voxel of ``R`` has none."""
    _isb, cands = boundary_candidates(R, nbrs, sym)
    return {int(v): list(cands.get(int(v), [])) for v in voxels}


# ---- the repair ----------------------------------------------------------------------------

FitFn = Callable[[int, np.ndarray, np.random.Generator], Tuple[np.ndarray, float, float, int]]
ScoreFn = Callable[[int, np.ndarray], Tuple[float, float]]


@dataclass
class RepairResult:
    """What one repair did. ``fits`` has one dict per fit (voxel, start voxel or None for the own
    end, cost, hit, evals) in the order made."""

    R: np.ndarray
    cost: np.ndarray
    hit: np.ndarray
    flagged: List[int] = field(default_factory=list)
    changed: List[int] = field(default_factory=list)
    fits: List[Dict[str, Any]] = field(default_factory=list)
    n_evals_fits: int = 0
    n_evals_score: int = 0

    @property
    def n_evals(self) -> int:
        return self.n_evals_fits + self.n_evals_score


def repair(
    R: np.ndarray,
    cost: np.ndarray,
    hit: np.ndarray,
    flagged: Sequence[int],
    cands: Dict[int, List[int]],
    fit: FitFn,
    score: Optional[ScoreFn] = None,
    seed: int = REPAIR_SEED,
) -> RepairResult:
    """The repair of the plan on a map after the pass. ``R``, ``cost``, ``hit``: the state after
    the pass (not changed; copies are returned). ``flagged``: the contested voxels (taken in index
    order); ``cands``: their candidate voxels read from the snapshot ``R``. ``fit(v, R0, rng)``: one
    MC local optimization of voxel ``v`` from ``R0`` -> (end, cost, hit, evaluations). ``score(v,
    R)`` is called, and one evaluation counted, for a flagged voxel whose stored hit or cost is
    NaN."""
    snap = np.asarray(R, dtype=np.float64).copy()
    out = RepairResult(
        R=snap.copy(),
        cost=np.asarray(cost, dtype=float).copy(),
        hit=np.asarray(hit, dtype=float).copy(),
        flagged=sorted(int(v) for v in flagged),
    )
    rng = np.random.default_rng(seed)  # its own generator: no other draw moves
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
            end, ec, eh, ev = fit(v, np.asarray(R0, dtype=np.float64), rng)
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


# ---- scoring (truth is used here, and only here) --------------------------------------------


def wrong_mask(R: np.ndarray, R_true: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """(n,) True where the end is more than ``WRONG_DEG`` from the truth (symmetry-reduced)."""
    return np.asarray(misorientation_deg(R, R_true, sym), dtype=float) > WRONG_DEG


def fixed_broken(wrong_old: np.ndarray, wrong_new: np.ndarray) -> Tuple[int, int]:
    """(fixed, broken) of the newer arm against the older: wrong in the older and right in the
    newer / the reverse, counted by voxel."""
    a, b = np.asarray(wrong_old, dtype=bool), np.asarray(wrong_new, dtype=bool)
    return int(np.sum(a & ~b)), int(np.sum(~a & b))


def gate_g12(sample: str, fixed: int, broken: int, share_pct: float) -> Dict[str, Any]:
    """The gate's test for one sample (g1 on mid, g2 on the others) as a dict with ``ok``."""
    cap = max_share(sample)
    if sample == "mid":
        ok_net = (fixed - broken) >= int(GATE["g1_min_net_fixed"])
        rule = "g1"
    else:
        ok_net = (broken - fixed) <= int(GATE["g2_max_net_broken"])
        rule = "g2"
    ok_cost = share_pct <= cap
    return dict(
        rule=rule,
        sample=sample,
        fixed=fixed,
        broken=broken,
        share_pct=share_pct,
        cap_pct=cap,
        ok_net=bool(ok_net),
        ok_cost=bool(ok_cost),
        ok=bool(ok_net and ok_cost),
    )
