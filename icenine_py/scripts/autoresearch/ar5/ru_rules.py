"""AR5 step 29 (task ar5ad), phase A: rule ``RU``, the gates and the predictions, as constants and
pure functions (no study imports besides ``px_rules``; no file is opened here).

``ru_freeze`` writes the constants into ``ar5/ru/frozen.json`` before any run. The code the live
runs add is ``make_index_seed`` (a wrapper of the seed step of the pick, ``SeedPick._index_seed``,
which the library's ``BFSReconstruction._index_seed`` calls) on top of step 28's wrappers.

Rule ``RU`` = before every seed step, try the orientations of the 8 nearest voxels accepted so far
(one strict cost each at the voxel, then one CMA fit from the lowest-cost one, first step 0.5
degree, 200 evaluations, then a second fit from the lowest-cost start over 0.5 degree away from
the first); a fit that ends at hit ratio 0.95 or more is the seed result, else the step-28 (``PXF``)
seed step runs unchanged. See ``RULE``.
"""

import time
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

import px_rules as PX

K_NEAREST = 8
SIGMA0_DEG = 0.5
MAX_EVALS = 200
SEP_DEG = 0.5  # the second start must be more than this from the first
ACCEPT = 0.95  # hit ratio of a fit end that is taken as the seed result
OLD_TOTALS = {"mid": 6_822_192, "held": 16_580_548}  # step 28's stored PXF J5SD runs
OLD_TOTAL = sum(OLD_TOTALS.values())  # 23,402,740
OLD_WALL_S = 4_545  # step 28's PXF runs, pooled
FULL_SAMPLES = ("mid", "held")
NEW_SAMPLE = "held8"
GENERATOR_SEED = 8  # held8 = held with generator seed 8: a new sample, not an extra random seed

RULE = (
    "Rule RU. (1) Before every seed step (seed or reseed) at voxel v, take the voxels accepted so "
    "far in the run (seed test passed or fitted in the BFS, not rejected since: state FITTED) "
    "and keep the 8 nearest to v by centre distance (ties by voxel index; any distance; fewer if "
    "fewer exist; none: go to (4)). (2) Take the strict cost of each one's orientation at v (one "
    "evaluation each). From the lowest-cost one (ties: nearer first), run one CMA local "
    "optimization: first step 0.5 degree, 200 evaluations, a fresh default_rng(0). If the end's "
    "hit ratio (one overlap evaluation, not counted, as in the library) is 0.95 or more, the end "
    "is the seed result. Then there is no indexing and no seed step, and the BFS goes on as after "
    "a passing seed step. (3) Else, run the same fit from the lowest-cost start more than 0.5 "
    "degree from the first (if any), with the same acceptance. (4) Else the PXF seed step runs "
    "unchanged. All try evaluations count as the seed's. Nothing else changes. No truth, no "
    "threshold read from data after this freeze, nothing crystal-specific (misorientation under "
    "the sample's symmetry group; any group works). Script layer only: a wrapper of the seed "
    "step of the pick that the library's BFSReconstruction._index_seed calls, on top of PXF's."
)

GATES: Dict[str, Any] = dict(
    r1=(
        "evaluations: new J5SD (PXF and RU) at least 15% fewer than the stored PXF runs "
        "(6,822,192 + 16,580,548 = 23,402,740), summed over mid and held"
    ),
    r2=(
        "accuracy against stored PXF (16 and 40 wrong; per sample): broken minus fixed at most 3; "
        "median error at most 0.001 degree higher; seed steps reaching the grid search at most 2 "
        "more"
    ),
    r3=(
        "wall: pooled wall of the new runs at least 15% lower than step 28's PXF runs (4,545 s), "
        "judged at any load and labelled contended; preflight verdicts and sampled load reported"
    ),
    h1=(
        "house rule on held8: new J5SD against S8, at least 20% fewer evaluations, broken minus "
        "fixed at most 5, median error at most 0.005 degree higher"
    ),
    h2=(
        "on held8 against the baseline: new against PXF J5SD, at least 10% fewer evaluations, "
        "broken minus fixed at most 5, median error at most 0.001 degree higher; walls of all "
        "three arms reported, contended"
    ),
)
PREDICTIONS = dict(
    seed_steps_ended_by_a_try="mid 72% to 85%, held 62% to 76%",
    r1="18% to 32% fewer (centre 28%)",
    r2="broken minus fixed -3 to +3; grid searches mid 0 to 1, held 0 to 3",
    r3="22% to 38% lower",
    h1="passes",
    h2="15% to 30% fewer",
)
STOP_RULES = (
    "A test fails. A try accepts an end under 0.95, or uses a voxel not accepted at that moment "
    "(checked from the log). An index call at a voxel that a try accepted. r1, r2 or r3 fails: "
    "no held8. held2 opened or listed, fine7 read, held8 built before the gates pass. A gate is "
    "never loosened after a run; every refused preflight is kept."
)
IN_SAMPLE = (
    "mid and held are in-sample for this rule (first step and second start chosen on the probe)"
)
TIMING_RULE = (
    "r3 is judged with any load. timing_preflight is run before each timed job and its verdict "
    "KEPT (a refusal is recorded, not waited on); the load is sampled every 30 s; every second is "
    "labelled with the verdict. One run at a time."
)


def fewer(new: float, old: float) -> float:
    return PX.fewer(new, old)


def judge_r1(new_total: int, old_total: int = OLD_TOTAL) -> Dict[str, Any]:
    f = fewer(new_total, old_total)
    return dict(new=int(new_total), old=int(old_total), fewer=float(f), ok=bool(f >= 0.15))


def judge_r2(broken_minus_fixed: int, median_higher: float, grid_extra: int) -> Dict[str, Any]:
    return PX.judge_c3(broken_minus_fixed, median_higher, grid_extra)


def judge_r3(new_wall: float, old_wall: float = float(OLD_WALL_S)) -> Dict[str, Any]:
    return PX.judge_t(fewer(new_wall, old_wall), 0.15)


def judge_h1(fewer_: float, broken_minus_fixed: int, median_higher: float) -> Dict[str, Any]:
    return PX.judge_h1(fewer_, broken_minus_fixed, median_higher)


def judge_h2(fewer_: float, broken_minus_fixed: int, median_higher: float) -> Dict[str, Any]:
    return dict(
        fewer=float(fewer_),
        net=int(broken_minus_fixed),
        median_higher=float(median_higher),
        ok=bool(fewer_ >= 0.10 and broken_minus_fixed <= 5 and median_higher <= 0.001),
    )


# ---- the nearest accepted voxels


def nearest_accepted(
    tree: Any,
    pos: np.ndarray,
    v: int,
    is_accepted: Callable[[int], bool],
    k: int = K_NEAREST,
) -> List[int]:
    """The ``k`` voxels nearest to voxel ``v`` (centre distance; ties by voxel index) among those
    for which ``is_accepted`` is true, ``v`` itself left out; fewer if fewer exist. ``tree`` is a
    ``scipy.spatial.cKDTree`` of ``pos`` (shape (N, d)). The query grows until the k-th accepted
    distance is strictly inside the queried radius, so that no tie at the cut is missed."""
    n = len(pos)
    kq = min(n, max(64, 4 * k))
    while True:
        _, ii = tree.query(pos[v], k=kq)
        ii = np.atleast_1d(ii).astype(int)
        d = np.linalg.norm(pos[ii] - pos[v], axis=1)
        keep = np.array([int(j) != v and bool(is_accepted(int(j))) for j in ii], dtype=bool)
        dk, ik = d[keep], ii[keep]
        order = np.lexsort((ik, np.round(dk, 9)))
        out = ik[order][:k]
        if kq >= n:
            return [int(j) for j in out]
        if len(dk) >= k:
            kth = np.round(dk[order][k - 1], 9)
            if kth < np.round(d.max(), 9):
                return [int(j) for j in out]
        kq = min(n, kq * 4)


# ---- the try


def second_start(
    order: Sequence[int],
    R: Sequence[np.ndarray],
    mis_deg: Callable[[np.ndarray, np.ndarray], float],
    sep: float = SEP_DEG,
) -> Optional[int]:
    """Position (in ``R``) of the lowest-cost start (``order`` is the cost order) that is more than
    ``sep`` degrees from the lowest-cost one; None if there is none."""
    best = int(order[0])
    for j in order[1:]:
        if mis_deg(np.asarray(R[int(j)]), np.asarray(R[best])) > sep:
            return int(j)
    return None


def run_try(
    start_ids: Sequence[int],
    orient_of: Callable[[int], np.ndarray],
    strict_cost: Callable[[np.ndarray], float],
    fit: Callable[[np.ndarray], Tuple[np.ndarray, int, float]],
    hit_of: Callable[[np.ndarray], float],
    mis_deg: Callable[[np.ndarray, np.ndarray], float],
    accept: float = ACCEPT,
) -> Tuple[Optional[np.ndarray], Dict[str, Any]]:
    """Rule RU around one seed step. ``start_ids``: the nearest accepted voxels (distance order).
    ``fit(R0)`` runs the CMA fit from ``R0`` and returns (end orientation, evaluations, end cost);
    ``hit_of(R)`` is the library's hit ratio of ``R`` at the voxel (not counted).

    Returns (the end that is the seed result, or None) and a row: starts, costs, each fit's
    evaluations and hit ratio and seconds, ``evals`` (costs plus fits), ``success``."""
    row: Dict[str, Any] = dict(
        tried=False,
        starts=[int(j) for j in start_ids],
        costs=[],
        fits=[],
        evals=0,
        success=False,
        second=False,
    )
    if len(start_ids) == 0:
        return None, row
    row["tried"] = True
    R = [np.asarray(orient_of(int(j)), np.float64) for j in start_ids]
    costs = [float(strict_cost(r)) for r in R]
    row["costs"] = costs
    row["evals"] = len(R)
    order = [int(j) for j in np.argsort(np.asarray(costs), kind="stable")]
    best = order[0]
    seq: List[Tuple[int, bool]] = [(best, False)]
    sec = second_start(order, R, mis_deg)
    for j, is_second in seq + ([(sec, True)] if sec is not None else []):
        t0 = time.perf_counter()
        r_end, ev, c_end = fit(R[j])
        dt = time.perf_counter() - t0
        hit = float(hit_of(r_end))
        row["evals"] += int(ev)
        row["fits"].append(
            dict(
                start=int(start_ids[j]),
                evals=int(ev),
                hit=hit,
                cost=float(c_end),
                s=float(dt),
                second=is_second,
            )
        )
        row["second"] = bool(is_second)
        if hit >= accept:
            row["success"] = True
            return np.asarray(r_end, np.float64), row
    return None, row


def make_index_seed(
    orig: Callable[..., Any],
    log: List[Dict[str, Any]],
    ctx: Callable[[Any, Any, int], Dict[str, Any]],
    result_cls: Callable[..., Any],
) -> Callable[..., Any]:
    """``SeedPick._index_seed`` with rule RU in front. ``ctx(self, vertices, phase_index)`` returns
    the pieces the rule needs for the voxel in progress: ``voxel``, ``reseed``, ``near`` (the
    nearest accepted voxel ids, checked by the caller), ``orient_of``, ``strict_cost``, ``fit``,
    ``hit_of``, ``mis_deg``.
    A success returns ``result_cls`` (a ``SeedResult``) with the fit's end, no indexing and no seed
    step; a failure runs ``orig`` unchanged and adds the try's evaluations and seconds to its
    result (they count as the seed's)."""

    def index_seed(self: Any, vertices: Any, phase_index: int, *a: Any, **k: Any) -> Any:
        c = ctx(self, vertices, phase_index)
        t0 = time.perf_counter()
        end, row = run_try(
            c["near"], c["orient_of"], c["strict_cost"], c["fit"], c["hit_of"], c["mis_deg"]
        )
        dt = time.perf_counter() - t0
        row.update(voxel=int(c["voxel"]), reseed=bool(c["reseed"]), try_s=float(dt))
        if end is not None:
            fin = row["fits"][-1]
            assert fin["hit"] >= ACCEPT, "a try accepted an end under 0.95"
            row.update(indexing_calls=0, seed_step=False)
            log.append(row)
            return result_cls(
                orientation=end,
                cost=float(fin["cost"]),
                taken=-1,
                n_cand=0,
                n_att=0,
                evals=int(row["evals"]),
                index_s=0.0,
                total_s=float(dt),
            )
        res = orig(self, vertices, phase_index, *a, **k)
        row.update(indexing_calls=1, seed_step=True)
        log.append(row)
        if row["tried"]:
            res.evals = int(res.evals) + int(row["evals"])
            res.total_s = float(res.total_s) + float(dt)
        return res

    return index_seed
