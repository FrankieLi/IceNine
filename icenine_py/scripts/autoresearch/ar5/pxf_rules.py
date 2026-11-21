"""AR5 step 28 (task ar5ab): rule ``PXF``, the gates and the predictions, as constants and pure
functions (no study imports besides ``px_rules``; no file is opened here).

``pxf_freeze`` writes the constants into ``ar5/pxf/frozen.json`` before any run. The only code the
live runs add is ``make_e1_step`` (a wrapper of the seed step of the pick, ``DeferPick.e1_step``)
and ``make_index_voxel`` (the indexer wrapper that passes ``max_pairs`` and logs the work).

Rule ``PXF`` = step 27's rule ``PX`` (indexer at ``max_pairs`` 50,000) plus a second try: a seed
step that fails the seed test after an indexing that kept a subsample is run again with the full
sample (600,000 pairs). See ``RULE``.
"""

import dataclasses
import time
from typing import Any, Callable, Dict, List, Sequence, Tuple

import px_rules as PX

FIRST_MAX_PAIRS = PX.NEW_MAX_PAIRS  # 50,000
FULL_MAX_PAIRS = PX.OLD_MAX_PAIRS  # 600,000
STORED_TOTALS = PX.STORED_TOTALS  # stored J5SD (step 22): mid 7,383,368; held 17,114,589
FULL_SAMPLES = ("mid", "held")
OLD_WALL_S = 8_219  # step 27's old runs, pooled (2,020 s mid + 6,199 s held)
TEST = 0.80  # the seed test (min_acceleration_threshold of the config; read from it at run time)

RULE = (
    "Rule PXF. (1) A seed step (seed or reseed) indexes with max_pairs 50,000, as PX. (2) If it "
    "fails the seed test (hit ratio under 0.80, or no candidate) and that indexing kept a "
    "subsample (the voxel had more than 50,000 ray pairs), the same voxel is indexed again with "
    "600,000 pairs and the whole seed step runs again on the new candidates, same pick, same "
    "generators, no candidate skipped. (3) If that step passes, it is the result. Else G1, G2b "
    "and the grid search run as now, from the second step's kept ends. (4) Both steps count as "
    "the seed's evaluations; both indexings count in c1. Nothing else changes (pair generator 0, "
    "n_seeds 150, top_k 20, tolerances, BFS, screen, repair). No truth, no threshold read from "
    "data, nothing crystal-specific. Script layer only: a wrapper of the seed step of the pick. "
    "A subsample is detected as n_pairs == 50,000 after the first indexing (a voxel with exactly "
    "50,000 pairs repeats an identical step; its evaluations are counted)."
)

GATES: Dict[str, Any] = dict(
    c1=(
        "indexer work: summed ray pairs and votes of the new live J5SD runs, the fallback "
        "indexings included, at least 50% lower than those of the stored J5SD seed voxels indexed "
        "again with base settings (step 27's t1 rows; mid and held pooled)"
    ),
    c2="evaluations: new J5SD at most 2% more than stored (7,383,368 + 17,114,589), summed",
    c3=(
        "accuracy against stored J5SD (23 and 45 wrong; mid and held): broken minus fixed at "
        "most 3 per sample; median error at most 0.001 degree higher; seed steps that reach the "
        "grid search (grid_fallbacks_left) at most 2 more per sample"
    ),
    t2=(
        "whole run: pooled wall of the new J5SD runs at least 25% lower than step 27's old runs "
        "(8,219 s), judged with any load; preflight verdicts and sampled load reported"
    ),
    h1=(
        "house rule on held7 (a non-exact change): new J5SD against S8, at least 20% fewer "
        "evaluations, broken minus fixed at most 5, median error at most 0.005 degree higher; "
        "wall seconds of both reported"
    ),
)
PREDICTIONS = dict(
    c1="68% to 76% fewer pairs",
    c2="-5% to +1% (centre -2.5%)",
    fallbacks="mid 7 to 14, held 50 to 70",
    to_rescue="mid 2 to 7, held 15 to 30",
    c3="broken minus fixed -6 to +3; grid searches mid 0 to 1, held 1 to 4",
    t2="35% to 50% lower (step 27's PX: 44.9%)",
    h1="passes",
)
EXTRA_RUN = (
    "Not a gate: plain J5SD (600,000 pairs, no fallback) on held7 after S8 and the new J5SD, "
    "reported with evaluations, wrong voxels, grid searches and wall seconds beside the others"
)
TIMING_RULE = (
    "t2 is judged with any load. timing_preflight is run before each timed job and its verdict "
    "KEPT (a refusal is recorded, not waited on); the load is sampled every 30 s; every second is "
    "labelled with the verdict. One run at a time."
)
IN_SAMPLE = "mid and held are in-sample for this rule (it was designed after seeing step 27's runs)"


def fewer(new: float, old: float) -> float:
    return PX.fewer(new, old)


def judge_c1(new_pairs: int, new_votes: int, old_pairs: int, old_votes: int) -> Dict[str, Any]:
    fp, fv = fewer(new_pairs, old_pairs), fewer(new_votes, old_votes)
    return dict(pairs_fewer=fp, votes_fewer=fv, ok=bool(fp >= 0.50 and fv >= 0.50))


judge_c2 = PX.judge_c2
judge_c3 = PX.judge_c3
judge_h1 = PX.judge_h1


def judge_t2(new_wall: float, old_wall: float = float(OLD_WALL_S)) -> Dict[str, Any]:
    return PX.judge_t(fewer(new_wall, old_wall), 0.25)


class Controller:
    """The indexer setting of the step in progress and the log of every indexing."""

    def __init__(self) -> None:
        self.max_pairs = FIRST_MAX_PAIRS
        self.stage = "first"
        self.log: List[Dict[str, Any]] = []


def make_index_voxel(
    orig: Callable[..., Tuple[Any, Any]], ctl: Controller, mode: str
) -> Callable[..., Tuple[Any, Any]]:
    """``index_voxel`` with the controller's setting (mode ``pxf``) or unchanged (mode ``old``:
    logging only). Every call logs ``n_rays``, ``n_pairs``, ``n_votes``, its seconds, the stage."""
    if mode not in ("pxf", "old"):
        raise ValueError("mode must be 'pxf' or 'old'")

    def wrapped(sc: Any, voxel_vertices: Any, crystal: Any, p: Any, *a: Any, **k: Any) -> Any:
        if mode == "pxf":
            p = dataclasses.replace(p, max_pairs=ctl.max_pairs)
        t0 = time.perf_counter()
        res, att = orig(sc, voxel_vertices, crystal, p, *a, **k)
        dt = time.perf_counter() - t0
        ctl.log.append(
            dict(
                n_rays=int(res.n_rays),
                n_pairs=int(res.n_pairs),
                n_votes=int(res.n_votes),
                s=float(dt),
                max_pairs=int(p.max_pairs),
                stage=ctl.stage if mode == "pxf" else "old",
            )
        )
        return res, att

    return wrapped


def summed(log: Sequence[Dict[str, Any]]) -> Dict[str, float]:
    return PX.summed(log)


def run_pxf_step(
    step: Callable[[], Any],
    hit_of: Callable[[Any], float],
    ctl: Controller,
    threshold: float = TEST,
) -> Tuple[Any, Dict[str, Any]]:
    """Rule PXF around one seed step. ``step()`` runs the seed step once (the indexer reads
    ``ctl.max_pairs``) and returns a result with ``orientation`` (None: no candidate) and
    ``evals``; ``hit_of(out)`` is the hit ratio the library's seed test uses.

    Returns the result to hand on and a row. The first step runs at 50,000 pairs. If it fails
    the seed test and its indexing kept a subsample, the second step runs at 600,000 pairs and
    is the result (its ``evals`` and ``total_s`` include the first step's)."""
    n0 = len(ctl.log)
    ctl.max_pairs, ctl.stage = FIRST_MAX_PAIRS, "first"
    out1 = step()
    rows1 = ctl.log[n0:]
    hit1 = None if out1.orientation is None else float(hit_of(out1))
    failed = hit1 is None or hit1 < threshold
    subsampled = bool(rows1) and rows1[-1]["n_pairs"] >= FIRST_MAX_PAIRS
    row: Dict[str, Any] = dict(
        failed_first=bool(failed),
        subsampled=subsampled,
        fallback=False,
        first=dict(
            evals=int(out1.evals),
            hit=hit1,
            n_pairs=int(rows1[-1]["n_pairs"]) if rows1 else 0,
            n_votes=int(rows1[-1]["n_votes"]) if rows1 else 0,
            index_s=float(sum(r["s"] for r in rows1)),
        ),
    )
    if not (failed and subsampled):
        return out1, row
    n1 = len(ctl.log)
    ctl.max_pairs, ctl.stage = FULL_MAX_PAIRS, "fallback"
    try:
        out2 = step()
    finally:
        ctl.max_pairs, ctl.stage = FIRST_MAX_PAIRS, "first"
    rows2 = ctl.log[n1:]
    hit2 = None if out2.orientation is None else float(hit_of(out2))
    ev1, ev2 = int(out1.evals), int(out2.evals)
    row.update(
        fallback=True,
        passed_second=bool(hit2 is not None and hit2 >= threshold),
        second=dict(
            evals=ev2,
            hit=hit2,
            n_pairs=int(rows2[-1]["n_pairs"]) if rows2 else 0,
            n_votes=int(rows2[-1]["n_votes"]) if rows2 else 0,
            index_s=float(sum(r["s"] for r in rows2)),
        ),
    )
    out2.evals = ev1 + ev2  # both steps count as the seed's evaluations
    if hasattr(out2, "total_s") and hasattr(out1, "total_s"):
        out2.total_s = float(out1.total_s) + float(out2.total_s)
    return out2, row


def make_e1_step(
    orig: Callable[..., Any],
    ctl: Controller,
    log: List[Dict[str, Any]],
    hit_of: Callable[[Any, Any, Any, Any], float],
    threshold_of: Callable[[Any], float],
) -> Callable[..., Any]:
    """The pick's ``e1_step`` with rule PXF around it. ``hit_of(self, out, vertices, phase)`` is
    the library's seed-test hit ratio; ``threshold_of(self)`` its threshold. Each call adds one
    row to ``log`` and, when the pick keeps an ``e1_log``, tags its new rows with the stage."""

    def e1_step(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2, tau: float = -1
    ) -> Any:
        e1 = getattr(self, "e1_log", None)
        n_before = len(e1) if e1 is not None else 0

        def step() -> Any:
            return orig(self, vertices, phase_index, n_cost=n_cost, n_refine=n_refine, tau=tau)

        out, row = run_pxf_step(
            step, lambda o: hit_of(self, o, vertices, phase_index), ctl, threshold_of(self)
        )
        row.update(voxel=getattr(self, "cur", None), reseed=getattr(self, "reseed", None))
        log.append(row)
        if e1 is not None:
            for i, x in enumerate(e1[n_before:]):
                x["pxf_stage"] = "first" if i == 0 else "fallback"
        return out

    return e1_step
