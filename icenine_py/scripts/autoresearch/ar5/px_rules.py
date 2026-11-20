"""AR5 step 27 (task ar5aa): rule ``PX``, the gates and the predictions, as constants and pure
functions (no study imports, no file is opened here).

``px_freeze`` writes the constants into ``ar5/px/frozen.json`` before any run. The arm wrapper
``make_index_voxel`` is the only code the live runs add: it passes ``max_pairs`` through
``IndexParams`` (arm ``new``) and logs ``n_pairs``, ``n_votes`` and the seconds of every indexed
voxel (both arms; the indexer function itself is the old ``index_voxel`` in both).
"""

import dataclasses
import time
from typing import Any, Callable, Dict, List, Sequence, Tuple

NEW_MAX_PAIRS = 50_000
OLD_MAX_PAIRS = 600_000
SEED_GENERATOR = 27  # 300 stored seed steps per sample (q1)
N_PER_SAMPLE = 300
FULL_SAMPLES = ("mid", "held")  # t1: every stored seed voxel
Q1_SAMPLES = ("mid", "held", "held3", "held4", "held5", "held6", "grad", "big")
BLOCK = 100

RULE = (
    "Rule PX: the indexer runs with max_pairs 50,000 in place of 600,000. Nothing else changes: "
    "the pair generator (seed 0), tolerances, n_seeds 150, top_k 20, the seed step, rescue, BFS, "
    "screen and repair. Both arms use the same indexer function (index_voxel). No truth, no "
    "threshold read from data after this freeze, nothing crystal-specific. Chosen on 400 voxels "
    "of mid and held: in-sample, so held7 is the clean test."
)

GATES: Dict[str, Any] = dict(
    q1=(
        "candidates, no evaluations: over the 2,400 voxels (300 stored seed steps of each of mid, "
        "held, held3, held4, held5, held6, grad, big; generator 27), the finally accepted "
        "orientation leaves the top 20 in at most 0.5% of the voxels where it was in it under "
        "base (gained ones reported apart)"
    ),
    c1=(
        "indexer work: summed ray pairs and votes of the new live J5SD runs at least 60% lower "
        "than those of the stored J5SD seed voxels indexed again with base settings (mid and "
        "held pooled; the old live runs' own logs reported beside)"
    ),
    c2="evaluations: new J5SD at most 2% more than stored (7,383,368 + 17,114,589), summed",
    c3=(
        "accuracy against stored J5SD (23 and 45 wrong; mid and held): broken minus fixed at "
        "most 3 per sample; median error at most 0.001 degree higher; seed steps that reach the "
        "grid search (grid_fallbacks_left) at most 2 more per sample"
    ),
    t1=(
        "paired, index only: old and new settings on the 5,311 stored seed voxels of mid and "
        "held, alternating blocks of 100, one session: summed indexing seconds at least 50% lower"
    ),
    t2=(
        "paired, whole run: J5SD old and new in the order old mid, new mid, new held, old held, "
        "one session, preflight before each, load sampled every 30 s: pooled wall at least 25% "
        "lower. t1 and t2 are judged with any load and reported with the preflight verdicts"
    ),
    h1=(
        "house rule on held7 (a non-exact change): new J5SD against S8, at least 20% fewer "
        "evaluations, broken minus fixed at most 5, median error at most 0.005 degree higher"
    ),
)
PREDICTIONS = dict(
    q1="0 to 6 lost",
    c1="pairs 70% to 80% lower",
    c2="within 1%",
    c3="-3 to +3",
    t1="55% to 80% lower",
    t2="mid 25% to 40%, held 40% to 60% less",
    h1="passes",
)
STORED_TOTALS = dict(mid=7_383_368, held=17_114_589)  # stored J5SD (step 22)
TIMING_RULE = (
    "t1 and t2 are judged with any load. timing_preflight is run before each timed job and its "
    "verdict KEPT (a refusal is recorded, not waited on); the load is sampled every 30 s; every "
    "second is labelled with the verdict"
)


def fewer(new: float, old: float) -> float:
    """Share fewer (positive: ``new`` is smaller)."""
    return 1.0 - float(new) / float(old) if old else 0.0


def judge_q1(lost: int, was_in: int, gained: int) -> Dict[str, Any]:
    return dict(lost=lost, was_in_top20=was_in, gained=gained, ok=bool(lost <= 0.005 * was_in))


def judge_c1(new_pairs: int, new_votes: int, old_pairs: int, old_votes: int) -> Dict[str, Any]:
    fp, fv = fewer(new_pairs, old_pairs), fewer(new_votes, old_votes)
    return dict(pairs_fewer=fp, votes_fewer=fv, ok=bool(fp >= 0.60 and fv >= 0.60))


def judge_c2(new_total: int) -> Dict[str, Any]:
    old = sum(STORED_TOTALS.values())
    return dict(
        new=int(new_total),
        stored=old,
        more=new_total / old - 1.0,
        ok=bool(new_total <= 1.02 * old),
    )


def judge_c3(broken_minus_fixed: int, median_higher: float, grid_extra: int) -> Dict[str, Any]:
    return dict(
        net=int(broken_minus_fixed),
        median_higher=float(median_higher),
        grid_extra=int(grid_extra),
        ok=bool(broken_minus_fixed <= 3 and median_higher <= 0.001 and grid_extra <= 2),
    )


def judge_t(fewer_seconds: float, need: float) -> Dict[str, Any]:
    return dict(fewer=float(fewer_seconds), need=need, ok=bool(fewer_seconds >= need))


def judge_h1(fewer_: float, broken_minus_fixed: int, median_higher: float) -> Dict[str, Any]:
    return dict(
        fewer=float(fewer_),
        net=int(broken_minus_fixed),
        median_higher=float(median_higher),
        ok=bool(fewer_ >= 0.20 and broken_minus_fixed <= 5 and median_higher <= 0.005),
    )


def make_index_voxel(
    orig: Callable[..., Tuple[Any, Any]],
    arm: str,
    log: List[Dict[str, Any]],
    max_pairs: int = NEW_MAX_PAIRS,
) -> Callable[..., Tuple[Any, Any]]:
    """``index_voxel`` with the arm's setting and a log row per call.

    ``arm`` ``new``: the ``IndexParams`` passed on carry ``max_pairs``; ``old``: unchanged. Both
    log ``n_rays``, ``n_pairs``, ``n_votes`` and the seconds of the call. The wrapped function
    is called with the same arguments otherwise."""
    if arm not in ("old", "new"):
        raise ValueError("arm must be 'old' or 'new'")

    def wrapped(sc: Any, voxel_vertices: Any, crystal: Any, p: Any, *a: Any, **k: Any) -> Any:
        if arm == "new":
            p = dataclasses.replace(p, max_pairs=max_pairs)
        t0 = time.perf_counter()
        res, att = orig(sc, voxel_vertices, crystal, p, *a, **k)
        dt = time.perf_counter() - t0
        log.append(
            dict(
                n_rays=int(res.n_rays),
                n_pairs=int(res.n_pairs),
                n_votes=int(res.n_votes),
                s=float(dt),
                max_pairs=int(p.max_pairs),
            )
        )
        return res, att

    return wrapped


def summed(log: Sequence[Dict[str, Any]]) -> Dict[str, float]:
    return dict(
        calls=len(log),
        pairs=int(sum(r["n_pairs"] for r in log)),
        votes=int(sum(r["n_votes"] for r in log)),
        seconds=float(sum(r["s"] for r in log)),
    )
