"""AR5 step 10 (task ar5j): the level-0 cut of the grid search (pure; no study imports).

The grid search's level 0 lists about 10,200 candidates, sorted by the discrete score, and runs a
quick MC refinement (10 steps, 5 restarts) on every one. The cut runs the quick MC only on the
first ``N`` of that list (N = 1,000); the quarter kept for level 1 is counted from the cut list
(the library takes ``int(len(list) * 0.25)`` of the list it is handed). Levels 1 to 3, the final
local optimization, rule ``G1``, 0.80 and 0.95 are unchanged. The cut is applied here by wrapping
the discrete stage the reconstructor calls (``run_discrete_search_spaced``) so that, on level 0
only, the list it returns is cut; the wrapper draws no random number and passes every argument
through, so the aligned draws of the grid call (voxel, "grid", 0) are those of arm ``G``.
Level 0 is the call whose start set is the reconstructor's full fundamental-zone set (the same
object; later levels build a new array).
"""

from typing import Any, Callable, Dict, List

TOP_N = 1000  # the quick MC of level 0 runs on this many listed candidates (frozen)
KEEP_FRACTION = 0.25  # the library's keep fraction (counted from the list handed on)

RULE = (
    "Level-0 cut (arm C = arm G + cut): in every grid search, the quick MC refinement of level 0 "
    "runs only on the first N = 1,000 candidates of the level-0 list (sorted by the discrete "
    "score); the quarter kept for level 1 is counted from the cut list. Levels 1 to 3, the final "
    "local optimization, rule G1, 0.80 and 0.95 are unchanged; no draw is added or moved."
)


def cut_list(candidates: List[Any], n: int) -> List[Any]:
    """The first ``n`` of the sorted list (the list itself when ``n`` is 0 or not smaller)."""
    if n <= 0 or len(candidates) <= n:
        return candidates
    return list(candidates[:n])


def n_keep_after_cut(n_listed: int, n: int) -> int:
    """How many candidates level 0 hands to level 1: the quarter of the (cut) list."""
    m = n_listed if n <= 0 else min(n_listed, n)
    return max(1, int(m * KEEP_FRACTION))


def make_cut_search(
    inner: Callable[..., Any],
    n: int,
    full_fz: Callable[[], Any],
    log: List[Dict[str, Any]],
) -> Callable[..., Any]:
    """``inner`` (the discrete stage) with the cut on level 0. ``full_fz()`` returns the
    reconstructor's full start set; a call is level 0 when its ``fz_orientations`` is that object.
    One row per call is appended to ``log``: ``level0``, ``n_listed``, ``n_run``."""

    def wrapped(**kw: Any) -> Any:
        out = inner(**kw)
        level0 = kw.get("fz_orientations") is full_fz()
        run = cut_list(out, n) if level0 else out
        log.append(dict(level0=bool(level0), n_listed=len(out), n_run=len(run)))
        return run

    return wrapped
