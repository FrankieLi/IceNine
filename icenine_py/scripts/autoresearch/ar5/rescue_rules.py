"""AR5 step 9 (task ar5i): rule ``G1``, the rescue of a failed seed step (pure; no study imports).

After a failed seed test the kept ends of the seed step are tried in order of cost: a CMA local
optimization from each (first step 0.25 degree, 1,000 evaluations, ``default_rng(0)``), at most 5;
the first end with hit ratio at least 0.95 is taken. If none reaches it, the grid fallback runs as
today. This file holds the constants, the choice rule, and the wrapper that puts the rule in front
of the grid fallback (``make_wrapper``); the tests use it with fakes. The wrapper issues no
random draw of the aligned-draw book (a rescue has its own ``default_rng(0)``), and the grid call
that follows a failed rescue is made unchanged (same arguments), so its aligned draws are those
of arm ``A``.
"""

from typing import Any, Callable, Dict, List, Optional

TAU = 0.95  # accept a rescued end at this hit ratio (the library's SeedIndexStopHitRatio)
SIGMA0_DEG = 0.25  # initial step of the rescue CMA, degrees
EVALS = 1000  # evaluations of one rescue CMA
MAX_TRIES = 5  # at most this many kept ends are tried (the seed step keeps at most 5)
SEED_TEST = 0.80  # the seed test (hit ratio); a step below it is "failed"

RULE = (
    "Rule G1: after a failed seed test (taken end hit ratio < 0.80), try the kept ends of the "
    "seed step in order of cost (stable; at most 5): a CMA local optimization from the end "
    "(first step 0.25 degree, 1,000 evaluations, default_rng(0)); take the first end with hit "
    "ratio >= 0.95 (one overlap evaluation, not counted) and stop; if none, run the grid "
    "fallback as today (aligned draws unchanged). Arm G = arm A (interleaved order, aligned "
    "draws, run S11) with G1."
)


def rescue_order(costs: List[float], max_tries: int = MAX_TRIES) -> List[int]:
    """Indices of the kept ends in the order they are tried: by cost, ties by position."""
    idx = sorted(range(len(costs)), key=lambda i: (float(costs[i]), i))
    return idx[:max_tries]


def rescue_choose(
    costs: List[float],
    try_fn: Callable[[int], Dict[str, Any]],
    tau: float = TAU,
    max_tries: int = MAX_TRIES,
) -> Dict[str, Any]:
    """Run the rule. ``try_fn(k)`` rescues kept end ``k`` and returns a dict with ``hit`` and
    ``evals`` (and anything else, kept in the try row). Returns ``rescued`` (bool), ``taken``
    (position in the tries, -1 if none), ``tries`` (rows with ``k``) and ``evals`` (all tries)."""
    tries: List[Dict[str, Any]] = []
    taken = -1
    for k in rescue_order(costs, max_tries):
        t = dict(try_fn(k))
        t["k"] = int(k)
        tries.append(t)
        if float(t["hit"]) >= tau:
            taken = len(tries) - 1
            break
    return dict(
        rescued=taken >= 0,
        taken=taken,
        tries=tries,
        evals=int(sum(int(t["evals"]) for t in tries)),
    )


def make_wrapper(
    grid_fn: Callable[..., Any],
    get_step: Callable[[], Optional[Dict[str, Any]]],
    key_fn: Callable[[Any], int],
    try_end: Callable[[Any, int, Any], Dict[str, Any]],
    rec: Any,
    log: List[Dict[str, Any]],
) -> Callable[..., Any]:
    """The grid fallback with rule G1 in front of it.

    ``get_step()`` returns the last seed step: ``key`` (voxel key), ``voxel``, ``kept`` (list of
    dicts ``R``, ``cost``, ``hit``) and ``used`` (set here); ``try_end(R, vertices, phase)``
    returns a dict with ``hit``, ``cost``, ``evals``, ``result`` (the local optimization's
    return). ``rec._last_global_evals`` / ``_last_local_evals`` are the counters behind
    ``last_eval_counts`` (the BFS reads the total after the call: it must include the rescues).
    One row per call is appended to ``log``.
    """

    def wrapped(voxel_vertices: Any, phase_index: int, rng: Any = None, **kw: Any) -> Any:
        step = get_step()
        key = key_fn(voxel_vertices)
        if step is None or step["used"] or step["key"] != key or not step["kept"]:
            out = grid_fn(voxel_vertices=voxel_vertices, phase_index=phase_index, rng=rng, **kw)
            log.append(
                dict(
                    voxel=None if step is None else step["voxel"],
                    attempted=False,
                    rescued=False,
                    tries=[],
                    rescue_evals=0,
                    grid_evals=int(rec.last_eval_counts[2]),
                )
            )
            return out
        step["used"] = True
        kept = step["kept"]

        def one(k: int) -> Dict[str, Any]:
            return try_end(kept[k], voxel_vertices, phase_index)

        r = rescue_choose([e["cost"] for e in kept], one)
        # log rows: no return object; an end orientation as a flat list (the log is json)
        rows = [
            {
                a: (b.ravel().tolist() if hasattr(b, "ravel") else b)
                for a, b in t.items()
                if a != "result"
            }
            for t in r["tries"]
        ]
        if r["rescued"]:
            res = r["tries"][r["taken"]]["result"]
            rec._last_global_evals = 0
            rec._last_local_evals = int(r["evals"])
            log.append(
                dict(
                    voxel=step["voxel"],
                    attempted=True,
                    rescued=True,
                    tries=rows,
                    rescue_evals=int(r["evals"]),
                    grid_evals=0,
                )
            )
            return res
        out = grid_fn(voxel_vertices=voxel_vertices, phase_index=phase_index, rng=rng, **kw)
        g = int(rec.last_eval_counts[2])
        rec._last_local_evals = int(rec._last_local_evals) + int(r["evals"])
        log.append(
            dict(
                voxel=step["voxel"],
                attempted=True,
                rescued=False,
                tries=rows,
                rescue_evals=int(r["evals"]),
                grid_evals=g,
            )
        )
        return out

    return wrapped
