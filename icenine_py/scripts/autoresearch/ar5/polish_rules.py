"""AR5 step 14 (task ar5n): rule ``G2``, a second local optimization from a failed rescue's ends,
before the grid search (pure; numpy only, no study imports).

After a failed seed test, rule ``G1`` (step 9) tries a CMA from each kept end (first step 0.25
degree, 1,000 evaluations) and takes the first end at hit ratio 0.95 or more. If none reaches it,
the grid search runs. Rule ``G2`` sits between: from every **rescue end** (the end of one ``G1``
try), in ``G1``'s order, a second, different local optimization; the first end at 0.95 or more is
taken; if none, the grid search runs unchanged (same arguments, so its aligned draws are those of
arm ``H``; ``G2`` takes no draw from the book and has its own generator ``default_rng(1)``).

Two forms, both frozen before any run:

* ``G2a`` (more of the same): one CMA from the rescue end, first step 0.5 degree, 3,000
  evaluations. At most 3,001 evaluations a try (the CMA's final overlap evaluation, which gives the
  hit ratio, is its last evaluation and is counted), 15,005 a case (5 ends).
* ``G2b`` (short starts nearby): the library's fixed local grid about the rescue end
  (``get_structured_local_grid``, side 1.5 degrees, level 2: 64 points about 0.375 degree apart;
  a candidate is ``grid[i] @ R``), each scored once with the strict cost (64 evaluations), then a
  CMA (first step 0.1 degree, 500 evaluations) from its 2 lowest-cost points (ties by position),
  in order of cost. At most 64 + 2 x 501 = 1,066 evaluations a try, 5,330 a case. The grid points
  themselves are never taken: only a CMA end can be accepted.

Every evaluation counts. A CMA's end hit ratio is that of its final overlap evaluation (in the
count), not a further one.
"""

import math
from typing import Any, Callable, Dict, List, Optional

import numpy as np

TAU = 0.95  # accept a G2 end at this hit ratio (the library's SeedIndexStopHitRatio)
RNG_SEED = 1  # G2's own generator, fresh for each CMA (G1 uses default_rng(0))
MAX_TRIES = 5  # at most 5 rescue ends (the seed step keeps at most 5)

FORMS: Dict[str, Dict[str, Any]] = {
    "G2a": dict(kind="cma", sigma0_deg=0.5, evals=3000),
    "G2b": dict(
        kind="grid_cma",
        side_deg=1.5,
        level=2,
        n_grid=64,
        n_starts=2,
        sigma0_deg=0.1,
        evals=500,
    ),
}

RULE = (
    "Rule G2: after a failed rescue (rule G1: no rescue end at hit ratio 0.95), try the rescue "
    "ends in G1's order (at most 5): G2a, a CMA from the end (first step 0.5 degree, 3,000 "
    "evaluations); or G2b, the library's local grid about the end (side 1.5 degrees, level 2, 64 "
    "points, candidate = grid[i] @ R, each scored once with the strict cost), then a CMA (first "
    "step 0.1 degree, 500 evaluations) from its 2 lowest-cost points. Own generator "
    "default_rng(1) for each CMA. Take the first CMA end with hit ratio >= 0.95 (the hit ratio of "
    "the CMA's final overlap evaluation, which is counted); if none, run the grid search "
    "unchanged (aligned draws as arm H). Every evaluation counts. Arm J = arm H with G2."
)


def try_bound(form: str) -> int:
    """Most evaluations of one try (one rescue end)."""
    f = FORMS[form]
    if f["kind"] == "cma":
        return int(f["evals"]) + 1
    return int(f["n_grid"]) + int(f["n_starts"]) * (int(f["evals"]) + 1)


def case_bound(form: str) -> int:
    """Most evaluations of one case (a failed rescue with 5 ends)."""
    return MAX_TRIES * try_bound(form)


def grid_matrices(form: str) -> np.ndarray:
    """(64, 3, 3) perturbation matrices of ``G2b`` (the library's local grid, side in radians)."""
    from icenine.sampling import generate_local_grid

    f = FORMS[form]
    return np.asarray(
        generate_local_grid(math.radians(float(f["side_deg"])), int(f["level"])), np.float64
    )


def candidates(grid: np.ndarray, R: np.ndarray) -> np.ndarray:
    """The grid about ``R``: ``grid[i] @ R`` (as the library's discrete search composes)."""
    return np.matmul(grid, np.asarray(R, np.float64).reshape(3, 3)[None])


def pick_starts(costs: List[float], n: int) -> List[int]:
    """Indices of the ``n`` lowest-cost grid points, ties by position."""
    idx = sorted(range(len(costs)), key=lambda i: (float(costs[i]), i))
    return idx[:n]


def polish_order(n_ends: int, max_tries: int = MAX_TRIES) -> List[int]:
    """The order of the tries: the rescue ends as ``G1`` tried them (the stored order)."""
    return list(range(min(n_ends, max_tries)))


def polish_choose(
    n_ends: int,
    try_fn: Callable[[int], Dict[str, Any]],
    tau: float = TAU,
    max_tries: int = MAX_TRIES,
) -> Dict[str, Any]:
    """Run the rule over the rescue ends. ``try_fn(j)`` runs ``G2`` from rescue end ``j`` and
    returns a dict with ``hit`` and ``evals`` (anything else is kept in the try row). Returns
    ``accepted`` (bool), ``taken`` (position in the tries, -1 if none), ``tries`` (rows with
    ``j``) and ``evals`` (all tries)."""
    tries: List[Dict[str, Any]] = []
    taken = -1
    for j in polish_order(n_ends, max_tries):
        t = dict(try_fn(j))
        t["j"] = int(j)
        tries.append(t)
        if float(t["hit"]) >= tau:
            taken = len(tries) - 1
            break
    return dict(
        accepted=taken >= 0,
        taken=taken,
        tries=tries,
        evals=int(sum(int(t["evals"]) for t in tries)),
    )


def one_try(
    form: str,
    R: np.ndarray,
    score: Callable[[np.ndarray], float],
    cma: Callable[[np.ndarray, float, int, Any], Dict[str, Any]],
) -> Dict[str, Any]:
    """One ``G2`` try from rescue end ``R``.

    ``score(R_candidate)`` returns the strict cost of one orientation (one counted evaluation);
    ``cma(R_start, sigma0_deg, max_evals, rng)`` runs one CMA and returns a dict with ``hit``,
    ``cost``, ``evals`` (counted, final overlap evaluation included), ``R`` (end) and ``result``.
    Returns ``hit``, ``cost``, ``evals``, ``R``, ``result`` of the best CMA end of the try (``G2a``:
    the one end; ``G2b``: the first of its 2 CMA ends at ``TAU``, else the one with the highest hit
    ratio, ties by cost), plus ``ends`` (every CMA end of the try) and ``grid_evals``.
    """
    f = FORMS[form]
    if f["kind"] == "cma":
        r = cma(np.asarray(R, np.float64), float(f["sigma0_deg"]), int(f["evals"]), None)
        r = dict(r)
        r["ends"] = [{k: v for k, v in r.items() if k != "result"}]
        r["grid_evals"] = 0
        return r
    grid = grid_matrices(form)
    cand = candidates(grid, R)
    costs = [float(score(cand[i])) for i in range(len(cand))]
    starts = pick_starts(costs, int(f["n_starts"]))
    ends: List[Dict[str, Any]] = []
    total = len(cand)
    pick: Optional[Dict[str, Any]] = None
    for i in starts:
        e = dict(cma(cand[i], float(f["sigma0_deg"]), int(f["evals"]), None))
        e["start"] = int(i)
        e["start_cost"] = costs[i]
        ends.append(e)
        total += int(e["evals"])
        if float(e["hit"]) >= TAU:
            pick = e
            break
    if pick is None:
        pick = sorted(ends, key=lambda e: (-float(e["hit"]), float(e["cost"])))[0]
    out = dict(pick)
    out["evals"] = int(total)
    out["grid_evals"] = len(cand)
    out["ends"] = [{k: v for k, v in e.items() if k != "result"} for e in ends]
    return out


def g2_try(
    form: str,
    rec: Any,
    orig_lo: Any,
    ratios: Any,
    R: np.ndarray,
    vertices: Any,
    phase: int,
) -> Dict[str, Any]:
    """One ``G2`` try from rescue end ``R`` (``orig_lo`` is the reconstructor's unwrapped
    ``local_optimization``; ``ratios(info)`` gives (hit ratio, confidence)). Every CMA has a fresh
    ``default_rng(1)``; the hit ratio is that of the CMA's final overlap evaluation (counted in
    its ``evals``), so no evaluation is added or hidden. ``G2b``'s grid scores are one strict
    evaluation each, counted."""

    def score(Rc: np.ndarray) -> float:
        return float(rec.evaluate_overlap(np.asarray(Rc, np.float64), vertices, phase).cost)

    def cma(Rs: np.ndarray, sigma: float, n: int, _rng: Any) -> Dict[str, Any]:
        r = orig_lo(
            voxel_vertices=vertices,
            phase_index=phase,
            initial_orientation=np.asarray(Rs, np.float64),
            rng=np.random.default_rng(RNG_SEED),
            cma_sigma0_deg=sigma,
            cma_max_evals=n,
            optimizer="cma",
        )
        ev = int(rec.last_local_optimization_evals)
        hit = float(ratios(r.overlap_info)[0])
        return dict(
            hit=hit, cost=float(r.cost), evals=ev, R=np.asarray(r.orientation, np.float64), result=r
        )

    return one_try(form, R, score, cma)


def make_wrapper(
    form: str,
    grid_fn: Callable[..., Any],
    get_step: Callable[[], Optional[Dict[str, Any]]],
    key_fn: Callable[[Any], int],
    g1_end: Callable[[Any, Any, Any], Dict[str, Any]],
    g2_end: Callable[[Any, Any, Any], Dict[str, Any]],
    rec: Any,
    log: List[Dict[str, Any]],
    g1_choose: Callable[..., Dict[str, Any]],
) -> Callable[..., Any]:
    """The grid fallback with rule ``G1`` and then rule ``G2`` in front of it.

    As step 9's wrapper (``rescue_rules.make_wrapper``), with the extra stage: ``g1_choose(costs,
    one)`` is ``rescue_rules.rescue_choose``; ``g1_end(R, vertices, phase)`` and ``g2_end(R,
    vertices, phase)`` return dicts with ``hit``, ``cost``, ``evals``, ``R``, ``result``. When
    ``G1`` fails, ``G2`` runs from the rescue ends in ``G1``'s order; an accepted end is returned
    as the fallback's result; if none, the grid call is made unchanged (same arguments), after
    which the counters hold the grid evaluations plus those of ``G1`` and ``G2``. ``G2`` draws
    nothing from ``rng`` (the aligned-draw book); its CMA fits have their own generators.
    One row per call is appended to ``log``.
    """

    def flat(t: Dict[str, Any]) -> Dict[str, Any]:
        return {
            a: (b.ravel().tolist() if hasattr(b, "ravel") else b)
            for a, b in t.items()
            if a not in ("result", "ends")
        }

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
                    polished=False,
                    tries=[],
                    polish_tries=[],
                    rescue_evals=0,
                    polish_evals=0,
                    grid_evals=int(rec.last_eval_counts[2]),
                )
            )
            return out
        step["used"] = True
        kept = step["kept"]

        def one(k: int) -> Dict[str, Any]:
            return g1_end(kept[k], voxel_vertices, phase_index)

        r = g1_choose([e["cost"] for e in kept], one)
        rows = [flat(t) for t in r["tries"]]
        if r["rescued"]:
            res = r["tries"][r["taken"]]["result"]
            rec._last_global_evals = 0
            rec._last_local_evals = int(r["evals"])
            log.append(
                dict(
                    voxel=step["voxel"],
                    attempted=True,
                    rescued=True,
                    polished=False,
                    tries=rows,
                    polish_tries=[],
                    rescue_evals=int(r["evals"]),
                    polish_evals=0,
                    grid_evals=0,
                )
            )
            return res
        ends = [t["R"] for t in r["tries"]]  # G1's rescue ends, in G1's order

        def two(j: int) -> Dict[str, Any]:
            return g2_end(ends[j], voxel_vertices, phase_index)

        p = polish_choose(len(ends), two)
        prows = [flat(t) for t in p["tries"]]
        if p["accepted"]:
            res = p["tries"][p["taken"]]["result"]
            rec._last_global_evals = 0
            rec._last_local_evals = int(r["evals"]) + int(p["evals"])
            log.append(
                dict(
                    voxel=step["voxel"],
                    attempted=True,
                    rescued=False,
                    polished=True,
                    tries=rows,
                    polish_tries=prows,
                    rescue_evals=int(r["evals"]),
                    polish_evals=int(p["evals"]),
                    grid_evals=0,
                )
            )
            return res
        out = grid_fn(voxel_vertices=voxel_vertices, phase_index=phase_index, rng=rng, **kw)
        g = int(rec.last_eval_counts[2])
        rec._last_local_evals = int(rec._last_local_evals) + int(r["evals"]) + int(p["evals"])
        log.append(
            dict(
                voxel=step["voxel"],
                attempted=True,
                rescued=False,
                polished=False,
                tries=rows,
                polish_tries=prows,
                rescue_evals=int(r["evals"]),
                polish_evals=int(p["evals"]),
                grid_evals=g,
            )
        )
        return out

    return wrapped
