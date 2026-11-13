"""AR5 step 8 (task ar5h): pure readers of the seed-step logs (``e1_log.json.gz``).

No sample, cache or truth is imported here, so the tests can use these functions on synthetic
steps. ``order_stats.py`` and ``order_trace.py`` import them (review fix of 2026-11-12: the
functions that decide which end a step took were untested; they were moved here unchanged).

A log row (one seed step) holds ``voxel``, ``stop`` (the index, in the order of the refined
candidate positions, of the end that reached hit ratio 0.95, or -1) and ``fits``: rows
``[pos, j, kind, hit, cost, evals, R]`` with ``kind`` "C" (the CMA fit) or "M" (the MC refinement).
"""

from typing import Any, Callable, Dict, List, Tuple

TAU = 0.95  # the early stop of the seed step (fit_rules.TAU; repeated so this file stays pure)


def kept_ends(e: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The kept end of each refined candidate position (lower cost of its fits; CMA on a tie)."""
    by: Dict[int, Dict[str, Any]] = {}
    for pos, j, kind, hit, cost, ev, R in e["fits"]:
        d = dict(pos=pos, j=j, kind=kind, hit=hit, cost=cost, evals=ev, R=R)
        if pos not in by or (kind == "M" and cost < by[pos]["cost"]):
            by[pos] = d
    return [by[p] for p in sorted(by)]


def taken_end(e: Dict[str, Any]) -> Dict[str, Any]:
    """The end the step took: the kept end at ``stop``, else the lowest cost (first on a tie)."""
    ends = kept_ends(e)
    if e["stop"] >= 0:
        return ends[e["stop"]]
    costs = [x["cost"] for x in ends]
    return ends[costs.index(min(costs))]


def log_counts(log: List[Dict[str, Any]], order: str, tau: float = TAU) -> Dict[str, Any]:
    """Fits, evaluations and the validity counts of one arm's seed log (``order`` "A" or "Y")."""
    c = dict(
        seed_steps=len(log),
        steps_with_fits=0,
        cma_fits=0,
        mc_fits=0,
        cma_evals=0,
        mc_evals=0,
        steps_with_mc=0,
        bad_mc_after_cma_pass=0,
        bad_cma_after_cma_pass=0,
        steps_taken_cma_end=0,
        steps_taken_mc_end=0,
    )
    for e in log:
        if not e["fits"]:
            continue
        c["steps_with_fits"] += 1
        cs = [f for f in e["fits"] if f[2] == "C"]
        ms = [f for f in e["fits"] if f[2] == "M"]
        c["cma_fits"] += len(cs)
        c["mc_fits"] += len(ms)
        c["cma_evals"] += sum(f[5] for f in cs)
        c["mc_evals"] += sum(f[5] for f in ms)
        c["steps_with_mc"] += int(bool(ms))
        passed = [f[0] for f in cs if f[3] >= tau]
        if passed:  # a CMA end reached 0.95 at position p
            p = min(passed)
            if order == "Y":
                c["bad_mc_after_cma_pass"] += len(ms)  # Y: no MC at all
                c["bad_cma_after_cma_pass"] += sum(1 for f in cs if f[0] > p)
            else:
                c["bad_mc_after_cma_pass"] += sum(1 for f in ms if f[0] >= p)
                c["bad_cma_after_cma_pass"] += sum(1 for f in cs if f[0] > p)
        t = taken_end(e)
        c["steps_taken_mc_end" if t["kind"] == "M" else "steps_taken_cma_end"] += 1
    return c


def index(log: List[Dict[str, Any]]) -> Dict[Tuple[int, int], Dict[str, Any]]:
    """Seed steps keyed by (voxel, occurrence at that voxel), in log order."""
    seen: Dict[int, int] = {}
    out: Dict[Tuple[int, int], Dict[str, Any]] = {}
    for e in log:
        v = int(e["voxel"])
        k = seen.get(v, 0)
        seen[v] = k + 1
        out[(v, k)] = e
    return out


def compare_logs(
    la: List[Dict[str, Any]],
    lb: List[Dict[str, Any]],
    wrong: Callable[[int, List[float]], bool],
) -> Dict[str, Any]:
    """Matched seed steps of two arms: taken ends that differ (and how many are wrong) and the fits
    both ran at the same candidate and kind (compared exactly). ``wrong(voxel, R)`` scores an end.
    """
    ia, ib = index(la), index(lb)
    keys = sorted(set(ia) & set(ib))
    differs = wrong_a = wrong_b = 0
    fits_cmp = fits_diff = 0
    for k in keys:
        ea, eb = ia[k], ib[k]
        if not ea["fits"] or not eb["fits"]:
            continue
        fa = {(f[1], f[2]): f for f in ea["fits"]}
        fb = {(f[1], f[2]): f for f in eb["fits"]}
        for key in set(fa) & set(fb):
            fits_cmp += 1
            x, y = fa[key], fb[key]
            if (x[3], x[4], x[5], x[6]) != (y[3], y[4], y[5], y[6]):
                fits_diff += 1
        ta, tb = taken_end(ea), taken_end(eb)
        if (ta["j"], ta["kind"]) != (tb["j"], tb["kind"]):
            differs += 1
            wrong_a += int(wrong(k[0], ta["R"]))
            wrong_b += int(wrong(k[0], tb["R"]))
    return dict(
        steps_a=len(la),
        steps_b=len(lb),
        matched=len(keys),
        only_a=len(set(ia) - set(ib)),
        only_b=len(set(ib) - set(ia)),
        taken_differs=differs,
        differing_wrong_first=wrong_a,
        differing_wrong_second=wrong_b,
        shared_fits_compared=fits_cmp,
        shared_fits_different=fits_diff,
    )


def same_voxel_prefix(la: List[Dict[str, Any]], lb: List[Dict[str, Any]]) -> int:
    """How many seed steps, from the start, visit the same voxel in both logs."""
    n = 0
    for ea, eb in zip(la, lb):
        if int(ea["voxel"]) != int(eb["voxel"]):
            break
        n += 1
    return n


def first_taken_difference(la: List[Dict[str, Any]], lb: List[Dict[str, Any]]) -> int:
    """The 0-based index of the first step, in both logs at once, at the same voxel with fits in
    both, whose taken end (candidate, kind) differs; -1 if there is none."""
    for i, (ea, eb) in enumerate(zip(la, lb)):
        if int(ea["voxel"]) != int(eb["voxel"]) or not ea["fits"] or not eb["fits"]:
            continue
        ta, tb = taken_end(ea), taken_end(eb)
        if (ta["j"], ta["kind"]) != (tb["j"], tb["kind"]):
            return i
    return -1


def flip_origin(wa_before: bool, wy_before: bool, wa_after: bool, wy_after: bool) -> str:
    """Where a voxel's flip after the pass (wrong in exactly one arm) was made: "before the pass"
    if the same arm was already the only wrong one before the pass, else "in the pass"."""
    if wa_after == wy_after:
        raise ValueError("not a flip after the pass")
    same = wa_before != wy_before and wa_before == wa_after
    return "before the pass" if same else "in the pass"
