"""AR5 step 12 (task ar5l): rule ``PH(k)``, stop a stalled seed-step MC refinement only when its
best already reaches hit ratio 0.95 (pure; no study imports).

``PH(k)`` acts on one MC refinement of the seed step (the stock MC local optimization from one
candidate: the MC stage then the variance stage, one run). A **stall** is ``k`` strict-cost
evaluations with no strict gain of the best cost, counted over the whole refinement. At a stall
there is one **stall check**: the hit ratio of the running best (the first lowest-cost orientation)
is evaluated once. If it is at least ``TAU`` (0.95) the refinement ends and the check is its final
evaluation of the kept orientation (counted, as the uncut refinement's final evaluation). Else the
refinement goes on, the check is counted as one evaluation, and the next check is ``k``
evaluations after the next strict gain (a gain re-arms the counter; with no gain there is no second
check). The refinement has its own generator, so a stop moves no other random draw.

This file holds the constants, the counter (``StallHitCounter``), the live wrapper (``Meter`` and
``patch_evaluate``), the stall events of a stored trace (``stall_events``), the walk of one stored
refinement with known hit ratios (``refine_walk``) and the walk of a whole seed step (``walk_step``,
arithmetic only). A trace is the list ``[n, best cost]`` at each strict gain, ``n`` counting
evaluations from 1, as ``order_run.Tracer`` records it.
"""

from typing import Any, Callable, Dict, List, Optional

import numpy as np

K = 3000  # frozen: a stall is this many evaluations with no strict gain
TAU = 0.95  # the hit ratio that ends a seed step, and the stall check's threshold
SEED_TEST = 0.80  # the seed test (unchanged)

RULE = (
    "Rule PH(3000): in the seed step, the MC refinement of a candidate (MC then variance stage) "
    "has a stall check after 3,000 strict-cost evaluations with no strict gain of its best cost: "
    "the hit ratio of the running best (the first lowest-cost orientation) is evaluated once; at "
    "0.95 or more the refinement stops and the check is its final evaluation; else it goes on (the "
    "check counted as one evaluation) and checks again 3,000 evaluations after the next gain. CMA "
    "fits, the 0.95 stop, the five candidates, the seed test 0.80, rule G1 and the grid search are "
    "unchanged. Arm H = arm G (step 9) with PH(3000)."
)


class StallHitCounter:
    """Counts the evaluations of one refinement. ``feed(cost)`` returns True when a stall check is
    due (``k`` evaluations since the last strict gain, counter armed). ``go_on()`` (the check did
    not stop the refinement) disarms it until the next strict gain. Draws no random number."""

    def __init__(self, k: int) -> None:
        self.k = int(k)
        self.n = 0  # evaluations seen (stall checks that went on are not in it)
        self.best = float("inf")
        self.best_n = 0  # evaluation number of the first lowest cost (0: none yet)
        self.since = 0  # evaluations since the last strict gain
        self.armed = True

    def feed(self, cost: float) -> bool:
        self.n += 1
        if cost < self.best:
            self.best = float(cost)
            self.best_n = self.n
            self.since = 0
            self.armed = True
        else:
            self.since += 1
        return self.k > 0 and self.armed and self.since >= self.k

    def go_on(self) -> None:
        self.armed = False


class StallHitStop(Exception):
    """Raised by the wrapped cost function when a stall check stops the refinement."""


class Meter:
    """State of the wrapped cost function over one refinement: ``on``, evaluations ``n``, the
    running best cost with its trace (``[n, cost]`` at each strict gain), ``best_R`` (the first
    lowest-cost orientation), ``checks`` (one ``[n, hit, stopped]`` per stall check: the hit ratio
    at every check is stored) and ``continued`` (checks that went on, each one counted evaluation).
    ``hit_fn(R)`` gives the hit ratio of an orientation (an overlap evaluation, not through the
    wrapped cost function). ``begin(k)`` starts a refinement (``k`` 0: only records)."""

    def __init__(self) -> None:
        self.on = False
        self.n = 0
        self.best = float("inf")
        self.trace: List[List[float]] = []
        self.k = 0
        self.counter = StallHitCounter(0)
        self.best_R: Any = None
        self.checks: List[List[Any]] = []
        self.continued = 0
        self.hit_fn: Optional[Callable[[Any], float]] = None

    def begin(self, k: int) -> None:
        self.on, self.n, self.best, self.trace = True, 0, float("inf"), []
        self.k = int(k)
        self.counter = StallHitCounter(self.k)
        self.best_R = None
        self.checks = []
        self.continued = 0


def patch_evaluate(orig: Callable[..., Any], meter: Meter, tau: float = TAU) -> Callable[..., Any]:
    """``orig`` (a cost function's ``evaluate(self, orientation, ...)``) wrapped: while
    ``meter.on`` every call is counted. At a stall the hit ratio of the kept best is evaluated once
    (``meter.hit_fn``, the wrapper off meanwhile); at ``tau`` or more the meter switches itself off
    and ``StallHitStop`` is raised after the evaluation; else the counter is disarmed and the
    refinement goes on. Changes no value, draws no random number; the orientation of the first
    lowest cost is copied."""

    def patched(self_: Any, orientation: Any, *args: Any, **kw: Any) -> Any:
        out = orig(self_, orientation, *args, **kw)
        if meter.on:
            c = float(out.cost)
            gain = c < meter.counter.best
            fire = meter.counter.feed(c)
            meter.n = meter.counter.n
            if gain:
                meter.best = c
                meter.trace.append([meter.n, c])
                meter.best_R = np.array(orientation, dtype=np.float64, copy=True)
            if fire:
                assert meter.hit_fn is not None and meter.best_R is not None
                meter.on = False  # the check is not itself counted by the wrapper
                hit = float(meter.hit_fn(meter.best_R))
                stop = hit >= tau
                meter.checks.append([meter.n, hit, bool(stop)])
                if stop:
                    raise StallHitStop()
                meter.on = True
                meter.counter.go_on()
                meter.continued += 1
        return out

    return patched


def stall_events(trace: List[List[float]], n_total: int, k: int) -> List[Dict[str, Any]]:
    """The stall checks an uncut refinement would reach, in order, from its stored trace and its
    evaluation count ``n_total`` (final evaluation included). One event per strict gain whose gap to
    the next gain exceeds ``k`` (``last`` False; the check is at evaluation ``n + k``) and, for the
    last gain, one if the check falls inside the refinement (``last`` True: ``n_total > n + k``, the
    final evaluation not counted as a loop evaluation). Each
    event: ``i`` (index in the trace), ``n`` (evaluation of the gain), ``at`` (``n + k``), ``cost``
    (best cost, the kept cost if it stops) and ``last``. Step 11's ``cut_point`` is the first event
    without the hit condition."""
    out: List[Dict[str, Any]] = []
    if k <= 0 or not trace:
        return out
    for i, (n, c) in enumerate(trace):
        n = int(n)
        if i + 1 < len(trace):
            if int(trace[i + 1][0]) - n > k:
                out.append(dict(i=i, n=n, at=n + k, cost=float(c), last=False))
        elif int(n_total) > n + k:
            out.append(dict(i=i, n=n, at=n + k, cost=float(c), last=True))
    return out


def refine_walk(
    trace: List[List[float]],
    n_total: int,
    k: int,
    hit_of: Callable[[int], Optional[float]],
    tau: float = TAU,
    unknown: str = "go",
) -> Dict[str, Any]:
    """Where ``PH(k)`` stops a stored refinement, given the hit ratio of the best at each gain.

    ``hit_of(i)`` is the hit ratio of the best after the ``i``-th gain, or None when unknown.
    An unknown check is resolved by ``unknown``: ``"go"`` (the refinement goes on: the plain
    reading, no saving, no change) or ``"stop"`` (it stops there); the walk lists them in
    ``unknown_events``. Returns ``kind`` (``none``: never stopped; ``after_last``: stopped after
    its last gain, the end unchanged; ``moved``: stopped before a later gain), ``evals`` (stop point
    plus the final evaluation, or ``n_total`` plus one per check that went on; never more than
    ``n_total`` + continued), ``cost`` (the kept cost), ``continued`` (checks that went on),
    ``stop_i`` (index of the gain whose best is kept, or None) and ``unknown_events``."""
    final = float(trace[-1][1]) if trace else float("inf")
    cont = 0
    unk: List[int] = []
    for ev in stall_events(trace, n_total, k):
        h = hit_of(int(ev["i"]))
        if h is None:
            unk.append(int(ev["i"]))
            stop = unknown == "stop"
        else:
            stop = float(h) >= tau
        if stop:
            return dict(
                kind="after_last" if ev["last"] else "moved",
                evals=min(int(n_total), int(ev["n"]) + k + 1) + cont,
                cost=float(ev["cost"]),
                continued=cont,
                stop_i=int(ev["i"]),
                unknown_events=unk,
            )
        cont += 1
    return dict(
        kind="none",
        evals=int(n_total) + cont,
        cost=final,
        continued=cont,
        stop_i=None,
        unknown_events=unk,
    )


def walk_step(
    fits: List[Dict[str, Any]],
    cuts: Dict[int, Dict[str, Any]],
    tau: float = TAU,
    n_avail: int = 5,
) -> Dict[str, Any]:
    """The seed step (interleaved order) with cuts, from its stored fits, by arithmetic.

    ``fits`` are the step's stored fits in run order: dicts ``pos``, ``kind`` (``C`` or ``M``),
    ``hit``, ``cost``, ``evals``. ``cuts[pos]`` is ``refine_walk``'s answer for the MC fit at that
    position, with ``hit`` (the hit ratio of the kept end when ``kind`` is ``moved``; a stop under
    this rule happens only at ``tau`` or more, so it defaults to the stored hit class: at least
    ``tau``). A candidate's kept end is its CMA end, or its MC end if strictly cheaper. The step
    stops at the first kept end with hit ratio at least ``tau``; else the lowest-cost kept end is
    taken (first of equals). Returns ``evals`` (the walked step's evaluations, with ``extra`` for
    candidates the stored step never ran, priced at the step's mean per-candidate cost), ``taken``
    (position), ``kept`` (per position: ``cost``, ``hit``, ``which``, ``moved``), ``stopped``,
    ``lost`` (a stored stop that no longer happens), ``extra_candidates`` and ``extra_evals``."""
    n_pos = 1 + max(int(f["pos"]) for f in fits)
    by: Dict[int, Dict[str, Dict[str, Any]]] = {p: {} for p in range(n_pos)}
    for f in fits:
        by[int(f["pos"])][f["kind"]] = f
    kept: List[Dict[str, Any]] = []
    evals = 0
    stopped_at: Optional[int] = None
    stored_stop: Optional[int] = None
    for p in range(n_pos):
        c = by[p]["C"]
        e = dict(cost=float(c["cost"]), hit=float(c["hit"]), which="C", moved=False)
        evals += int(c["evals"])
        stored_kept = e
        if "M" in by[p]:
            m = by[p]["M"]
            cut = cuts.get(p)
            evals += int(m["evals"]) if cut is None else int(cut["evals"])
            m_cost, m_hit, moved = float(m["cost"]), float(m["hit"]), False
            if cut is not None and cut["kind"] == "moved":
                moved = True
                m_cost = float(cut["cost"])
                m_hit = float(cut.get("hit", tau))
            if float(m["cost"]) < float(c["cost"]):
                stored_kept = dict(cost=float(m["cost"]), hit=float(m["hit"]), which="M")
            if m_cost < float(c["cost"]):
                e = dict(cost=m_cost, hit=m_hit, which="M", moved=moved)
        kept.append(e)
        if stored_stop is None and stored_kept["hit"] >= tau:
            stored_stop = p
        if e["hit"] >= tau:
            stopped_at = p
            break
    lost = stored_stop is not None and stopped_at is None
    extra = extra_evals = 0
    if lost:
        per = float(sum(int(f["evals"]) for f in fits)) / max(1, n_pos)
        extra = max(0, n_avail - n_pos)
        extra_evals = int(round(per * extra))
        evals += extra_evals
    taken = (
        stopped_at
        if stopped_at is not None
        else min(range(len(kept)), key=lambda i: (kept[i]["cost"], i))
    )
    return dict(
        evals=int(evals),
        taken=int(taken),
        kept=kept,
        stopped=stopped_at is not None,
        lost=bool(lost),
        extra_candidates=int(extra),
        extra_evals=int(extra_evals),
    )
