"""AR5 step 11 (task ar5k): rule ``P(k)``, patience of the seed step's MC refinement (pure).

``P(k)`` stops one MC refinement of the seed step (the stock MC local optimization from one
candidate: the MC stage, then the variance stage, one run) after ``k`` strict-cost evaluations with
no strict gain of the best cost, counted over the whole refinement. It keeps the first
lowest-cost orientation seen (a later evaluation of equal cost does not replace it) and makes one
final evaluation of that orientation, as the uncut refinement does (the library ends a refinement
with one overlap evaluation of its kept orientation, counted). The refinement has its own
generator (``default_rng(0)``), so a stop moves no other random draw.

This file holds the constants, the counter the live wrapper uses (``PatienceCounter``), the stop
point on a stored trace (``cut_point``: a trace is the list ``[n, best cost]`` at each strict
gain, ``n`` counting evaluations from 1, as ``order_run.Tracer`` records it) and the walk of a
whole seed step with the cut (``walk_step``, arithmetic only). No study imports.
"""

from typing import Any, Dict, List, Optional

K = 3000  # frozen: stop after this many evaluations with no strict gain
TAU = 0.95  # the hit ratio that ends a seed step (unchanged)
SEED_TEST = 0.80  # the seed test (unchanged)

RULE = (
    "Rule P(3000): in the seed step, the MC refinement of a candidate (MC then variance stage) "
    "stops after 3,000 strict-cost evaluations with no strict gain of its best cost, keeps the "
    "first lowest-cost orientation, and makes one final evaluation of it. CMA fits, the 0.95 stop, "
    "the "
    "five candidates, the seed test 0.80, rule G1 and the grid search are unchanged. Arm S = arm G "
    "(step 9) with P(3000)."
)


class PatienceCounter:
    """Counts evaluations of one refinement. ``feed(cost)`` returns True when the refinement must
    end now (``k`` evaluations since the last strict gain). Draws no random number."""

    def __init__(self, k: int) -> None:
        self.k = int(k)
        self.n = 0  # evaluations seen
        self.best = float("inf")
        self.best_n = 0  # evaluation number of the first lowest cost (0: none yet)
        self.since = 0  # evaluations since the last strict gain

    def feed(self, cost: float) -> bool:
        self.n += 1
        if cost < self.best:
            self.best = float(cost)
            self.best_n = self.n
            self.since = 0
        else:
            self.since += 1
        return self.k > 0 and self.since >= self.k


def run_costs(costs: List[float], k: int) -> Dict[str, Any]:
    """Drive a counter with the cost sequence of one refinement (the last cost is the final
    evaluation of the kept orientation). Returns ``stopped`` (the counter fired before the sequence
    ended), ``best_n`` (1-based evaluation of the kept orientation), ``evals`` (evaluations
    including the final one) and ``cost`` (the kept cost)."""
    c = PatienceCounter(k)
    for i, x in enumerate(costs):
        if c.feed(x) and i + 1 + 1 < len(costs):
            return dict(stopped=True, best_n=c.best_n, evals=c.n + 1, cost=c.best)
    return dict(stopped=False, best_n=c.best_n, evals=len(costs), cost=c.best)


def cut_point(
    trace: List[List[float]], n_total: int, k: int, worse_than: Optional[float] = None
) -> Dict[str, Any]:
    """Where ``P(k)`` stops a stored refinement. ``trace`` is ``[[n, best cost], ...]`` (one row
    per strict gain), ``n_total`` the evaluations of the uncut refinement (final one included).
    Returns ``kind`` ``none`` (never stopped earlier), ``after_last`` (stopped after its last
    improvement: the end is unchanged) or ``moved`` (stopped before a later improvement: the end is
    the best at the stop), ``evals`` (new count: stop point plus the final evaluation, never more
    than ``n_total``) and ``cost`` (the kept cost; the stored end cost unless ``moved``).
    ``worse_than`` (what-if tables only, not the rule): stop only while the best cost is above it.
    """
    final = float(trace[-1][1])
    if k <= 0 or not trace:
        return dict(kind="none", evals=int(n_total), cost=final)
    for i, (n, c) in enumerate(trace):
        n = int(n)
        if worse_than is not None and not float(c) > worse_than:
            continue
        if i + 1 < len(trace):
            if int(trace[i + 1][0]) - n > k:
                return dict(kind="moved", evals=min(int(n_total), n + k + 1), cost=float(c))
        elif int(n_total) > n + k and n + k + 1 < int(n_total):
            return dict(kind="after_last", evals=n + k + 1, cost=float(c))
    return dict(kind="none", evals=int(n_total), cost=final)


def walk_step(
    fits: List[Dict[str, Any]],
    cuts: Dict[int, Dict[str, Any]],
    tau: float = TAU,
    pessimistic: bool = False,
    n_avail: int = 5,
) -> Dict[str, Any]:
    """The seed step (interleaved order) with cuts, from its stored fits, by arithmetic.

    ``fits`` are the step's stored fits in run order: dicts ``pos``, ``kind`` (``C`` or ``M``),
    ``hit``, ``cost``, ``evals``. ``cuts[pos]`` is the ``cut_point`` of the MC fit at that position.
    A candidate's kept end is its CMA end, or its MC end if strictly cheaper. The step stops at the
    first kept end with hit ratio at least ``tau``; else the lowest-cost kept end is taken (first of
    equals). A cut that moved an MC end leaves its hit ratio unknown: optimistic (default) keeps
    the stored hit class; ``pessimistic`` treats an end that reached ``tau`` as lost. When the MC
    end becomes dearer than the CMA end, the CMA end (hit known, under ``tau``) is kept.
    Returns ``evals`` (finisher evaluations of the walked step, with ``extra`` for candidates the
    stored step never ran, estimated as the step's mean per-candidate cost), ``taken`` (position),
    ``kept`` (per position: ``cost``, ``hit``, ``which`` C or M, ``moved``), ``stopped``, ``lost``
    (a stored stop that no longer happens) and ``extra_candidates``."""
    n_pos = 1 + max(int(f["pos"]) for f in fits)
    by: Dict[int, Dict[str, Dict[str, Any]]] = {p: {} for p in range(n_pos)}
    for f in fits:
        by[int(f["pos"])][f["kind"]] = f
    kept: List[Dict[str, Any]] = []
    evals = 0
    stopped_at: Optional[int] = None
    lost = False
    stored_stop = None
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
                if pessimistic and m_hit >= tau:
                    m_hit = 0.0
            if float(m["cost"]) < float(c["cost"]):
                stored_kept = dict(cost=float(m["cost"]), hit=float(m["hit"]), which="M")
            if m_cost < float(c["cost"]):
                e = dict(cost=m_cost, hit=m_hit, which="M", moved=moved)
            # else: the CMA end stays (its hit ratio is under tau, since the MC fit ran)
        kept.append(e)
        if stored_stop is None and stored_kept["hit"] >= tau:
            stored_stop = p
        if e["hit"] >= tau:
            stopped_at = p
            break
    stored_stopped = stored_stop is not None
    # a stored stop that no longer happens: the step goes on to the candidates it never ran
    extra = 0
    extra_evals = 0
    if stored_stopped and stopped_at is None:
        lost = True
        per = float(sum(int(f["evals"]) for f in fits)) / max(1, n_pos)
        # the candidates the stored step never ran (up to ``n_avail``) have no stored fit
        extra = max(0, n_avail - n_pos)
        extra_evals = int(round(per * extra))
        evals += extra_evals
    if stopped_at is not None:
        taken = stopped_at
    else:
        taken = min(range(len(kept)), key=lambda i: (kept[i]["cost"], i))
    return dict(
        evals=int(evals),
        taken=int(taken),
        kept=kept,
        stopped=stopped_at is not None,
        lost=lost,
        extra_candidates=int(extra),
        extra_evals=int(extra_evals),
    )
