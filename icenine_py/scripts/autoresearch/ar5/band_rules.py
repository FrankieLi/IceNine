"""AR5 step 15 (task ar5o): rule ``PB(k)``, a two-sided stall check in the seed step's MC
refinement (pure; numpy and ``stallhit_rules`` only).

``PB(k)`` acts on one MC refinement of the seed step (the stock MC local optimization from one
candidate: the MC stage then the variance stage, one run), as ``PH(k)`` of step 12 does, with the
same counter. A **stall** is ``k`` strict-cost evaluations with no strict gain of the best cost,
counted over the whole refinement. At a stall there is one **stall check**: the hit ratio of the
running best (the first lowest-cost orientation) is evaluated once. The refinement stops if that
hit ratio is at least ``TAU`` (0.95, the seed step's stop value) or under ``LOW`` (0.80, the seed
test value); the check is then the final evaluation of the kept orientation (counted). Between the
two it goes on, the check counted as one evaluation, and the next check is ``k`` evaluations after
the next strict gain (a gain re-arms the counter; with no gain there is no second check). The
refinement has its own generator, so a stop moves no other random draw. ``PB(k)`` replaces
``PH(3000)`` (there is no 3,000 check under it). Two forms: ``PB(1000)`` and ``PB(500)``.

This file holds the constants, the decision (``decide``), the meter that records or applies the
rule (``BandMeter`` and ``patch_evaluate``) and the walk of one refinement from recorded checks
(``walk_events``, arithmetic only). In the recording replay the meter runs ``PH(3000)`` exactly as
arm ``J`` does (``low`` None) and, besides, records the hit ratio of the running best at every 500
and every 1000 stall (auxiliary counters, never stopping): ``events``.
"""

from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np
import stallhit_rules as SR

TAU = SR.TAU  # 0.95: the seed step's stop value
LOW = SR.SEED_TEST  # 0.80: the seed test value
FORMS: Dict[str, int] = {"PB1000": 1000, "PB500": 500}
RECORD_KS = (500, 1000)  # the stall sizes the recording replay records
PH_K = SR.K  # 3,000: the stall size of arm J's own check (rule PH(3000)), kept in the recording

RULE = (
    "Rule PB(k): in the seed step, the MC refinement of a candidate (MC then variance stage) has "
    "a stall check after k strict-cost evaluations with no strict gain of its best cost: the hit "
    "ratio of the running best (the first lowest-cost orientation) is evaluated once; at 0.95 or "
    "more or under 0.80 the refinement stops and the check is its final evaluation (counted); "
    "else it goes on (the check counted as one evaluation) and checks again k evaluations after "
    "the next gain. It replaces PH(3000). Forms PB(1000) and PB(500). CMA fits, the five "
    "candidates, the 0.95 stop, the seed test 0.80, rules G1 and G2b and the grid search are "
    "unchanged. Arm JB = arm J (step 14) with PB(k)."
)


def decide(hit: float, tau: float = TAU, low: Optional[float] = LOW) -> str:
    """``"high"`` (stop, the best is at ``tau`` or more), ``"low"`` (stop, under ``low``; ``low``
    None never stops low) or ``"go"``."""
    if float(hit) >= tau:
        return "high"
    if low is not None and float(hit) < low:
        return "low"
    return "go"


class BandMeter(SR.Meter):
    """``stallhit_rules.Meter`` with a low side and auxiliary recording counters.

    ``begin(k, low, aux)``: ``k`` the stall size of the rule that acts (0: only records), ``low``
    the low threshold (None: ``PH`` as step 12, stops high only), ``aux`` stall sizes that are
    recorded and never act. ``events[a]`` holds, per stall of size ``a``, ``[n, cost, hit]`` (loop
    evaluation number, the running best's cost and hit ratio); each auxiliary counter is disarmed
    after a check until the next strict gain, as the acting one is."""

    def __init__(self) -> None:
        super().__init__()
        self.low: Optional[float] = None
        self.aux: Dict[int, SR.StallHitCounter] = {}
        self.events: Dict[int, List[List[float]]] = {}
        self.stop_kind = ""

    def begin(self, k: int, low: Optional[float] = None, aux: Sequence[int] = ()) -> None:
        super().begin(k)
        self.low = low
        self.aux = {int(a): SR.StallHitCounter(int(a)) for a in aux}
        self.events = {int(a): [] for a in aux}
        self.stop_kind = ""


def patch_evaluate(
    orig: Callable[..., Any], meter: BandMeter, tau: float = TAU
) -> Callable[..., Any]:
    """``orig`` (a cost function's ``evaluate``) wrapped as ``stallhit_rules.patch_evaluate`` does,
    with the two-sided decision and the auxiliary recording. While ``meter.on`` every call is
    counted; the hit ratio of a check is one overlap evaluation made with the wrapper off (not
    counted by it; the live run adds it to the evaluation count). Changes no value and draws no
    random number; the orientation of the first lowest cost is copied."""

    def patched(self_: Any, orientation: Any, *args: Any, **kw: Any) -> Any:
        out = orig(self_, orientation, *args, **kw)
        if meter.on:
            c = float(out.cost)
            gain = c < meter.counter.best
            fired_aux = [a for a, ctr in meter.aux.items() if ctr.feed(c)]
            fire = meter.counter.feed(c)
            meter.n = meter.counter.n
            if gain:
                meter.best = c
                meter.trace.append([meter.n, c])
                meter.best_R = np.array(orientation, dtype=np.float64, copy=True)
            if fired_aux or fire:
                assert meter.hit_fn is not None and meter.best_R is not None
                meter.on = False  # the checks are not themselves counted by the wrapper
                hit = float(meter.hit_fn(meter.best_R))
                for a in fired_aux:
                    meter.events[a].append([meter.n, meter.best, hit])
                    meter.aux[a].go_on()
                if fire:
                    d = decide(hit, tau, meter.low)
                    meter.checks.append([meter.n, hit, d != "go"])
                    if d != "go":
                        meter.stop_kind = d
                        raise SR.StallHitStop()
                    meter.counter.go_on()
                    meter.continued += 1
                meter.on = True
        return out

    return patched


def walk_events(
    events: List[List[float]],
    n_base: int,
    tau: float = TAU,
    low: Optional[float] = LOW,
) -> Dict[str, Any]:
    """Where ``PB(k)`` stops a refinement, from its recorded stall events (``[n, cost, hit]``, in
    order, for that ``k``) and ``n_base``, the evaluations of the refinement as ``PB`` would run it
    to its own end with no check at all (arm J's count minus J's own ``PH`` checks that went on).

    Returns ``kind`` (``high``, ``low`` or ``none``), ``evals`` (a stop: ``n + 1`` plus one per
    earlier check that went on; no stop: ``n_base`` plus one per check that went on), ``cost`` (the
    kept cost: the stop event's best; None when it never stops), ``continued``, ``i`` (the stop
    event's index, or None)."""
    cont = 0
    for i, ev in enumerate(events):
        d = decide(ev[2], tau, low)
        if d != "go":
            return dict(
                kind=d,
                evals=int(ev[0]) + 1 + cont,
                cost=float(ev[1]),
                hit=float(ev[2]),
                continued=cont,
                i=i,
            )
        cont += 1
    return dict(kind="none", evals=int(n_base) + cont, cost=None, hit=None, continued=cont, i=None)


def cut_entry(w: Dict[str, Any], final_cost: float, tol: float = 1e-12) -> Dict[str, Any]:
    """``walk_events``'s answer as an entry of ``stallhit_rules.walk_step``'s ``cuts``: ``kind``
    ``moved`` when it stops before a later gain (its kept cost is above the refinement's final best
    cost, so the end differs from J's), ``after_last`` when it stops with J's own end, ``none``
    when it never stops (``evals`` then still differ from J's: ``n_base`` plus the checks)."""
    if w["kind"] == "none":
        return dict(kind="none", evals=int(w["evals"]), cost=float(final_cost))
    moved = float(w["cost"]) > float(final_cost) + tol
    return dict(
        kind="moved" if moved else "after_last",
        evals=int(w["evals"]),
        cost=float(w["cost"]),
        hit=float(w["hit"]),
    )
