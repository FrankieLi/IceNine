"""AR5 step 15 (task ar5o): the seed step with rule ``PB(k)`` (script layer; library unchanged).

``BandPick`` is step 12's ``StallHitPick`` (arm ``H``: interleaved order, 20 costed, 5 refined) with
the MC refinement's stall check replaced by ``band_rules`` (two-sided) and, in the recording mode,
auxiliary recording of the hit ratio of the running best at every 500 and 1000 stall. Three modes,
chosen by the dictionary ``BAND`` before a step is run:

* ``dict(k=3000, low=None, aux=(500, 1000))``: the recording replay. The acting rule is arm ``J``'s
  own ``PH(3000)`` (stops high only), so every fit equals the stored one; the auxiliary counters
  only record (``events`` in the trace of each refinement).
* ``dict(k=1000, low=0.80, aux=())`` or ``k=500``: ``PB(k)``, live (the changed paths and the
  whole-layer runs of arm ``JB``).
* ``k=0``: off (as arm ``G`` of step 9).

A refinement has its own generator (``default_rng(0)`` in the stock MC local optimization), so a
stop moves no random draw. The stall check's hit ratio is one overlap evaluation made outside the
wrapped cost function; when the refinement goes on it is counted as one evaluation, when it stops
the check is the final evaluation of the kept orientation (counted).
"""

from typing import Any, Dict, Optional

import band_common as BC  # noqa: E402,F401  (first: guards, cache, freeze, held3)

import mix_common as MC  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import rescue_run as RRUN  # noqa: E402
import band_rules as BR  # noqa: E402
import stallhit_rules as SR  # noqa: E402

BAND: Dict[str, Any] = dict(k=0, low=None, aux=())
STATS: Dict[str, int] = dict(
    mc_refinements=0,
    checks_total=0,
    checks_went_on=0,
    stops_high=0,
    stops_low=0,
    evals_in_stopped=0,
)


class BandTracer(BR.BandMeter, ORUN.Tracer):
    """``order_run.Tracer`` with the band rule (``band_rules.BandMeter`` and ``patch_evaluate``).
    Draws no random number and changes no value. Only one is made (it patches
    ``VoxelCostFunction.evaluate`` for the process)."""

    def __init__(self) -> None:
        from icenine.cost_functions import VoxelCostFunction

        BR.BandMeter.__init__(self)
        VoxelCostFunction.evaluate = BR.patch_evaluate(VoxelCostFunction.evaluate, self)


class BandPick(RRUN.RescuePick):
    """``RescuePick`` whose MC fits obey ``BAND`` (see the module text)."""

    ptracer: Optional[BandTracer] = None

    def _fit(
        self, kind: str, rec: Any, j: int, Rj: np.ndarray, cj: float, vertices: Any, phase: int
    ) -> Dict[str, Any]:
        if BandPick.ptracer is None:
            BandPick.ptracer = BandTracer()
            # the base ``_fit`` finds a tracer set (it reads the module's ``OrderPick`` name, which
            # ``main`` rebinds to this class; the replay does not rebind it)
            ORUN.OrderPick.tracer = BandPick.ptracer
            BandPick.tracer = BandPick.ptracer
        if kind != "M":
            return super()._fit(kind, rec, j, Rj, cj, vertices, phase)
        tr = BandPick.ptracer
        tr.begin(int(BAND["k"]), BAND["low"], tuple(BAND["aux"]))

        def hit_fn(R: np.ndarray) -> float:
            info = rec.evaluate_overlap(np.asarray(R, np.float64), vertices, phase)
            return float(self.bfs._ratios(info)[0])

        tr.hit_fn = hit_fn
        cut = False
        try:
            Rend, c, ev = MC._orig_refine_arm("M", rec, self.orig_lo, Rj, cj, vertices, phase)
            ev = int(ev) + tr.continued  # each stall check that went on is one counted evaluation
        except SR.StallHitStop:
            cut = True
            assert tr.best_R is not None
            Rend = tr.best_R
            # the stall check is the final evaluation of the kept orientation (counted)
            final = rec._make_local_cost_fn().evaluate(
                orientation=Rend, voxel_vertices=vertices, phase_index=phase
            )
            c = float(final.cost)
            ev = tr.n + 1 + tr.continued
        tr.on = False
        info = rec.evaluate_overlap(np.asarray(Rend, np.float64), vertices, phase)
        hit = float(self.bfs._ratios(info)[0])
        ORUN.COUNTS["mc_fits"] += 1
        STATS["mc_refinements"] += 1
        STATS["checks_total"] += len(tr.checks)
        STATS["checks_went_on"] += tr.continued
        if cut:
            STATS["stops_high" if tr.stop_kind == "high" else "stops_low"] += 1
            STATS["evals_in_stopped"] += int(ev)
        ORUN.TRACES.append(
            dict(
                voxel=self.cur,
                step=ORUN.COUNTS["seed_steps"],
                j=int(j),
                evals=int(ev),
                n_calls=int(ev),
                trace=tr.trace,
                cut=cut,
                stop_kind=tr.stop_kind,
                checks=[list(x) for x in tr.checks],
                continued=int(tr.continued),
                events={str(a): [list(x) for x in v] for a, v in tr.events.items()},
            )
        )
        return dict(j=int(j), kind="M", hit=hit, cost=float(c), evals=int(ev), R=np.asarray(Rend))


def summary() -> Dict[str, Any]:
    return dict(k=BAND["k"], low=BAND["low"], aux=list(BAND["aux"]), rule=BR.RULE, **STATS)
