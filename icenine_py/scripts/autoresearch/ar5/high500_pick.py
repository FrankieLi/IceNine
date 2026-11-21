"""AR5 step 16 (task ar5p): the seed step with rule ``PH(500)`` (script layer; library unchanged).

``High500Pick`` is step 15's ``BandPick`` (itself step 12's ``StallHitPick``: arm ``H``,
interleaved order, 20 costed, 5 refined) with the MC refinement's stall check taken from
``band_rules`` with NO low side (``low`` None always; the rule is ``high500_rules``). Modes, chosen
by the dictionary ``H5`` before a step is run:

* ``dict(k=500, aux=())``: ``PH(500)``, live (the replay of the moved steps and the whole-layer
  runs of arm ``J5``).
* ``dict(k=3000, aux=())``: arm ``J``'s own ``PH(3000)`` (every fit equals the stored one).
* ``k=0``: off (as arm ``G`` of step 9).

A refinement has its own generator (``default_rng(0)`` in the stock MC local optimization), so a
stop moves no random draw. The stall check's hit ratio is one overlap evaluation made outside the
wrapped cost function; when the refinement goes on it is counted as one evaluation, when it stops
the check is the final evaluation of the kept orientation (counted).
"""

from typing import Any, Dict, Optional

import high500_common as HC  # noqa: E402,F401  (first: guards, cache, freeze, held3)

import mix_common as MC  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import rescue_run as RRUN  # noqa: E402
import band_rules as BR  # noqa: E402
import high500_rules as HR  # noqa: E402
import high500_rules as HR  # noqa: E402
import stallhit_rules as SR  # noqa: E402

H5: Dict[str, Any] = dict(k=0, aux=())
STATS: Dict[str, int] = dict(
    mc_refinements=0,
    checks_total=0,
    checks_went_on=0,
    stops_high=0,
    stops_low=0,
    evals_in_stopped=0,
)


class High500Tracer(BR.BandMeter, ORUN.Tracer):
    """``order_run.Tracer`` with the high-stop rule (``band_rules.BandMeter``, no low side).
    Draws no random number and changes no value. Only one is made (it patches
    ``VoxelCostFunction.evaluate`` for the process)."""

    def __init__(self) -> None:
        from icenine.cost_functions import VoxelCostFunction

        BR.BandMeter.__init__(self)
        VoxelCostFunction.evaluate = BR.patch_evaluate(VoxelCostFunction.evaluate, self)


class High500Pick(RRUN.RescuePick):
    """``RescuePick`` whose MC fits obey ``H5`` (see the module text)."""

    ptracer: Optional[High500Tracer] = None

    def _fit(
        self, kind: str, rec: Any, j: int, Rj: np.ndarray, cj: float, vertices: Any, phase: int
    ) -> Dict[str, Any]:
        if High500Pick.ptracer is None:
            High500Pick.ptracer = High500Tracer()
            # the base ``_fit`` finds a tracer set (it reads the module's ``OrderPick`` name, which
            # ``main`` rebinds to this class; the replay does not rebind it)
            ORUN.OrderPick.tracer = High500Pick.ptracer
            High500Pick.tracer = High500Pick.ptracer
        if kind != "M":
            return super()._fit(kind, rec, j, Rj, cj, vertices, phase)
        tr = High500Pick.ptracer
        tr.begin(int(H5["k"]), None, tuple(H5["aux"]))

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
    return dict(k=H5["k"], low=None, aux=list(H5["aux"]), rule=HR.RULE, **STATS)
