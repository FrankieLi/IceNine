"""AR5 step 11, phase A step 5: one whole-layer run of arm ``S`` (script layer; library unchanged).

  RESEED_SAMPLE=<big|grad|mid|held> uv run python scripts/autoresearch/ar5/patience_run.py \
      [--cap-min 240] [--smoke N] [--guards-only]

Arm ``S`` is arm ``G`` of step 9 (``rescue_run.py``: arm ``A`` of step 8 with rule ``G1`` in front
of the grid fallback, interleaved order, aligned draws) with rule ``P(3000)`` (``patience_rules``):
every MC refinement of a seed step stops after 3,000 strict-cost evaluations with no strict gain,
keeps the first lowest-cost orientation and makes one final evaluation of it (counted). It is
``rescue_run.py`` with the run named ``S``, the pick subclassed (``PatiencePick``) and the cost
function's ``evaluate`` wrapped (``PatienceTracer``, a drop-in for ``order_run.Tracer`` that also
ends a refinement by an exception caught in ``_fit``). A refinement has its own generator
(``default_rng(0)``), so a stop moves no random draw. Besides the files of ``rescue_run.py`` it
writes ``patience.json`` (its counts). The boundary pass is made after by ``order_pass.py --name
S``. ``held`` is refused until ``ar5/patience/frozen.json`` is committed and clean and
``patience_freeze`` agrees (``patience_common``), which replaces step 9's and step 8's freeze
checks (step 9's phase B changed the library after those freezes). ``--k N`` (default 3000; the
replay uses 0 as a check) is refused on ``held`` unless it is 3000.
"""

import json
import sys
from typing import Any, Dict, List, Optional

import patience_common as PC  # noqa: E402  (first: sets the cache, guards held)

import fit_common as FC  # noqa: E402
import mix_common as MC  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import patience_rules as PR  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_run as RRUN  # noqa: E402


PatienceStop = PR.PatienceStop


class PatienceTracer(PR.Meter, ORUN.Tracer):
    """``order_run.Tracer`` (best cost at each improvement of the strict cost calls while ``on``)
    with the patience counter (``patience_rules.Meter`` and ``patch_evaluate``): when ``k`` > 0 and
    ``k`` evaluations pass with no strict gain it switches itself off and raises ``PatienceStop``.
    It keeps the first lowest-cost orientation (``best_R``). Draws no random number and changes no
    value. Only one is made (it patches ``VoxelCostFunction.evaluate`` for the process)."""

    def __init__(self) -> None:
        from icenine.cost_functions import VoxelCostFunction

        PR.Meter.__init__(self)
        VoxelCostFunction.evaluate = PR.patch_evaluate(VoxelCostFunction.evaluate, self)


PATIENCE: Dict[str, Any] = {"k": PR.K}
STATS: Dict[str, int] = dict(mc_refinements=0, cut=0, evals_in_cut=0)


class PatiencePick(RRUN.RescuePick):
    """``RescuePick`` whose MC fits obey rule ``P(k)`` (``PATIENCE['k']``; 0: off, as arm G)."""

    ptracer: Optional[PatienceTracer] = None

    def _fit(
        self, kind: str, rec: Any, j: int, Rj: np.ndarray, cj: float, vertices: Any, phase: int
    ) -> Dict[str, Any]:
        if PatiencePick.ptracer is None:
            PatiencePick.ptracer = PatienceTracer()
            # the base ``_fit`` finds a tracer set (it reads the module's ``OrderPick`` name, which
            # ``main`` rebinds to this class; the replay does not rebind it)
            ORUN.OrderPick.tracer = PatiencePick.ptracer
            PatiencePick.tracer = PatiencePick.ptracer
        if kind != "M":
            return super()._fit(kind, rec, j, Rj, cj, vertices, phase)
        tr = PatiencePick.ptracer
        tr.begin(int(PATIENCE["k"]))
        cut = False
        try:
            Rend, c, ev = MC._orig_refine_arm("M", rec, self.orig_lo, Rj, cj, vertices, phase)
        except PatienceStop:
            cut = True
            assert tr.best_R is not None
            Rend = tr.best_R
            # the final evaluation of the kept orientation (counted, as the uncut refinement's)
            final = rec._make_local_cost_fn().evaluate(
                orientation=Rend, voxel_vertices=vertices, phase_index=phase
            )
            c = float(final.cost)
            ev = tr.n + 1
        tr.on = False
        info = rec.evaluate_overlap(np.asarray(Rend, np.float64), vertices, phase)
        hit = float(self.bfs._ratios(info)[0])
        ORUN.COUNTS["mc_fits"] += 1
        STATS["mc_refinements"] += 1
        if cut:
            STATS["cut"] += 1
            STATS["evals_in_cut"] += int(ev)
        ORUN.TRACES.append(
            dict(
                voxel=self.cur,
                step=ORUN.COUNTS["seed_steps"],
                j=int(j),
                evals=int(ev),
                n_calls=int(ev),
                trace=tr.trace,
                cut=cut,
            )
        )
        return dict(j=int(j), kind="M", hit=hit, cost=float(c), evals=int(ev), R=np.asarray(Rend))


def summary() -> Dict[str, Any]:
    return dict(k=PATIENCE["k"], rule=PR.RULE, **STATS)


def main() -> None:
    args = sys.argv[1:]
    k = PR.K
    if "--k" in args:
        i = args.index("--k")
        k = int(args[i + 1])
        del args[i : i + 2]
    if FC.SAMPLE == "held" and k != PR.K:
        raise SystemExit("patience rule: held runs only with the frozen k")
    PATIENCE["k"] = k
    sys.argv = ["order_run.py", "--arm", "A", "--name", "S"] + args
    ORUN.OrderPick = PatiencePick  # type: ignore[misc]
    ORUN.OrderWrap = RRUN.RescueWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"S_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(RRUN.LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(RRUN.summary(), indent=1) + "\n")
    (out / "patience.json").write_text(json.dumps(summary(), indent=1) + "\n")


if __name__ == "__main__":
    main()
