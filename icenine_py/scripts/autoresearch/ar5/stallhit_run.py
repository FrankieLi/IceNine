"""AR5 step 12, phase A step 6 (run): one whole-layer run of arm ``H`` (script layer; library
unchanged).

  RESEED_SAMPLE=<big|grad|mid|held|held2> uv run python scripts/autoresearch/ar5/stallhit_run.py \
      [--cap-min 300] [--smoke N] [--guards-only] [--k 3000]

Arm ``H`` is arm ``G`` of step 9 (``rescue_run.py``: arm ``A`` of step 8 with rule ``G1`` in front
of the grid fallback, interleaved order, aligned draws) with rule ``PH(3000)`` (``stallhit_rules``):
in every MC refinement of a seed step, a stall (3,000 strict-cost evaluations with no strict gain)
has one stall check, the hit ratio of the running best (the first lowest-cost orientation); at 0.95
or more the refinement stops and the check is its final evaluation (counted), else it goes on (the
check counted as one evaluation) and checks again 3,000 evaluations after the next gain. It is
``rescue_run.py`` with the run named ``H``, the pick subclassed (``StallHitPick``) and the cost
function's ``evaluate`` wrapped (``StallHitTracer``). A refinement has its own generator
(``default_rng(0)``), so a stop moves no random draw. Besides the files of ``rescue_run.py`` it
writes ``stallhit.json`` (its counts) and puts the hit ratio of every stall check into the traces.
The boundary pass is made after by ``order_pass.py --name H``. ``held`` and ``held2`` are refused
until ``ar5/stallhit/frozen.json`` is committed and clean and ``stallhit_freeze`` agrees
(``stallhit_common``). ``--k N`` (default 3000; the replay uses 0 as a check) is refused on both
held samples unless it is 3000.
"""

import json
import sys
from typing import Any, Dict, Optional

import stallhit_common as SC  # noqa: E402,F401  (first: sets the cache, guards held and held2)

import fit_common as FC  # noqa: E402
import mix_common as MC  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_run as RRUN  # noqa: E402
import stallhit_rules as SR  # noqa: E402

StallHitStop = SR.StallHitStop


class StallHitTracer(SR.Meter, ORUN.Tracer):
    """``order_run.Tracer`` with the stall-hit rule (``stallhit_rules.Meter`` and
    ``patch_evaluate``). Draws no random number and changes no value. Only one is made (it patches
    ``VoxelCostFunction.evaluate`` for the process)."""

    def __init__(self) -> None:
        from icenine.cost_functions import VoxelCostFunction

        SR.Meter.__init__(self)
        VoxelCostFunction.evaluate = SR.patch_evaluate(VoxelCostFunction.evaluate, self)


STALLHIT: Dict[str, Any] = {"k": SR.K}
STATS: Dict[str, int] = dict(
    mc_refinements=0, stalls=0, cut=0, evals_in_cut=0, checks_went_on=0, checks_total=0
)


class StallHitPick(RRUN.RescuePick):
    """``RescuePick`` whose MC fits obey rule ``PH(k)`` (``STALLHIT['k']``; 0: off, as arm G)."""

    ptracer: Optional[StallHitTracer] = None

    def _fit(
        self, kind: str, rec: Any, j: int, Rj: np.ndarray, cj: float, vertices: Any, phase: int
    ) -> Dict[str, Any]:
        if StallHitPick.ptracer is None:
            StallHitPick.ptracer = StallHitTracer()
            # the base ``_fit`` finds a tracer set (it reads the module's ``OrderPick`` name, which
            # ``main`` rebinds to this class; the replay does not rebind it)
            ORUN.OrderPick.tracer = StallHitPick.ptracer
            StallHitPick.tracer = StallHitPick.ptracer
        if kind != "M":
            return super()._fit(kind, rec, j, Rj, cj, vertices, phase)
        tr = StallHitPick.ptracer
        tr.begin(int(STALLHIT["k"]))

        def hit_fn(R: np.ndarray) -> float:
            info = rec.evaluate_overlap(np.asarray(R, np.float64), vertices, phase)
            return float(self.bfs._ratios(info)[0])

        tr.hit_fn = hit_fn
        cut = False
        try:
            Rend, c, ev = MC._orig_refine_arm("M", rec, self.orig_lo, Rj, cj, vertices, phase)
            ev = int(ev) + tr.continued  # each stall check that went on is one counted evaluation
        except StallHitStop:
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
        STATS["stalls"] += int(bool(tr.checks))
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
                checks=[list(x) for x in tr.checks],
                continued=int(tr.continued),
            )
        )
        return dict(j=int(j), kind="M", hit=hit, cost=float(c), evals=int(ev), R=np.asarray(Rend))


def summary() -> Dict[str, Any]:
    return dict(k=STALLHIT["k"], rule=SR.RULE, **STATS)


def main() -> None:
    args = sys.argv[1:]
    k = SR.K
    if "--k" in args:
        i = args.index("--k")
        k = int(args[i + 1])
        del args[i : i + 2]
    if FC.SAMPLE in ("held", "held2") and k != SR.K:
        raise SystemExit("stallhit rule: held and held2 run only with the frozen k")
    STALLHIT["k"] = k
    sys.argv = ["order_run.py", "--arm", "A", "--name", "H"] + args
    ORUN.OrderPick = StallHitPick  # type: ignore[misc]
    ORUN.OrderWrap = RRUN.RescueWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"H_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(RRUN.LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(RRUN.summary(), indent=1) + "\n")
    (out / "stallhit.json").write_text(json.dumps(summary(), indent=1) + "\n")


if __name__ == "__main__":
    main()
