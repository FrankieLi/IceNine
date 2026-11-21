"""AR5 step 9, phase A: one whole-layer run of arm ``G`` (script layer; the library is unchanged).

  RESEED_SAMPLE=<big|grad|mid|held> uv run python scripts/autoresearch/ar5/rescue_run.py \
      [--cap-min 240] [--smoke N] [--guards-only]

Arm ``G`` is arm ``A`` of step 8 (run ``S11``: ``SeedIndexCMAFirst 200``,
``BFSNeighborCMANoFallback 1``, interleaved order, aligned draws) with rule ``G1`` in front of the
grid fallback (``rescue_rules``): after a failed seed test the kept ends of the seed step are tried
in order of cost with a CMA local optimization (0.25 degree, 1,000 evaluations, ``default_rng(0)``),
the first end with hit ratio at least 0.95 is taken, else the grid fallback runs as in ``A`` (same
aligned draws). It is ``order_run.py --arm A --name G`` with the pick and the wrapper subclassed.
The rescue evaluations are added to the fallback call's evaluation count (so ``run.json`` totals
add up); ``rescue_log.json.gz`` holds one row per fallback call (tries, evaluations, the grid
evaluations left) and ``rescue.json`` its counts. The boundary pass is made after by
``order_pass.py --name G``. ``held`` is refused until ``ar5/rescue/frozen.json`` is committed
(``rescue_common``) and ``order_freeze.py --check`` agrees (``order_run``).
"""

import json
import sys
from typing import Any, Dict, List

import rescue_common as RCM  # noqa: E402  (first: sets the cache, guards held)

import fit_common as FC  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402
import reseed_s8 as S8  # noqa: E402
import rescue_rules as RR  # noqa: E402

LOG: List[Dict[str, Any]] = []


class RescuePick(ORUN.OrderPick):
    """``OrderPick`` that remembers the kept ends of its last seed step (for the wrapper)."""

    last_step: Any = None

    def e1_step(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2, tau: float = -1
    ) -> Any:
        out = super().e1_step(vertices, phase_index, n_cost, n_refine, tau)
        entry = self.e1_log[-1]
        kept = [
            dict(R=R, cost=float(c), hit=float(h))
            for R, c, h in zip(out.finisher_R, out.finisher_costs, entry["fin_hit"])
        ]
        self.last_step = dict(
            key=OD.voxel_key(vertices), voxel=int(self.cur), kept=kept, used=False
        )
        return out


class RescueWrap(ORUN.OrderWrap):
    """``OrderWrap`` with rule ``G1`` in front of the (aligned) grid fallback."""

    def __init__(self, bfs: Any) -> None:
        super().__init__(bfs)

        def pick() -> Any:
            return S8.state["pick"]

        def try_end(e: Dict[str, Any], vertices: Any, phase: int) -> Dict[str, Any]:
            p = pick()
            return RCM.cma_rescue(self.rec, p.orig_lo, p.bfs._ratios, e["R"], vertices, phase)

        self.rec.reconstruct_voxel = RR.make_wrapper(
            self.rec.reconstruct_voxel,
            lambda: pick().last_step,
            OD.voxel_key,
            try_end,
            self.rec,
            LOG,
        )


def summary() -> Dict[str, Any]:
    att = [r for r in LOG if r["attempted"]]
    return dict(
        fallback_calls=len(LOG),
        rescues_attempted=len(att),
        rescued=sum(1 for r in att if r["rescued"]),
        tries=sum(len(r["tries"]) for r in att),
        rescue_evals=int(sum(r["rescue_evals"] for r in LOG)),
        grid_fallbacks_left=sum(1 for r in LOG if r["grid_evals"] > 0),
        grid_evals=int(sum(r["grid_evals"] for r in LOG)),
        rule=RR.RULE,
    )


def main() -> None:
    args = sys.argv[1:]
    sys.argv = ["order_run.py", "--arm", "A", "--name", "G"] + args
    ORUN.OrderPick = RescuePick  # type: ignore[misc]
    ORUN.OrderWrap = RescueWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"G_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(summary(), indent=1) + "\n")


if __name__ == "__main__":
    main()
