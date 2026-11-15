"""AR5 step 16, phase A step 5 (run): one whole-layer run of arm ``J5`` (script layer; library
unchanged).

  RESEED_SAMPLE=<big|grad|mid|held|held3> HIGH500_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/high500_run.py [--cap-min 300] [--smoke N]

Arm ``J5`` is arm ``J`` of step 14 (``polish_run.py``: arm ``H`` of step 12 with rule ``G2b``
between a failed ``G1`` rescue and the grid search; interleaved order, aligned draws) with rule
``PH(500)`` (``high500_rules``) in place of ``PH(3000)`` in the seed step's MC refinements. It is
``polish_run.py`` with the run named ``J5``, the pick replaced (``high500_pick.High500Pick``) and
``G2b`` read from step 14's frozen file. ``G2b`` takes no draw from the aligned book and the grid
call after a failed ``G2b`` is made unchanged, so its draws are those of ``J``. Writes, besides the
files of ``order_run.py``: ``polish_log.json.gz`` and ``polish.json`` (as ``J``), ``high500.json``
(the stall checks: stops, checks that went on). The boundary pass is made after by
``high500_pass.py --name J5``.

``held`` and ``held3`` are refused unless ``ar5/high500/frozen.json`` is committed and agrees
(``high500_common``).
"""

import json
import sys
from typing import Any, Dict, List

import high500_common as BC  # noqa: E402  (first: sets the cache, guards held2, held and held3)

import fit_common as FC  # noqa: E402
import high500_pick as BP  # noqa: E402
import high500_rules as HR  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import polish_rules as PR  # noqa: E402
import reseed_common as RC  # noqa: E402
import reseed_s8 as S8  # noqa: E402
import rescue_common as RCM  # noqa: E402
import rescue_rules as RR  # noqa: E402

LOG: List[Dict[str, Any]] = []
G2 = BC.g2_form()


class High500Wrap(ORUN.OrderWrap):
    """``OrderWrap`` with rules ``G1`` then ``G2b`` in front of the (aligned) grid fallback (the
    wrapper of ``polish_run.PolishWrap``)."""

    def __init__(self, bfs: Any) -> None:
        super().__init__(bfs)

        def pick() -> Any:
            return S8.state["pick"]

        def g1_end(e: Dict[str, Any], vertices: Any, phase: int) -> Dict[str, Any]:
            p = pick()
            return RCM.cma_rescue(self.rec, p.orig_lo, p.bfs._ratios, e["R"], vertices, phase)

        def g2_end(R: Any, vertices: Any, phase: int) -> Dict[str, Any]:
            p = pick()
            return PR.g2_try(G2, self.rec, p.orig_lo, p.bfs._ratios, R, vertices, phase)

        self.rec.reconstruct_voxel = PR.make_wrapper(
            G2,
            self.rec.reconstruct_voxel,
            lambda: pick().last_step,
            OD.voxel_key,
            g1_end,
            g2_end,
            self.rec,
            LOG,
            RR.rescue_choose,
        )


def summary() -> Dict[str, Any]:
    att = [r for r in LOG if r["attempted"]]
    return dict(
        form=G2,
        fallback_calls=len(LOG),
        rescues_attempted=len(att),
        rescued=sum(1 for r in att if r["rescued"]),
        polished=sum(1 for r in att if r["polished"]),
        g1_tries=sum(len(r["tries"]) for r in att),
        g2_tries=sum(len(r["polish_tries"]) for r in att),
        rescue_evals=int(sum(r["rescue_evals"] for r in LOG)),
        polish_evals=int(sum(r["polish_evals"] for r in LOG)),
        grid_fallbacks_left=sum(1 for r in LOG if r["grid_evals"] > 0),
        grid_evals=int(sum(r["grid_evals"] for r in LOG)),
        rule=PR.RULE,
    )


def main() -> None:
    args = sys.argv[1:]
    BP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "A", "--name", "J5"] + args
    ORUN.OrderPick = BP.High500Pick  # type: ignore[misc]
    ORUN.OrderWrap = High500Wrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"J5_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(LOG, out / "polish_log.json.gz")
    (out / "polish.json").write_text(json.dumps(summary(), indent=1) + "\n")
    (out / "high500.json").write_text(json.dumps(BP.summary(), indent=1) + "\n")


if __name__ == "__main__":
    main()
