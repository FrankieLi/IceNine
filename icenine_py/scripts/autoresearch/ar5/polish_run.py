"""AR5 step 14, phase A step 6 (run): one whole-layer run of arm ``J`` (script layer; library
unchanged).

  RESEED_SAMPLE=<big|grad|mid|held|held2> POLISH_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/polish_run.py [--cap-min 300] [--smoke N] [--guards-only]

Arm ``J`` is arm ``H`` of step 12 (``stallhit_run.py``: arm ``G`` of step 9 with rule ``PH(3000)``;
interleaved order, aligned draws, rule ``G1`` in front of the grid fallback) with rule ``G2``
(``polish_rules``, the form frozen in ``ar5/polish/frozen.json``) between a failed ``G1`` rescue and
the grid search. It is ``stallhit_run.py`` with the run named ``J`` and the wrapper replaced
(``PolishWrap``). ``G2`` takes no draw from the aligned book (its CMA fits have their own
``default_rng(1)``), and the grid call after a failed ``G2`` is made unchanged, so its draws are
those of ``H``. Writes, besides the files of ``order_run.py``: ``polish_log.json.gz`` (one row per
fallback call: the ``G1`` tries, the ``G2`` tries, evaluations, the grid evaluations left),
``polish.json`` (counts) and ``stallhit.json``. The boundary pass is made after by
``polish_pass.py --name J``. Every sample needs the freeze check (``polish_freeze.py --check``);
``held2`` is read only after it (``polish_common``).
"""

import json
import os
import sys
from typing import Any, Dict, List

import polish_common as PC  # noqa: E402  (first: sets the cache, guards held2)

import polish_freeze as PF  # noqa: E402

PF.check_or_exit()
os.environ["STALLHIT_CACHE"] = str(PC.PCACHE)  # stallhit_common must not point back at step 12's

import stallhit_common as SC  # noqa: E402,F401  (its freeze guard on held and held2 stays on)

import fit_common as FC  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import polish_rules as PR  # noqa: E402
import reseed_common as RC  # noqa: E402
import reseed_s8 as S8  # noqa: E402
import rescue_common as RCM  # noqa: E402
import rescue_rules as RR  # noqa: E402
import stallhit_run as SRUN  # noqa: E402

LOG: List[Dict[str, Any]] = []
FORM = PF.frozen_form()


class PolishWrap(ORUN.OrderWrap):
    """``OrderWrap`` with rules ``G1`` then ``G2`` in front of the (aligned) grid fallback."""

    def __init__(self, bfs: Any) -> None:
        super().__init__(bfs)

        def pick() -> Any:
            return S8.state["pick"]

        def g1_end(e: Dict[str, Any], vertices: Any, phase: int) -> Dict[str, Any]:
            p = pick()
            return RCM.cma_rescue(self.rec, p.orig_lo, p.bfs._ratios, e["R"], vertices, phase)

        def g2_end(R: Any, vertices: Any, phase: int) -> Dict[str, Any]:
            p = pick()
            return PC.g2_try(FORM, self.rec, p.orig_lo, p.bfs._ratios, R, vertices, phase)

        self.rec.reconstruct_voxel = PR.make_wrapper(
            FORM,
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
        form=FORM,
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
    SRUN.STALLHIT["k"] = SRUN.SR.K
    sys.argv = ["order_run.py", "--arm", "A", "--name", "J"] + args
    ORUN.OrderPick = SRUN.StallHitPick  # type: ignore[misc]
    ORUN.OrderWrap = PolishWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"J_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(LOG, out / "polish_log.json.gz")
    (out / "polish.json").write_text(json.dumps(summary(), indent=1) + "\n")
    (out / "stallhit.json").write_text(json.dumps(SRUN.summary(), indent=1) + "\n")


if __name__ == "__main__":
    main()
