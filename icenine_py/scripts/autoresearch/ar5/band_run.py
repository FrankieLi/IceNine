"""AR5 step 15, phase A step 7 (run): one whole-layer run of arm ``JB`` (script layer; library
unchanged).

  RESEED_SAMPLE=<big|grad|mid|held|held3> BAND_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/band_run.py --form PB1000|PB500 [--cap-min 300] [--smoke N]

Arm ``JB`` is arm ``J`` of step 14 (``polish_run.py``: arm ``H`` of step 12 with rule ``G2b``
between a failed ``G1`` rescue and the grid search; interleaved order, aligned draws) with rule
``PB(k)`` (``band_rules``) in place of ``PH(3000)`` in the seed step's MC refinements. It is
``polish_run.py`` with the run named ``JB``, the pick replaced (``band_pick.BandPick``) and ``G2b``
read from step 14's frozen file. ``G2b`` takes no draw from the aligned book and the grid call
after a failed ``G2b`` is made unchanged, so its draws are those of ``J``. Writes, besides the files
of ``order_run.py``: ``polish_log.json.gz`` and ``polish.json`` (as ``J``), ``band.json`` (the
stall checks: stops high and low, checks that went on). The boundary pass is made after by
``band_pass.py --name JB``.

``held`` and ``held3`` are refused unless ``ar5/band/frozen.json`` is committed and agrees
(``band_common``) and the form is the one the gate chose (``ar5/band/chosen.json``, committed).
"""

import json
import sys
from typing import Any, Dict, List

import band_common as BC  # noqa: E402  (first: sets the cache, guards held2, held and held3)

import band_pick as BP  # noqa: E402
import band_rules as BR  # noqa: E402
import fit_common as FC  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import polish_rules as PR  # noqa: E402
import reseed_common as RC  # noqa: E402
import reseed_s8 as S8  # noqa: E402
import rescue_common as RCM  # noqa: E402
import rescue_rules as RR  # noqa: E402

LOG: List[Dict[str, Any]] = []
G2 = BC.g2_form()


class BandWrap(ORUN.OrderWrap):
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


def chosen_form() -> str:
    """The form the gate chose (``ar5/band/chosen.json``, committed), or ``SystemExit(2)``."""
    import subprocess

    p = BC.BD / "chosen.json"
    if not p.exists():
        BC.BG.refuse("band rule: ar5/band/chosen.json does not exist (the gate has not chosen)")
    tracked = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(p)],
        cwd=str(BC.ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not tracked:
        BC.BG.refuse("band rule: ar5/band/chosen.json is not committed")
    return str(json.loads(p.read_text())["chosen"])


def main() -> None:
    args = sys.argv[1:]
    if "--form" not in args:
        raise SystemExit("band_run.py needs --form PB1000|PB500")
    i = args.index("--form")
    form = args[i + 1]
    del args[i : i + 2]
    if form not in BR.FORMS:
        raise SystemExit(f"unknown form {form}")
    if RC.SAMPLE in ("held", "held3") and form != chosen_form():
        BC.BG.refuse(f"band rule: held and held3 run only with the chosen form, not {form}")
    BP.BAND.update(k=BR.FORMS[form], low=BR.LOW, aux=())
    sys.argv = ["order_run.py", "--arm", "A", "--name", "JB"] + args
    ORUN.OrderPick = BP.BandPick  # type: ignore[misc]
    ORUN.OrderWrap = BandWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"JB_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(LOG, out / "polish_log.json.gz")
    (out / "polish.json").write_text(json.dumps(summary(), indent=1) + "\n")
    (out / "band.json").write_text(json.dumps(dict(form=form, **BP.summary()), indent=1) + "\n")


if __name__ == "__main__":
    main()
