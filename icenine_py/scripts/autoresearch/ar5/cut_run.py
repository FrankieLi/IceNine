"""AR5 step 10, phase A step 5: one whole-layer run of arm ``C`` (script layer; library unchanged).

  RESEED_SAMPLE=<big|grad|mid|held> uv run python scripts/autoresearch/ar5/cut_run.py \
      [--cap-min 240] [--smoke N] [--guards-only]

Arm ``C`` is arm ``G`` of step 9 (``rescue_run.py``: arm ``A`` of step 8 with rule ``G1`` in front
of the grid fallback, interleaved order, aligned draws) with the level-0 cut of ``cut_rules``: every
grid search runs the quick MC of level 0 only on the first N = 1,000 listed candidates. It is
``rescue_run.py`` with the run named ``C`` and the wrapper subclassed to install the cut. Each grid
search keeps its aligned generator (voxel, "grid", 0), so the cut moves no other draw. Besides the
files of ``rescue_run.py`` it writes ``cut_log.json.gz`` (one row per discrete-stage call: level 0
or not, the list length and the number of candidates run) and ``cut.json`` (its counts). On
``big`` no grid search runs, so the run must equal arm ``G`` there. The boundary pass is made after
by ``order_pass.py --name C``. ``held`` is refused until ``ar5/cut/frozen.json`` is committed and
clean and ``cut_freeze`` agrees (``cut_common``). ``cut_freeze`` replaces step 9's and step 8's
freeze checks, which ``cut_common`` switches off in this process (step 9's phase B changed the
library after those freezes); ``order_common.require_frozen()`` still applies.
"""

import json
import sys
from typing import Any, Dict, List

import cut_common as CC  # noqa: E402  (first: sets the cache, guards held)

import cut_rules as CR  # noqa: E402
import fit_common as FC  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_run as RRUN  # noqa: E402

CUT_LOG: List[Dict[str, Any]] = []


class CutWrap(RRUN.RescueWrap):
    """``RescueWrap`` with the level-0 cut installed in front of the discrete stage."""

    def __init__(self, bfs: Any) -> None:
        super().__init__(bfs)
        CC.install_cut(self.rec, CR.TOP_N, CUT_LOG)


def summary() -> Dict[str, Any]:
    l0 = [r for r in CUT_LOG if r["level0"]]
    return dict(
        top_n=CR.TOP_N,
        rule=CR.RULE,
        discrete_calls=len(CUT_LOG),
        level0_calls=len(l0),
        level0_listed=[int(r["n_listed"]) for r in l0],
        level0_run=[int(r["n_run"]) for r in l0],
        level0_cut_applied=sum(1 for r in l0 if r["n_run"] < r["n_listed"]),
    )


def main() -> None:
    args = sys.argv[1:]
    sys.argv = ["order_run.py", "--arm", "A", "--name", "C"] + args
    ORUN.OrderPick = RRUN.RescuePick  # type: ignore[misc]
    ORUN.OrderWrap = CutWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"C_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(RRUN.LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(RRUN.summary(), indent=1) + "\n")
    RC.write_gz(CUT_LOG, out / "cut_log.json.gz")
    (out / "cut.json").write_text(json.dumps(summary(), indent=1) + "\n")


if __name__ == "__main__":
    main()
