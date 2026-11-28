"""AR5 step 22, phase B: one whole-layer run of the script arm ``J5Dl`` (arm ``J5l`` of step 16 with
rule ``DEF``: arm ``J5D`` without ``G2b``, with the library's own draws).

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/defer_hl.py [--cap-min 300]

``J5l`` is step 12's arm ``Hl`` with rule ``PH(500)``: the interleaved order with the library's own
draws (one generator for the whole run) plus rule ``G1`` in front of the grid fallback. ``J5Dl`` is
that arm with the seed step run by ``defer_pick.DeferPick`` (rule ``DEF``: a stalled refinement is
paused while untried candidates remain). It is the script arm that the library run with
``SeedIndexMCStallHit 500`` and ``SeedIndexMCDefer 1`` must equal before the boundary pass: the
library has neither ``G2b`` nor aligned draws, so ``J5D`` itself (which has both) is not the
comparison, as ``J5`` was not in step 16 (``J5l`` was). Writes ``<cache>/J5Dl_<sample>/`` (the files
of ``order_run.py``, ``rescue_log.json.gz``, ``rescue.json``, ``defer.json``,
``defer_steps.json.gz``, ``stallhit.json``) and ``J5Dl.done``. The boundary pass is not made.
``mid`` only. Goes through ``defer_common`` (this step's freeze, which must still agree: run it
BEFORE the library key is added).
"""

import json
import sys

import defer_common as DC  # noqa: E402,F401  (first: guards, freeze, cache)

import defer_pick as DP  # noqa: E402
import fit_common as FC  # noqa: E402
import high500_pick as HP  # noqa: E402
import high500_rules as HR  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_run as RRUN  # noqa: E402


def main() -> None:
    args = sys.argv[1:]
    assert "--smoke" not in args, "no smoke run: a truncated map does not end quickly on mid"
    assert FC.SAMPLE == "mid", "the library check is made on mid"
    DP.DEFER["on"] = True
    HP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "Al", "--name", "J5Dl"] + args
    ORUN.OrderPick = DP.DeferPick  # type: ignore[misc]
    ORUN.OrderWrap = RRUN.RescueWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    out = FC.CACHE / f"J5Dl_{FC.SAMPLE}"
    RC.write_gz(RRUN.LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(RRUN.summary(), indent=1) + "\n")
    (out / "stallhit.json").write_text(json.dumps(HP.summary(), indent=1) + "\n")
    (out / "defer.json").write_text(json.dumps(DP.summary(), indent=1) + "\n")
    RC.write_gz(DP.LOG, out / "defer_steps.json.gz")
    (out / "J5Dl.done").write_text("done\n")


if __name__ == "__main__":
    main()
