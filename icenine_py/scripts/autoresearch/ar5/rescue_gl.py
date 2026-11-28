"""AR5 step 9, phase B: one whole-layer run of the script arm ``Gl`` (arm ``G`` with library draws).

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/rescue_gl.py [--cap-min 240]

Arm ``Gl`` is step 8's arm ``Al`` (run ``S11``: interleaved order, the library's own draws, one
generator for the whole run) with rule ``G1`` in front of the grid fallback, exactly as
``rescue_run.py`` puts it in front of arm ``A`` (same ``RescuePick``, ``RescueWrap``, rule and
log). It is the script arm that the library run with ``SeedRescueCMA 0.25 1000`` must equal
(plan, phase B). ``rescue_run.py`` itself is frozen (``rescue/frozen.json``) and is not edited.
Writes ``<cache>/Gl_<sample>/`` (the files of ``order_run.py`` plus ``rescue_log.json.gz``,
``rescue.json``) and ``Gl.done`` there at the end. The boundary pass is made after by
``order_pass.py --name Gl --aligned 0`` (the library's single generator).
"""

import json
import sys

import rescue_common as RCM  # noqa: E402,F401  (first: sets the cache, guards held)

import fit_common as FC  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_run as RRUN  # noqa: E402


def main() -> None:
    args = sys.argv[1:]
    assert "--smoke" not in args, "no smoke run: a truncated map does not end quickly on mid"
    sys.argv = ["order_run.py", "--arm", "Al", "--name", "Gl"] + args
    ORUN.OrderPick = RRUN.RescuePick  # type: ignore[misc]
    ORUN.OrderWrap = RRUN.RescueWrap  # type: ignore[misc]
    ORUN.main()
    out = FC.CACHE / f"Gl_{FC.SAMPLE}"
    RC.write_gz(RRUN.LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(RRUN.summary(), indent=1) + "\n")
    (out / "Gl.done").write_text("done\n")


if __name__ == "__main__":
    main()
