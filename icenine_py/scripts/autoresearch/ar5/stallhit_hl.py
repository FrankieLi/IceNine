"""AR5 step 12, phase B: one whole-layer run of the script arm ``Hl`` (arm ``H`` with library
draws).

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/stallhit_hl.py [--cap-min 240]

Arm ``Hl`` is step 9's script arm ``Gl`` (``rescue_gl.py``: step 8's arm ``Al``, the interleaved
order with the library's own draws, one generator for the whole run, plus rule ``G1`` in front of
the grid fallback) with rule ``PH(3000)`` in every seed-step MC refinement, exactly as
``stallhit_run.py`` puts it on arm ``G`` (same ``StallHitPick``, ``RescueWrap``, rule, counting and
logs). It is the script arm that the library run with ``SeedIndexMCStallHit 3000`` must equal
before the boundary pass (plan, phase B). ``stallhit_run.py`` is frozen (``stallhit/frozen.json``)
and is not edited. Writes ``<cache>/Hl_<sample>/`` (the files of ``order_run.py`` plus
``rescue_log.json.gz``, ``rescue.json``, ``stallhit.json``) and ``Hl.done`` there at the end. The
boundary pass is made after by ``order_pass.py --name Hl --aligned 0`` (the library's single
generator). mid, grad and big only.
"""

import json
import sys

import stallhit_common as SC  # noqa: E402,F401  (first: sets the cache)

import fit_common as FC  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_run as RRUN  # noqa: E402
import stallhit_rules as SR  # noqa: E402
import stallhit_run as SRUN  # noqa: E402


def main() -> None:
    args = sys.argv[1:]
    assert "--smoke" not in args, "no smoke run: a truncated map does not end quickly on mid"
    assert FC.SAMPLE in ("mid", "grad", "big"), "phase B runs mid, grad or big only"
    SRUN.STALLHIT["k"] = SR.K
    sys.argv = ["order_run.py", "--arm", "Al", "--name", "Hl"] + args
    ORUN.OrderPick = SRUN.StallHitPick  # type: ignore[misc]
    ORUN.OrderWrap = RRUN.RescueWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    out = FC.CACHE / f"Hl_{FC.SAMPLE}"
    RC.write_gz(RRUN.LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(RRUN.summary(), indent=1) + "\n")
    (out / "stallhit.json").write_text(json.dumps(SRUN.summary(), indent=1) + "\n")
    (out / "Hl.done").write_text("done\n")


if __name__ == "__main__":
    main()
