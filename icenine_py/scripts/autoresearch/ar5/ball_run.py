"""AR5 step 23 (task ar5w): one whole-layer run of arm ``J5DC`` (script layer; library unchanged).

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/ball_launch.py ball_run.py [--cap-min 300] [--smoke N] \
      [--ball-off] [--name J5DC]

Arm ``J5DC`` is arm ``J5D`` (``defer_run.py``: ``J5`` with rule ``DEF``) with rule ``CB`` in the
seed step's MC refinements (``ball_pick.BallPick``: the second stage is the converging ball).
``--ball-off`` runs the same layer with the library's stage (arm ``J5D``'s seed step; the run is
named ``J5DC0`` unless ``--name`` says otherwise). ``G2b`` takes no draw from the aligned book and
the grid call after a failed ``G2b`` is made unchanged, as in ``J5``. Writes, besides the files of
``order_run.py``: ``polish_log.json.gz`` and ``polish.json`` (as ``J5``), ``high500.json`` (the
stall checks: all zero under CB), ``defer.json`` (pauses and resumes: all zero under CB),
``defer_steps.json.gz``, ``ball.json`` (counts of the refinements' second stages) and
``ball_steps.json.gz`` (one row per refinement). The boundary pass is made after by
``high500_pass.py --name J5DC`` (through ``ball_launch.py``).
"""

import json
import sys
from typing import Any, Dict, List

import ball_common as BC  # noqa: E402,F401  (first: guards, freeze, cache)

import ball_pick as BP  # noqa: E402
import defer_pick as DP  # noqa: E402
import fit_common as FC  # noqa: E402
import high500_pick as HP  # noqa: E402
import high500_rules as HR  # noqa: E402
import high500_run as HRUN  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402


def _pop(args: List[str], flag: str, has_value: bool) -> Any:
    if flag not in args:
        return None
    i = args.index(flag)
    v = args[i + 1] if has_value else True
    del args[i : i + (2 if has_value else 1)]
    return v


def main() -> None:
    args = sys.argv[1:]
    off = bool(_pop(args, "--ball-off", False))
    name = _pop(args, "--name", True) or ("J5DC0" if off else "J5DC")
    BP.BALL["on"] = not off
    DP.DEFER["on"] = True
    HP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "A", "--name", name] + args
    ORUN.OrderPick = BP.BallPick  # type: ignore[misc]
    ORUN.OrderWrap = HRUN.High500Wrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    smoke = "--smoke" in args and int(args[args.index("--smoke") + 1]) > 0
    out = FC.CACHE / f"{name}_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(HRUN.LOG, out / "polish_log.json.gz")
    (out / "polish.json").write_text(json.dumps(HRUN.summary(), indent=1) + "\n")
    (out / "high500.json").write_text(json.dumps(HP.summary(), indent=1) + "\n")
    d: Dict[str, Any] = DP.summary()
    d["name"] = name
    (out / "defer.json").write_text(json.dumps(d, indent=1) + "\n")
    RC.write_gz(DP.LOG, out / "defer_steps.json.gz")
    b: Dict[str, Any] = BP.summary()
    b["name"] = name
    (out / "ball.json").write_text(json.dumps(b, indent=1) + "\n")
    RC.write_gz(BP.LOG, out / "ball_steps.json.gz")
    if not off and (BP.STATS["checks"] or BP.STATS["pauses"] or HP.STATS["checks_total"]):
        raise SystemExit("STOP RULE: a check or a pause fired under CB")


if __name__ == "__main__":
    main()
