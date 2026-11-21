"""AR5 step 22 (task ar5v): one whole-layer run of arm ``J5D`` (script layer; library unchanged).

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/defer_launch.py defer_run.py [--cap-min 300] [--smoke N] \
      [--defer-off] [--name J5D]

Arm ``J5D`` is arm ``J5`` (``high500_run.py``: arm ``H`` of step 12 with rule ``G2b`` between a
failed ``G1`` rescue and the grid search; interleaved order, aligned draws; rule ``PH(500)`` in the
seed step's MC refinements) with rule ``DEF`` (``defer_pick.DeferPick``: a stalled refinement is
paused while untried candidates remain). ``--defer-off`` runs the same layer with nothing paused
(the seed step is ``J5``'s; gate g0), under the name given (default ``J5D0``). ``G2b`` takes no draw
from the aligned book and the grid call after a failed ``G2b`` is made unchanged, as in ``J5``.
Writes, besides the files of ``order_run.py``: ``polish_log.json.gz`` and ``polish.json`` (as
``J5``), ``high500.json`` (the stall checks), ``defer.json`` (pauses, resumes, steps by outcome,
the re-run's evaluations, listed apart and not counted) and ``defer_steps.json.gz`` (one row per
seed step: the order the fits ran, paused and resumed positions, where the end was found). The
boundary pass is made after by ``high500_pass.py --name J5D`` (through ``defer_launch.py``).
"""

import json
import sys
from typing import Any, Dict, List

import defer_common as DC  # noqa: E402,F401  (first: guards, freeze, cache)

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
    off = bool(_pop(args, "--defer-off", False))
    name = _pop(args, "--name", True) or ("J5D0" if off else "J5D")
    DP.DEFER["on"] = not off
    HP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "A", "--name", name] + args
    ORUN.OrderPick = DP.DeferPick  # type: ignore[misc]
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


if __name__ == "__main__":
    main()
