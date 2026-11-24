"""AR5 step 24 (task ar5x): one whole-layer run of arm ``J5DB`` (script layer; library unchanged).

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/bf_launch.py bf_run.py [--cap-min 300] [--smoke N] \
      [--bf-off] [--name J5DB]

Arm ``J5DB`` is arm ``J5D`` (``defer_run.py``: ``J5`` with rule ``DEF``) with rule ``BF``
(``bf_pick.BFPick``: the converging ball first on every candidate; J5D's stock second stage only
when no ball end reaches 0.95). ``--bf-off`` runs the same layer with the library's stage (arm
``J5D``'s seed step; the run is named ``J5DB0`` unless ``--name`` says otherwise). ``G2b`` takes no
draw from the aligned book and the grid call after a failed ``G2b`` is made unchanged, as in
``J5``. Writes, besides the files of ``order_run.py``: ``polish_log.json.gz`` and ``polish.json``
(as ``J5``), ``high500.json`` (the stall checks of the fallbacks), ``defer.json`` (pauses and
resumes of the fallbacks), ``defer_steps.json.gz``, ``ball.json`` (counts of the ball stages) and
``ball_steps.json.gz`` (one row per ball refinement), ``bf.json`` and ``bf_steps.json.gz`` (one row
per seed step: ball pass only or fallback, evaluations) and ``bf_features.json.gz`` (one feature
row per fallback step, no extra evaluation). The boundary pass is made after by
``high500_pass.py --name J5DB`` (through ``bf_launch.py``).
"""

import json
import sys
from typing import Any, Dict, List

import bf_common as BC  # noqa: E402,F401  (first: guards, freeze, cache)

import ball_pick as BP  # noqa: E402
import bf_pick as BFP  # noqa: E402
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
    off = bool(_pop(args, "--bf-off", False))
    name = _pop(args, "--name", True) or ("J5DB0" if off else "J5DB")
    BFP.BF["on"] = not off
    BP.BALL["on"] = not off
    DP.DEFER["on"] = True
    HP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "A", "--name", name] + args
    ORUN.OrderPick = BFP.BFPick  # type: ignore[misc]
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
    f: Dict[str, Any] = BFP.summary()
    f["name"] = name
    (out / "bf.json").write_text(json.dumps(f, indent=1) + "\n")
    RC.write_gz(BFP.LOG, out / "bf_steps.json.gz")
    RC.write_gz(BFP.FEATURES, out / "bf_features.json.gz")
    if not off and (BP.STATS["checks"] or BP.STATS["pauses"] or BFP.STATS["stage1_mismatch"]):
        raise SystemExit("STOP RULE: a check or a pause fired in the ball pass")


if __name__ == "__main__":
    main()
