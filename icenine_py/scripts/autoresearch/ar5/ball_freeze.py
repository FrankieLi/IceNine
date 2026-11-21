"""AR5 step 23, phase A step 3: freeze rule ``CB``, the gate, the decision rules, the predictions
and the scripts that define and set up arms ``J5DC`` and ``J5SDC``, before any ``J5DC`` run (live,
on any sample), before any replay and before any sample is read by this step's scripts.

  uv run python scripts/autoresearch/ar5/ball_freeze.py --write  # ar5/ball/frozen.json
  uv run python scripts/autoresearch/ar5/ball_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text and constants (``ball_rules``), the gate, the decision rules,
the predictions, the held7 spec, and the sha256 of the scripts that define and set up the arms
(this step's ``OWN`` list), of every script step 22's freeze lists (its own scripts and every older
step's) and of EVERY LIBRARY FILE (``icenine/*.py``, ``pyproject.toml``, ``uv.lock``): the check
pins file hashes, never names only. ``--check`` also needs the file committed and unchanged. No
hash of any ``held7`` or ``held2`` data file is taken (that would read it).

NOT in the list, on purpose: the scripts that apply the frozen rule to stored data or judge runs
(``ball_g0``, ``ball_replay``, ``ball_gate``, ``ball_summary``, ``ball_build_held7``,
``ball_run_all.sh``). Everything that decides a voxel (the rule, the live pick, the run) is in
``ball_rules``, ``ball_pick`` and ``ball_run``.

Why a layer of its own: step 22's freeze (``defer_freeze``) no longer agrees at this commit (the
library's key phases changed ``icenine/config_file.py``), and so do the older steps' checks. This
step's own check runs FIRST (``ball_common``); only when it passes are the older checks switched
off, in the script layer of this process only (the precedent: ``defer_common``). No hash of any
older step is re-recorded and no older script and no older frozen file is edited.
"""

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ball"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import ball_rules as BRU  # noqa: E402
import defer_freeze as DFZ  # noqa: E402  (its hashing lists only; pure)

OWN = [
    "ball_rules.py",
    "ball_guard.py",
    "ball_samples.py",
    "ball_common.py",
    "ball_pick.py",
    "ball_run.py",
    "ball_launch.py",
    "ball_link_cache.sh",
]
CONSTANTS = dict(
    rule=BRU.RULE,
    arms=[BRU.ARM_DC, BRU.ARM_SDC],
    fail_halve=BRU.FAIL_HALVE,
    floor_deg=BRU.FLOOR_DEG,
    cap_proposals=BRU.CAP_PROPOSALS,
    start_radius="the box (search_box_side passed to the variance stage, 0.329 degree)",
    tau=BRU.TAU,
    held7="build of held: 250 Voronoi grains on the 24,570-voxel grid, gradients 0 to 4 degrees, "
    "generator seed 7 (a new sample for a clean test, not an extra random seed), clean images, "
    "Q-max 8",
    samples_disclosed=list(BRU.DISCLOSED),
    live_samples=list(BRU.LIVE_SAMPLES),
    caps_min=BRU.CAPS_MIN,
    lost_end_denominator="all replayed seed steps of the sample (also reported: steps whose J5D "
    "run had an end at 0.95)",
    freeze_note="the freeze check pins the sha256 of every script and library file; step 22's own "
    "check is switched off in this process after this one passes",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    cur = DFZ.current()  # step 22's list (with steps 21 to 14 and the AR3 scripts) and the library
    files = dict(cur["files"])
    for s in sorted(OWN):
        files[f"scripts/autoresearch/ar5/{s}"] = sha(HERE / s)
    return dict(files=files, library=dict(cur["library"]))


def write() -> None:
    cur = current()
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(
        json.dumps(
            dict(
                constants=CONSTANTS,
                predictions=BRU.PREDICTIONS,
                gate=BRU.GATE,
                judge=BRU.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no J5DC run (live, on any sample), no replay and no g0 run was made by the "
                    "tracked scripts, and nothing of held7 was built or read, before this file was "
                    "committed; step 22's stored records were only checked for their counts "
                    "(read only); nothing of held2 or fine7 is used by this step; held7 does not "
                    "exist yet"
                ),
                **cur,
            ),
            indent=1,
        )
        + "\n"
    )
    print("wrote", FROZEN)


def problems() -> str:
    """'' when frozen.json is committed, clean and agrees with the files; else the reason."""
    if not FROZEN.exists():
        return "ball frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "ball frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "ball frozen.json differs from its commit"
    saved = json.loads(FROZEN.read_text())
    cur = current()
    for part in ("files", "library"):
        for k, v in cur[part].items():
            if saved[part].get(k) != v:
                return f"{k} differs from the frozen hash"
        for k in saved[part]:
            if k not in cur[part]:
                return f"{k} is frozen but missing"
    return ""


def check_or_exit() -> None:
    p = problems()
    if p:
        print("freeze check failed:", p, file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    if "--write" in sys.argv:
        write()
    elif "--check" in sys.argv:
        check_or_exit()
        print("freeze check passed")
    else:
        raise SystemExit("use --write or --check")
