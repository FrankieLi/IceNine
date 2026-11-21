"""AR5 step 22, phase A step 3: freeze rule ``DEF``, the gate, the decision rules, the predictions,
the plan's table and the scripts that define and set up arms ``J5D`` and ``J5SD``, before any
``J5D`` run (live, on any sample) and before any sample is read by this step's scripts.

  uv run python scripts/autoresearch/ar5/defer_freeze.py --write  # ar5/defer/frozen.json
  uv run python scripts/autoresearch/ar5/defer_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text and constants (``defer_rules``), the gate, the decision rules,
the predictions, the plan's table, the held6 spec, and the sha256 of the scripts that define and set
up the arms (this step's ``OWN`` list), of every script that step 21's freeze lists (step 21's own
scripts, step 20's and step 19's repair, step 18's list with steps 14 to 17's and the AR3 scripts)
and of EVERY LIBRARY FILE (``icenine/*.py``, ``pyproject.toml``, ``uv.lock``): the check pins file
hashes, it never matches by name only. ``--check`` also needs the file committed and unchanged. No
hash of any ``held6`` or ``held2`` data file is taken (that would read it).

NOT in the list, on purpose: the scripts that apply the frozen rule to stored data or judge runs
(``defer_replay``, ``defer_g0``, ``defer_gate``, ``defer_summary``, ``defer_build_held6``,
``defer_run_all.sh``). Everything that decides a voxel (the rule, the seed step's state machine,
the live pick, the run) is in ``defer_rules``, ``defer_pick`` and ``defer_run``.

Why a layer of its own: step 21's freeze (``screen_freeze``) no longer agrees at this commit
(``icenine/boundary_pass.py`` changed with step 21's key phase), and so do the older steps' checks.
This step's own check runs FIRST (``defer_common``); only when it passes are the older checks
switched off, in the script layer of this process only (the precedent of ``screen_common``,
``cmarepair_common`` and ``stallhit_common``). No hash of any older step is re-recorded and no
older script and no older frozen file is edited.
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "defer"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import defer_rules as DR  # noqa: E402
import screen_freeze as SFZ  # noqa: E402  (its hashing lists only; pure)

OWN = [
    "defer_rules.py",
    "defer_guard.py",
    "defer_samples.py",
    "defer_common.py",
    "defer_pick.py",
    "defer_run.py",
    "defer_launch.py",
    "defer_link_cache.sh",
]
CONSTANTS = dict(
    rule=DR.RULE,
    arms=[DR.ARM_D, DR.ARM_SD],
    tau=DR.TAU,
    stall_k=DR.K,
    n_candidates=DR.N_CAND,
    plan_table=DR.PLAN_TABLE,
    held6="build of held: 250 Voronoi grains on the 24,570-voxel grid, gradients 0 to 4 degrees, "
    "generator seed 6 (a new sample for a clean test, not an extra random seed), clean images, "
    "Q-max 8",
    samples_disclosed=list(DR.DISCLOSED),
    caps_min=DR.CAPS_MIN,
    resume_in_phase_a="a paused refinement is resumed by an exact re-run of the paused part: the "
    "re-run's evaluations up to the pause are checked identical to the first run (trace, check "
    "evaluation and hit ratio) and are NOT counted; they are listed apart (rerun_evals). Every "
    "other evaluation is counted, the checks included.",
    freeze_note="the freeze check pins the sha256 of every script and library file; step 21's "
    "own check is switched off in this process after this one passes",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    cur = SFZ.current()  # step 21's list (with step 20, 19, 18 and 14 to 17) and the library
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
                predictions=DR.PREDICTIONS,
                gate=DR.GATE,
                judge=DR.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no J5D run (live, on any sample) was made by the tracked scripts and nothing "
                    "of held6 was built or read before this file was committed; the replay of "
                    "J5's stored seed steps (arithmetic, read only; defer_replay.py) was run "
                    "before it and reproduced the plan's table; nothing of held2 is used by this "
                    "step; held6 does not exist yet"
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
        return "defer frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "defer frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "defer frozen.json differs from its commit"
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
