"""AR5 step 11, phase A step 3: freeze rule P(3000) and the files before ``held`` is read.

  uv run python scripts/autoresearch/ar5/patience_freeze.py --write   # writes patience/frozen.json
  uv run python scripts/autoresearch/ar5/patience_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text, k, and the sha256 of the scripts that define and run the rule
(and of the step 8 and step 9 scripts it builds on, the same list as step 9's freeze) and of every
library file. ``--check`` also needs the file committed and unchanged: no seed log, trace, record,
rescue log or map of ``held`` is opened by this step's scripts before the freeze commit. No hash
of any ``held`` file is taken (that would read it). Scripts that never touch ``held`` (the what-if,
the replay, the gate, the launcher) are not in the list, so they may change after the freeze.

Why a new freeze (as step 10's ``cut_freeze``): step 9's phase B changed library files, so the
step 8 and step 9 freeze checks no longer agree. This file freezes the same script list and the
library as it is now; the older scripts keep their own check files, which ``patience_common``
switches off in this step's processes after this check has passed. No import of the study's other
modules except ``patience_rules`` (pure) and ``rescue_freeze``'s list.
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
FROZEN = ICE / "benchmarks" / "autoresearch" / "ar5" / "patience" / "frozen.json"
sys.path.insert(0, str(HERE))
import patience_rules as PR  # noqa: E402
import rescue_freeze as RF  # noqa: E402  (pure: hashing only)

SCRIPTS = ["patience_rules.py", "patience_common.py", "patience_run.py"] + list(RF.SCRIPTS)
CONSTANTS = dict(
    k=PR.K,
    stage="seed step MC refinement only (MC stage then variance stage, one counter)",
    tau=PR.TAU,
    seed_test=PR.SEED_TEST,
    candidates=5,
    rule_g1="unchanged (rescue_rules.py)",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    files = {f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in SCRIPTS}
    lib = {str(p.relative_to(ICE)): sha(p) for p in sorted((ICE / "icenine").glob("*.py"))}
    lib["pyproject.toml"] = sha(ICE / "pyproject.toml")
    lib["uv.lock"] = sha(ICE / "uv.lock")
    return dict(files=files, library=lib)


def write() -> None:
    cur = current()
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(
        json.dumps(
            dict(
                rule=PR.RULE,
                constants=CONSTANTS,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no held file was hashed or opened before this file was committed; steps 8 "
                    "and 9 froze the library before step 9's phase B (key SeedRescueCMA, off by "
                    "default), so their frozen.json no longer agree; this file freezes the same "
                    "script list and the library as it is now (as step 10's cut_freeze)"
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
        return "patience frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "patience frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "patience frozen.json differs from its commit"
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
