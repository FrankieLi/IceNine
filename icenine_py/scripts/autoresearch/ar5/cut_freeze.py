"""AR5 step 10, phase A step 3: freeze the level-0 cut and the files before ``held`` is read.

  uv run python scripts/autoresearch/ar5/cut_freeze.py --write   # writes cut/frozen.json
  uv run python scripts/autoresearch/ar5/cut_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text, N, and the sha256 of the scripts that define and run the cut
(and of the step 8 and step 9 scripts it builds on, the same list as step 9's freeze) and of every
library file. ``--check`` also needs the file committed and unchanged: no seed log, trace, rescue
log, map or level record of ``held`` is opened by this step's scripts before the freeze commit. No
hash of any ``held`` file is taken (that would read it). Scripts that never touch ``held`` (the
replay, the gate, the launcher) are not in the list, so they may change after the freeze. No
import of the study's other modules except ``cut_rules`` (pure) and ``rescue_freeze``'s list.
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
FROZEN = ICE / "benchmarks" / "autoresearch" / "ar5" / "cut" / "frozen.json"
sys.path.insert(0, str(HERE))
import cut_rules as CR  # noqa: E402
import rescue_freeze as RF  # noqa: E402  (pure: hashing only)

SCRIPTS = ["cut_rules.py", "cut_common.py", "cut_run.py"] + list(RF.SCRIPTS)
CONSTANTS = dict(
    top_n=CR.TOP_N,
    level=0,
    keep_fraction=CR.KEEP_FRACTION,
    tau=0.95,
    seed_test=0.80,
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
                rule=CR.RULE,
                constants=CONSTANTS,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no held file was hashed or opened before this file was committed; steps 8 "
                    "and 9 froze the library before step 9's phase B (key SeedRescueCMA, off by "
                    "default), so their frozen.json no longer agrees; this file freezes the "
                    "same script list and the library as it is now"
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
        return "cut frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "cut frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "cut frozen.json differs from its commit"
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
