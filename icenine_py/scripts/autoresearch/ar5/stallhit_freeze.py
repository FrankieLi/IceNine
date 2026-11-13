"""AR5 step 12, phase A step 3: freeze rule PH(3000) and the files before ``held2`` is built and
before ``held`` is read again.

  uv run python scripts/autoresearch/ar5/stallhit_freeze.py --write   # writes stallhit/frozen.json
  uv run python scripts/autoresearch/ar5/stallhit_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text, k, and the sha256 of the scripts that define, run and replay
the rule (this step's, step 11's list that it builds on, and the AR3 scripts that build and render
``held2``: ``held2`` is ``held`` with generator seed 2) and of every library file. ``--check`` also
needs the file committed and unchanged: no ``held2`` file exists before the freeze commit, and no
log, trace, record or map of ``held`` is opened by this step before it. No hash of any ``held`` or
``held2`` data file is taken (that would read it). Scripts that never touch the held samples (the
gate, the link and launcher scripts) are not in the list, so they may change after the freeze (the
what-if is in it: the replay imports it). Nothing is tuned after a run.

No import of the study's other modules except ``stallhit_rules`` (pure) and ``rescue_freeze``'s
list.
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
FROZEN = ICE / "benchmarks" / "autoresearch" / "ar5" / "stallhit" / "frozen.json"
sys.path.insert(0, str(HERE))
import rescue_freeze as RF  # noqa: E402  (pure: hashing only)
import stallhit_rules as SR  # noqa: E402

SCRIPTS = [
    "stallhit_rules.py",
    "stallhit_common.py",
    "stallhit_run.py",
    "stallhit_replay.py",
    "stallhit_whatif.py",
    "stallhit_make_held2.py",
] + list(RF.SCRIPTS)
AR3_SCRIPTS = [
    "ar3_common.py",
    "gradient_common.py",
    "reseed_common.py",
    "held_make_sample.py",
    "held_sample_lib.py",
    "held_stack_note.py",
]
CONSTANTS = dict(
    k=SR.K,
    stage="seed step MC refinement only (MC stage then variance stage, one counter)",
    tau=SR.TAU,
    seed_test=SR.SEED_TEST,
    candidates=5,
    rule_g1="unchanged (rescue_rules.py)",
    stall_check="one hit-ratio evaluation of the running best; at 0.95 or more stop, else go on",
    held2="build of held: 250 Voronoi grains, gradients 0 to 4 degrees, generator seed 2",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    files = {f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in SCRIPTS}
    files.update(
        {f"scripts/autoresearch/ar3/{s}": sha(HERE.parent / "ar3" / s) for s in AR3_SCRIPTS}
    )
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
                rule=SR.RULE,
                constants=CONSTANTS,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no held or held2 file was hashed or opened before this file was committed; "
                    "held2 did not exist; steps 8, 9 and 11 froze the library before step 9's "
                    "phase B (key SeedRescueCMA, off by default), so their frozen.json no longer "
                    "agree; this file freezes their script list and the library as it is now"
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
        return "stallhit frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "stallhit frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "stallhit frozen.json differs from its commit"
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
