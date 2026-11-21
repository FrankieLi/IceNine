"""AR5 step 28, phase A: freeze rule ``PXF``, the gates, the predictions and the hashes of this
step's scripts and of the library, before any live run.

  uv run python scripts/autoresearch/ar5/pxf_freeze.py --write  # ar5/pxf/frozen.json
  uv run python scripts/autoresearch/ar5/pxf_freeze.py --check  # exit 2 unless it agrees

``frozen.json`` holds the rule text, the gates, the predictions and the sha256 of this step's
scripts that define the arms and the checks (``OWN``), of every script step 27's freeze lists (and
so of every older step; ``px_freeze.current``) and of EVERY LIBRARY FILE (``icenine/*.py``,
``pyproject.toml``, ``uv.lock``): phase A changes no library file. ``--check`` also needs the file
committed and unchanged. No older hash is re-recorded and no older file is edited.

Not in the list, on purpose: scripts that only summarise or judge (``pxf_score``,
``pxf_summary``) and ``px_link_cache.sh`` (step 27's, listed by its freeze). Nothing of ``held2``,
``fine7`` or ``held7`` is hashed, opened or listed. ``held7``'s own hashes go into
``ar5/pxf/held7_frozen.json``, committed after its build and before any method runs on it.
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "pxf"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import px_freeze as PXFZ  # noqa: E402  (step 27's list, with every older step's; pure)
import pxf_rules as PR  # noqa: E402

OWN = [
    "pxf_rules.py",
    "pxf_common.py",
    "pxf_launch.py",
    "pxf_run_all.sh",
    "pxf_build_held7.py",
]
CONSTANTS = dict(
    rule=PR.RULE,
    first_max_pairs=PR.FIRST_MAX_PAIRS,
    full_max_pairs=PR.FULL_MAX_PAIRS,
    seed_test=PR.TEST,
    timing_rule=PR.TIMING_RULE,
    stored_totals=PR.STORED_TOTALS,
    old_wall_s=PR.OLD_WALL_S,
    in_sample=PR.IN_SAMPLE,
    extra_run=PR.EXTRA_RUN,
    indexer="index_voxel (the old indexer); index_voxel_fast is not used",
    held7="built only after c1, c2 and c3 pass; its hashes are committed in held7_frozen.json "
    "before any method runs on it; no truth of it is read before S8, the new J5SD and the plain "
    "J5SD have end markers",
    caps_minutes=dict(S8=480, others=300),
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    cur = PXFZ.current()  # step 27's list (with every older step's) and the library
    files = dict(cur["files"])
    for s in sorted(OWN):
        files[f"scripts/autoresearch/ar5/{s}"] = sha(HERE / s)
    files["tests/test_ar5_pxf.py"] = sha(ICE / "tests" / "test_ar5_pxf.py")
    return dict(files=files, library=dict(cur["library"]))


def write() -> None:
    cur = current()
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(
        json.dumps(
            dict(
                constants=CONSTANTS,
                predictions=PR.PREDICTIONS,
                gates=PR.GATES,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "this is the latest freeze; the first one (gates, rule, predictions) is in "
                    "freeze_history.json with every re-freeze, its reason and what had already "
                    "run; the sizing (stored rows of step 27's mid and held runs) was made "
                    "before the first freeze; nothing of held2 or fine7 is used or read"
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
        return "pxf frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "pxf frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "pxf frozen.json differs from its commit"
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
