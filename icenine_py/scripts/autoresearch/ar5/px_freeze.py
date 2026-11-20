"""AR5 step 27, phase A: freeze rule ``PX``, the gates, the predictions, the voxel lists and the
hashes of this step's scripts and of the library, before any q1, t1 or live run.

  uv run python scripts/autoresearch/ar5/px_freeze.py --write  # ar5/px/frozen.json
  uv run python scripts/autoresearch/ar5/px_freeze.py --check  # exit 2 unless it agrees

``frozen.json`` holds the rule text, the gates, the predictions, the sha256 of ``voxels.json``
(the q1 and t1 lists), of this step's scripts that define the arms and the checks (``OWN``), of
every script step 24's freeze lists (and so of every older step) and of EVERY LIBRARY FILE
(``icenine/*.py``, ``pyproject.toml``, ``uv.lock``): phase A changes no library file. ``--check``
also needs the file committed and unchanged. Step 24's own check (``bf_freeze``) is switched off in
the processes of this step after this one passes (``px_common``); no older hash is re-recorded and
no older file is edited.

Not in the list, on purpose: scripts that only summarise or judge (``px_score``, ``px_summary``).
Nothing of ``held2``, ``fine7`` or ``held7`` is hashed, opened or listed. ``held7``'s own hashes go
into ``ar5/px/held7_frozen.json``, committed after its build and before any method runs on it.
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "px"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import bf_freeze as BFF  # noqa: E402  (its hashing lists only; pure)
import px_rules as PR  # noqa: E402

OWN = [
    "px_rules.py",
    "px_common.py",
    "px_launch.py",
    "px_voxels.py",
    "px_q1.py",
    "px_t1.py",
    "px_run_all.sh",
    "px_link_cache.sh",
    "px_build_held7.py",
]
CONSTANTS = dict(
    rule=PR.RULE,
    new_max_pairs=PR.NEW_MAX_PAIRS,
    old_max_pairs=PR.OLD_MAX_PAIRS,
    seed_generator=PR.SEED_GENERATOR,
    n_per_sample=PR.N_PER_SAMPLE,
    q1_samples=list(PR.Q1_SAMPLES),
    t1_samples=list(PR.FULL_SAMPLES),
    block=PR.BLOCK,
    timing_rule=PR.TIMING_RULE,
    stored_totals=PR.STORED_TOTALS,
    indexer="index_voxel (the old indexer) in both arms; index_voxel_fast is not used",
    held7="built only after q1, c1, c2 and c3 pass; its hashes are committed in held7_frozen.json "
    "before any method runs on it; no truth of it is read before S8 and new J5SD have end markers",
    q1_orientation="the best fit by hit ratio of the stored seed step (as the sizing scripts)",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    cur = BFF.current()  # step 24's list (with every older step's) and the library
    files = dict(cur["files"])
    for s in sorted(OWN):
        files[f"scripts/autoresearch/ar5/{s}"] = sha(HERE / s)
    files["benchmarks/autoresearch/ar5/px/voxels.json"] = sha(BD / "voxels.json")
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
                    "no q1, t1 or live run was made by the tracked scripts before this file was "
                    "committed; the sizing measurements (400 stored voxels of mid and held) were "
                    "made before; nothing of held2, fine7 or held7 is used or read"
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
        return "px frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "px frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "px frozen.json differs from its commit"
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
