"""AR5 step 26, phase A: freeze rule ``IX``, the gates, the predictions, the exactness voxel lists
and the hashes of this step's scripts and of the library, before any exactness, timing or live run.

  uv run python scripts/autoresearch/ar5/ix_freeze.py --write  # ar5/ix/frozen.json
  uv run python scripts/autoresearch/ar5/ix_freeze.py --check  # exit 2 unless it agrees

``frozen.json`` holds the rule text, the gates, the predictions, the sha256 of ``voxels.json``
(the exactness lists), of this step's scripts that define the arms and the checks (``OWN``), of
every script step 24's freeze lists (and so of every older step) and of EVERY LIBRARY FILE
(``icenine/*.py``, ``pyproject.toml``, ``uv.lock``). ``--check`` also needs the file committed and
unchanged. The library file ``ray_indexing.py`` carries the new function by then: the freeze pins
the new function as it is judged. Step 24's own check (``bf_freeze``) is switched off in the
processes of this step after this one passes (``ix_common``); no older hash is re-recorded and no
older file is edited.

Not in the list, on purpose: scripts that only summarise or judge (``ix_summary``, ``ix_gates``).
Nothing of ``held2``, ``fine7`` or ``held7`` is hashed, opened or listed.
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ix"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import bf_freeze as BFF  # noqa: E402  (its hashing lists only; pure)
import ix_rules as IR  # noqa: E402

OWN = [
    "ix_rules.py",
    "ix_common.py",
    "ix_launch.py",
    "ix_voxels.py",
    "ix_exact.py",
    "ix_run_all.sh",
]
CONSTANTS = dict(
    rule=IR.RULE,
    seed_generator=IR.SEED_GENERATOR,
    n_per_sample=IR.N_PER_SAMPLE,
    full_samples=list(IR.FULL_SAMPLES),
    sampled=list(IR.SAMPLED),
    timing_samples=list(IR.TIMING_SAMPLES),
    block=IR.BLOCK,
    timing_rule=IR.TIMING_RULE,
    held7="not built, not read: tier 1 makes no accuracy claim",
    freeze_note="the freeze check pins the sha256 of every script and library file; step 24's "
    "own check is switched off in this process after this one passes",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    cur = BFF.current()  # step 24's list (with every older step's) and the library
    files = dict(cur["files"])
    for s in sorted(OWN):
        files[f"scripts/autoresearch/ar5/{s}"] = sha(HERE / s)
    files["benchmarks/autoresearch/ar5/ix/voxels.json"] = sha(BD / "voxels.json")
    return dict(files=files, library=dict(cur["library"]))


def write() -> None:
    cur = current()
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(
        json.dumps(
            dict(
                constants=CONSTANTS,
                predictions=IR.PREDICTIONS,
                gates=IR.GATES,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no exactness (e1), live (e2) or timed (w1, w2) run was made by the tracked "
                    "scripts before this file was committed; the phase 0 profile (200 mid "
                    "voxels, old indexer only, contended) and design measurements on mid "
                    "voxels were made before; gate z0 was NOT judged at the freeze (the "
                    "machine was not quiet; see the findings); nothing of held2, fine7 or "
                    "held7 is used"
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
        return "ix frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "ix frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "ix frozen.json differs from its commit"
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
