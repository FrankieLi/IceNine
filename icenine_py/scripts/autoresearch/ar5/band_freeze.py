"""AR5 step 15, phase A step 4: freeze rule ``PB(k)`` (both forms), the gate and the scripts that
define them before any stall hit ratio is read, and before ``held`` or ``held3`` is read.

  uv run python scripts/autoresearch/ar5/band_freeze.py --write   # writes ar5/band/frozen.json
  uv run python scripts/autoresearch/ar5/band_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text, both forms, the plan's gate and pass rules (``band_rules.GATE``
and ``JUDGE``), and the sha256 of the scripts that define, record and judge the rule (this step's
list below, step 12's and step 14's lists that it builds on, the AR3 scripts) and of every library
file. ``--check`` also needs the file committed and unchanged. No hash of any ``held``, ``held3``
or ``held2`` data file is taken (that would read it). The whole-layer run and summary scripts
(``band_run*``, ``band_summary``) are written after this commit and are not in the list: they
apply ``band_rules`` and ``band_pick`` as frozen. Nothing is tuned after a run.

This step's own freeze is checked FIRST for ``held`` and ``held3`` (``band_common``); only when it
passes are the older steps' checks switched off there (they cannot agree any more: step 12's
freeze hashed the library before its phase B added the opt-in key ``SeedIndexMCStallHit``). No
hash of any older step is re-recorded and no older frozen file is edited.

No import of the study's other modules except ``band_rules`` (pure), ``stallhit_freeze`` (hashing
lists) and ``polish_freeze``'s file list (``OWN``).
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "band"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import band_rules as BR  # noqa: E402
import polish_freeze as PF  # noqa: E402  (pure: hashing only)
import stallhit_freeze as SF  # noqa: E402  (pure: hashing only)

OWN = [
    "band_rules.py",
    "band_guard.py",
    "band_samples.py",
    "band_common.py",
    "band_pick.py",
    "band_record.py",
    "band_paths.py",
    "band_gate.py",
    "band_link_cache.sh",
]
CONSTANTS = dict(
    tau=BR.TAU,
    low=BR.LOW,
    forms=BR.FORMS,
    record_ks=list(BR.RECORD_KS),
    ph_k=BR.PH_K,
    g2_form="G2b (step 14's frozen form; polish_rules.py)",
    candidates=5,
    seed_test=BR.LOW,
    held3="build of held: 250 Voronoi grains, gradients 0 to 4 degrees, generator seed 3",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    files = {f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in OWN}
    prior = sorted(set(SF.SCRIPTS) | set(PF.OWN) | {"polish_rules.py"})
    files.update({f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in prior})
    files.update(
        {f"scripts/autoresearch/ar3/{s}": sha(HERE.parent / "ar3" / s) for s in SF.AR3_SCRIPTS}
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
                rule=BR.RULE,
                constants=CONSTANTS,
                gate=BR.GATE,
                judge=BR.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no stall hit ratio of any recording was read, and no log of held or held3 "
                    "was opened by this step, before this file was committed; nothing of held2 "
                    "is used by this step"
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
        return "band frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "band frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "band frozen.json differs from its commit"
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
