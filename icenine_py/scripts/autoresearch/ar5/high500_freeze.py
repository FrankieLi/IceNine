"""AR5 step 16, phase A step 3: freeze rule ``PH(500)``, the gate and the scripts that define them
before any script runs on the stored rows, and before ``held`` or ``held3`` is read.

  uv run python scripts/autoresearch/ar5/high500_freeze.py --write  # writes ar5/high500/frozen.json
  uv run python scripts/autoresearch/ar5/high500_freeze.py --check  # exit 2 unless it agrees

``frozen.json`` holds the rule text, the form, the plan's gate and pass rules
(``high500_rules.GATE`` and ``JUDGE``), and the sha256 of the scripts that define, price and judge
the rule (this step's list below, step 15's rule and sample scripts and step 12's and step 14's
lists that they build on, the AR3 scripts) and of every library file. ``--check`` also needs the
file committed and unchanged. No hash of any ``held``, ``held3`` or ``held2`` data file is taken
(that would read it). The whole-layer run, launch and summary scripts (``high500_run*``,
``high500_make_held3``, ``high500_pass``, ``high500_launch``) are written after this commit and are
not in the list: they apply ``high500_rules`` and ``high500_pick`` as frozen. Nothing is tuned
after a run.

This step's own freeze is checked FIRST for ``held`` and ``held3`` (``high500_common``); only
then are older steps' stale checks (steps 8, 9 and 12: they hash the library before step 12's
phase B added the opt-in key ``SeedIndexMCStallHit``) replaced, in this layer only. Step 14's
freeze check stays active where it still passes. No older hash is re-recorded and no older frozen
file is edited.

No import of the study's other modules except ``high500_rules`` (pure), ``stallhit_freeze`` and
``polish_freeze`` (hashing lists).
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "high500"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import high500_rules as HR  # noqa: E402
import polish_freeze as PF  # noqa: E402  (pure: hashing only)
import stallhit_freeze as SF  # noqa: E402  (pure: hashing only)

OWN = [
    "high500_rules.py",
    "high500_guard.py",
    "high500_stored.py",
    "high500_common.py",
    "high500_pick.py",
    "high500_paths.py",
    "high500_gate.py",
    "high500_link_cache.sh",
]
# step 15's pure rule and sample scripts that this step's layer imports (not edited, only hashed)
BAND_USED = ["band_rules.py", "band_samples.py", "band_guard.py"]
CONSTANTS = dict(
    tau=HR.TAU,
    k=HR.K,
    ph_k_replaced=HR.PH_K,
    seed_test=HR.SEED_TEST,
    g2_form="G2b (step 14's frozen form; polish_rules.py)",
    candidates=5,
    low_side="none",
    held3="build of held: 250 Voronoi grains, gradients 0 to 4 degrees, generator seed 3",
    j_totals={s: v["total"] for s, v in HR.EXPECT.items()},
    j_mc={s: v["mc"] for s, v in HR.EXPECT.items()},
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    files = {f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in OWN + BAND_USED}
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
                rule=HR.RULE,
                constants=CONSTANTS,
                predictions_1=HR.PRED1,
                gate=HR.GATE,
                judge=HR.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no stall hit ratio of any recording was read by this step, and no log of "
                    "held or held3 was opened by it, before this file was committed; nothing of "
                    "held2 is used by this step"
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
        return "high500 frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "high500 frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "high500 frozen.json differs from its commit"
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
