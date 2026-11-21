"""AR5 step 17, phase A step 3: freeze the rule family, the selection, the gate and the scripts that
define and record them before the recording is read, and before ``held``, ``held3`` (or ``held4``)
is read.

  uv run python scripts/autoresearch/ar5/fitskip_freeze.py --write  # ar5/fitskip/frozen.json
  uv run python scripts/autoresearch/ar5/fitskip_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text, the 27 rules, the selection and its tie rule, the gate, the
pass rules, the predictions, the stored counts, and the sha256 of the scripts that define, record
and judge the rule (this step's list below, step 16's, step 15's and the older lists they build on,
the AR3 scripts) and of every library file. ``--check`` also needs the file committed and
unchanged. No hash of any ``held``, ``held3``, ``held4`` or ``held2`` data file is taken (that
would read it). The summary, replay and selection scripts are written after this commit and are
not in the list: they apply ``fitskip_rules`` as frozen. Nothing is tuned after a recording.

This step's own freeze is checked FIRST for ``held``, ``held3`` and ``held4`` (``fitskip_common``);
only when it passes are the older steps' checks (8, 9, 12, 14, 15 and 16: they no longer agree, a
docstring line of the library moved their hashes) switched off, in the script layer of this process
only. No hash of any older step is re-recorded and no older frozen file is edited.

No import of the study's other modules except ``fitskip_rules`` (pure), ``high500_freeze``,
``polish_freeze`` and ``stallhit_freeze`` (hashing lists).
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "fitskip"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import fitskip_rules as FR  # noqa: E402
import high500_freeze as H5F  # noqa: E402  (pure: hashing only)
import polish_freeze as PF  # noqa: E402  (pure: hashing only)
import stallhit_freeze as SF  # noqa: E402  (pure: hashing only)

OWN = [
    "fitskip_rules.py",
    "fitskip_guard.py",
    "fitskip_hook.py",
    "fitskip_common.py",
    "fitskip_run.py",
    "fitskip_pass.py",
    "fitskip_p1.py",
    "fitskip_final.py",
    "fitskip_link_cache.sh",
]
# step 16's run, launch and pass scripts that define arm J5 (hashed, not edited)
H5_RUN = ["high500_run.py", "high500_launch.py", "high500_pass.py", "high500_pick.py"]
CONSTANTS = dict(
    ts=list(FR.TS),
    scopes=list(FR.SCOPES),
    ds=list(FR.DS),
    check_at=FR.CHECK_AT,
    accept=FR.ACCEPT,
    wrong_skip_max_share="1/1000 of the sample's fits",
    final_changed_max=FR.FINAL_CHANGED_MAX,
    rules=[r.name for r in FR.family()],
    tie_completion=FR.TIE_COMPLETION,
    j5_stored=FR.EXPECT,
    held4="build of held: 250 Voronoi grains on the 24,570-voxel grid, gradients 0 to 4 degrees,"
    " generator seed 4 (a new sample for a carry-over test, not an extra random seed)",
    live_forms="S and D live; C replayed only (the library cannot stop a CMA run after 20 evals)",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    names = OWN + H5_RUN + H5F.OWN + H5F.BAND_USED
    files = {f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in sorted(set(names))}
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
                rule=FR.RULE,
                selection=FR.SELECTION,
                constants=CONSTANTS,
                predictions=FR.PREDICTIONS,
                gate=FR.GATE,
                judge=FR.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no row of the recording was read by this step, and no log of held, held3 or "
                    "held4 was opened by it, before this file was committed; nothing of held2 is "
                    "used by this step; held4 does not exist yet"
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
        return "fitskip frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "fitskip frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "fitskip frozen.json differs from its commit"
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
