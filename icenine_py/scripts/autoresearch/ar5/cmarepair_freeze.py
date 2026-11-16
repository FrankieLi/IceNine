"""AR5 step 20, phase A step 3: freeze the repair rule, the gate, the decision rules, the stored
counts and the scripts that define and set up the repair, before any repair fit runs from the
tracked scripts and before ``held``, ``held3`` (or ``held4``) is read.

  uv run python scripts/autoresearch/ar5/cmarepair_freeze.py --write  # ar5/cmarepair/frozen.json
  uv run python scripts/autoresearch/ar5/cmarepair_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text and constants (``cmarepair_rules``), the gate, the decision
rules, the predictions, the stored counts, the held4 spec, and the sha256 of the scripts that
define and set up the repair (this step's ``OWN`` list, and step 19's ``repair_rules.py``, whose
flag and keep rule this step imports), of every older script that a run of this step imports
through ``learnskip_freeze`` (step 18's list, which holds steps 14 to 17's and the AR3 scripts) and
of every library file. ``--check`` also needs the file committed and unchanged. No hash of any
``held``, ``held3``, ``held4`` or ``held2`` data file is taken (that would read it).

NOT in the list, on purpose: the scripts that apply the frozen rule (``cmarepair_run`` builds the
reconstructor and calls the frozen ``cmarepair_rules.repair`` with the CMA fit; ``cmarepair_flags``,
``cmarepair_gate``, ``cmarepair_build_held4``, ``cmarepair_run_all.sh``, the summary). Everything
that decides a voxel (the flag, the candidates, the order, the keep rule, the generators, the
optimizer settings, the counting) is in ``cmarepair_rules`` and the step 19 module it imports,
which are in the list.

This step's own freeze is checked FIRST for ``held``, ``held3`` and ``held4``
(``cmarepair_common``); only when it passes are the older steps' checks (8, 9, 12, 14, 15, 16 and
18: they no longer agree, docstring lines of the library moved their hashes) switched off, in the
script layer of this process only: the precedent of ``stallhit_common.py`` and
``polish_replace_freeze.py``. No hash of any older step is re-recorded and no older frozen file is
edited.
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "cmarepair"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import cmarepair_rules as CR  # noqa: E402
import learnskip_freeze as LFZ  # noqa: E402  (its hashing lists only; pure)

OWN = [
    "cmarepair_rules.py",
    "cmarepair_guard.py",
    "cmarepair_data.py",
    "cmarepair_common.py",
    "cmarepair_launch.py",
    "cmarepair_link_cache.sh",
    "repair_rules.py",  # step 19's frozen flag, imported (not copied)
]
CONSTANTS = dict(
    rule=CR.RULE,
    ratio_max=CR.RATIO_MAX,
    margin_max=CR.MARGIN_MAX,
    far_deg=CR.FAR_DEG,
    cand_deg=CR.CAND_DEG,
    wrong_deg=CR.WRONG_DEG,
    repair_seed=CR.REPAIR_SEED,
    repair_optimizer=CR.REPAIR_OPTIMIZER,
    cma_sigma0_deg=CR.CMA_SIGMA0_DEG,
    cma_max_evals=CR.CMA_MAX_EVALS,
    evals_per_fit=CR.EVALS_PER_FIT,
    per_fit_generator="default_rng([0, v, c + 1]), c = -1 for the voxel's own end",
    j5_total=CR.J5_TOTAL,
    j5_wrong=CR.J5_WRONG,
    flagged=CR.FLAGGED,
    fits=CR.FITS,
    j5l=CR.J5L,
    p1_evals=CR.P1_EVALS,
    p1_fixed=CR.P1_FIXED,
    p1_broken=CR.P1_BROKEN,
    p1_j5l=CR.P1_J5L,
    neighbours="radius 2 side lengths (mic.get_neighbors), the pass's own relation",
    shares="repair evaluations over the sample's J5 run with the pass (stats_<s>.json of step 16,"
    " stages.total)",
    held4="build of held: 250 Voronoi grains on the 24,570-voxel grid, gradients 0 to 4 degrees,"
    " generator seed 4 (a new sample for a carry-over test, not an extra random seed), clean"
    " images, Q-max 8",
    samples_disclosed=list(CR.SAMPLES),
    caps_min=dict(S8=480, J5=300, build=300),
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    cur = LFZ.current()  # step 18's script list (with steps 14 to 17's) and the library
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
                predictions=CR.PREDICTIONS,
                gate=CR.GATE,
                judge=CR.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no repair fit was run from the tracked scripts, no flag was computed by "
                    "this step's scripts, and no log of held, held3 or held4 was opened by it, "
                    "before this file was committed; nothing of held2 is used by this step; held4 "
                    "does not exist yet"
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
        return "cmarepair frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "cmarepair frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "cmarepair frozen.json differs from its commit"
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
