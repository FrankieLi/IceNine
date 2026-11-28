"""AR5 step 21, phase A step 3: freeze the screen rule, the gate, the decision rules, the stored
counts and the scripts that define and set up arm J5S, before any J5S pass or repair runs from the
tracked scripts and before any sample is read by them.

  uv run python scripts/autoresearch/ar5/screen_freeze.py --write  # ar5/screen/frozen.json
  uv run python scripts/autoresearch/ar5/screen_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text and constants (``screen_rules``), the gate, the decision rules,
the predictions, the stored counts, the held5 spec, and the sha256 of the scripts that define and
set up the arm (this step's ``OWN`` list, step 20's ``cmarepair_rules.py`` and step 19's
``repair_rules.py``, whose repair and flag this step imports), of every older script that a run of
this step imports through ``learnskip_freeze`` (step 18's list, which holds steps 14 to 17's and the
AR3 scripts) and of every library file. ``--check`` also needs the file committed and unchanged. No
hash of any ``held5`` or ``held2`` data file is taken (that would read it).

NOT in the list, on purpose: the scripts that apply the frozen rule (``screen_pass``,
``screen_run``, ``screen_gate``, ``screen_build_held5``, ``screen_run_all.sh``, the summary).
Everything that decides a voxel (the screen try, the wide flag, the candidates, the order, the keep
rule, the generators, the optimizer settings, the counting) is in ``screen_rules`` and the two
modules it imports, which are in the list.

This step's own freeze is checked FIRST (``screen_common``); only when it passes are the older
steps' checks (8, 9, 12, 14, 15 and 16: they no longer agree, the library moved their hashes)
switched off, in the script layer of this process only: the precedent of ``cmarepair_common``,
``stallhit_common.py`` and ``polish_replace_freeze.py``. Step 20's own check
(``cmarepair_freeze``) is not used at all: the library changed after its freeze. No hash of any
older step is re-recorded and no older script and no older frozen file is edited.
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "screen"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import learnskip_freeze as LFZ  # noqa: E402  (its hashing lists only; pure)
import screen_rules as SR  # noqa: E402

OWN = [
    "screen_rules.py",
    "screen_guard.py",
    "screen_samples.py",
    "screen_common.py",
    "screen_data.py",
    "screen_launch.py",
    "screen_link_cache.sh",
    "cmarepair_rules.py",  # step 20's repair, imported (not copied)
    "repair_rules.py",  # step 19's flag, imported by step 20's
]
CONSTANTS = dict(
    rule=SR.RULE,
    arm=SR.ARM,
    screen_max_evals=SR.SCREEN_MAX_EVALS,
    evals_per_screen=SR.EVALS_PER_SCREEN,
    pass_sigma0_deg=SR.PASS_SIGMA0_DEG,
    full_max_evals=SR.FULL_MAX_EVALS,
    ratio_max=SR.RATIO_MAX,
    wide_deg=SR.WIDE_DEG,
    margin_max=SR.CR.MARGIN_MAX,
    far_deg=SR.CR.FAR_DEG,
    wrong_deg=SR.WRONG_DEG,
    repair_optimizer=SR.CR.REPAIR_OPTIMIZER,
    repair_cma_sigma0_deg=SR.CR.CMA_SIGMA0_DEG,
    repair_cma_max_evals=SR.REPAIR_MAX_EVALS,
    repair_evals_per_fit=SR.EVALS_PER_FIT,
    repair_per_fit_generator="default_rng([0, v, c + 1]), c = -1 for the voxel's own end",
    pass_generator="order_draws.make_rng(v, 'pass', k), k = tries made at the voxel so far",
    neighbours="radius 2 side lengths (mic.get_neighbors), the pass's own relation",
    totals_j5_j5c_s8=SR.TOTALS,
    g0_j5c_repair_evals=SR.G0_J5C_EVALS,
    p1=dict(
        pass_evals=SR.P1_PASS_EVALS,
        flagged=SR.P1_FLAGGED,
        repair_evals=SR.P1_REPAIR_EVALS,
        total=SR.P1_TOTAL,
        fixed=SR.P1_FIXED,
        broken=SR.P1_BROKEN,
    ),
    library_note="Reconstructor.local_optimization refuses cma_max_evals < 2; the script layer "
    "sets the budget 1 on the optimizer the library makes (CMAOptimizer takes 1)",
    held5="build of held: 250 Voronoi grains on the 24,570-voxel grid, gradients 0 to 4 degrees, "
    "generator seed 5 (a new sample for a clean test, not an extra random seed), clean images, "
    "Q-max 8",
    samples_disclosed=list(SR.SAMPLES),
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
                predictions=SR.PREDICTIONS,
                gate=SR.GATE,
                judge=SR.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no J5S pass or repair was run from the tracked scripts, no flag was "
                    "computed by this step's scripts, and nothing of held5 was read by it, before "
                    "this file was committed; nothing of held2 is used by this step; held5 does "
                    "not exist yet"
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
        return "screen frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "screen frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "screen frozen.json differs from its commit"
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
