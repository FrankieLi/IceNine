"""AR5 step 18, phase A step 3: freeze the model, the features, the threshold procedure, the gate,
the pass rules and the scripts that define and run them, before any model is trained and before
``held`` , ``held3`` (or ``held4``) is read.

  uv run python scripts/autoresearch/ar5/learnskip_freeze.py --write  # ar5/learnskip/frozen.json
  uv run python scripts/autoresearch/ar5/learnskip_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the model text, the 16 features, the threshold procedure, the gate, the pass
rules, the predictions, the stored counts, and the sha256 of the scripts that define and run the
rule (this step's list below, step 17's, step 16's, step 15's and the older lists they build on,
the AR3 scripts) and of every library file. ``--check`` also needs the file committed and
unchanged. No hash of any ``held``, ``held3``, ``held4`` or ``held2`` data file is taken (that would
read it). The replay, gate, export and summary scripts are written after this commit and are not in
the list: they apply the frozen modules (``learnskip_rules``, ``_features``, ``_model``). Nothing is
tuned after a recording is read.

This step's own freeze is checked FIRST for ``held``, ``held3`` and ``held4``
(``learnskip_common``); only when it passes are the older steps' checks (8, 9, 12, 14, 15 and 16:
they no longer agree, docstring lines of the library moved their hashes) switched off, in the
script layer of this process only: the precedent of ``stallhit_common.py`` and
``polish_replace_freeze.py``. Step 17's own freeze (``ar5/fitskip/frozen.json``) is stale for the
same reason and is not used by this step. No hash of any older step is re-recorded and no older
frozen file is edited.

No import of the study's other modules except ``learnskip_rules`` and ``learnskip_features``
(pure), ``fitskip_freeze`` (its lists), ``high500_freeze``, ``polish_freeze`` and
``stallhit_freeze`` (hashing lists).
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "learnskip"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import fitskip_freeze as FFZ  # noqa: E402  (its lists only)
import high500_freeze as H5F  # noqa: E402  (pure: hashing only)
import learnskip_features as LFE  # noqa: E402
import learnskip_rules as LR  # noqa: E402
import polish_freeze as PF  # noqa: E402  (pure: hashing only)
import stallhit_freeze as SF  # noqa: E402  (pure: hashing only)

OWN = [
    "learnskip_rules.py",
    "learnskip_features.py",
    "learnskip_model.py",
    "learnskip_guard.py",
    "learnskip_samples.py",
    "learnskip_common.py",
    "learnskip_hook.py",
    "learnskip_run.py",
    "learnskip_pass.py",
    "learnskip_launch.py",
]
CONSTANTS = dict(
    model=LR.MODEL_SPEC,
    features=list(LFE.FEATURES),
    n_features=LFE.N_FEATURES,
    distance_cap_deg=LFE.DIST_CAP,
    threshold_procedure="quantiles at 0, 0.1%, ..., 80% of all out-of-fold scores; the largest sum"
    " of net shares over the training samples among those where every training sample has wrong"
    " skips at most 0.1% of its fits and changed final answers at most 2; none: no skip; a tie"
    " takes the smaller threshold; skip when the score is under the threshold",
    wrong_skip_max_share="1/1000 of the sample's fits",
    changed_final_max=LR.FC_MAX,
    j5_stored=LR.EXPECT,
    fits=LR.FITS,
    plan_table=LR.PLAN_TABLE,
    plan_gate_model_mid_threshold=LR.SIZING_GATE_MODEL_MID_THRESHOLD,
    plan_final_threshold=LR.SIZING_FINAL_THRESHOLD,
    plan_final_rounds=LR.SIZING_FINAL_ROUNDS,
    held4="build of held: 250 Voronoi grains on the 24,570-voxel grid, gradients 0 to 4 degrees,"
    " generator seed 4 (a new sample for a carry-over test, not an extra random seed), clean"
    " images, Q-max 8",
    live_form="one counted check of the start per fit; a skipped fit is rejected, takes and"
    " discards its aligned draw, has no fit end; a model call counts zero evaluations",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    names = OWN + FFZ.OWN + FFZ.H5_RUN + H5F.OWN + H5F.BAND_USED
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
                constants=CONSTANTS,
                predictions=LR.PREDICTIONS,
                gate=LR.GATE,
                judge=LR.JUDGE,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no model was trained by this step, no row of the recording was read by it, "
                    "and no log of held, held3 or held4 was opened by it, before this file was "
                    "committed; nothing of held2 is used by this step; held4 does not exist yet"
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
        return "learnskip frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "learnskip frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "learnskip frozen.json differs from its commit"
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
