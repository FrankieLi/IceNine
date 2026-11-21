"""AR5 step 29, phase A: freeze rule ``RU``, the gates, the predictions and the hashes of this
step's scripts and of the library, before any run.

  uv run python scripts/autoresearch/ar5/ru_freeze.py --write  # ar5/ru/frozen.json
  uv run python scripts/autoresearch/ar5/ru_freeze.py --check  # exit 2 unless it agrees

``frozen.json`` holds the rule text, the gates, the predictions, the stop rules and the sha256 of
this step's scripts that define the arms and the checks (``OWN``), of every script step 28's freeze
lists (and so of every older step; ``pxf_freeze.current``) and of EVERY LIBRARY FILE
(``icenine/*.py``, ``pyproject.toml``, ``uv.lock``): phase A changes no library file. ``--check``
also needs the file committed and unchanged. No older hash is re-recorded and no older file is
edited.

Not in the list, on purpose: scripts that only summarise or judge (``ru_score``, ``ru_summary``).
Nothing of ``held2``, ``fine7`` or ``held8`` is hashed, opened or listed. ``held8``'s own hashes go
into ``ar5/ru/held8_frozen.json``, committed after its build and before any method runs on it.
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
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ru"
FROZEN = BD / "frozen.json"
sys.path.insert(0, str(HERE))
import pxf_freeze as PXFZ  # noqa: E402  (step 28's list, with every older step's; pure)
import ru_rules as RR  # noqa: E402

OWN = [
    "ru_rules.py",
    "ru_common.py",
    "ru_launch.py",
    "ru_samples.py",
    "ru_run_all.sh",
    "ru_build_held8.py",
]
CONSTANTS = dict(
    rule=RR.RULE,
    k_nearest=RR.K_NEAREST,
    sigma0_deg=RR.SIGMA0_DEG,
    max_evals=RR.MAX_EVALS,
    second_start_min_sep_deg=RR.SEP_DEG,
    accept_hit_ratio=RR.ACCEPT,
    timing_rule=RR.TIMING_RULE,
    old_totals=RR.OLD_TOTALS,
    old_total=RR.OLD_TOTAL,
    old_wall_s=RR.OLD_WALL_S,
    in_sample=RR.IN_SAMPLE,
    stop_rules=RR.STOP_RULES,
    held8="built only after r1, r2 and r3 pass; its hashes are committed in held8_frozen.json "
    "before any method runs on it; no truth of it is read before S8, the PXF J5SD and the new "
    "J5SD have end markers",
    run_order_held8="S8, then PXF J5SD, then new J5SD (one at a time)",
    caps_minutes=dict(S8=480, others=300),
    sample_data_rebuilt=dict(
        note="mid and held were missing from the main checkout and were put back before any run: "
        "mid's two truth mics copied from the 20 Nov rebuild (worktree agent-a2ce49fff826b9a52), "
        "held rebuilt with the stored builders (held_make_sample.py, render_full.py --tag held, "
        "held_stack_note.py --build); both checked against the stored records",
        mid_gradient_truth_sha256=(
            "3767c827bffc8eb966543ad7227a65dcd249aa10ad07cefae4da7775a00e4188"
        ),
        held_gradient_truth_sha256=(
            "613dd43d3fd9278d6db638346a41e36d4232969e4faef337784b8fdb276f4f99"
        ),
        held_stack_nonzero=10897893,
        held_render_spots_total=179718,
        held_tracked_files_unchanged="grainmap.npy, sample.json, held_check.json (git diff clean)",
    ),
    old_arm="step 28's PXF runs on mid and held, stored in the main checkout's ar5/pxf/cache "
    "(live_pxf/J5D_<s>, live_pxf/j5sd/J5S_<s>), read only, not rerun",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    cur = PXFZ.current()  # step 28's list (with every older step's) and the library
    files = dict(cur["files"])
    for s in sorted(OWN):
        files[f"scripts/autoresearch/ar5/{s}"] = sha(HERE / s)
    files["tests/test_ar5_ru.py"] = sha(ICE / "tests" / "test_ar5_ru.py")
    return dict(files=files, library=dict(cur["library"]))


def write() -> None:
    cur = current()
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(
        json.dumps(
            dict(
                constants=CONSTANTS,
                predictions=RR.PREDICTIONS,
                gates=RR.GATES,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "this is the latest freeze; the first one (gates, rule, predictions) is in "
                    "freeze_history.json with every re-freeze, its reason and what had already "
                    "run; the sizing probe (stored rows of step 28's mid and held runs) was made "
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
        return "ru frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "ru frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "ru frozen.json differs from its commit"
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
