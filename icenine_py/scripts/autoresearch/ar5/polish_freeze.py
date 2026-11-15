"""AR5 step 14, phase A step 4: freeze rule ``G2`` (the chosen form and values) before ``held2``'s
rescue logs are opened.

  uv run python scripts/autoresearch/ar5/polish_freeze.py --write FORM   # writes polish/frozen.json
  uv run python scripts/autoresearch/ar5/polish_freeze.py --check         # exit 2 unless it agrees

``frozen.json`` holds the rule text, both forms' values, the chosen form (chosen by the fixed rule
of the plan from ``polish/gate1.json``, whose sha256 is in it), the sha256 of the scripts that
define, run and judge the rule (this step's, and the lists of step 12 and of steps 8, 9, 11 that
it builds on) and of every library file. ``--check`` also needs the file committed and unchanged.
No hash of any ``held2`` data file is taken (that would read it). The whole-layer run and summary
scripts (``polish_run*``, ``polish_summary``) are written after this commit and are not in the
list: they apply ``polish_rules`` and ``polish_common`` as frozen. Nothing is tuned after a run.

No import of the study's other modules except ``polish_rules`` (pure), ``stallhit_freeze`` and
``rescue_freeze`` (pure: hashing only).
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
PO = ICE / "benchmarks" / "autoresearch" / "ar5" / "polish"
FROZEN = PO / "frozen.json"
GATE1 = PO / "gate1.json"
sys.path.insert(0, str(HERE))
import polish_rules as PR  # noqa: E402
import stallhit_freeze as SF  # noqa: E402

OWN = [
    "polish_rules.py",
    "polish_guard.py",
    "polish_common.py",
    "polish_cases.py",
    "polish_probe.py",
    "polish_gate.py",
    "polish_link_cache.sh",
]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    files = {f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in OWN}
    prior = sorted(set(SF.SCRIPTS))
    files.update({f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in prior})
    files.update(
        {f"scripts/autoresearch/ar3/{s}": sha(HERE.parent / "ar3" / s) for s in SF.AR3_SCRIPTS}
    )
    lib = {str(p.relative_to(ICE)): sha(p) for p in sorted((ICE / "icenine").glob("*.py"))}
    lib["pyproject.toml"] = sha(ICE / "pyproject.toml")
    lib["uv.lock"] = sha(ICE / "uv.lock")
    return dict(files=files, library=lib)


def write(form: str) -> None:
    if form not in PR.FORMS:
        raise SystemExit(f"unknown form {form}")
    g = json.loads(GATE1.read_text())
    if g.get("chosen") != form:
        raise SystemExit(f"gate1.json chose {g.get('chosen')!r}, not {form!r}: refusing")
    cur = current()
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(
        json.dumps(
            dict(
                rule=PR.RULE,
                chosen=form,
                forms=PR.FORMS,
                constants=dict(
                    tau=PR.TAU,
                    rng=f"default_rng({PR.RNG_SEED})",
                    max_tries=PR.MAX_TRIES,
                    case_bound=PR.case_bound(form),
                    seed_test=0.80,
                ),
                gate1_sha256=sha(GATE1),
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note=(
                    "no held2 rescue log, record or map was opened before this file was "
                    "committed; nothing under any passes_held2 folder was opened, listed or "
                    "statted; the form was chosen by the plan's fixed rule from gate1.json"
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
        return "polish frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "polish frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "polish frozen.json differs from its commit"
    saved = json.loads(FROZEN.read_text())
    cur = current()
    for part in ("files", "library"):
        for k, v in cur[part].items():
            if saved[part].get(k) != v:
                return f"{k} differs from the frozen hash"
        for k in saved[part]:
            if k not in cur[part]:
                return f"{k} is frozen but missing"
    if GATE1.exists() and saved.get("gate1_sha256") != sha(GATE1):
        return "gate1.json differs from the frozen hash"
    return ""


def check_or_exit() -> None:
    p = problems()
    if p:
        print("freeze check failed:", p, file=sys.stderr)
        raise SystemExit(2)


def frozen_form() -> str:
    return str(json.loads(FROZEN.read_text())["chosen"])


if __name__ == "__main__":
    if "--write" in sys.argv:
        write(sys.argv[sys.argv.index("--write") + 1])
    elif "--check" in sys.argv:
        check_or_exit()
        print("freeze check passed")
    else:
        raise SystemExit("use --write FORM or --check")
