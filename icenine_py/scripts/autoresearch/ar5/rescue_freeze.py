"""AR5 step 9, phase A step 3: freeze rule G1 and the files before ``held`` is read.

  uv run python scripts/autoresearch/ar5/rescue_freeze.py --write   # writes rescue/frozen.json
  uv run python scripts/autoresearch/ar5/rescue_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text, the constants and the sha256 of the scripts that define and run
rule G1 (and of the step 8 scripts it builds on) and of every library file. ``--check`` also needs
the file committed and unchanged (the plan: no ``held`` seed log, fallback end, map or count is
opened before the freeze commit). No hash of any ``held`` file is taken (that would read it). No
import of the study's other modules: this file is safe to import from the guards.
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
FROZEN = ICE / "benchmarks" / "autoresearch" / "ar5" / "rescue" / "frozen.json"
SCRIPTS = [
    "rescue_rules.py",
    "rescue_common.py",
    "rescue_run.py",
    "order_rules.py",
    "order_draws.py",
    "order_common.py",
    "order_run.py",
    "order_pass.py",
    "order_freeze.py",
    "order_logs.py",
    "mix_rules.py",
    "mix_common.py",
    "fit_rules.py",
    "fit_common.py",
    "fit_run.py",
    "fit_pass.py",
    "fit_part2_stats.py",
]
CONSTANTS = dict(
    tau=0.95, sigma0_deg=0.25, evals=1000, max_tries=5, seed_test=0.80, rescue_rng="default_rng(0)"
)


def _rule() -> str:
    sys.path.insert(0, str(HERE))
    import rescue_rules as RR

    return RR.RULE


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current() -> Dict[str, Dict[str, str]]:
    files = {f"scripts/autoresearch/ar5/{s}": sha(HERE / s) for s in SCRIPTS}
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
                rule=_rule(),
                constants=CONSTANTS,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note="no held file was hashed or opened before this file was committed",
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
        return "frozen.json does not exist"
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log:
        return "frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "frozen.json differs from its commit"
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
