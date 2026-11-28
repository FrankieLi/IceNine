"""AR5 step 8, Phase A step 3: freeze the rule and the files before ``held`` is read.

  uv run python scripts/autoresearch/ar5/order_freeze.py --write   # writes order/frozen.json
  uv run python scripts/autoresearch/ar5/order_freeze.py --check   # exit 2 unless it agrees

``frozen.json`` holds the rule text, the constants, and the sha256 of the scripts that define and
run the two orders and of every file of the library (``icenine/*.py``, ``pyproject.toml``,
``uv.lock``). ``--check`` also needs the file committed and unchanged (the plan's order rule: no
``held`` seed log, replay or count is read before the freeze commit). No import of the study's
other modules: this file is safe to import from the guards.
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
FROZEN = ICE / "benchmarks" / "autoresearch" / "ar5" / "order" / "frozen.json"
SCRIPTS = [
    "order_rules.py",
    "order_draws.py",
    "order_common.py",
    "order_run.py",
    "order_pass.py",
    "order_whatif.py",
    "mix_rules.py",
    "mix_common.py",
    "fit_rules.py",
    "fit_common.py",
    "fit_run.py",
    "fit_pass.py",
    "fit_part1_lib.py",
]
RULE = (
    "Interleaved order A: per candidate (strict-cost order, at most 5) a CMA fit (first step 0.05 "
    "degree, 200 evaluations), then the stock MC refinement if its hit ratio is under 0.95; the "
    "kept end of a candidate is the lower-cost of the two (a tie keeps the CMA end); the step "
    "stops "
    "at the first kept end with hit ratio >= 0.95; else the lowest cost of all kept ends. "
    "Sweep order Y: CMA fits on the candidates in order until a CMA end reaches 0.95 (taken, step "
    "stops); if none does, the MC refinement on the candidates in order, each kept end the lower "
    "cost of its two, stop at the first kept end >= 0.95; else the lowest cost of all kept ends. "
    "Both arms: run S11 (SeedIndexCMAFirst 200, BFSNeighborCMANoFallback 1), no N1 rule, aligned "
    "draws (voxel, kind of call, count at the voxel), pass with rule hit."
)
CONSTANTS = dict(tau=0.95, cma_evals=200, cma_sigma0_deg=0.05, candidates_refined=5, seed_test=0.80)


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
                rule=RULE,
                constants=CONSTANTS,
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
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
