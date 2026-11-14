"""AR5 step 2 (task ar5b): shared set-up of the fit study.

Import this module first in every ``fit_*.py`` script, with the environment variable
``RESEED_SAMPLE`` set to ``big``, ``grad`` or ``mid`` (the launcher ``fit_bg.sh`` does it). It
adds the AR3 script folder to the path, sends the cache to ``benchmarks/autoresearch/ar5/fit/
cache`` (untracked) and then imports step 8's ``reseed_common``.

Stored runs (untracked, in sibling worktrees; read, never written): run ``S8`` and its pass on the
three samples. ``mid`` is read only by scripts run after the pick file is committed (the order
rule of the plan); ``fit_guard.require_pick`` enforces it in code.

The four seed-refinement arms of Part 1 (``refine_arm``): ``M`` (the stock MC refinement of the
library, as run ``S8`` calls it), ``C60`` and ``C200`` (one CMA fit, first step 0.05 degree, at
most 60 or 200 evaluations, final evaluation included) and ``C0`` (no refinement: the candidate as
costed).
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
AR3 = HERE.parent / "ar3"
sys.path.insert(0, str(AR3))
sys.path.insert(0, str(HERE))

ICE = HERE.parents[2]  # icenine_py
FIT = ICE / "benchmarks" / "autoresearch" / "ar5" / "fit"  # tracked results (small)
CACHE_DIR = FIT / "cache"  # untracked
os.environ.setdefault("RESEED_CACHE", str(CACHE_DIR))
CACHE_DIR.mkdir(parents=True, exist_ok=True)
(CACHE_DIR / "logs").mkdir(parents=True, exist_ok=True)

PICK_FILE = FIT / "pick.json"


def pick_commit() -> str:
    """The hash of the commit that holds the pick file ('' when it is not committed)."""
    import subprocess

    out = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(PICK_FILE)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


# The order rule of the plan: no script reads mid before the pick file is committed.
if os.environ.get("RESEED_SAMPLE") == "mid" and not pick_commit():
    raise SystemExit("order rule: mid is not read before pick.json is committed")

import numpy as np  # noqa: E402
import reseed_common as RC  # noqa: E402

SAMPLE = RC.SAMPLE
CACHE = RC.CACHE
ARMS = ("M", "C60", "C200", "C0")
CMA_SIGMA0_DEG = 0.05
CMA_EVALS = {"C60": 60, "C200": 200}
N_COST = 20  # candidates costed by S8's seed step
N_REFINE = 5  # candidates refined at most

_WT = RC.WORKTREES
_RESEED = _WT / "agent-a635b2f904279bdc9/icenine_py/benchmarks/autoresearch/ar3/reseed/cache"
_START = _WT / "agent-a9c88b610fb10e1b6/icenine_py/benchmarks/autoresearch/ar3/start/cache"
STORED_S8: Dict[str, Path] = {
    "big": _RESEED / "S8_big",
    "grad": _RESEED / "S8_grad",
    "mid": _START / "S8_mid",
}
STORED_PASS: Dict[str, Path] = {
    "big": _RESEED / "passes_big",
    "grad": _RESEED / "passes_grad",
    "mid": _START / "passes_mid",
}


def refine_arm(
    arm: str,
    rec: Any,
    orig_lo: Any,
    R: "np.ndarray",
    cost: float,
    vertices: Any,
    phase: int,
) -> Tuple["np.ndarray", float, int]:
    """One seed refinement of the candidate ``R`` (strict cost ``cost``). Returns (end orientation,
    end cost, evaluations of the refinement). ``orig_lo`` is the reconstructor's unwrapped
    ``local_optimization`` (the run's ``Wrap`` replaces the instance attribute)."""
    import numpy as np
    from icenine.orientation_search import SearchCandidate

    if arm == "M":
        r = rec.refine_from_candidates(
            [SearchCandidate(orientation=R.copy(), cost=float(cost))],
            vertices,
            phase,
            rng=np.random.default_rng(0),
        )
        return np.asarray(r.orientation, np.float64), float(r.cost), int(rec.last_eval_counts[1])
    if arm == "C0":
        return np.asarray(R, np.float64).copy(), float(cost), 0
    r = orig_lo(
        voxel_vertices=vertices,
        phase_index=phase,
        initial_orientation=np.asarray(R, np.float64),
        rng=np.random.default_rng(0),
        cma_sigma0_deg=CMA_SIGMA0_DEG,
        cma_max_evals=CMA_EVALS[arm],
        optimizer="cma",
    )
    return (
        np.asarray(r.orientation, np.float64),
        float(r.cost),
        int(rec.last_local_optimization_evals),
    )


def stored_e1(sample: str) -> List[Dict[str, Any]]:
    """The seed-step log of stored S8 (one row per step)."""
    return RC.read_gz(STORED_S8[sample] / "e1_log.json.gz")
