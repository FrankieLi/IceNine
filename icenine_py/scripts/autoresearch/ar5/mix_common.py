"""AR5 step 3 (task ar5c): the mixed seed refinement in the script layer.

Import this module first in every ``mix_*.py`` script that runs seed steps, with the environment
variable ``RESEED_SAMPLE`` set (the launcher ``mix_bg.sh`` does it). It sends the cache to
``benchmarks/autoresearch/ar5/mix/cache`` (untracked; it holds links to step 2's replay) and then
imports step 2's ``fit_common``, whose arm function ``refine_arm`` it extends with two arms,
without changing any step 2 file:

* ``X60`` and ``X200`` (``mixed_refine``): the CMA fit of ``C60`` or ``C200``; if its end hit ratio
  is below 0.95, the stock MC refinement (arm ``M``) from the same candidate; the end kept is the
  one with the lower cost (a tie keeps the CMA end); the evaluations of both are counted. The seed
  step's early stop (``fit_run.FitPick.e1_step``) then tests the hit ratio of the kept end. The
  rule itself is ``mix_rules.mixed_candidate`` (tested); this is its running version.

Nothing under ``icenine/`` is touched. The hit ratio of the CMA end is one extra overlap
evaluation that is not counted, as the hit ratio of every end is not counted in S8.
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
MIX = ICE / "benchmarks" / "autoresearch" / "ar5" / "mix"  # tracked results (small)
CACHE_DIR = MIX / "cache"  # untracked
os.environ["RESEED_CACHE"] = str(CACHE_DIR)
CACHE_DIR.mkdir(parents=True, exist_ok=True)
(CACHE_DIR / "logs").mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))

import fit_common as FC  # noqa: E402
import fit_rules as FR  # noqa: E402
import numpy as np  # noqa: E402

PICK_FILE = MIX / "pick.json"
X_ARMS = {"X60": "C60", "X200": "C200"}
HOLDER: Dict[str, Any] = {"bfs": None}  # the check script sets it; a run reads the run's pick
LOG: List[Dict[str, Any]] = []  # one row per candidate refined with a mixed arm (analysis only)
_orig_refine_arm = FC.refine_arm


def _bfs() -> Any:
    if HOLDER["bfs"] is not None:
        return HOLDER["bfs"]
    import reseed_s8 as S8

    return S8.state["pick"].bfs


def hit_of(rec: Any, R: "np.ndarray", vertices: Any, phase: int) -> float:
    info = rec.evaluate_overlap(np.asarray(R, np.float64), vertices, phase)
    return float(_bfs()._ratios(info)[0])


def mixed_refine(
    arm: str,
    rec: Any,
    orig_lo: Any,
    R: "np.ndarray",
    cost: float,
    vertices: Any,
    phase: int,
) -> Tuple["np.ndarray", float, int]:
    """One mixed seed refinement of the candidate ``R`` (strict cost ``cost``). Returns (end
    orientation, end cost, evaluations of both refinements), as ``fit_common.refine_arm``."""
    rc, cc, ec = _orig_refine_arm(X_ARMS[arm], rec, orig_lo, R, cost, vertices, phase)
    hc = hit_of(rec, rc, vertices, phase)
    if hc >= FR.TAU:
        LOG.append(dict(c_hit=hc, mc_ran=False, evals=int(ec)))
        return rc, cc, ec
    rm, cm, em = _orig_refine_arm("M", rec, orig_lo, R, cost, vertices, phase)
    LOG.append(dict(c_hit=hc, mc_ran=True, evals=int(ec + em)))
    if cm < cc:
        return rm, cm, int(ec + em)
    return rc, cc, int(ec + em)


def refine_arm_mixed(
    arm: str,
    rec: Any,
    orig_lo: Any,
    R: "np.ndarray",
    cost: float,
    vertices: Any,
    phase: int,
) -> Tuple["np.ndarray", float, int]:
    if arm in X_ARMS:
        return mixed_refine(arm, rec, orig_lo, R, cost, vertices, phase)
    return _orig_refine_arm(arm, rec, orig_lo, R, cost, vertices, phase)


def install() -> None:
    """Add the arms ``X60`` and ``X200`` to step 2's script layer (idempotent)."""
    FC.ARMS = tuple(a for a in FC.ARMS if a not in X_ARMS) + tuple(X_ARMS)
    FC.refine_arm = refine_arm_mixed


install()
