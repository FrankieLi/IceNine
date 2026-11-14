"""AR5 step 9 (task ar5i): shared set-up of the rescue study (script layer; no library change).

Import this module first in every ``rescue_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid`` or ``held``). It

* points step 8's cache at ``benchmarks/autoresearch/ar5/rescue/cache`` (untracked) before
  importing ``order_common`` (which reads ``ORDER_CACHE``); arm ``A`` of step 8 is linked there as
  ``A_<sample>`` and ``passes_<sample>/{D4_A,before_A.npy,after_A.npy}`` (read only);
* refuses ``held`` until ``ar5/rescue/frozen.json`` is committed and agrees
  (``rescue_freeze.problems``): no held seed log, fallback end, map or count is read before;
* provides ``cma_rescue`` (one rescue CMA from one kept end) used by the run and the replay.
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RESCUE = ICE / "benchmarks" / "autoresearch" / "ar5" / "rescue"  # tracked results (small)
RCACHE = Path(os.environ.get("RESCUE_CACHE", os.environ.get("ORDER_CACHE", str(RESCUE / "cache"))))
os.environ["ORDER_CACHE"] = str(RCACHE)
sys.path.insert(0, str(HERE))

import rescue_freeze as RF  # noqa: E402

if os.environ.get("RESEED_SAMPLE") == "held":
    _p = RF.problems()
    if _p:
        raise SystemExit(f"rescue rule: held is not read before the freeze commit ({_p})")

import order_common as OC  # noqa: E402,F401  (sets the cache; installs the arms)
import rescue_rules as RR  # noqa: E402

import numpy as np  # noqa: E402


def cma_rescue(
    rec: Any, orig_lo: Any, ratios: Any, R: np.ndarray, vertices: Any, phase: int
) -> Dict[str, Any]:
    """One rescue: CMA from ``R`` (first step 0.25 degree, 1,000 evaluations, ``default_rng(0)``).
    ``orig_lo`` is the reconstructor's unwrapped ``local_optimization``; ``ratios(info)`` gives
    (hit ratio, confidence) as ``BFSReconstruction._ratios``. The hit ratio is one overlap
    evaluation, not counted. Returns hit, cost, evals, R (end) and ``result`` (the return value)."""
    r = orig_lo(
        voxel_vertices=vertices,
        phase_index=phase,
        initial_orientation=np.asarray(R, np.float64),
        rng=np.random.default_rng(0),
        cma_sigma0_deg=RR.SIGMA0_DEG,
        cma_max_evals=RR.EVALS,
        optimizer="cma",
    )
    ev = int(rec.last_local_optimization_evals)
    Rend = np.asarray(r.orientation, np.float64)
    info = rec.evaluate_overlap(Rend, vertices, phase)
    return dict(hit=float(ratios(info)[0]), cost=float(r.cost), evals=ev, R=Rend, result=r)
