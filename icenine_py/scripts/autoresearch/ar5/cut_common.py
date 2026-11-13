"""AR5 step 10 (task ar5j): shared set-up of the level-0 cut study (script layer only).

Import this module first in every ``cut_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid`` or ``held``). It

* points step 9's cache variable at ``benchmarks/autoresearch/ar5/cut/cache`` (untracked; the stored
  arm ``A`` and ``G`` runs, passes and the per-level record are linked there read only by
  ``cut_link_cache.sh``) before importing ``rescue_common`` (which reads ``RESCUE_CACHE``);
* refuses ``held`` until ``ar5/cut/frozen.json`` is committed and agrees (``cut_freeze.problems``),
  before anything of ``held`` is read; ``rescue_common`` adds step 9's own freeze check;
* provides ``install_cut`` (put the level-0 cut in front of the reconstructor's discrete stage).
"""

import os
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
CUT = ICE / "benchmarks" / "autoresearch" / "ar5" / "cut"  # tracked results (small)
CCACHE = Path(os.environ.get("CUT_CACHE", str(CUT / "cache")))  # untracked
os.environ["RESCUE_CACHE"] = str(CCACHE)
os.environ["ORDER_CACHE"] = str(CCACHE)
sys.path.insert(0, str(HERE))

import cut_freeze as CF  # noqa: E402

if os.environ.get("RESEED_SAMPLE") == "held":
    _p = CF.problems()
    if _p:
        raise SystemExit(f"cut rule: held is not read before the freeze commit ({_p})")

import rescue_common as RCM  # noqa: E402,F401  (sets the cache; step 9 guard for held)
import cut_rules as CR  # noqa: E402

_ORIG: Optional[Callable[..., Any]] = None


def install_cut(rec: Any, n: int, log: List[Dict[str, Any]]) -> None:
    """Wrap ``icenine.reconstructor.run_discrete_search_spaced`` with the level-0 cut of ``n``
    (``n`` 0 puts the plain function back). ``rec`` is the reconstructor whose full start set
    marks level 0."""
    global _ORIG
    import icenine.reconstructor as R

    if _ORIG is None:
        _ORIG = R.run_discrete_search_spaced
    if n <= 0:
        R.run_discrete_search_spaced = _ORIG
    else:
        R.run_discrete_search_spaced = CR.make_cut_search(
            _ORIG, n, lambda: rec.setup.fz_orientations, log
        )
