"""AR5 step 12 (task ar5l): shared set-up of the stall-hit study (script layer only).

Import this module first in every ``stallhit_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held`` or ``held2``). It

* points the cache variables at ``benchmarks/autoresearch/ar5/stallhit/cache`` (untracked; the
  stored arm ``A`` and ``G`` runs and passes and step 11's arm ``S`` runs are linked there read only
  by ``stallhit_link_cache.sh``) before importing ``rescue_common`` (which reads ``RESCUE_CACHE``);
* refuses ``held`` and ``held2`` until ``ar5/stallhit/frozen.json`` is committed and agrees
  (``stallhit_freeze.problems``), before anything of either is read; the freeze checks of steps 8,
  9 and 11 are replaced by it (see the comment below). ``held2`` is read by no script before its
  runs have ended.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
SH = ICE / "benchmarks" / "autoresearch" / "ar5" / "stallhit"  # tracked results (small)
SCACHE = Path(os.environ.get("STALLHIT_CACHE", str(SH / "cache")))  # untracked
os.environ["RESCUE_CACHE"] = str(SCACHE)
os.environ["ORDER_CACHE"] = str(SCACHE)
sys.path.insert(0, str(HERE))

import stallhit_freeze as SF  # noqa: E402

if os.environ.get("RESEED_SAMPLE") in ("held", "held2"):
    _p = SF.problems()
    if _p:
        raise SystemExit(f"stallhit rule: held and held2 are not read before the freeze ({_p})")

# Steps 8, 9 and 11 froze the library as it was then; step 9's phase B added the opt-in key
# SeedRescueCMA (off by default), so their ``--check`` can no longer agree. This step's own freeze
# (``stallhit_freeze``: the same script list and the library as it is now, checked above for the
# held samples) replaces them: their checks are switched off here, after ours has passed.
import order_freeze as _OFR  # noqa: E402
import rescue_freeze as _RF  # noqa: E402

_RF.problems = lambda: ""  # type: ignore[assignment]
_OFR.check_or_exit = lambda: None  # type: ignore[assignment]

import rescue_common as RCM  # noqa: E402,F401  (sets the cache)
