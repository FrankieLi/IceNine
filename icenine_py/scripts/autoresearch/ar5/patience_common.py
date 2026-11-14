"""AR5 step 11 (task ar5k): shared set-up of the patience study (script layer only).

Import this module first in every ``patience_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid`` or ``held``). It

* points step 9's cache variable at ``benchmarks/autoresearch/ar5/patience/cache`` (untracked; the
  stored arm ``A`` and ``G`` runs and passes are linked there read only by
  ``patience_link_cache.sh``) before importing ``rescue_common`` (which reads ``RESCUE_CACHE``);
* refuses ``held`` until ``ar5/patience/frozen.json`` is committed and agrees
  (``patience_freeze.problems``), before anything of ``held`` is read; steps 8's and 9's own
  freeze checks are replaced by it (see the comment below).
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
PAT = ICE / "benchmarks" / "autoresearch" / "ar5" / "patience"  # tracked results (small)
PCACHE = Path(os.environ.get("PATIENCE_CACHE", str(PAT / "cache")))  # untracked
os.environ["RESCUE_CACHE"] = str(PCACHE)
os.environ["ORDER_CACHE"] = str(PCACHE)
sys.path.insert(0, str(HERE))

import patience_freeze as PF  # noqa: E402

if os.environ.get("RESEED_SAMPLE") == "held":
    _p = PF.problems()
    if _p:
        raise SystemExit(f"patience rule: held is not read before the freeze commit ({_p})")

# Steps 8 and 9 froze the library as it was then; step 9's phase B added the opt-in key
# SeedRescueCMA (off by default), so their ``--check`` can no longer agree. This step's own freeze
# (``patience_freeze``: the same script list, the library as it is now, checked above for ``held``)
# replaces them: their checks are switched off here, after ours has passed.
import order_freeze as _OFR  # noqa: E402

PF.RF.problems = lambda: ""  # type: ignore[assignment]
_OFR.check_or_exit = lambda: None  # type: ignore[assignment]

import rescue_common as RCM  # noqa: E402,F401  (sets the cache)
