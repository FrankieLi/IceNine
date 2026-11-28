"""AR5 step 14 (task ar5n): shared set-up of the polish study (script layer; no library change).

Import this module first in every ``polish_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held`` or ``held2``). It

* points the cache variables at ``benchmarks/autoresearch/ar5/polish/cache`` (untracked; step 12's
  stored arm-``H`` files are linked there per entry, read only, by ``polish_link_cache.sh``)
  before importing ``rescue_common`` (which reads ``RESCUE_CACHE``);
* refuses ``held2`` until ``ar5/polish/frozen.json`` is committed and agrees
  (``polish_freeze.problems``), before anything of it is read; ``held`` is a disclosed sample;
* switches off the freeze checks of steps 8, 9 and 12 (this step's own freeze replaces them; the
  library is unchanged by this step);
* provides ``load_gz`` / ``load_json``: every read of a stored file goes through
  ``polish_guard.guarded`` (no ``passes_held2``; no ``*_held2`` run folder before the freeze);
* provides ``g2_try``: one ``G2`` try from one rescue end with the running reconstructor.
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
PO = ICE / "benchmarks" / "autoresearch" / "ar5" / "polish"  # tracked results (small)
PCACHE = Path(os.environ.get("POLISH_CACHE", str(PO / "cache")))  # untracked
os.environ["RESCUE_CACHE"] = str(PCACHE)
os.environ["ORDER_CACHE"] = str(PCACHE)
sys.path.insert(0, str(HERE))

import polish_freeze as PF  # noqa: E402
import polish_guard as PG  # noqa: E402

if os.environ.get("RESEED_SAMPLE") == "held2":
    _p = PF.problems()
    if _p:
        raise SystemExit(f"polish rule: held2 is not read before the freeze ({_p})")

import order_freeze as _OFR  # noqa: E402
import rescue_freeze as _RF  # noqa: E402

_RF.problems = lambda: ""  # type: ignore[assignment]
_OFR.check_or_exit = lambda: None  # type: ignore[assignment]

import rescue_common as RCM  # noqa: E402,F401  (sets the cache)
import polish_rules as PR  # noqa: E402

import numpy as np  # noqa: E402


def load_gz(path: Path) -> Any:
    with gzip.open(PG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(PG.guarded(path).read_text())


g2_try = PR.g2_try
