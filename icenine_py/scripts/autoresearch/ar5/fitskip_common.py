"""AR5 step 17 (task ar5q): shared set-up of the fit-skip study (script layer; no library change).

Import this module first in every ``fitskip_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held``, ``held3``; ``held4`` only after the
gate, and only once its build exists). In this order it

1. refuses ``RESEED_SAMPLE=held2`` (and any other name that is not a sample of this step) before
   anything is imported (``fitskip_guard``); ``held4`` is refused until ``ar5/fitskip/gate.json``
   exists;
2. for ``held``, ``held3`` and ``held4``: refuses unless this step's own freeze
   (``ar5/fitskip/frozen.json``, ``fitskip_freeze.problems``) is committed and agrees, before
   anything of them is read;
3. points every cache variable at ``benchmarks/autoresearch/ar5/fitskip/cache`` (untracked; step
   16's stored arm-``J5`` files are linked there per entry, read only, by
   ``fitskip_link_cache.sh``) before step 16's set-up is imported (it reads ``HIGH500_CACHE``);
4. only after 2 has passed (or for a sample that needs no freeze): replaces the stale freeze
   checks of steps 8, 9, 12, 14, 15 and 16 in this process only (see the comment below);
5. imports step 16's set-up (``high500_common``: ``held2`` guard, the held3 boot of step 15,
   ``rescue_common`` and the cache). Its own freeze check is the replaced one.

``load_gz`` and ``load_json`` read stored files through ``fitskip_guard.guarded`` (no path that
mentions ``held2``; none that mentions ``held4`` before its runs end).
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
FS = ICE / "benchmarks" / "autoresearch" / "ar5" / "fitskip"  # tracked results (small)
FCACHE = Path(os.environ.get("FITSKIP_CACHE", str(FS / "cache")))  # untracked
sys.path.insert(0, str(HERE))

import fitskip_guard as FG  # noqa: E402

SAMPLE = os.environ.get("RESEED_SAMPLE", "")
FG.refuse_sample(SAMPLE)  # 1. held2 (or any foreign name) is refused before anything else

import fitskip_freeze as FF  # noqa: E402

if SAMPLE in FG.HELD_SAMPLES:  # 2. our own freeze first
    _p = FF.problems()
    if _p:
        FG.refuse(f"fitskip rule: held, held3 and held4 are not read before the freeze ({_p})")

for _v in (
    "RESCUE_CACHE",
    "ORDER_CACHE",
    "STALLHIT_CACHE",
    "POLISH_CACHE",
    "RESEED_CACHE",
    "BAND_CACHE",
    "HIGH500_CACHE",
    "FITSKIP_CACHE",
):
    os.environ[_v] = str(FCACHE)  # 3. one cache for everything this step writes

import high500_freeze as _H5F  # noqa: E402
import order_freeze as _OFR  # noqa: E402
import polish_freeze as _PF  # noqa: E402
import rescue_freeze as _RF  # noqa: E402
import stallhit_freeze as _SF  # noqa: E402

# 4. Steps 8, 9, 12, 14, 15 and 16 froze the library as it was then. A docstring line of the library
# (step 16's phase B, and step 12's key before it) moved their hashes, so their ``--check`` can no
# longer agree; ``high500_common`` / ``order_run`` / ``polish_run`` call them for the held samples.
# This step's own freeze (the same script lists and the library as it is now, checked above for
# the held samples) replaces them in this process only: the precedent of ``stallhit_common.py``
# and ``polish_replace_freeze.py``. No hash of any older step is re-recorded, no older script and
# no older frozen file is edited.
_RF.problems = lambda: ""  # type: ignore[assignment]
_OFR.check_or_exit = lambda: None  # type: ignore[assignment]
_SF.problems = lambda: ""  # type: ignore[assignment]
_PF.problems = lambda: ""  # type: ignore[assignment]
_PF.check_or_exit = lambda: None  # type: ignore[assignment]
_H5F.problems = lambda: ""  # type: ignore[assignment]
_H5F.check_or_exit = lambda: None  # type: ignore[assignment]

import high500_common as BC  # noqa: E402,F401  (5.: sets the cache, boots held3)
import rescue_common as RCM  # noqa: E402,F401


def load_gz(path: Path) -> Any:
    with gzip.open(FG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(FG.guarded(path).read_text())
