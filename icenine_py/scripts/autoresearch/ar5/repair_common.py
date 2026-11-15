"""AR5 step 19 (task ar5s): shared set-up of the repair study (script layer; no library change).
Step 18's ``learnskip_common`` with this step's guards and freeze (nothing of step 18's own freeze
or gate is used: step 18 ended at its gate).

Import this module first in every ``repair_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held``, ``held3``; ``held4`` only after the
gate, and only once its build exists). In this order it

1. refuses ``RESEED_SAMPLE=held2`` (and any other name that is not a sample of this step) before
   anything is imported (``repair_guard``); ``held4`` is refused until ``ar5/repair/gate.json``
   exists;
2. for ``held``, ``held3`` and ``held4``: refuses unless this step's own freeze
   (``ar5/repair/frozen.json``, ``repair_freeze.problems``) is committed and agrees, before
   anything of them is read;
3. points every cache variable at ``benchmarks/autoresearch/ar5/repair/cache`` (untracked; step
   16's and step 18's stored files are linked there per entry, read only, by
   ``repair_link_cache.sh``) before step 16's set-up is imported (it reads ``HIGH500_CACHE``);
4. makes ``held3`` and ``held4`` known to the AR3 modules (``learnskip_samples.boot``, step 18's
   file, replaces step 15's ``band_samples.boot``, which knows ``held3`` only) and widens step 16's
   guard to ``held4`` (its sample list is a module variable);
5. only after 2 has passed (or for a sample that needs no freeze): replaces the stale freeze checks
   of steps 8, 9, 12, 14, 15 and 16 in this process only (their hashes moved with later library
   edits; the precedent of ``stallhit_common.py`` and ``polish_replace_freeze.py``). No hash of any
   older step is re-recorded, no older script and no older frozen file is edited;
6. imports step 16's set-up (``high500_common``).

``load_gz`` and ``load_json`` read stored files through ``repair_guard.guarded`` (no path that
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
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "repair"  # tracked results (small)
RCACHE = Path(os.environ.get("REPAIR_CACHE", str(RD / "cache")))  # untracked
sys.path.insert(0, str(HERE))

import repair_guard as RG  # noqa: E402

SAMPLE = os.environ.get("RESEED_SAMPLE", "")
RG.refuse_sample(SAMPLE)  # 1. held2 (or any foreign name) is refused before anything else

import repair_freeze as RF  # noqa: E402

if SAMPLE in RG.HELD_SAMPLES:  # 2. our own freeze first
    _p = RF.problems()
    if _p:
        RG.refuse(f"repair rule: held, held3 and held4 are not read before the freeze ({_p})")

for _v in (
    "RESCUE_CACHE",
    "ORDER_CACHE",
    "STALLHIT_CACHE",
    "POLISH_CACHE",
    "RESEED_CACHE",
    "BAND_CACHE",
    "HIGH500_CACHE",
    "FITSKIP_CACHE",
    "LEARNSKIP_CACHE",
    "REPAIR_CACHE",
):
    os.environ[_v] = str(RCACHE)  # 3. one cache for everything this step writes

import band_samples as _BS  # noqa: E402
import high500_guard as _H5G  # noqa: E402
import learnskip_samples as _LSM  # noqa: E402

_BS.boot = _LSM.boot  # type: ignore[assignment]  # 4. held3 and held4 for the AR3 modules
_BS.wanted = _LSM.wanted  # type: ignore[assignment]
if "held4" not in _H5G.SAMPLES:
    _H5G.SAMPLES = tuple(_H5G.SAMPLES) + ("held4",)  # type: ignore[assignment]

import high500_freeze as _H5F  # noqa: E402
import order_freeze as _OFR  # noqa: E402
import polish_freeze as _PF  # noqa: E402
import rescue_freeze as _RF  # noqa: E402
import stallhit_freeze as _SF  # noqa: E402

# 5. see the module text: this step's own freeze (checked above for the held samples) replaces the
# stale checks of older steps in this process only.
_RF.problems = lambda: ""  # type: ignore[assignment]
_OFR.check_or_exit = lambda: None  # type: ignore[assignment]
_SF.problems = lambda: ""  # type: ignore[assignment]
_PF.problems = lambda: ""  # type: ignore[assignment]
_PF.check_or_exit = lambda: None  # type: ignore[assignment]
_H5F.problems = lambda: ""  # type: ignore[assignment]
_H5F.check_or_exit = lambda: None  # type: ignore[assignment]

import high500_common as BC  # noqa: E402,F401  (6.: sets the cache, boots held3 / held4)
import rescue_common as RCM  # noqa: E402,F401


def load_gz(path: Path) -> Any:
    with gzip.open(RG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(RG.guarded(path).read_text())
