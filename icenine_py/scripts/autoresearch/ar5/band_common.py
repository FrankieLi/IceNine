"""AR5 step 15 (task ar5o): shared set-up of the band study (script layer; no library change).

Import this module first in every ``band_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held`` or ``held3``). In this order it

1. refuses ``RESEED_SAMPLE=held2`` (and any other name that is not a sample of this step) before
   anything is imported (``band_guard``); nothing of ``held2`` is used by step 15;
2. for ``held`` and ``held3``: refuses unless this step's own freeze (``ar5/band/frozen.json``,
   ``band_freeze.problems``) is committed and agrees, before anything of either is read;
3. makes ``held3`` known to the AR3 modules at run time (``band_samples.boot``; no older script is
   edited);
4. points every cache variable at ``benchmarks/autoresearch/ar5/band/cache`` (untracked; step 14's
   stored arm-``J`` files are linked there per entry, read only, by ``band_link_cache.sh``) before
   ``rescue_common`` is imported (it reads ``RESCUE_CACHE``);
5. only for ``held`` and ``held3``, and only after 2 has passed: switches off the older freeze
   checks (steps 8, 9, 12 and 14; see the comment below), in this process only.

``load_gz`` and ``load_json`` read stored files through ``band_guard.guarded`` (no path that
mentions ``held2``).
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
BD = ICE / "benchmarks" / "autoresearch" / "ar5" / "band"  # tracked results (small)
BCACHE = Path(os.environ.get("BAND_CACHE", str(BD / "cache")))  # untracked
sys.path.insert(0, str(HERE))

import band_guard as BG  # noqa: E402

SAMPLE = os.environ.get("RESEED_SAMPLE", "")
BG.refuse_sample(SAMPLE)  # 1. held2 (or any foreign name) is refused before anything else

import band_freeze as BF  # noqa: E402

if SAMPLE in BG.HELD_SAMPLES:  # 2. our own freeze first
    _p = BF.problems()
    if _p:
        BG.refuse(f"band rule: held and held3 are not read before the freeze ({_p})")

import band_samples as BS  # noqa: E402

BS.boot()  # 3. held3 for the AR3 modules (a no-op for the other samples)

for _v in (
    "RESCUE_CACHE",
    "ORDER_CACHE",
    "STALLHIT_CACHE",
    "POLISH_CACHE",
    "RESEED_CACHE",
    "BAND_CACHE",
):
    os.environ[_v] = str(BCACHE)  # 4. one cache for everything this step writes

import order_freeze as _OFR  # noqa: E402
import polish_freeze as _PF  # noqa: E402
import rescue_freeze as _RF  # noqa: E402
import stallhit_freeze as _SF  # noqa: E402

# 5. Steps 8, 9 and 12 froze the library as it was then; step 12's phase B added the opt-in key
# SeedIndexMCStallHit (off by default), so their ``--check`` can no longer agree, and ``stallhit_
# common`` / ``order_run`` / ``polish_run`` call them for ``held``. Step 14 missed this on ``held``
# because its smoke ran on ``big`` only. This step's own freeze (the same script lists and the
# library as it is now, checked above for the held samples) replaces them there. No hash of any
# older step is re-recorded, no older script and no older frozen file is edited.
if SAMPLE in BG.HELD_SAMPLES:
    _RF.problems = lambda: ""  # type: ignore[assignment]
    _OFR.check_or_exit = lambda: None  # type: ignore[assignment]
    _SF.problems = lambda: ""  # type: ignore[assignment]
    _PF.problems = lambda: ""  # type: ignore[assignment]
    _PF.check_or_exit = lambda: None  # type: ignore[assignment]

import rescue_common as RCM  # noqa: E402,F401  (sets the cache)


def load_gz(path: Path) -> Any:
    with gzip.open(BG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(BG.guarded(path).read_text())


def g2_form() -> str:
    """The form of step 14's rule ``G2`` (frozen there); arm ``J`` and arm ``JB`` use it."""
    chosen = str(
        json.loads((ICE / "benchmarks/autoresearch/ar5/polish/frozen.json").read_text())["chosen"]
    )
    assert chosen == "G2b", chosen
    return chosen
