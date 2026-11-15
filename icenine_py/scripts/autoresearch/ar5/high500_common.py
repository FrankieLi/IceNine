"""AR5 step 16 (task ar5p): shared set-up of the high-stop study (script layer; no library change).

Import this module first in every ``high500_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held`` or ``held3``). In this order it

1. refuses ``RESEED_SAMPLE=held2`` (and any other name that is not a sample of this step) before
   anything is imported (``high500_guard``); nothing of ``held2`` is used by step 16;
2. for ``held`` and ``held3``: refuses unless this step's own freeze (``ar5/high500/frozen.json``,
   ``high500_freeze.problems``) is committed and agrees, before anything of either is read;
3. makes ``held3`` known to the AR3 modules at run time (``band_samples.boot``, step 15's file,
   hashed in this step's freeze; no older script is edited);
4. points every cache variable at ``benchmarks/autoresearch/ar5/high500/cache`` (untracked; step
   14's stored arm-``J`` files and step 15's recording rows are linked there per entry, read only,
   by ``high500_link_cache.sh``) before ``rescue_common`` is imported (it reads ``RESCUE_CACHE``);
5. only for ``held`` and ``held3``, and only after 2 has passed: replaces the stale freeze checks
   of steps 8, 9 and 12 (see the comment below), in this process only. Step 14's check
   (``polish_freeze``) is left active: it still agrees.

``load_gz`` and ``load_json`` read stored files through ``high500_guard.guarded`` (no path that
mentions ``held2``; none that mentions ``held3`` before its runs end).
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
HD = ICE / "benchmarks" / "autoresearch" / "ar5" / "high500"  # tracked results (small)
HCACHE = Path(os.environ.get("HIGH500_CACHE", str(HD / "cache")))  # untracked
sys.path.insert(0, str(HERE))

import high500_guard as HG  # noqa: E402

SAMPLE = os.environ.get("RESEED_SAMPLE", "")
HG.refuse_sample(SAMPLE)  # 1. held2 (or any foreign name) is refused before anything else

import high500_freeze as HF  # noqa: E402

if SAMPLE in HG.HELD_SAMPLES:  # 2. our own freeze first
    _p = HF.problems()
    if _p:
        HG.refuse(f"high500 rule: held and held3 are not read before the freeze ({_p})")

import band_samples as BS  # noqa: E402

BS.boot()  # 3. held3 for the AR3 modules (a no-op for the other samples)

for _v in (
    "RESCUE_CACHE",
    "ORDER_CACHE",
    "STALLHIT_CACHE",
    "POLISH_CACHE",
    "RESEED_CACHE",
    "BAND_CACHE",
    "HIGH500_CACHE",
):
    os.environ[_v] = str(HCACHE)  # 4. one cache for everything this step writes

import order_freeze as _OFR  # noqa: E402
import rescue_freeze as _RF  # noqa: E402
import stallhit_freeze as _SF  # noqa: E402

# 5. Steps 8, 9 and 12 froze the library as it was then; step 12's phase B added the opt-in key
# SeedIndexMCStallHit (off by default), so their ``--check`` can no longer agree, and
# ``stallhit_common`` / ``order_run`` call them for ``held``. This step's own freeze (the same
# script lists and the library as it is now, checked above for the held samples) replaces them
# there. Step 14's check is not replaced (it agrees today). No hash of any older step is
# re-recorded, no older script and no older frozen file is edited.
if SAMPLE in HG.HELD_SAMPLES:
    _RF.problems = lambda: ""  # type: ignore[assignment]
    _OFR.check_or_exit = lambda: None  # type: ignore[assignment]
    _SF.problems = lambda: ""  # type: ignore[assignment]

import rescue_common as RCM  # noqa: E402,F401  (sets the cache)


def load_gz(path: Path) -> Any:
    with gzip.open(HG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(HG.guarded(path).read_text())


def g2_form() -> str:
    """The form of step 14's rule ``G2`` (frozen there); arm ``J`` and arm ``J5`` use it."""
    chosen = str(
        json.loads((ICE / "benchmarks/autoresearch/ar5/polish/frozen.json").read_text())["chosen"]
    )
    assert chosen == "G2b", chosen
    return chosen
