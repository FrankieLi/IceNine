"""AR5 step 24 (task ar5x): shared set-up of the ball-first study (script layer; no library change).
Step 23's ``ball_common`` with this step's guards, freeze and cache.

Import this module first in every script of this step that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held``, ``held3``, ``held4``, ``held5``,
``held6``; ``held7`` only after the gate is judged). In this order it

1. refuses ``RESEED_SAMPLE`` with ``held2`` or ``fine7`` in it (and any other name that is not a
   sample of this step) before anything is imported (``bf_guard``); ``held7`` is refused until
   ``<cache>/logs/gate_passed`` exists;
2. for EVERY sample: refuses unless this step's own freeze (``ar5/bf/frozen.json``,
   ``bf_freeze.problems``) is committed and agrees (every script and library file pinned by its
   sha256), before anything is read: nothing of this step runs before the freeze is committed;
3. only after 2 has passed: puts this step's guards in the place of step 23's (in the module
   ``ball_guard``, this process only) and switches off step 23's freeze check. No hash of any
   older step is re-recorded, no older script and no older frozen file is edited. The precedent:
   ``ball_common``, ``defer_common``;
4. imports step 23's set-up (``ball_common``: step 22's, step 21's, ... ``held7`` known) with the
   cache variables pointed at ``BF_CACHE`` (default ``benchmarks/autoresearch/ar5/bf/cache``,
   untracked; a link to a folder outside the worktree).

``load_gz`` and ``load_json`` read stored files through ``bf_guard.guarded``.
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "bf"  # tracked results (small)
sys.path.insert(0, str(HERE))
os.environ.setdefault("SCREEN_CACHE", str(RD / "cache"))  # untracked
os.environ.setdefault("DEFER_CACHE", os.environ["SCREEN_CACHE"])
os.environ.setdefault("BALL_CACHE", os.environ["DEFER_CACHE"])
os.environ.setdefault("BF_CACHE", os.environ["BALL_CACHE"])

import bf_guard as BG  # noqa: E402

SAMPLE = os.environ.get("RESEED_SAMPLE", "")
BG.refuse_sample(SAMPLE)  # 1. held2 / fine7 (or any foreign name, or held7 too early) first

import bf_freeze as BFZ  # noqa: E402

_p = BFZ.problems()  # 2. our own freeze first, for every sample
if _p:
    BG.refuse(f"bf rule: nothing runs before the freeze is committed and agrees ({_p})")

import ball_freeze as _BALLFZ  # noqa: E402
import ball_guard as _BALLG  # noqa: E402

_BALLFZ.problems = lambda: ""  # type: ignore[assignment]  # 3. step 23's check is off (ours passed)
for _n in (
    "refuse",
    "refuse_sample",
    "mentions_held2",
    "mentions_fine7",
    "mentions_held7",
    "held7_open",
    "producer",
    "guarded",
    "install",
    "base_cache",
    "SAMPLES",
    "HELD_SAMPLES",
    "DISCLOSED",
    "MARKER",
    "GATE_PASSED",
    "PRODUCER_OK",
):
    setattr(_BALLG, _n, getattr(BG, _n))

import ball_common as BC  # noqa: E402,F401  (4.: step 23's, which sets up step 22's)

DC = BC.DC
RCACHE = DC.RCACHE
DG = BG


def load_gz(path: Path) -> Any:
    with gzip.open(BG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(BG.guarded(path).read_text())
