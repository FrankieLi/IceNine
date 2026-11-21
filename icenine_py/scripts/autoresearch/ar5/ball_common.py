"""AR5 step 23 (task ar5w): shared set-up of the ball study (script layer; no library change).
Step 22's ``defer_common`` with this step's guards, freeze and cache, and ``held7`` made known.

Import this module first in every script of this step that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held``, ``held3``, ``held4``, ``held5``,
``held6``; ``held7`` only after the gate is judged). In this order it

1. refuses ``RESEED_SAMPLE`` with ``held2`` or ``fine7`` in it (and any other name that is not a
   sample of this step) before anything is imported (``ball_guard``); ``held7`` is refused until
   ``<cache>/logs/gate_passed`` exists;
2. for EVERY sample: refuses unless this step's own freeze (``ar5/ball/frozen.json``,
   ``ball_freeze.problems``) is committed and agrees (every script and library file pinned by its
   sha256), before anything is read: nothing of this step runs before the freeze is committed;
3. only after 2 has passed: puts this step's guards in the place of step 22's (in the module
   ``defer_guard``, this process only) and switches off step 22's freeze check (its hashes moved
   with the library's key phases). No hash of any older step is re-recorded, no older script and no
   older frozen file is edited. The precedent: ``defer_common``, ``screen_common``;
4. imports step 22's set-up (``defer_common``: step 21's, step 16's, ``held6`` known) with the
   cache variables pointed at ``BALL_CACHE`` (default ``benchmarks/autoresearch/ar5/ball/cache``,
   untracked; a link to a folder outside the worktree);
5. makes ``held7`` known to the AR3 modules (``ball_samples``).

``load_gz`` and ``load_json`` read stored files through ``ball_guard.guarded``.
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ball"  # tracked results (small)
sys.path.insert(0, str(HERE))
os.environ.setdefault("SCREEN_CACHE", str(RD / "cache"))  # untracked
os.environ.setdefault("DEFER_CACHE", os.environ["SCREEN_CACHE"])
os.environ.setdefault("BALL_CACHE", os.environ["DEFER_CACHE"])

import ball_guard as BG  # noqa: E402

SAMPLE = os.environ.get("RESEED_SAMPLE", "")
BG.refuse_sample(SAMPLE)  # 1. held2 / fine7 (or any foreign name, or held7 too early) first

import ball_freeze as BF  # noqa: E402

_p = BF.problems()  # 2. our own freeze first, for every sample
if _p:
    BG.refuse(f"ball rule: nothing runs before the freeze is committed and agrees ({_p})")

import defer_freeze as _DFZ  # noqa: E402
import defer_guard as _DG  # noqa: E402

_DFZ.problems = lambda: ""  # type: ignore[assignment]  # 3. step 22's check is off (ours passed)
_DFZ.check_or_exit = lambda: None  # type: ignore[assignment]
for _n in ("refuse_sample", "guarded", "mentions_held2", "producer", "refuse"):
    setattr(_DG, _n, getattr(BG, _n))
_DG.SAMPLES = BG.SAMPLES  # type: ignore[assignment]
_DG.HELD_SAMPLES = BG.HELD_SAMPLES  # type: ignore[assignment]
_DG.DISCLOSED = BG.DISCLOSED  # type: ignore[assignment]
_DG.base_cache = BG.base_cache  # type: ignore[assignment]
_DG.MARKER = BG.MARKER
_DG.PRODUCER_OK = BG.PRODUCER_OK

import defer_common as DC  # noqa: E402,F401  (4.)
import ball_samples as _BSM  # noqa: E402,F401  (5. adds seed 7 to the shared table)

RCACHE = DC.RCACHE
DG = BG


def load_gz(path: Path) -> Any:
    with gzip.open(BG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(BG.guarded(path).read_text())
