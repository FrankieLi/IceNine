"""AR5 step 22 (task ar5v): shared set-up of the defer study (script layer; no library change).
Step 21's ``screen_common`` with this step's guards, freeze and cache, and ``held6`` made known.

Import this module first in every script of this step that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held``, ``held3``, ``held4``, ``held5``;
``held6`` only after the gate is judged). In this order it

1. refuses ``RESEED_SAMPLE=held2`` (and any other name that is not a sample of this step) before
   anything is imported (``defer_guard``); ``held6`` is refused until ``<cache>/logs/gate_passed``
   exists;
2. for EVERY sample: refuses unless this step's own freeze (``ar5/defer/frozen.json``,
   ``defer_freeze.problems``) is committed and agrees (every script and library file pinned by its
   sha256), before anything is read: nothing of this step runs before the freeze is committed;
3. only after 2 has passed: puts this step's guards in the place of step 21's (``defer_guard.
   install``) and switches off the freeze checks of step 21 and of steps 8, 9, 12, 14, 15 and 16
   in this process only (their hashes moved with later library edits; step 21's moved with its key
   phase: ``icenine/boundary_pass.py``). No hash of any older step is re-recorded, no older script
   and no older frozen file is edited. The precedent: ``screen_common``, ``cmarepair_common``;
4. makes ``held6`` known to the AR3 modules (``defer_samples``: step 21's table plus seed 6) and
   widens step 16's guard to ``held4``, ``held5`` and ``held6`` (its sample list is a module
   variable); ``held3`` is disclosed here, its logs are readable;
5. imports step 21's set-up (``screen_common``: cache variables pointed at ``SCREEN_CACHE``,
   default ``benchmarks/autoresearch/ar5/defer/cache``, untracked) and through it step 16's.

``load_gz`` and ``load_json`` read stored files through ``defer_guard.guarded`` (no path that
mentions ``held2``; none that mentions ``held6`` before its runs end).
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "defer"  # tracked results (small)
sys.path.insert(0, str(HERE))
os.environ.setdefault("SCREEN_CACHE", str(RD / "cache"))  # untracked
os.environ.setdefault("DEFER_CACHE", os.environ["SCREEN_CACHE"])

import defer_guard as DG  # noqa: E402

SAMPLE = os.environ.get("RESEED_SAMPLE", "")
DG.refuse_sample(SAMPLE)  # 1. held2 (or any foreign name, or held6 too early) is refused first

import defer_freeze as DF  # noqa: E402

_p = DF.problems()  # 2. our own freeze first, for every sample
if _p:
    DG.refuse(f"defer rule: nothing runs before the freeze is committed and agrees ({_p})")

import screen_freeze as _SCF  # noqa: E402
import screen_guard as _SG  # noqa: E402

DG.install(_SG)  # 3. this step's guards in the place of step 21's
_SCF.problems = lambda: ""  # type: ignore[assignment]
_SCF.check_or_exit = lambda: None  # type: ignore[assignment]

import high500_guard as _H5G  # noqa: E402

for _name in ("held4", "held5", "held6"):  # 4.
    if _name not in _H5G.SAMPLES:
        _H5G.SAMPLES = tuple(_H5G.SAMPLES) + (_name,)  # type: ignore[assignment]
_H5G.held3_open = lambda cache=None: True  # type: ignore[assignment]  # held3 is disclosed

import defer_samples as _DSM  # noqa: E402,F401  (adds seed 6 to the shared table)
import screen_common as SC  # noqa: E402,F401  (5.)

RCACHE = SC.RCACHE


def load_gz(path: Path) -> Any:
    with gzip.open(DG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(DG.guarded(path).read_text())
