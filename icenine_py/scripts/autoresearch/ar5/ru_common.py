"""AR5 step 29 (task ar5ad): shared set-up of the ``RU`` study (script layer; no library change).
Step 28's ``pxf_common`` with this step's cache and freeze, and ``held8`` made known to the guards
(under the same gate as ``held7`` was: ``<cache>/logs/gate_passed`` must exist, and its logs and
maps stay closed until ``<cache>/logs/bf_runs.done``).

Import this module first in every script of this step that touches a sample, with
``RESEED_SAMPLE`` set. In order it

1. points every cache variable at ``RU_CACHE`` (default: the run cache in the main checkout,
   ``ar5/ru/cache``, untracked; inside the worktree when it is not a worktree);
2. unless ``RU_PRE_FREEZE=1`` (only scripts that write the freeze), refuses unless this step's
   freeze (``ar5/ru/frozen.json``) is committed and agrees (``ru_freeze.problems``: every script of
   this step, of step 28 and of every older step, and every library file, by sha256);
3. teaches step 24's guard (``bf_guard``) and step 16's guard (``high500_guard``) the name
   ``held8`` (same treatment as ``held7`` in step 28: the name only; step 28's guard fixes), puts
   seed 8 into the shared sample table (``ru_samples``), and refuses ``held2``/``fine7`` and
   ``held8`` before ``gate_passed`` exists;
4. imports step 28's set-up (``pxf_common``), whose own freeze check also has to pass.
"""

import os
import sys
from pathlib import Path
from typing import Any, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ru"  # tracked results (small)


def default_cache() -> Path:
    """The run cache: in the main checkout when this is one of its worktrees, else beside RD."""
    parts = ICE.parts
    if ".claude" in parts and "worktrees" in parts:
        main = Path(*parts[: parts.index(".claude")])
        return main / "icenine_py" / "benchmarks" / "autoresearch" / "ar5" / "ru" / "cache"
    return RD / "cache"


sys.path.insert(0, str(HERE))
os.environ.setdefault("RU_CACHE", str(default_cache()))
for _k in ("PXF_CACHE", "SCREEN_CACHE", "DEFER_CACHE", "BALL_CACHE", "BF_CACHE", "IX_CACHE"):
    os.environ.setdefault(_k, os.environ["RU_CACHE"])

SAMPLE = os.environ.get("RESEED_SAMPLE", "")


def refuse(msg: str) -> None:
    print(msg, file=sys.stderr)
    raise SystemExit(2)


PRE_FREEZE = os.environ.get("RU_PRE_FREEZE") == "1"
if not PRE_FREEZE:
    import ru_freeze as RF  # noqa: E402

    _p = RF.problems()
    if _p:
        refuse(f"ru rule: nothing runs before the freeze is committed and agrees ({_p})")

import bf_guard as BG  # noqa: E402

LATE = ("held7", "held8")
BG.SAMPLES = tuple(BG.SAMPLES) + ("held8",)  # type: ignore[assignment]
BG.HELD_SAMPLES = tuple(BG.HELD_SAMPLES) + ("held8",)  # type: ignore[assignment]
BG.PRODUCER_OK = tuple(BG.PRODUCER_OK) + ("held8_clean.npy",)  # type: ignore[assignment]
_orig_refuse_sample = BG.refuse_sample


def _refuse_sample(name: str, cache: Optional[Path] = None) -> None:
    """Step 24's refusal (``held2``, ``fine7``, foreign names, ``held7`` before the gate), plus
    ``held8`` only once ``<cache>/logs/gate_passed`` exists."""
    _orig_refuse_sample(name, cache)
    if name == "held8" and not ((cache or BG.base_cache()) / BG.GATE_PASSED).exists():
        BG.refuse("ru rule: held8 is not built or run before the gate is judged (gate_passed)")


def _mentions_late(p: Any, cache: Optional[Path] = None) -> bool:
    """A part of the path (below the cache folder, when it is under it) mentions ``held7`` or
    ``held8``."""
    root = cache or BG.base_cache()
    q = Path(p)
    for base in (root, root.resolve()):
        try:
            q = q.relative_to(base)
            break
        except ValueError:
            continue
    return any(any(n in x for n in LATE) for x in BG._parts(q))


BG.refuse_sample = _refuse_sample  # type: ignore[assignment]
BG.mentions_held7 = _mentions_late  # type: ignore[assignment]
BG.refuse_sample(SAMPLE)

import high500_guard as _H5G  # noqa: E402

for _name in ("held7", "held8"):
    if _name not in _H5G.SAMPLES:
        _H5G.SAMPLES = tuple(_H5G.SAMPLES) + (_name,)  # type: ignore[assignment]

import ru_samples as _RSM  # noqa: E402,F401  (seed 8 into the shared table before the boots)

# Step 28's own freeze (ar5/pxf/frozen.json) pins the library as it was in step 28's phase A; its
# phase B (merged as fa31ffac) added opt-in keys to icenine/config_file.py and reconstructor.py, so
# that check can no longer agree. This step's freeze (checked above, before anything else) pins the
# current library and every script of step 28 and of every older step by sha256, so step 28's check
# is switched off for this process (PXF_PRE_FREEZE=1 skips only that check). No older hash is
# re-recorded and no older file is edited.
os.environ["PXF_PRE_FREEZE"] = "1"

import pxf_common as PC  # noqa: E402  (step 28's set-up; its freeze check is the one above)

DC = PC.DC
RCACHE = PC.RCACHE
CACHE = Path(os.environ["RU_CACHE"])
HERE_RU = HERE
load_gz = PC.load_gz
load_json = PC.load_json
setup = PC.setup
