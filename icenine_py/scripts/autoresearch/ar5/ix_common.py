"""AR5 step 26 (task ar5z): shared set-up of the faster-indexer study (script layer).
Step 24's ``bf_common`` with this step's cache and freeze; nothing older is edited.

Import this module first in every script of this step that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid``, ``held``, ``held3`` to ``held6``). In order it

1. refuses ``RESEED_SAMPLE`` with ``held2``, ``fine7`` or ``held7`` in it (step 26 uses none of
   them: ``held7`` stays unbuilt under tier 1) and any name that is not a disclosed sample;
2. unless ``IX_PRE_FREEZE=1`` (phase 0 only: the profile and gate z0, which call the old indexer
   and change nothing), refuses unless this step's freeze (``ar5/ix/frozen.json``) is committed and
   agrees (``ix_freeze.problems``: every script of this step and every library file by sha256);
3. only then switches off step 24's own freeze check (its library hashes moved with this step's
   new function; this process only; no hash of an older step is re-recorded) and imports step
   24's set-up, with every cache variable pointed at ``IX_CACHE`` (default
   ``benchmarks/autoresearch/ar5/ix/cache``, untracked; a link to a folder outside the worktree).
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ix"  # tracked results (small)
sys.path.insert(0, str(HERE))
os.environ.setdefault("IX_CACHE", str(RD / "cache"))  # untracked
for _k in ("SCREEN_CACHE", "DEFER_CACHE", "BALL_CACHE", "BF_CACHE"):
    os.environ.setdefault(_k, os.environ["IX_CACHE"])

DISCLOSED = ("big", "grad", "mid", "held", "held3", "held4", "held5", "held6")
SAMPLE = os.environ.get("RESEED_SAMPLE", "")


def refuse(msg: str) -> None:
    print(msg, file=sys.stderr)
    raise SystemExit(2)


if any(x in SAMPLE for x in ("held2", "fine7", "held7")) or (SAMPLE and SAMPLE not in DISCLOSED):
    refuse(f"ix rule: sample {SAMPLE!r} is not a disclosed sample of step 26")

PRE_FREEZE = os.environ.get("IX_PRE_FREEZE") == "1"
if not PRE_FREEZE:
    import ix_freeze as IXF  # noqa: E402

    _p = IXF.problems()
    if _p:
        refuse(f"ix rule: nothing runs before the freeze is committed and agrees ({_p})")

import bf_freeze as _BFF  # noqa: E402

_BFF.problems = lambda: ""  # type: ignore[assignment]  # step 24's check is off (ours passed)

import bf_common as BCM  # noqa: E402,F401

DC = BCM.DC
RCACHE = BCM.RCACHE
CACHE = Path(os.environ["IX_CACHE"])


def load_gz(path: Path) -> Any:
    with gzip.open(BCM.BG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(BCM.BG.guarded(path).read_text())


def stored_steps(sample: str) -> List[Tuple[int, bool]]:
    """(voxel, reseed) of every stored seed step of ``sample``'s J5D run, in run order."""
    e = load_gz(CACHE / f"J5D_{sample}" / "e1_log.json.gz")
    return [(int(x["voxel"]), bool(x["reseed"])) for x in e]


def setup() -> Dict[str, Any]:
    """The environment of ``RESEED_SAMPLE``: ``env`` and ``rec``."""
    import reseed_common as RC

    env = RC.make_env()
    return dict(env=env, rec=env.hooks.rec)
