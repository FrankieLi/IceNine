"""AR5 step 27 (task ar5aa): shared set-up of the ``max_pairs`` study (script layer; no library
change). Step 26's ``ix_common`` with this step's cache and freeze; ``held7`` allowed under step
24's guard (``bf_guard``) once ``<cache>/logs/gate_passed`` exists.

Import this module first in every script of this step that touches a sample, with
``RESEED_SAMPLE`` set. In order it

1. points every cache variable at ``PX_CACHE`` (default: the run cache in the main checkout,
   ``ar5/px/cache``, untracked; inside the worktree when it is not a worktree);
2. unless ``PX_PRE_FREEZE=1`` (only scripts that write the freeze), refuses unless this step's
   freeze (``ar5/px/frozen.json``) is committed and agrees (``px_freeze.problems``: every script of
   this step and every library file by sha256);
3. refuses ``held2`` and ``fine7`` in the sample name and ``held7`` until the gate file exists (step
   24's ``bf_guard``, which also keeps ``held7``'s logs and maps closed until
   ``<cache>/logs/bf_runs.done``);
4. only then switches off step 24's own freeze check (its library hashes are the older ones; this
   process only; no older hash is re-recorded) and imports step 24's set-up.
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "px"  # tracked results (small)


def default_cache() -> Path:
    """The run cache: in the main checkout when this is one of its worktrees, else beside RD."""
    parts = ICE.parts
    if ".claude" in parts and "worktrees" in parts:
        main = Path(*parts[: parts.index(".claude")])
        return main / "icenine_py" / "benchmarks" / "autoresearch" / "ar5" / "px" / "cache"
    return RD / "cache"


sys.path.insert(0, str(HERE))
os.environ.setdefault("PX_CACHE", str(default_cache()))
for _k in ("SCREEN_CACHE", "DEFER_CACHE", "BALL_CACHE", "BF_CACHE", "IX_CACHE"):
    os.environ.setdefault(_k, os.environ["PX_CACHE"])

SAMPLE = os.environ.get("RESEED_SAMPLE", "")


def refuse(msg: str) -> None:
    print(msg, file=sys.stderr)
    raise SystemExit(2)


PRE_FREEZE = os.environ.get("PX_PRE_FREEZE") == "1"
if not PRE_FREEZE:
    import px_freeze as PXF  # noqa: E402

    _p = PXF.problems()
    if _p:
        refuse(f"px rule: nothing runs before the freeze is committed and agrees ({_p})")

import bf_guard as BG  # noqa: E402

BG.refuse_sample(SAMPLE)  # held2, fine7, foreign names; held7 until gate_passed

import bf_freeze as _BFF  # noqa: E402

_BFF.problems = lambda: ""  # type: ignore[assignment]  # step 24's check is off (ours passed)

import bf_common as BCM  # noqa: E402,F401

DC = BCM.DC
RCACHE = BCM.RCACHE
CACHE = Path(os.environ["PX_CACHE"])


def load_gz(path: Path) -> Any:
    with gzip.open(BG.guarded(path), "rt") as f:
        return json.load(f)


def load_json(path: Path) -> Any:
    return json.loads(BG.guarded(path).read_text())


def stored_steps(sample: str) -> List[Tuple[int, bool]]:
    """(voxel, reseed) of every stored seed step of ``sample``'s J5D run, in run order."""
    e = load_gz(CACHE / f"J5D_{sample}" / "e1_log.json.gz")
    return [(int(x["voxel"]), bool(x["reseed"])) for x in e]


def setup() -> Dict[str, Any]:
    """The environment of ``RESEED_SAMPLE``: ``env`` and ``rec``."""
    import reseed_common as RC

    env = RC.make_env()
    return dict(env=env, rec=env.hooks.rec)
