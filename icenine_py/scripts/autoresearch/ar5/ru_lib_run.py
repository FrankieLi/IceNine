"""AR5 step 29, phase B: a library run with the key ``SeedReuseKnown 1`` on top of ``LibPXF``'s keys
(``LibRU``).

  RESEED_SAMPLE=mid RU_CACHE=<cache> uv run python scripts/autoresearch/ar5/ru_lib_run.py \
      [--cap-min 300]

Step 28's ``pxf_lib_run.py`` (``LibPXF``) with one more key line. Everything is the library's own:
the keys are written to the config; no method is replaced and nothing is counted by patching. The
library form of the script arm ``RUl`` (``ru_hl.py``). This script does not go through
``ru_common``: that freeze pins the library files and so refuses once the key exists (phase B
changes ``icenine/``); the script arm ``RUl`` ran before, under the freeze. Writes
``<cache>/LibRU_<sample>/`` (the files of ``high500_lib_run.py``) and ``run.done``. ``mid`` only.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]


def _default_cache() -> Path:
    """The run cache: in the main checkout when this is one of its worktrees."""
    parts = ICE.parts
    if ".claude" in parts and "worktrees" in parts:
        main = Path(*parts[: parts.index(".claude")])
        return main / "icenine_py" / "benchmarks" / "autoresearch" / "ar5" / "ru" / "cache"
    return ICE / "benchmarks" / "autoresearch" / "ar5" / "ru" / "cache"


_cache = os.environ.get("RU_CACHE", str(_default_cache()))
os.environ["HIGH500_CACHE"] = _cache  # step 16's script reads this one
sys.path.insert(0, str(HERE))

import high500_lib_run as H  # noqa: E402

H.KEYS = H.KEYS + [
    "SeedIndexMCDefer 1",
    "RayIndexMaxPairs 50000",
    "RayIndexFallbackFull 1",
    "SeedReuseKnown 1",
]
H.NAME = "LibRU"

if __name__ == "__main__":
    H.main()
