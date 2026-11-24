"""AR5 step 28, phase B: a library run with the keys ``SeedIndexMCDefer 1``,
``RayIndexMaxPairs 50000`` and ``RayIndexFallbackFull 1`` on (``LibPXF``).

  RESEED_SAMPLE=mid PXF_CACHE=<cache> uv run python scripts/autoresearch/ar5/pxf_lib_run.py \
      [--cap-min 300]

Step 22's ``defer_lib_run.py`` (``LibJ5D``) with two more key lines. Everything is the library's
own: the keys are written to the config; no method is replaced and nothing is counted by
patching. The library form of the script arm ``PXFl`` (``pxf_hl.py``). This script does not go
through ``pxf_common``: that freeze pins the library files and so refuses once the keys exist
(phase B changes ``icenine/``); the script arm ``PXFl`` ran before, under the freeze. Writes
``<cache>/LibPXF_<sample>/`` (the files of ``high500_lib_run.py``) and ``run.done``. ``mid`` only.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
_cache = os.environ.get(
    "PXF_CACHE", str(ICE / "benchmarks" / "autoresearch" / "ar5" / "pxf" / "cache")
)
os.environ["HIGH500_CACHE"] = _cache  # step 16's script reads this one
sys.path.insert(0, str(HERE))

import high500_lib_run as H  # noqa: E402

H.KEYS = H.KEYS + [
    "SeedIndexMCDefer 1",
    "RayIndexMaxPairs 50000",
    "RayIndexFallbackFull 1",
]
H.NAME = "LibPXF"

if __name__ == "__main__":
    H.main()
