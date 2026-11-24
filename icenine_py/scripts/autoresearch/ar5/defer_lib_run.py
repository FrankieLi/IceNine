"""AR5 step 22, phase B: a library run with the keys ``SeedIndexMCStallHit 500`` and
``SeedIndexMCDefer 1`` on (``LibJ5D``).

  RESEED_SAMPLE=mid DEFER_CACHE=<cache> uv run python \
      scripts/autoresearch/ar5/defer_lib_run.py [--cap-min 300]

Step 16's ``high500_lib_run.py`` (``LibJ5``: step 9's ``LibG``, the keys of run ``S11L`` of step 5
plus ``SeedRescueCMA 0.25 1000``, plus ``SeedIndexMCStallHit 500``) with one more key line,
``SeedIndexMCDefer 1``. Everything is the library's own: the keys are written to the config; no
method is replaced and nothing is counted by patching. The library form of the script arm
``J5Dl`` (``defer_hl.py``). This script does not go through ``defer_common``: this step's freeze
pins the library files and so refuses once the key exists (phase B changes ``icenine/``); the
script arm ``J5Dl`` ran before, under the freeze. Writes ``<cache>/LibJ5D_<sample>/`` (the files of
``high500_lib_run.py``) and ``run.done`` at the end. ``mid``, ``grad`` and ``big`` only.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
_cache = os.environ.get(
    "DEFER_CACHE", str(ICE / "benchmarks" / "autoresearch" / "ar5" / "defer" / "cache")
)
os.environ["HIGH500_CACHE"] = _cache  # step 16's script reads this one
sys.path.insert(0, str(HERE))

import high500_lib_run as H  # noqa: E402

H.KEYS = H.KEYS + ["SeedIndexMCDefer 1"]
H.NAME = "LibJ5D"

if __name__ == "__main__":
    H.main()
