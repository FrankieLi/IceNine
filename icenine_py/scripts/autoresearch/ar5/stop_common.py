"""AR5 step 1: shared set-up. Import this module first in every ``stop_*.py`` that needs a sample.

The sample is chosen with the environment variable ``RESEED_SAMPLE`` (``big``, ``grad`` or ``mid``;
as in the scripts of AR3 steps 8 and 9, whose ``reseed_common`` this module imports). The stored
runs of S8 are read through links in ``benchmarks/autoresearch/ar5/cache/stored/<sample>/``
(untracked; they point into the worktrees of steps 8 and 9) and are never written.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]  # icenine_py
AR5 = ICE / "benchmarks" / "autoresearch" / "ar5"
CACHE = AR5 / "cache"
SAMPLE = os.environ.get("RESEED_SAMPLE", "grad")
os.environ.setdefault("RESEED_CACHE", str(CACHE / "stored" / SAMPLE))
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar3"))
sys.path.insert(0, str(HERE))

import reseed_common as RC  # noqa: E402,F401  (sets AR3_SAMPLE from RESEED_SAMPLE)

STORED = CACHE / "stored" / SAMPLE
S8_DIR = STORED / f"S8_{SAMPLE}"
PASS_DIR = STORED / f"passes_{SAMPLE}"
TRIES_DIR = PASS_DIR / "D4_S8"
N = 24570
WRONG = 1.0  # degrees: a voxel is wrong above this error
SAMPLES = ("big", "grad", "mid")


def errors_path(sample: str) -> Path:
    base = ICE / "benchmarks" / "autoresearch" / "ar3"
    if sample == "mid":
        return base / "start" / "errors_mid.npz"
    return base / "reseed" / "slim" / f"errors_{sample}.npz"
