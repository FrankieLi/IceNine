"""fine1: the paths of the sample ``fine7`` and the constants of its freeze."""

import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]  # icenine_py
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
BENCH = ICE / "benchmarks" / "autoresearch" / "fine1"  # tracked: records, freeze
CACHE = BENCH / "cache"  # gitignored, a link to the main checkout's folder
SRC_GRID = EX / "SimInput" / "rand_500grains_1mm_neworient_s0_grid.mic"  # generation 6, 24,570
TAG = "fine7"
SEED = 7
N_CELLS = 250
MIN_SEP_DEG = 5.0
LEVELS = (0.0, 0.2, 0.5, 1.0, 2.0, 4.0)
SLICE_N = 2000
SIDE_FINE_MM = 0.6 / 128  # 0.0046875
N_FINE = 98280

GRID = CACHE / "grid" / "fine7_grid.mic"  # generation 7, 98,280 voxels, orientations zero
SLICE = CACHE / "grid" / "fine7_slice2000.mic"
SAMPLE_DIR = CACHE / "sample"
TRUTH = SAMPLE_DIR / "gradient_truth.mic"
IMAGES = CACHE / "images"  # <IMAGES>/fine7/clean/500Grains.sim*.d0/.d1
FRAMES = IMAGES / TAG / "clean"
STACK = CACHE / "stacks" / "fine7_clean.npy"
RUNS = CACHE / "runs"
LOGS = CACHE / "logs"

# tracked records
GRID_CHECK = BENCH / "grid_check.json"
SAMPLE_JSON = BENCH / "sample" / "sample.json"
GRAINMAP = BENCH / "sample" / "grainmap.npz"
CAPS_JSON = BENCH / "sample" / "caps.json"
GRADIENT_CHECK = BENCH / "gradient_check.json"
SLICE_JSON = BENCH / "slice.json"
STACK_JSON = BENCH / "stack_check.json"
FROZEN = BENCH / "frozen.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()
