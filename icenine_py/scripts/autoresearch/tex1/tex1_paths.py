"""tex1: the paths of the sample ``tex24`` and the constants of its freeze.

``tex24`` lives on the grid of ``fine7`` (98,280 voxels, side 0.0046875 mm), copied unchanged into
``<cache>/grid/`` and checked against the sha256 recorded in ``fine1/frozen.json``.
"""

import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]  # icenine_py
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
BENCH = ICE / "benchmarks" / "autoresearch" / "tex1"  # tracked: records, freeze
CACHE = BENCH / "cache"  # gitignored, a link to the main checkout's folder
FINE1_FROZEN = ICE / "benchmarks" / "autoresearch" / "fine1" / "frozen.json"
TAG = "tex24"
SEED = 24
N_CELLS = 250
MIN_SEP_DEG = 5.0
SLICE_N = 2000
SIDE_FINE_MM = 0.6 / 128  # 0.0046875
N_FINE = 98280
# classes (owner task of 2026-11-24)
NOISE_FRACTION = 0.20  # of the grains: voxel-to-voxel noise
GRADIENT_FRACTION = 0.50  # of the grains: smooth orientation gradient
NOISE_RNG = [SEED, 1]  # default_rng([24, 1]) picks the noise grains
GRADIENT_RNG = [SEED, 2]  # default_rng([24, 2]) picks the gradient grains
NOISE_BALL_DEG = 1.0  # per-voxel noise: rotation vector uniform in this ball (pairs within 2)
GRADIENT_CAP_DEG = 0.5  # largest edge-pair misorientation of a gradient grain
GRADIENT_MIN_DEG = 0.1  # strength t_g uniform on [0.1, 0.5]
CURVATURE = 0.25

GRID = CACHE / "grid" / "fine7_grid.mic"  # generation 7, 98,280 voxels, orientations zero
SLICE = CACHE / "grid" / "tex24_slice2000.mic"  # the 2,000 voxels nearest the centre
SAMPLE_DIR = CACHE / "sample"
TRUTH = SAMPLE_DIR / "tex24_truth.mic"
IMAGES = CACHE / "images"  # <IMAGES>/tex24/clean/500Grains.sim*.d0/.d1
FRAMES = IMAGES / TAG / "clean"
STACK = CACHE / "stacks" / "tex24_clean.npy"
RUNS = CACHE / "runs"
LOGS = CACHE / "logs"
SCORE = CACHE / "score"

# tracked records
SAMPLE_JSON = BENCH / "sample" / "sample.json"
GRAINMAP = BENCH / "sample" / "grainmap.npz"
SLICE_JSON = BENCH / "slice.json"
STACK_JSON = BENCH / "stack_check.json"
FROZEN = BENCH / "frozen.json"
RENDER_RECORD = ICE / "scripts" / "phase_d" / "results" / "render_tex24_nodrop.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()
