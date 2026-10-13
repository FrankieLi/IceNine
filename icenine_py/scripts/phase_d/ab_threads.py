"""Contended A/B of torch's thread count: 30 `local_optimization` fits (from 0.3 deg off the truth)
and one no-start `reconstruct_voxel` (region voxel 18), clean images, mc config.

Usage (from icenine_py/): uv run python scripts/phase_d/ab_threads.py default|1
Run the two arms staggered (not at the same time) for a cleaner comparison; the 2026-10-08 result
(benchmarks/phase_d_pilot/ab_threads.json) had both arms running together with the 5 pilot runs.
"""

import os
import sys

thr = sys.argv[1]
if thr == "1":
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[k] = "1"
import json  # noqa: E402
import re  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402
from scipy.spatial.transform import Rotation  # noqa: E402

if thr == "1":
    torch.set_num_threads(1)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EX = ROOT.parent / "Examples" / "Example2.ManyGrains"
REGION = ROOT / "benchmarks" / "phase_d_pilot"
text = (HERE / "configs" / "ReconstructPhaseD_mc_clean.config").read_text()
text = re.sub(r"(?m)^SampleFilename\s+\S+", f"SampleFilename {REGION / 'region_grid.mic'}", text)
cache = HERE / "cache" / "pilot"
cache.mkdir(parents=True, exist_ok=True)
cfg_path = cache / f"ab_{thr}.config"
cfg_path.write_text(text)
os.chdir(EX)
from icenine.config_file import ConfigFile  # noqa: E402
from icenine.experimental_data import ExperimentalData  # noqa: E402
from icenine.mic_file import MicFile  # noqa: E402
from icenine.reconstructor import (  # noqa: E402
    AdaptiveVoxelReconstructor,
    _get_voxel_vertices,
    setup_reconstruction,
)

cfg = ConfigFile.from_file(str(cfg_path))
stack = EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
setup = setup_reconstruction(cfg, exp_data=ExperimentalData.from_binary_memmap(stack, 180, 2))
rec = AdaptiveVoxelReconstructor(setup)
truth = MicFile.read("SimInput/rand_500grains_1mm_neworient_s0.mic")
idx = np.load(REGION / "region_index.npy")
rng = np.random.default_rng(3)
ts, nev = [], 0
for k in range(30):
    tv = truth.voxels[int(idx[rng.integers(0, len(idx))])]
    R = np.asarray(tv.orientation, dtype=np.float64)
    ax = rng.normal(size=3)
    ax /= np.linalg.norm(ax)
    start = (Rotation.from_rotvec(np.radians(0.3) * ax).as_matrix() @ R).astype(np.float32)
    t0 = time.time()
    rec.local_optimization(_get_voxel_vertices(tv), tv.phase, start, rng)
    ts.append(time.time() - t0)
    nev += rec.last_local_optimization_evals
v = truth.voxels[int(idx[18])]
t0 = time.time()
rec.reconstruct_voxel(
    voxel_vertices=_get_voxel_vertices(v), phase_index=v.phase, rng=np.random.default_rng(0)
)
print(
    json.dumps(
        {
            "threads": thr,
            "torch_threads": torch.get_num_threads(),
            "neighbour_median_s": float(np.median(ts)),
            "neighbour_mean_evals": nev / 30,
            "seed_s": time.time() - t0,
        }
    )
)
