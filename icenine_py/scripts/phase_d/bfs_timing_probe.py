"""Per-voxel cost probe on the Phase D clean images (memmap loader, mc config): wall time of
a no-start `reconstruct_voxel` (BFS seed) and of the classic neighbour step
`local_optimization` from a start 0.3 deg off the truth. Single process, other jobs may be
running ("contended"). Checks the todo doc's per-seed / per-neighbour time assumptions.

Usage (from icenine_py/): uv run python scripts/phase_d/bfs_timing_probe.py [--n-seed 4 --n-nb 20]
"""

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-seed", type=int, default=4)
    ap.add_argument("--n-nb", type=int, default=20)
    ap.add_argument("--variant", default="clean")
    a = ap.parse_args()
    os.chdir(EX)
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import (
        AdaptiveVoxelReconstructor,
        _get_voxel_vertices,
        setup_reconstruction,
    )

    cfg = ConfigFile.from_file(str(HERE / "configs" / f"ReconstructPhaseD_mc_{a.variant}.config"))
    stack = EX / "ScatteringData_PhaseD" / "stacks" / f"{a.variant}.npy"
    setup = setup_reconstruction(cfg, exp_data=ExperimentalData.from_binary_memmap(stack, 180, 2))
    rec = AdaptiveVoxelReconstructor(setup)
    mic = setup.sample.get_mic()
    truth = MicFile.read("SimInput/rand_500grains_1mm_neworient_s0.mic")
    rng = np.random.default_rng(3)
    idx = rng.choice(len(mic.voxels), a.n_seed + a.n_nb, replace=False)
    seeds, nbs = [], []
    for k, i in enumerate(idx):
        v, R = mic.voxels[int(i)], np.asarray(truth.voxels[int(i)].orientation, dtype=np.float64)
        verts = _get_voxel_vertices(v)
        t0 = time.time()
        if k < a.n_seed:
            res = rec.reconstruct_voxel(voxel_vertices=verts, phase_index=v.phase, rng=rng)
            R_out = res.orientation
        else:
            ax = rng.normal(size=3)
            ax /= np.linalg.norm(ax)
            start = Rotation.from_rotvec(np.radians(0.3) * ax).as_matrix() @ R
            R_out = rec.local_optimization(
                verts, v.phase, start.astype(np.float32), rng
            ).orientation
        dt = time.time() - t0
        err = np.degrees(np.linalg.norm(Rotation.from_matrix(np.asarray(R_out) @ R.T).as_rotvec()))
        (seeds if k < a.n_seed else nbs).append({"voxel": int(i), "s": dt, "err_deg": float(err)})
        print(k, int(i), f"{dt:.1f}s", f"{err:.3f} deg", flush=True)
    out = {
        "variant": a.variant,
        "label": "single process, contended (other jobs on the machine)",
        "seed_s": [x["s"] for x in seeds],
        "seed_err_deg": [x["err_deg"] for x in seeds],
        "neighbour_median_s": float(np.median([x["s"] for x in nbs])),
        "neighbour_mean_s": float(np.mean([x["s"] for x in nbs])),
        "neighbour_err_deg_median": float(np.median([x["err_deg"] for x in nbs])),
        "n_neighbours": len(nbs),
    }
    (HERE / "results" / f"bfs_timing_probe_{a.variant}.json").write_text(
        json.dumps(out, indent=2) + "\n"
    )
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
