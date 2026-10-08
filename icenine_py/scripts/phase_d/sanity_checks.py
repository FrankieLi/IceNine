"""Phase D sanity checks on the rendered full-sample images.

(a) cost-function check: for a handful of voxels the true (new) orientation vs orientations
    0.25/0.5/1.0 deg away, scored with the reconstruction's hard VoxelCostFunction (Q-max 8),
    loading the images through the Phase D reconstruction config (so this also tests the loader).
(b) serial-path containment: a few voxels re-rendered alone with the serial forward model; the
    fraction of their lit pixels that are lit in the full clean images (expected ~100%).
(c) spot/pixel statistics per frame from the render summary.

Usage (from icenine_py/): uv run python scripts/phase_d/sanity_checks.py --variant clean
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
EX = ROOT / "Examples" / "Example2.ManyGrains"


def cost_check(variant: str, n_vox: int, seed: int) -> dict:
    from icenine.config_file import ConfigFile
    from icenine.reconstructor import AdaptiveVoxelReconstructor, _get_voxel_vertices
    from icenine.reconstructor import setup_reconstruction

    cfg = ConfigFile.from_file(str(HERE / "configs" / f"ReconstructPhaseD_mc_{variant}.config"))
    t0 = time.time()
    setup = setup_reconstruction(cfg)
    t_load = time.time() - t0
    rec = AdaptiveVoxelReconstructor(setup)
    mic = setup.sample.get_mic()
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(mic.voxels), n_vox, replace=False)
    rows = []
    for i in idx:
        v = mic.voxels[int(i)]
        R = np.asarray(v.orientation, dtype=np.float64)
        verts = _get_voxel_vertices(v)

        def q(Rx: np.ndarray) -> float:
            return float(rec.evaluate_overlap(Rx.astype(np.float32), verts, v.phase).quality)

        row = {"voxel": int(i), "q_true": q(R)}
        for ang in (0.25, 0.5, 1.0):
            vals = []
            for _ in range(3):
                ax = rng.normal(size=3)
                ax /= np.linalg.norm(ax)
                dR = Rotation.from_rotvec(np.radians(ang) * ax).as_matrix()
                vals.append(q(dR @ R))
            row[f"q_{ang}deg_mean"] = float(np.mean(vals))
            row[f"q_{ang}deg_max"] = float(np.max(vals))
        rows.append(row)
        print(row, flush=True)
    return {"variant": variant, "load_s": t_load, "rows": rows}


def containment(n_vox: int, seed: int, tag: str) -> dict:
    from icenine.config_file import ConfigFile
    from icenine.forward_simulation import ForwardSimulation
    from icenine.image_data import ImageData
    from icenine.sample import Sample
    from icenine.simulation import Simulation
    from noise import read_ascii_image

    cfg = ConfigFile.from_file(str(EX / "ConfigFiles" / "Example2.Simulation.config"))
    cfg.sample_filename = "SimInput/rand_500grains_1mm_neworient_s0.mic"
    sim = ForwardSimulation(cfg)
    sim.exp_setup.initialize_experiment()
    sim.simulator = Simulation(sim.exp_setup)
    dets = sim.exp_setup.get_detector_list()
    sample = Sample()
    sim.exp_setup.initialize_sample(sample, dets[0])
    mic = sample.get_mic()
    rng = np.random.default_rng(seed)
    pick = rng.choice(len(mic.voxels), n_vox, replace=False)
    out = []
    root = EX / "ScatteringData_PhaseD" / tag / "clean"
    for i in pick:
        mic_all = mic.voxels
        mic.voxels = [mic_all[int(i)]]
        imgs = [
            [ImageData(d.num_rows, d.num_cols) for d in dets]
            for _ in sim.exp_setup.get_omega_range_list()
        ]
        sim._simulate_peaks(imgs, dets, sample, sim.exp_setup.get_range_to_index_map())
        mic.voxels = mic_all
        tot = hit = 0
        for oi, row in enumerate(imgs):
            for di, im in enumerate(row):
                a = im._pixels_dense.numpy() > 0
                if not a.any():
                    continue
                full = (
                    read_ascii_image(
                        str(root / f"500Grains.sim{oi:05d}.d{di}"), a.shape[0], a.shape[1]
                    )
                    > 0
                )
                tot += int(a.sum())
                hit += int((a & full).sum())
        out.append({"voxel": int(i), "lit": tot, "in_full": hit, "frac": hit / max(tot, 1)})
        print(out[-1], flush=True)
    return {"rows": out}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="clean")
    ap.add_argument("--n-vox", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="full")
    ap.add_argument("--skip-cost", action="store_true")
    ap.add_argument("--skip-containment", action="store_true")
    a = ap.parse_args()
    os.chdir(EX)
    res = {}
    if not a.skip_containment:
        res["containment"] = containment(4, a.seed, a.tag)
    if not a.skip_cost:
        res["cost"] = cost_check(a.variant, a.n_vox, a.seed)
    (HERE / "results").mkdir(exist_ok=True)
    name = f"sanity_{a.tag}_{a.variant}.json"
    (HERE / "results" / name).write_text(json.dumps(res, indent=2) + "\n")
