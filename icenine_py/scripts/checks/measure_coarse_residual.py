"""Measure the orientation error the adaptive coarse search hands to FindOptimal (Stage 1, decision D2).

Runs AdaptiveVoxelReconstructor.reconstruct_voxel on each Example2.ThreeVoxels voxel with
ReconstructQ8.config (Q_max = 8, 4 levels from a 5 deg grid radius) against the Python-simulated
data (ScatteringData_Python), several seeds per voxel. MCOptimizer.optimize is wrapped to record
every call: the calls with the full max_mc_steps are FindOptimal, and their starting orientations
(in candidate order, best first) are the coarse search's hand-off. Reports the symmetry-reduced
misorientation of the best candidate, and of the closest candidate among those tried, to the
ground-truth orientation in the .mic file, plus the orientation offset delta (rotation vector
applied in the sample frame) of the best candidate.

Usage: cd icenine_py && uv run python scripts/checks/measure_coarse_residual.py [--seeds 5]
"""

import argparse
import math
import os
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ICENINE_PY = Path(__file__).resolve().parents[2]
EXAMPLE = ICENINE_PY.parent / "Examples" / "Example2.ThreeVoxels"


def main():
    from icenine import orientation_search
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import (
        AdaptiveVoxelReconstructor,
        _get_voxel_vertices,
        setup_reconstruction,
    )
    from icenine.sampling import get_misorientation, matrix_to_quaternion
    from icenine.symmetry import create_cubic_symmetry

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument(
        "--out",
        default=str(ICENINE_PY / "benchmarks" / "toy_orientation_stage1" / "coarse_residual.npz"),
    )
    args = parser.parse_args()
    out_path = Path(args.out).resolve()

    os.chdir(EXAMPLE)
    config = ConfigFile.from_file("ConfigFiles/ReconstructQ8.config")
    config.out_file_basename = "3Grains.sim"
    exp_data = ExperimentalData.from_image_directory(
        directory="ScatteringData_Python",
        basename="3Grains.sim",
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
    )
    setup = setup_reconstruction(config, exp_data=exp_data)
    recon = AdaptiveVoxelReconstructor(setup)
    max_steps = recon.params.max_mc_steps
    mic = MicFile.read("SimInput/three_voxels.mic")
    sym = create_cubic_symmetry(4.0)
    sym_q = np.array(
        [
            matrix_to_quaternion(np.array(m))
            for m in sym.get_rotation_matrices()
            if np.linalg.det(m) > 0
        ]
    )

    calls = []
    original = orientation_search.MCOptimizer.optimize

    def recording(
        self,
        initial_orientation,
        angular_box_side,
        angular_step,
        max_mc_steps,
        max_restarts,
        *a,
        **k,
    ):
        calls.append(
            (np.array(initial_orientation, dtype=np.float64), angular_box_side, max_mc_steps)
        )
        return original(
            self,
            initial_orientation,
            angular_box_side,
            angular_step,
            max_mc_steps,
            max_restarts,
            *a,
            **k,
        )

    orientation_search.MCOptimizer.optimize = recording

    def misori_deg(R, R_true):
        return math.degrees(
            get_misorientation(matrix_to_quaternion(R), matrix_to_quaternion(R_true), sym_q)
        )

    def offset_deg(R, R_true):
        """Rotation-vector offset of R from R_true, using the symmetry-equivalent of R closest to R_true."""
        best = None
        for S in sym.get_rotation_matrices():
            S = np.array(S)
            if np.linalg.det(S) < 0:
                continue
            Rs = R @ S
            d = Rotation.from_matrix(Rs @ R_true.T).as_rotvec()
            if best is None or np.linalg.norm(d) < np.linalg.norm(best):
                best = d
        return np.degrees(best)

    rows = []
    for vi, voxel in enumerate(mic.voxels):
        for seed in range(args.seeds):
            calls.clear()
            t0 = time.time()
            result = recon.reconstruct_voxel(
                voxel_vertices=_get_voxel_vertices(voxel),
                phase_index=voxel.phase,
                rng=np.random.default_rng(seed),
            )
            final = [c for c in calls if c[2] == max_steps]
            if not final:
                print(f"voxel {vi} seed {seed}: no FindOptimal calls recorded")
                continue
            box_deg = math.degrees(final[0][1])
            errs = [misori_deg(c[0], voxel.orientation) for c in final]
            delta = offset_deg(final[0][0], voxel.orientation)
            rows.append(
                dict(
                    voxel=vi,
                    seed=seed,
                    best=errs[0],
                    closest=min(errs),
                    n=len(final),
                    box=box_deg,
                    delta=delta,
                    final=misori_deg(result.orientation, voxel.orientation),
                )
            )
            print(
                f"voxel {vi} seed {seed}: hand-off best {errs[0]:.3f} deg (delta {np.round(delta, 3)}), "
                f"closest of {len(final)} {min(errs):.3f} deg, FindOptimal box {box_deg:.3f} deg, "
                f"final {rows[-1]['final']:.3f} deg  ({time.time() - t0:.0f}s)"
            )

    orientation_search.MCOptimizer.optimize = original
    best = np.array([r["best"] for r in rows])
    closest = np.array([r["closest"] for r in rows])
    deltas = np.array([r["delta"] for r in rows])
    print(
        f"\n{len(rows)} runs. Hand-off error of the best candidate (deg): median {np.median(best):.3f}, "
        f"90th pct {np.percentile(best, 90):.3f}, max {best.max():.3f}"
    )
    print(
        f"Closest candidate among those passed to FindOptimal: median {np.median(closest):.3f}, max {closest.max():.3f}"
    )
    print(
        f"|delta| components of the best candidate, RMS: x {np.sqrt(np.mean(deltas[:, 0]**2)):.3f}, "
        f"y {np.sqrt(np.mean(deltas[:, 1]**2)):.3f}, z {np.sqrt(np.mean(deltas[:, 2]**2)):.3f} deg"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_path,
        best=best,
        closest=closest,
        deltas=deltas,
        voxel=np.array([r["voxel"] for r in rows]),
        seed=np.array([r["seed"] for r in rows]),
        final=np.array([r["final"] for r in rows]),
    )
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
