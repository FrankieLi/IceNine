#!/usr/bin/env python3
"""
Generate a windowed local-orientation-refinement dataset for the toy NN.

Picks one voxel from Example2.ThreeVoxels' ground-truth .mic, defines its ROI
peak set at the ground-truth orientation, then renders windowed detector
crops for many small random perturbations around that orientation. Caches the
result to a .pt file consumed by OrientationDataset / train_toy_orientation_nn.py.

Usage:
  cd icenine_py
  uv run python scripts/generate_toy_orientation_dataset.py --smoke-test
  uv run python scripts/generate_toy_orientation_dataset.py --n-samples 20000
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np
import torch

project_root = Path(__file__).parent.parent.parent


def setup_example(example_dir: Path, basename: str = "3Grains.sim"):
    """Load physics + ground-truth mic. Same pattern as benchmarks/bench_riemannian_optimization.py."""
    from icenine.config_file import ConfigFile
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.mic_file import MicFile
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    os.chdir(example_dir)

    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = basename

    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()
    detector_list = exp_setup.get_detector_list()
    range_map = exp_setup.get_range_to_index_map()

    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])
    simulator = Simulation(exp_setup)
    structure_list = sample.get_structure_list()

    mic_path = Path(config.sample_filename)
    if not mic_path.is_absolute():
        mic_path = example_dir / mic_path
    mic = MicFile.read(str(mic_path))

    return mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, _get_voxel_vertices


def select_voxel(mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices, voxel_index=None):
    """Pick a voxel whose ROI set is non-empty (peaks actually land on a detector)."""
    from icenine.orientation_nn import define_roi_set

    candidates = [voxel_index] if voxel_index is not None else range(len(mic.voxels))
    for idx in candidates:
        voxel = mic.voxels[idx]
        vertices = get_vertices(voxel)
        orientation = torch.from_numpy(voxel.orientation).float()
        roi_list = define_roi_set(
            orientation, vertices, sample, detector_list, range_map, exp_setup, structure_list,
            simulator, phase_index=voxel.phase,
        )
        if roi_list:
            return idx, voxel, vertices, orientation, roi_list
    raise RuntimeError("No voxel in the .mic produced a non-empty ROI set")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example", default="threevoxels", choices=["threevoxels"])
    parser.add_argument("--voxel-index", type=int, default=None, help="Force a specific voxel; default auto-selects")
    parser.add_argument("--n-samples", type=int, default=500)
    parser.add_argument("--max-angle-deg", type=float, default=2.0, help="Bound on perturbation angle")
    parser.add_argument("--window-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", default=None, help="Output .pt path (default: derived from --n-samples)")
    parser.add_argument("--smoke-test", action="store_true", help="Shorthand for --n-samples 300")
    args = parser.parse_args()

    if args.smoke_test:
        args.n_samples = 300

    example_dir = project_root / "Examples" / "Example2.ThreeVoxels"
    print(f"Loading {args.example} example from {example_dir} ...")
    mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices = setup_example(example_dir)
    print(f"  Loaded .mic with {len(mic.voxels)} voxels")

    idx, voxel, voxel_vertices, nominal_orientation, roi_list = select_voxel(
        mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices,
        voxel_index=args.voxel_index,
    )
    print(f"  Selected voxel {idx}: {len(roi_list)} ROI peaks")

    from icenine.orientation_nn import render_local_windows, sample_local_perturbations

    rng = np.random.default_rng(args.seed)
    matrices, quats = sample_local_perturbations(
        nominal_orientation, args.n_samples, args.max_angle_deg, rng
    )

    n_peaks = len(roi_list)
    windows = torch.zeros(args.n_samples, n_peaks, args.window_size, args.window_size, dtype=torch.float32)
    quaternions = torch.zeros(args.n_samples, 4, dtype=torch.float32)
    n_missing_total = 0

    t0 = time.time()
    for i, (mat, q) in enumerate(zip(matrices, quats)):
        perturbed_orientation = torch.from_numpy(mat).float()
        sample_windows, missing = render_local_windows(
            perturbed_orientation, roi_list, voxel_vertices, sample, detector_list,
            range_map, exp_setup, simulator, window_size=args.window_size,
        )
        windows[i] = sample_windows
        quaternions[i] = torch.from_numpy(q).float()
        n_missing_total += int(missing.sum().item())

        if (i + 1) % 100 == 0 or (i + 1) == args.n_samples:
            elapsed = time.time() - t0
            print(f"  {i + 1}/{args.n_samples} samples ({elapsed:.1f}s, {n_missing_total} missing peaks so far)")

    missing_rate = n_missing_total / (args.n_samples * n_peaks) if n_peaks else 0.0
    print(f"Missing-peak rate: {missing_rate:.4%} (peak dropped out of ROI set under perturbation)")

    output_path = args.output
    if output_path is None:
        output_dir = Path(__file__).parent
        output_path = output_dir / f"toy_orientation_dataset_{args.example}_v{idx}_n{args.n_samples}.pt"

    torch.save(
        {
            "windows": windows,
            "quaternions": quaternions,
            "n_peaks": n_peaks,
            "window_size": args.window_size,
            "voxel_index": idx,
            "nominal_orientation": nominal_orientation,
            "max_angle_deg": args.max_angle_deg,
            "example": args.example,
        },
        output_path,
    )
    print(f"Saved dataset to {output_path}")


if __name__ == "__main__":
    main()
