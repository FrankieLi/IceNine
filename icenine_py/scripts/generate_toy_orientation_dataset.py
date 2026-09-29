#!/usr/bin/env python3
"""
Generate the windowed local-orientation-refinement datasets for the toy NN.

Picks one voxel from Example2.ThreeVoxels' ground-truth .mic, defines its ROI
peak set at the ground-truth orientation, then renders thresholded windows
around the ROI peaks for many small orientation offsets.

  train: offsets drawn from the prior (uniform in a ball of radius --prior-radius)
  test:  --test-per-bin offsets at each fixed magnitude in --test-magnitudes,
         random directions (the bench_hp_sweep-style protocol)

Windows are stored as uint8 lit / not-lit values (thresholded), so the network
cannot read orientation information from exact simulated intensities.

Usage:
  cd icenine_py
  uv run python scripts/generate_toy_orientation_dataset.py --smoke-test
  uv run python scripts/generate_toy_orientation_dataset.py --n-train 1500 --test-per-bin 30
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np
import torch

project_root = Path(__file__).parent.parent.parent
DEFAULT_EXAMPLE = project_root / "Examples" / "Example2.ThreeVoxels"


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


def build_problem(example_dir: Path, voxel_index=None):
    """Set up physics and pick a voxel with a non-empty ROI set.

    Returns a dict with everything both this script and the Bayes baseline need.
    """
    from icenine.orientation_nn import define_roi_set

    mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices = setup_example(example_dir)
    candidates = [voxel_index] if voxel_index is not None else range(len(mic.voxels))
    for idx in candidates:
        voxel = mic.voxels[idx]
        vertices = get_vertices(voxel)
        R_nom = voxel.orientation.astype(np.float64)
        roi_list = define_roi_set(
            torch.from_numpy(R_nom).float(), vertices, sample, detector_list, range_map, exp_setup,
            structure_list, simulator, phase_index=voxel.phase,
        )
        if roi_list:
            return dict(
                voxel_index=idx, voxel=voxel, vertices=vertices, R_nom=R_nom, roi_list=roi_list, sample=sample,
                detector_list=detector_list, range_map=range_map, exp_setup=exp_setup, simulator=simulator,
            )
    raise RuntimeError("No voxel in the .mic produced a non-empty ROI set")


def render_dataset(problem, offsets_deg: np.ndarray, window_size: int, label: str) -> torch.Tensor:
    """Render thresholded uint8 windows (N, n_peaks, W, W) for the given offsets."""
    from icenine.orientation_eval import offsets_to_matrices
    from icenine.orientation_nn import render_local_windows

    mats = offsets_to_matrices(offsets_deg, problem["R_nom"])
    n_peaks = len(problem["roi_list"])
    windows = torch.zeros(len(offsets_deg), n_peaks, window_size, window_size, dtype=torch.uint8)
    t0 = time.time()
    for i, mat in enumerate(mats):
        w, _missing = render_local_windows(
            torch.from_numpy(mat).float(), problem["roi_list"], problem["vertices"], problem["sample"],
            problem["detector_list"], problem["range_map"], problem["exp_setup"], problem["simulator"],
            window_size=window_size,
        )
        windows[i] = (w > 0).to(torch.uint8)
        if (i + 1) % 100 == 0 or (i + 1) == len(mats):
            print(f"  [{label}] {i + 1}/{len(mats)} ({time.time() - t0:.0f}s)")
    return windows


def main():
    from icenine.orientation_eval import sample_fixed_magnitude_offsets, sample_prior_offsets

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--voxel-index", type=int, default=None, help="Force a specific voxel; default auto-selects")
    parser.add_argument("--n-train", type=int, default=1500)
    parser.add_argument("--test-per-bin", type=int, default=30)
    parser.add_argument("--test-magnitudes", type=float, nargs="+", default=[0.25, 0.5, 1.0, 2.0])
    parser.add_argument("--prior-radius", type=float, default=2.5, help="Training prior: uniform ball radius (deg)")
    parser.add_argument("--window-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", default=str(Path(__file__).parent))
    parser.add_argument("--tag", default="stage0")
    parser.add_argument("--smoke-test", action="store_true", help="Tiny run: 60 train, 5 per bin")
    args = parser.parse_args()
    if args.smoke_test:
        args.n_train, args.test_per_bin, args.tag = 60, 5, "smoke"

    outdir_abs = Path(args.outdir).resolve()  # build_problem() changes directory
    print(f"Setting up {DEFAULT_EXAMPLE} ...")
    problem = build_problem(DEFAULT_EXAMPLE, args.voxel_index)
    n_peaks = len(problem["roi_list"])
    print(f"  voxel {problem['voxel_index']}: {n_peaks} ROI peaks")

    rng = np.random.default_rng(args.seed)
    train_offsets = sample_prior_offsets(args.n_train, args.prior_radius, rng)
    test_offsets, test_mags = [], []
    for mag in args.test_magnitudes:
        test_offsets.append(sample_fixed_magnitude_offsets(args.test_per_bin, mag, rng))
        test_mags += [mag] * args.test_per_bin
    test_offsets = np.concatenate(test_offsets)

    meta = dict(
        window_size=args.window_size, n_peaks=n_peaks, voxel_index=problem["voxel_index"],
        R_nom=torch.from_numpy(problem["R_nom"]), prior_radius_deg=args.prior_radius, seed=args.seed,
        example="threevoxels",
    )
    outdir = outdir_abs
    for name, offsets, extra in (
        ("train", train_offsets, {}),
        ("test", test_offsets, {"magnitudes_deg": torch.tensor(test_mags)}),
    ):
        windows = render_dataset(problem, offsets, args.window_size, name)
        path = outdir / f"toy_orientation_{args.tag}_{name}.pt"
        torch.save({"windows": windows, "offsets_deg": torch.from_numpy(offsets).float(), **meta, **extra}, path)
        print(f"Saved {path}  ({windows.numel() / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
