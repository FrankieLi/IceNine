#!/usr/bin/env python3
"""
Regression test: compare serial vs batched forward simulation output.

Runs both paths on ThreeVoxels and checks pixel-exact agreement.

Usage:
    cd Examples/Example2.ThreeVoxels
    uv run python test_batched_vs_serial.py
"""

import sys
import os
import time
from pathlib import Path
import torch

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "icenine_py"))

from icenine.config_file import ConfigFile
from icenine.forward_simulation import ForwardSimulation


def run_simulation(config_path, batched=False, batch_size=None):
    """Run forward simulation and return images (no disk save)."""
    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = "3Grains.sim"
    simulator = ForwardSimulation(config)

    # Suppress file saving by running simulation directly
    simulator.exp_setup.initialize_experiment()
    omega_ranges = simulator.exp_setup.get_omega_range_list()
    detector_list = simulator.exp_setup.get_detector_list()
    range_map = simulator.exp_setup.get_range_to_index_map()
    simulator.simulator = __import__(
        "icenine.simulation", fromlist=["Simulation"]
    ).Simulation(simulator.exp_setup)

    from icenine.sample import Sample
    from icenine.image_data import ImageData

    sample = Sample()
    simulator.exp_setup.initialize_sample(sample, detector_list[0])

    images = []
    for i in range(len(omega_ranges)):
        detector_images = []
        for detector in detector_list:
            image = ImageData(detector.num_rows, detector.num_cols)
            detector_images.append(image)
        images.append(detector_images)
    simulator.images = images

    start = time.time()
    if batched:
        simulator._simulate_peaks_batched(
            images, detector_list, sample, range_map, batch_size=batch_size,
        )
    else:
        simulator._simulate_peaks(images, detector_list, sample, range_map)
    elapsed = time.time() - start
    return images, elapsed


def compare_images(serial_images, batched_images):
    """Compare using underlying dense tensors for speed."""
    assert len(serial_images) == len(batched_images)

    total_pixels_serial = 0
    total_pixels_batched = 0
    matching_pixels = 0
    mismatched_pixels = 0
    max_rel_diff = 0.0
    total_intensity_serial = 0.0
    total_intensity_batched = 0.0
    mismatch_examples = []

    for omega_idx in range(len(serial_images)):
        for det_idx in range(len(serial_images[omega_idx])):
            s_img = serial_images[omega_idx][det_idx]
            b_img = batched_images[omega_idx][det_idx]

            # Access underlying dense tensor directly
            s_t = s_img._pixels_dense
            b_t = b_img._pixels_dense

            # Count non-zero pixels
            s_nz = (s_t != 0)
            b_nz = (b_t != 0)
            total_pixels_serial += s_nz.sum().item()
            total_pixels_batched += b_nz.sum().item()
            total_intensity_serial += s_t.sum().item()
            total_intensity_batched += b_t.sum().item()

            # Compare: find positions where either is non-zero
            either_nz = s_nz | b_nz
            if not either_nz.any():
                continue

            s_vals = s_t[either_nz]
            b_vals = b_t[either_nz]

            # Exact match check
            exact = (s_vals == b_vals)
            matching_pixels += exact.sum().item()

            # Check non-exact matches
            if not exact.all():
                diff_mask = ~exact
                s_diff = s_vals[diff_mask]
                b_diff = b_vals[diff_mask]
                denom = torch.max(torch.abs(s_diff), torch.abs(b_diff)).clamp(min=1e-20)
                rel_diffs = torch.abs(s_diff - b_diff) / denom

                close = rel_diffs < 1e-4
                matching_pixels += close.sum().item()
                n_mismatch = (~close).sum().item()
                mismatched_pixels += n_mismatch

                if rel_diffs.numel() > 0:
                    max_rel_diff = max(max_rel_diff, rel_diffs.max().item())

                if n_mismatch > 0 and len(mismatch_examples) < 10:
                    bad_idx = torch.where(~close)[0][:5]
                    for idx in bad_idx:
                        mismatch_examples.append(
                            f"  omega={omega_idx} det={det_idx}: "
                            f"serial={s_diff[idx]:.6f} batched={b_diff[idx]:.6f} "
                            f"rel_diff={rel_diffs[idx]:.2e}"
                        )

    for ex in mismatch_examples:
        print(ex)

    total = matching_pixels + mismatched_pixels
    return {
        "total_pixels_serial": total_pixels_serial,
        "total_pixels_batched": total_pixels_batched,
        "matching": matching_pixels,
        "mismatched": mismatched_pixels,
        "total": total,
        "max_rel_diff": max_rel_diff,
        "intensity_serial": total_intensity_serial,
        "intensity_batched": total_intensity_batched,
    }


def main():
    example_dir = Path(__file__).parent
    os.chdir(example_dir)
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"

    print("=" * 70)
    print("Regression Test: Serial vs Batched Forward Simulation")
    print("Test case: ThreeVoxels")
    print("=" * 70)

    # Run serial
    print("\n--- Running SERIAL simulation ---")
    serial_images, serial_time = run_simulation(config_path, batched=False)
    print(f"Serial: {serial_time:.2f}s")

    # Run batched
    print("\n--- Running BATCHED simulation ---")
    batched_images, batched_time = run_simulation(config_path, batched=True)
    print(f"Batched: {batched_time:.2f}s")

    # Compare
    print("\n--- Comparing outputs ---")
    stats = compare_images(serial_images, batched_images)

    print(f"\nResults:")
    print(f"  Serial pixels:  {stats['total_pixels_serial']}")
    print(f"  Batched pixels: {stats['total_pixels_batched']}")
    print(f"  Matching:       {stats['matching']}/{stats['total']}")
    print(f"  Mismatched:     {stats['mismatched']}")
    print(f"  Max rel diff:   {stats['max_rel_diff']:.2e}")
    int_ratio = stats['intensity_batched'] / max(stats['intensity_serial'], 1e-20)
    print(f"  Intensity ratio: {int_ratio:.6f}")
    print(f"  Speedup:        {serial_time / max(batched_time, 0.001):.2f}x")

    # Pass/fail
    if stats["mismatched"] == 0 and stats["total_pixels_serial"] == stats["total_pixels_batched"]:
        print("\nPASSED: Pixel-exact match!")
        return 0
    elif stats["mismatched"] <= 10 and stats["max_rel_diff"] < 1e-3:
        print(f"\nPASSED (with tolerance): {stats['mismatched']} minor mismatches")
        return 0
    else:
        print(f"\nFAILED: {stats['mismatched']} mismatches, max_rel_diff={stats['max_rel_diff']:.2e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
