#!/usr/bin/env python3
"""
Regression test: serial vs batched on a subset of ManyGrains voxels.

Uses first N_TEST voxels to compare serial and batched output.
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
from icenine.sample import Sample
from icenine.image_data import ImageData
from icenine.simulation import Simulation

N_TEST = 100  # Number of voxels for comparison


def setup_simulation(config_path):
    """Set up simulation infrastructure (shared between runs)."""
    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = "500Grains.sim"
    simulator = ForwardSimulation(config)
    simulator.exp_setup.initialize_experiment()
    omega_ranges = simulator.exp_setup.get_omega_range_list()
    detector_list = simulator.exp_setup.get_detector_list()
    range_map = simulator.exp_setup.get_range_to_index_map()
    simulator.simulator = Simulation(simulator.exp_setup)
    sample = Sample()
    simulator.exp_setup.initialize_sample(sample, detector_list[0])
    return simulator, omega_ranges, detector_list, range_map, sample


def make_images(omega_ranges, detector_list):
    images = []
    for i in range(len(omega_ranges)):
        detector_images = []
        for det in detector_list:
            detector_images.append(ImageData(det.num_rows, det.num_cols))
        images.append(detector_images)
    return images


def compare_images(serial_images, batched_images):
    total_serial = 0
    total_batched = 0
    matching = 0
    mismatched = 0
    max_rel = 0.0

    for oi in range(len(serial_images)):
        for di in range(len(serial_images[oi])):
            s_t = serial_images[oi][di]._pixels_dense
            b_t = batched_images[oi][di]._pixels_dense
            s_nz = (s_t != 0)
            b_nz = (b_t != 0)
            total_serial += s_nz.sum().item()
            total_batched += b_nz.sum().item()
            either = s_nz | b_nz
            if not either.any():
                continue
            sv = s_t[either]
            bv = b_t[either]
            exact = (sv == bv)
            matching += exact.sum().item()
            if not exact.all():
                diff_mask = ~exact
                sv_d = sv[diff_mask]
                bv_d = bv[diff_mask]
                denom = torch.max(torch.abs(sv_d), torch.abs(bv_d)).clamp(min=1e-20)
                rd = torch.abs(sv_d - bv_d) / denom
                close = rd < 1e-4
                matching += close.sum().item()
                n_bad = (~close).sum().item()
                mismatched += n_bad
                if rd.numel() > 0:
                    max_rel = max(max_rel, rd.max().item())
                if n_bad > 0:
                    bad = torch.where(~close)[0][:3]
                    for b in bad:
                        print(f"    omega={oi} det={di}: s={sv_d[b]:.2f} b={bv_d[b]:.2f}")

    return total_serial, total_batched, matching, mismatched, max_rel


def main():
    example_dir = Path(__file__).parent
    os.chdir(example_dir)
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"

    print("=" * 70)
    print(f"Regression Test: Serial vs Batched (first {N_TEST} voxels)")
    print("Test case: ManyGrains")
    print("=" * 70)

    # Set up once
    sim, omega_ranges, det_list, range_map, sample = setup_simulation(config_path)
    mic = sample.get_mic()
    orig_voxels = mic.voxels[:]

    # Limit to N_TEST voxels
    mic.voxels = orig_voxels[:N_TEST]

    # Serial
    print(f"\n--- Serial ({N_TEST} voxels) ---")
    serial_images = make_images(omega_ranges, det_list)
    sim.images = serial_images
    t0 = time.time()
    sim._simulate_peaks(serial_images, det_list, sample, range_map)
    serial_time = time.time() - t0
    print(f"Serial: {serial_time:.2f}s")

    # Batched
    print(f"\n--- Batched ({N_TEST} voxels) ---")
    batched_images = make_images(omega_ranges, det_list)
    sim.images = batched_images
    t0 = time.time()
    sim._simulate_peaks_batched(batched_images, det_list, sample, range_map)
    batched_time = time.time() - t0
    print(f"Batched: {batched_time:.2f}s")

    # Compare
    print("\n--- Comparing ---")
    ts, tb, match, mismatch, max_rel = compare_images(serial_images, batched_images)
    total = match + mismatch
    print(f"  Serial pixels:  {ts}")
    print(f"  Batched pixels: {tb}")
    print(f"  Matching:       {match}/{total}")
    print(f"  Mismatched:     {mismatch}")
    print(f"  Max rel diff:   {max_rel:.2e}")
    print(f"  Speedup:        {serial_time / max(batched_time, 0.001):.2f}x")

    # Accept bin-boundary mismatches: same total pixel count, <0.1% mismatch rate
    mismatch_rate = mismatch / max(total, 1) * 100
    if ts == tb and mismatch_rate < 0.1:
        if mismatch == 0:
            print("\nPASSED: Pixel-exact match!")
        else:
            print(f"\nPASSED: {mismatch} bin-boundary mismatches ({mismatch_rate:.3f}%) — expected")

        # Now run full batched to verify no crash + benchmark
        mic.voxels = orig_voxels
        print(f"\n--- Full batched ({len(orig_voxels)} voxels) ---")
        full_images = make_images(omega_ranges, det_list)
        sim.images = full_images
        t0 = time.time()
        sim._simulate_peaks_batched(full_images, det_list, sample, range_map)
        full_time = time.time() - t0
        total_px = sum(
            (full_images[oi][di]._pixels_dense != 0).sum().item()
            for oi in range(len(full_images))
            for di in range(len(full_images[oi]))
        )
        print(f"  Time: {full_time:.2f}s ({full_time/60:.2f}min)")
        print(f"  Total pixels: {total_px}")
        print(f"  Per-voxel: {full_time/len(orig_voxels)*1000:.2f}ms")
        return 0
    else:
        print(f"\nFAILED: {mismatch} mismatches ({mismatch_rate:.3f}%)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
