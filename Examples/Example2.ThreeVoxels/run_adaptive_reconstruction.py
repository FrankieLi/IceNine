"""
Benchmark: Python ADAPTIVE reconstruction on ThreeVoxels example.

Uses AdaptiveVoxelReconstructor (the port of C++ DiscreteRefinement)
for identical algorithm comparison with C++ mode 'r'.

Usage:
    cd icenine_py
    uv run python -u ../Examples/Example2.ThreeVoxels/run_adaptive_reconstruction.py

Uses the same config and experimental data as the C++ benchmark:
    ConfigFiles/ReconstructBenchmark.config
    ScatteringData/3Grains.sim*.d{0,1}
"""

import os
import sys
import time
import math

import numpy as np

# Ensure we're in the example directory (config uses relative paths)
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

from icenine.config_file import ConfigFile
from icenine.experimental_data import ExperimentalData
from icenine.reconstructor import (
    setup_reconstruction,
    AdaptiveVoxelReconstructor,
    _get_voxel_vertices,
)


def matrix_to_euler_zxz(R):
    """Convert rotation matrix to Bunge ZXZ Euler angles (degrees)."""
    Phi = math.acos(max(-1, min(1, R[2, 2])))
    if abs(math.sin(Phi)) > 1e-6:
        phi1 = math.atan2(R[0, 2], -R[1, 2])
        phi2 = math.atan2(R[2, 0], R[2, 1])
    else:
        phi1 = math.atan2(R[0, 1], R[0, 0])
        phi2 = 0.0
    return (math.degrees(phi1) % 360, math.degrees(Phi) % 360, math.degrees(phi2) % 360)


def misorientation_angle(R1, R2):
    """Compute misorientation angle in degrees (ignoring symmetry)."""
    dR = R1.T @ R2
    trace = np.clip(np.trace(dR), -1, 3)
    return math.degrees(math.acos(max(-1, min(1, (trace - 1) / 2))))


def main():
    print("=" * 70, flush=True)
    print("Python ADAPTIVE Reconstruction Benchmark — ThreeVoxels (MaxQ=8)", flush=True)
    print("=" * 70, flush=True)

    t_wall_start = time.time()

    # --- Load config ---
    config = ConfigFile.from_file("ConfigFiles/ReconstructBenchmark.config")
    print(f"Config: MaxQ={config.max_q}, MaxMCSteps={config.max_mc_steps}, "
          f"MaxLocalRes={config.max_local_resolution}, "
          f"LocalGridRadius={math.degrees(config.local_orientation_grid_radius):.1f}°, "
          f"MCRadiusScale={config.mc_radius_scale_factor}", flush=True)

    # --- Load experimental data (same C++ ScatteringData/) ---
    t0 = time.time()
    exp_data = ExperimentalData.from_image_directory(
        directory="ScatteringData",
        basename="3Grains.sim",
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
    )
    t_load = time.time() - t0
    print(f"Experimental data loaded: {t_load:.1f}s "
          f"({exp_data.n_omega_intervals}x{exp_data.n_detectors}, "
          f"{exp_data.count_bright_pixels()} bright px)", flush=True)

    # --- Setup reconstruction ---
    t0 = time.time()
    setup = setup_reconstruction(config, exp_data=exp_data)
    t_setup = time.time() - t0
    print(f"Setup: {t_setup:.1f}s ({len(setup.fz_orientations)} FZ orientations)",
          flush=True)

    # --- Save ground truth before reconstruction overwrites orientations ---
    mic = setup.sample.get_mic()
    gt_orientations = [v.orientation.copy() for v in mic.voxels]
    gt_eulers = [matrix_to_euler_zxz(R) for R in gt_orientations]

    print(f"\nGround truth:", flush=True)
    for i, e in enumerate(gt_eulers):
        print(f"  Voxel {i}: phase={mic.voxels[i].phase} "
              f"Euler=({e[0]:.2f}, {e[1]:.2f}, {e[2]:.2f})", flush=True)

    # --- Reconstruct using AdaptiveVoxelReconstructor ---
    print(f"\n{'=' * 70}", flush=True)
    print("ADAPTIVE Reconstruction (DiscreteRefinement port)", flush=True)
    print(f"{'=' * 70}", flush=True)

    reconstructor = AdaptiveVoxelReconstructor(setup)
    rng = np.random.default_rng(42)

    results = []
    per_voxel_stats = []
    t_recon_start = time.time()

    for idx, voxel in enumerate(mic.voxels):
        t_voxel = time.time()
        elapsed = t_voxel - t_recon_start
        print(f"\n  Voxel {idx}/{len(mic.voxels)} (phase={voxel.phase}, "
              f"elapsed={elapsed:.1f}s)", flush=True)

        vertices = _get_voxel_vertices(voxel)

        result = reconstructor.reconstruct_voxel(
            voxel_vertices=vertices,
            phase_index=voxel.phase,
            rng=rng,
        )

        voxel_time = time.time() - t_voxel
        global_evals, local_evals, total_evals = reconstructor.last_eval_counts

        results.append(result)

        oi = result.overlap_info
        hit = (oi.pixel_overlap / oi.pixel_on_detector
               if oi and oi.pixel_on_detector > 0 else 0)
        conf = (oi.peak_overlap / oi.peak_on_detector
                if oi and oi.peak_on_detector > 0 else 0)
        us_per_eval = (voxel_time * 1e6 / total_evals) if total_evals > 0 else 0

        per_voxel_stats.append({
            'time': voxel_time,
            'global_evals': global_evals,
            'local_evals': local_evals,
            'total_evals': total_evals,
            'us_per_eval': us_per_eval,
        })

        print(f"  Voxel {idx} done: cost={result.cost:.4f}, hit_ratio={hit:.4f}, "
              f"confidence={conf:.4f}, time={voxel_time:.1f}s, "
              f"evals={total_evals} (global={global_evals}, local={local_evals}), "
              f"us/eval={us_per_eval:.1f}", flush=True)

        # Update voxel orientation
        voxel.orientation = result.orientation

    t_recon = time.time() - t_recon_start

    # --- Per-voxel results ---
    print(f"\n{'=' * 70}", flush=True)
    print("RESULTS", flush=True)
    print(f"{'=' * 70}", flush=True)

    for i, result in enumerate(results):
        euler = matrix_to_euler_zxz(result.orientation)
        misori = misorientation_angle(gt_orientations[i], result.orientation)
        oi = result.overlap_info
        hit = (oi.pixel_overlap / oi.pixel_on_detector
               if oi and oi.pixel_on_detector > 0 else 0)
        conf = (oi.peak_overlap / oi.peak_on_detector
                if oi and oi.peak_on_detector > 0 else 0)
        stats = per_voxel_stats[i]

        print(f"\nVoxel {i}:", flush=True)
        print(f"  Cost:        {result.cost:.6f}", flush=True)
        print(f"  Hit ratio:   {hit:.4f}", flush=True)
        print(f"  Confidence:  {conf:.4f}", flush=True)
        if oi:
            print(f"  Pixels:      {oi.pixel_overlap}/{oi.pixel_on_detector}", flush=True)
            print(f"  Peaks:       {oi.peak_overlap}/{oi.peak_on_detector}", flush=True)
            print(f"  Quality:     {oi.quality:.6f}", flush=True)
        print(f"  Euler (recon): ({euler[0]:.2f}, {euler[1]:.2f}, {euler[2]:.2f})", flush=True)
        print(f"  Euler (truth): ({gt_eulers[i][0]:.2f}, {gt_eulers[i][1]:.2f}, "
              f"{gt_eulers[i][2]:.2f})", flush=True)
        print(f"  Misori:      {misori:.2f} deg", flush=True)
        print(f"  Time:        {stats['time']:.1f}s", flush=True)
        print(f"  Evals:       {stats['total_evals']} "
              f"(global={stats['global_evals']}, local={stats['local_evals']})", flush=True)
        print(f"  us/eval:     {stats['us_per_eval']:.1f}", flush=True)

    # --- Timing summary ---
    t_wall = time.time() - t_wall_start
    total_evals = sum(s['total_evals'] for s in per_voxel_stats)
    avg_us = (t_recon * 1e6 / total_evals) if total_evals > 0 else 0

    print(f"\n{'=' * 70}", flush=True)
    print("TIMING SUMMARY", flush=True)
    print(f"{'=' * 70}", flush=True)
    print(f"  Data loading:   {t_load:.1f}s", flush=True)
    print(f"  Setup:          {t_setup:.1f}s", flush=True)
    print(f"  Reconstruction: {t_recon:.1f}s  ({t_recon/len(results):.1f}s/voxel)",
          flush=True)
    print(f"  Total wall:     {t_wall:.1f}s", flush=True)
    print(f"  Total evals:    {total_evals}", flush=True)
    print(f"  Avg us/eval:    {avg_us:.1f}", flush=True)


if __name__ == "__main__":
    main()
