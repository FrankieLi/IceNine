#!/usr/bin/env python3
"""
Cost Function Benchmark — per-operation profiling for C++/Python comparison.

Measures the same 6 operations as Src/CostFunctionBenchmark.cpp:
  1. Single GetScatteringOmegas call (+ batched version)
  2. GetObservablePeaks equivalent (omega + eta filtering)
  3. Single peak overlap: 0 detectors hit
  4. Single peak overlap: 1 detector hit
  5. Single peak overlap: 2 detectors hit
  6. Full VoxelCostFunction.evaluate() (single voxel)

Usage:
  cd icenine_py && uv run python benchmarks/benchmark_cost_function.py
"""

import math
import os
import time
from pathlib import Path

import numpy as np
import torch

# Change to example directory (configs use relative paths)
project_root = Path(__file__).parent.parent.parent
example_dir = project_root / "Examples" / "Example2.ThreeVoxels"
os.chdir(example_dir)

from icenine.config_file import ConfigFile
from icenine.cost_functions import (
    OverlapInfo,
    VoxelCostFunction,
    calculate_diffraction_overlap,
    calculate_diffraction_overlap_batched,
)
from icenine.diffraction_core import (
    get_scattering_omegas_torch,
    build_reflected_ray,
    get_illuminated_pixel,
)
from icenine.experiment_setup import XDMExperimentSetup
from icenine.experimental_data import ExperimentalData
from icenine.mic_file import MicFile
from icenine.reconstructor import _get_voxel_vertices
from icenine.sample import Sample
from icenine.simulation import Simulation


# =============================================================================
# Timer helper
# =============================================================================


def benchmark(name: str, n_iters: int, fn, warmup: int = 3):
    """Run fn n_iters times, report mean/stddev/min in microseconds."""
    for _ in range(min(warmup, n_iters)):
        fn()

    times = []
    for _ in range(n_iters):
        t0 = time.perf_counter()
        fn()
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1e6)  # to microseconds

    arr = np.array(times)
    mean = arr.mean()
    std = arr.std()
    mn = arr.min()

    print(
        f"  {name:<45s}  mean={mean:>10.2f} us"
        f"  stddev={std:>8.2f} us  min={mn:>10.2f} us  (n={n_iters})"
    )
    return {"mean_us": mean, "stddev_us": std, "min_us": mn, "n": n_iters}


# =============================================================================
# Setup
# =============================================================================

print("=== IceNine Python Cost Function Benchmark ===")
print()

print("--- Initialization ---")

config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
config = ConfigFile.from_file(str(config_path))
config.out_file_basename = "3Grains.sim"

exp_setup = XDMExperimentSetup(config)
exp_setup.initialize_experiment()
detector_list = exp_setup.get_detector_list()
range_map = exp_setup.get_range_to_index_map()

sample = Sample()
exp_setup.initialize_sample(sample, detector_list[0])
simulator = Simulation(exp_setup)
structure_list = sample.get_structure_list()

exp_data = ExperimentalData.from_image_directory(
    directory=str(example_dir / "ScatteringData_Python"),
    basename="3Grains.sim",
    ext="d",
    serial_length=5,
    n_omega=180,
    n_detectors=2,
    num_rows=2048,
    num_cols=2048,
)

# MaxQ=8 Å⁻¹ to match C++ ReconstructBenchmark.config
MAX_Q = 8.0

cost_fn = VoxelCostFunction(
    simulator=simulator,
    detector_list=detector_list,
    range_map=range_map,
    exp_data=exp_data,
    sample=sample,
    structure_list=structure_list,
    mode="hard",
    max_q=MAX_Q,
)

# Load ground truth voxel 2
mic = MicFile.read(str(example_dir / "SimInput" / "three_voxels.mic"))
voxel = mic.voxels[2]
vertices = _get_voxel_vertices(voxel)
orientation = voxel.orientation
phase_index = voxel.phase

print(f"  Config loaded, data files read.")
print(f"  Detectors: {len(detector_list)}")
print(f"  Using voxel 2, phase={phase_index}")
print(f"  Voxel center: ({voxel.position[0]:.6f}, {voxel.position[1]:.7f}, {voxel.position[2]})")
print(f"  Orientation matrix:")
for r in range(3):
    print(f"    [{orientation[r, 0]:>10.5f}{orientation[r, 1]:>10.5f}{orientation[r, 2]:>10.5f} ]")

# Get reciprocal vectors for this phase
g_hkl_batch, g_mag_batch = cost_fn._phase_recip_vecs[phase_index]
print(f"  Reciprocal vectors: {len(g_hkl_batch)}")

# Transform to sample frame
orientation_t = torch.from_numpy(orientation).float()
g_lab_batch = (orientation_t @ g_hkl_batch.T).T

print()

# =============================================================================
# Benchmark 1: Single GetScatteringOmegas
# =============================================================================

print("--- Benchmark 1: Single GetScatteringOmegas ---")

# Single reciprocal vector
g_single = g_lab_batch[0:1]
g_mag_single = g_mag_batch[0:1]

benchmark(
    "get_scattering_omegas_torch (single)",
    10000,
    lambda: get_scattering_omegas_torch(
        g_single,
        g_mag_single,
        simulator.beam_energy,
        simulator.beam_deflection_chi,
    ),
)

# All reciprocal vectors (batched)
benchmark(
    f"get_scattering_omegas_torch (all {len(g_hkl_batch)} recip vecs)",
    1000,
    lambda: get_scattering_omegas_torch(
        g_lab_batch,
        g_mag_batch,
        simulator.beam_energy,
        simulator.beam_deflection_chi,
    ),
)
print()

# =============================================================================
# Benchmark 2: Observable peaks (omega calc + eta filtering)
# =============================================================================

print("--- Benchmark 2: Observable peaks (omega + eta filtering) ---")


def compute_observable_peaks():
    """Same as VoxelCostFunction.evaluate lines 711-760: omega + eta filter."""
    bragg = get_scattering_omegas_torch(
        g_lab_batch,
        g_mag_batch,
        simulator.beam_energy,
        simulator.beam_deflection_chi,
    )
    peak_omegas = []
    peak_normals = []
    obs_mask = bragg.observable
    if not obs_mask.any():
        return peak_omegas, peak_normals
    obs_indices = torch.where(obs_mask)[0].tolist()
    beam_dir = simulator.beam_direction
    base_rot = sample.sample_to_lab_matrix[:3, :3]

    for idx in obs_indices:
        g_vec = g_lab_batch[idx]
        g_mag = torch.norm(g_vec)
        normal = g_vec / g_mag
        for omega_val in [bragg.omega1[idx].item(), bragg.omega2[idx].item()]:
            cos_w = math.cos(omega_val)
            sin_w = math.sin(omega_val)
            rz = torch.tensor(
                [[cos_w, -sin_w, 0], [sin_w, cos_w, 0], [0, 0, 1]],
                dtype=torch.float32,
            )
            rot = rz @ base_rot
            lab_normal = rot @ normal
            reflected = beam_dir - 2.0 * torch.dot(beam_dir, lab_normal) * lab_normal
            rd_norm = torch.norm(reflected)
            if rd_norm > 0:
                ry = abs(reflected[1].item()) / rd_norm.item()
                rz_val = abs(reflected[2].item()) / rd_norm.item()
                eta = math.atan2(ry, rz_val)
                if eta >= cost_fn.eta_limit:
                    continue
            peak_omegas.append(omega_val)
            peak_normals.append(normal)
    return peak_omegas, peak_normals


# Run once to get peak lists for subsequent benchmarks
peak_omegas, peak_normals = compute_observable_peaks()

benchmark("Observable peaks (omega + eta filter)", 100, compute_observable_peaks)
print(f"  Observable peaks generated: {len(peak_omegas)}")
print()

# =============================================================================
# Benchmark 3-5: Per-peak overlap by detector hit count
# =============================================================================

print("--- Benchmark 3-5: Per-peak overlap (categorized by detector hits) ---")

# Classify peaks by how many detectors they hit
orig_matrix = sample.sample_to_lab_matrix.clone()
beam_dir = simulator.beam_direction

n_hit_0, n_hit_1, n_hit_2 = 0, 0, 0
peak_0_det, peak_1_det, peak_2_det = -1, -1, -1

for p_idx in range(len(peak_omegas)):
    omega = peak_omegas[p_idx]
    normal = peak_normals[p_idx]
    omega_index = range_map.angle_to_wedge_index(omega)
    if omega_index is None:
        continue

    sample.sample_to_lab_matrix = orig_matrix.clone()
    sample.rotate_z(omega)
    rot_matrix = sample.sample_to_lab_matrix[:3, :3]
    lab_normal = rot_matrix @ normal
    reflected_dir = beam_dir - 2.0 * torch.dot(beam_dir, lab_normal) * lab_normal

    dets_hit = 0
    for det_idx in range(len(detector_list)):
        all_hit = True
        for vi in range(3):
            reflected_ray = build_reflected_ray(sample, vertices[vi], reflected_dir)
            hit, _, _ = get_illuminated_pixel(detector_list[det_idx], reflected_ray)
            if not hit.item():
                all_hit = False
                break
        if all_hit:
            dets_hit += 1

    if dets_hit == 0 and peak_0_det < 0:
        peak_0_det = p_idx
    if dets_hit == 1 and peak_1_det < 0:
        peak_1_det = p_idx
    if dets_hit == 2 and peak_2_det < 0:
        peak_2_det = p_idx
    if dets_hit == 0:
        n_hit_0 += 1
    elif dets_hit == 1:
        n_hit_1 += 1
    else:
        n_hit_2 += 1

sample.sample_to_lab_matrix = orig_matrix

print(f"  Peak classification: 0-det={n_hit_0} 1-det={n_hit_1} 2-det={n_hit_2}")


def bench_single_peak_overlap(label, peak_idx, n_iters):
    """Benchmark a single peak through calculate_diffraction_overlap."""
    if peak_idx < 0:
        print(f"  {label:<45s}  SKIPPED (no such peak)")
        return
    single_omega = [peak_omegas[peak_idx]]
    single_normal = [peak_normals[peak_idx]]
    benchmark(
        label,
        n_iters,
        lambda: calculate_diffraction_overlap(
            sample=sample,
            voxel_vertices=vertices,
            peak_omegas=single_omega,
            peak_normals=single_normal,
            detector_list=detector_list,
            range_map=range_map,
            exp_data=exp_data,
            beam_direction=beam_dir,
            mode="hard",
        ),
    )


bench_single_peak_overlap("Single peak: 0 detectors hit", peak_0_det, 1000)
bench_single_peak_overlap("Single peak: 1 detector hit", peak_1_det, 1000)
bench_single_peak_overlap("Single peak: 2 detectors hit", peak_2_det, 1000)

# All peaks (serial path)
benchmark(
    f"calc_diffraction_overlap (all {len(peak_omegas)} peaks, serial)",
    10,
    lambda: calculate_diffraction_overlap(
        sample=sample,
        voxel_vertices=vertices,
        peak_omegas=peak_omegas,
        peak_normals=peak_normals,
        detector_list=detector_list,
        range_map=range_map,
        exp_data=exp_data,
        beam_direction=beam_dir,
        mode="hard",
    ),
)

# All peaks (batched path)
benchmark(
    f"calc_diffraction_overlap_batched (all {len(peak_omegas)} peaks)",
    10,
    lambda: calculate_diffraction_overlap_batched(
        sample=sample,
        voxel_vertices=vertices,
        peak_omegas=peak_omegas,
        peak_normals=peak_normals,
        detector_list=detector_list,
        range_map=range_map,
        exp_data=exp_data,
        beam_direction=beam_dir,
        mode="hard",
        pixel_radius=0,
    ),
)
print()

# =============================================================================
# Benchmark 6: Full VoxelCostFunction.evaluate()
# =============================================================================

print("--- Benchmark 6: Full VoxelCostFunction (single voxel) ---")

result = [None]


def run_full_eval():
    result[0] = cost_fn.evaluate(
        orientation=orientation,
        voxel_vertices=vertices,
        phase_index=phase_index,
    )


benchmark("VoxelCostFunction.evaluate()", 10, run_full_eval)

info = result[0]
print(f"  Result: quality={info.quality:.2f}  cost={info.cost:.2f}  hit_ratio={info.hit_ratio:.2f}")
print(f"  Pixels: overlap={info.pixel_overlap}  on_det={info.pixel_on_detector}")
print(f"  Peaks: overlap={info.peak_overlap}  on_det={info.peak_on_detector}")
print()

print("=== Benchmark Complete ===")
