#!/usr/bin/env python3
"""
Cost Function Validation — apples-to-apples C++ vs Python comparison.

Uses identical config (ReconstructBenchmark.config) and data (ScatteringData/)
as the C++ CostFunctionBenchmark to validate quality metric agreement.

Key differences from benchmark_cost_function.py:
  - Uses ReconstructBenchmark.config (not Example2.Simulation.config)
  - Reads from ScatteringData/ (not ScatteringData_Python/)
  - Passes eta_limit from config (86° → radians)
  - Prints per-peak diagnostic output for comparison with C++

Usage:
  cd icenine_py && uv run python benchmarks/validate_cost_function.py
"""

import math
import os
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
    count_qualified_peaks,
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
# Setup — IDENTICAL to C++ CostFunctionBenchmark
# =============================================================================

print("=== IceNine Python Cost Function Validation ===")
print("  (Matching C++ CostFunctionBenchmark configuration)")
print()

print("--- Initialization ---")

# Use SAME config as C++ benchmark
config_path = example_dir / "ConfigFiles" / "ReconstructBenchmark.config"
config = ConfigFile.from_file(str(config_path))

# The config has InfileBasename = "ScatteringData/3Grains.sim"
# We need to read from the SAME directory as C++
config.out_file_basename = "3Grains.sim"

exp_setup = XDMExperimentSetup(config)
exp_setup.initialize_experiment()
detector_list = exp_setup.get_detector_list()
range_map = exp_setup.get_range_to_index_map()

sample = Sample()
exp_setup.initialize_sample(sample, detector_list[0])
simulator = Simulation(exp_setup)
structure_list = sample.get_structure_list()

# Read from SAME scattering data directory as C++
exp_data = ExperimentalData.from_image_directory(
    directory=str(example_dir / "ScatteringData"),
    basename="3Grains.sim",
    ext="d",
    serial_length=5,
    n_omega=180,
    n_detectors=2,
    num_rows=2048,
    num_cols=2048,
)

# MaxQ=8 from config, eta_limit=86° from config (converted to radians)
MAX_Q = 8.0
eta_limit = config.eta_limit  # parsed as radians from "EtaLimit 86"

print(f"  Config: {config_path.name}")
print(f"  Data dir: ScatteringData/ (same as C++)")
print(f"  MaxQ: {MAX_Q}")
print(f"  eta_limit: {eta_limit:.6f} rad ({math.degrees(eta_limit):.1f}°)")
print(f"  Detectors: {len(detector_list)}")

cost_fn = VoxelCostFunction(
    simulator=simulator,
    detector_list=detector_list,
    range_map=range_map,
    exp_data=exp_data,
    sample=sample,
    structure_list=structure_list,
    mode="hard",
    max_q=MAX_Q,
    eta_limit=eta_limit,
)

# Load ground truth voxel 2
mic = MicFile.read(str(example_dir / "SimInput" / "three_voxels.mic"))
voxel = mic.voxels[2]
vertices = _get_voxel_vertices(voxel)
orientation = voxel.orientation
phase_index = voxel.phase

print(f"  Using voxel 2, phase={phase_index}")
print(f"  Voxel center: ({voxel.position[0]:.6f}, {voxel.position[1]:.7f}, {voxel.position[2]})")
print(f"  Orientation matrix:")
for r in range(3):
    print(f"    [{orientation[r, 0]:>10.5f}{orientation[r, 1]:>10.5f}{orientation[r, 2]:>10.5f} ]")

g_hkl_batch, g_mag_batch = cost_fn._phase_recip_vecs[phase_index]
print(f"  Reciprocal vectors: {len(g_hkl_batch)}")
print()

# =============================================================================
# Full VoxelCostFunction.evaluate() — aggregate result
# =============================================================================

print("--- Full VoxelCostFunction.evaluate() ---")

result = cost_fn.evaluate(
    orientation=orientation,
    voxel_vertices=vertices,
    phase_index=phase_index,
)

print(f"  quality={result.quality:.6f}  cost={result.cost:.6f}  hit_ratio={result.hit_ratio:.6f}")
print(f"  Pixels: overlap={result.pixel_overlap}  on_det={result.pixel_on_detector}")
print(f"  Peaks: overlap={result.peak_overlap}  on_det={result.peak_on_detector}")
print(f"  n_quality_points={result.n_quality_points}")
print()

# =============================================================================
# Per-peak diagnostic — compute observable peaks then evaluate each individually
# =============================================================================

print("--- Per-Peak Diagnostics ---")

orientation_t = torch.from_numpy(orientation).float()
g_lab_batch = (orientation_t @ g_hkl_batch.T).T

bragg_result = get_scattering_omegas_torch(
    g_lab_batch, g_mag_batch,
    simulator.beam_energy, simulator.beam_deflection_chi,
)

# Collect observable peaks (same logic as VoxelCostFunction.evaluate)
peak_omegas = []
peak_normals = []

obs_mask = bragg_result.observable
obs_indices = torch.where(obs_mask)[0].tolist()
beam_dir = simulator.beam_direction
base_rot = sample.sample_to_lab_matrix[:3, :3]

for idx in obs_indices:
    g_vec = g_lab_batch[idx]
    g_mag = torch.norm(g_vec)
    normal = g_vec / g_mag

    for omega_val in [bragg_result.omega1[idx].item(), bragg_result.omega2[idx].item()]:
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
            if eta >= eta_limit:
                continue

        peak_omegas.append(omega_val)
        peak_normals.append(normal)

print(f"  Observable peaks (after eta filter): {len(peak_omegas)}")
print()

# Evaluate each peak individually
print(f"  {'Peak':>4s}  {'omega_deg':>10s}  {'pix_ovlp':>8s}  {'pix_det':>8s}  "
      f"{'pk_ovlp':>7s}  {'pk_det':>6s}  {'n_det_ovlp':>10s}  {'quality_i':>10s}")
print(f"  {'----':>4s}  {'----------':>10s}  {'--------':>8s}  {'--------':>8s}  "
      f"{'-------':>7s}  {'------':>6s}  {'----------':>10s}  {'----------':>10s}")

for p_idx in range(len(peak_omegas)):
    single_omega = [peak_omegas[p_idx]]
    single_normal = [peak_normals[p_idx]]

    info = calculate_diffraction_overlap(
        sample=sample,
        voxel_vertices=vertices,
        peak_omegas=single_omega,
        peak_normals=single_normal,
        detector_list=detector_list,
        range_map=range_map,
        exp_data=exp_data,
        beam_direction=beam_dir,
        mode="hard",
    )

    omega_deg = math.degrees(peak_omegas[p_idx])

    # Compute per-peak quality contribution
    if info.pixel_on_detector > 0 and info.n_quality_points > 0:
        quality_i = info.quality
    else:
        quality_i = -1.0  # not counted

    print(f"  {p_idx:>4d}  {omega_deg:>10.4f}  {info.pixel_overlap:>8d}  {info.pixel_on_detector:>8d}  "
          f"{info.peak_overlap:>7d}  {info.peak_on_detector:>6d}  {info.n_quality_points:>10d}  "
          f"{quality_i:>10.6f}")

print()

# Recompute aggregate quality from per-peak data to verify
print("--- Aggregate Quality Recomputation ---")
running_quality = 0.0
n_points = 0
for p_idx in range(len(peak_omegas)):
    single_omega = [peak_omegas[p_idx]]
    single_normal = [peak_normals[p_idx]]

    info = calculate_diffraction_overlap(
        sample=sample,
        voxel_vertices=vertices,
        peak_omegas=single_omega,
        peak_normals=single_normal,
        detector_list=detector_list,
        range_map=range_map,
        exp_data=exp_data,
        beam_direction=beam_dir,
        mode="hard",
    )

    if info.pixel_on_detector > 0:
        n_det = len(detector_list)
        if info.n_quality_points > 0:
            # The single-peak info.quality already has the correct per-peak quality
            cur_quality = info.quality
        else:
            cur_quality = 0.0
        running_quality += (cur_quality - running_quality) / (n_points + 1)
        n_points += 1

print(f"  Recomputed quality: {running_quality:.6f}  (n_points={n_points})")
print(f"  VoxelCostFunction quality: {result.quality:.6f}  (n_points={result.n_quality_points})")
print()
print("=== Validation Complete ===")
