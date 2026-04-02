#!/usr/bin/env python3
"""
Performance benchmark: ExperimentalImageStack vs sequential image access.

Measures:
  1. Stack construction time
  2. Sequential get_image() vs batch tensor lookup
  3. Sequential pixel lookup vs batch grid_sample
  4. End-to-end Stage D: hard overlap vs centroid sampling

Usage:
  cd icenine_py && uv run python benchmarks/bench_image_stack.py
"""

import os
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

# Change to example directory (configs use relative paths)
project_root = Path(__file__).parent.parent.parent
example_dir = project_root / "Examples" / "Example2.ThreeVoxels"
os.chdir(example_dir)

from icenine.config_file import ConfigFile
from icenine.cost_functions import VoxelCostFunction
from icenine.differentiable_cost import ExperimentalImageStack, MultiScaleImageStack
from icenine.experiment_setup import XDMExperimentSetup
from icenine.experimental_data import ExperimentalData
from icenine.mic_file import MicFile
from icenine.reconstructor import _get_voxel_vertices
from icenine.sample import Sample
from icenine.simulation import Simulation


def benchmark(name: str, n_iters: int, fn, warmup: int = 3):
    """Run fn n_iters times, report mean in microseconds."""
    for _ in range(min(warmup, n_iters)):
        fn()

    times = []
    for _ in range(n_iters):
        t0 = time.perf_counter()
        fn()
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1e6)

    arr = np.array(times)
    print(
        f"  {name:<50s}  mean={arr.mean():>10.1f} us"
        f"  min={arr.min():>10.1f} us  (n={n_iters})"
    )
    return arr.mean()


# =============================================================================
# Setup
# =============================================================================

print("=== ExperimentalImageStack Performance Benchmark ===")
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

print("  Loading experimental data...")
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
print(f"  Loaded: {exp_data.n_omega_intervals} omega x {exp_data.n_detectors} det")
print()

# =============================================================================
# Benchmark 1: Stack Construction Time
# =============================================================================

print("--- Benchmark 1: Stack Construction ---")

t0 = time.perf_counter()
image_stack = ExperimentalImageStack(exp_data, binary=True)
t_construct = time.perf_counter() - t0
print(f"  ExperimentalImageStack (binary=True): {t_construct:.3f}s")
print(f"  {image_stack}")

t0 = time.perf_counter()
multi_stack = MultiScaleImageStack(image_stack, downsample_factors=[1, 4, 8])
t_multi = time.perf_counter() - t0
print(f"  MultiScaleImageStack (3 scales, max_pool): {t_multi:.3f}s")
print(f"  {multi_stack}")
print()

# =============================================================================
# Benchmark 2: Sequential vs Batch Image Lookup
# =============================================================================

print("--- Benchmark 2: Image Lookup (M=200 random images) ---")

rng = np.random.RandomState(42)
M = 200
omega_indices = rng.randint(0, 180, size=M)
det_indices = rng.randint(0, 2, size=M)

# Sequential: get_image() + get_binary_numpy()
def seq_lookup():
    results = []
    for i in range(M):
        img = exp_data.get_image(int(omega_indices[i]), int(det_indices[i]))
        results.append(img.get_binary_numpy())
    return results

mean_seq = benchmark("Sequential get_image() + get_binary_numpy()", 20, seq_lookup)

# Batch: tensor indexing
flat_indices = torch.tensor(
    omega_indices * 2 + det_indices, dtype=torch.long
)

def batch_lookup():
    return image_stack.get_images_batch(flat_indices)

mean_batch = benchmark("Batch tensor indexing", 20, batch_lookup)

print(f"  Speedup: {mean_seq / mean_batch:.1f}x")
print()

# =============================================================================
# Benchmark 3: Sequential Pixel Lookup vs Batch grid_sample
# =============================================================================

print("--- Benchmark 3: Pixel Sampling (M=200 points) ---")

# Generate random pixel coordinates within image bounds
centroids_col = rng.uniform(10, 2038, size=M).astype(np.float32)
centroids_row = rng.uniform(10, 2038, size=M).astype(np.float32)

# Sequential: integer pixel lookup from binary numpy
binary_images = seq_lookup()  # pre-load

def seq_pixel_lookup():
    total = 0
    for i in range(M):
        cx = int(centroids_col[i])
        cy = int(centroids_row[i])
        total += binary_images[i][cy, cx]
    return total

mean_seq_px = benchmark("Sequential binary_img[cy, cx]", 50, seq_pixel_lookup)

# Batch: grid_sample
batch_images = image_stack.get_images_batch(flat_indices)  # (M, 1, H, W)
H, W = 2048, 2048
grid_x = torch.from_numpy(2.0 * centroids_col / W - 1.0).unsqueeze(0).unsqueeze(0)
grid_y = torch.from_numpy(2.0 * centroids_row / H - 1.0).unsqueeze(0).unsqueeze(0)

# For batch grid_sample, we need per-image grids: (M, 1, 1, 2)
grid_coords = torch.stack([
    torch.from_numpy(2.0 * centroids_col / W - 1.0),
    torch.from_numpy(2.0 * centroids_row / H - 1.0),
], dim=1).unsqueeze(1).unsqueeze(1)  # (M, 1, 1, 2)

def batch_grid_sample():
    return F.grid_sample(
        batch_images, grid_coords,
        mode='bilinear', padding_mode='zeros', align_corners=False
    )

mean_batch_px = benchmark("Batch F.grid_sample (bilinear)", 50, batch_grid_sample)

print(f"  Speedup: {mean_seq_px / mean_batch_px:.1f}x")
print()

# =============================================================================
# Benchmark 4: End-to-end — hard evaluate() vs centroid grid_sample
# =============================================================================

print("--- Benchmark 4: End-to-End Cost Function Eval ---")

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

# Hard evaluate
def hard_eval():
    return cost_fn.evaluate(orientation, vertices, phase_index)

info = hard_eval()
print(f"  Hard eval reference: quality={info.quality:.4f}, cost={info.cost:.4f}, "
      f"hit_ratio={info.hit_ratio:.4f}")

mean_hard = benchmark("VoxelCostFunction.evaluate() [hard]", 20, hard_eval)

# Centroid grid_sample equivalent (manual pipeline)
# We replicate what DifferentiableCostFunction will do:
# Stages A-C from cost_fn, then centroid sampling instead of Stage D
def centroid_eval():
    """Evaluate using centroid sampling against image stack."""
    # Stages A-C: same as VoxelCostFunction.evaluate() up to pixel coords
    g_hkl_batch, g_mag_batch = cost_fn._phase_recip_vecs[phase_index]
    orientation_t = torch.from_numpy(orientation).float()
    g_lab_batch = (orientation_t @ g_hkl_batch.T).T

    from icenine.diffraction_core import get_scattering_omegas_torch
    bragg_result = get_scattering_omegas_torch(
        g_lab_batch, g_mag_batch,
        cost_fn.simulator.beam_energy,
        cost_fn.simulator.beam_deflection_chi,
    )

    obs_mask = bragg_result.observable
    if not obs_mask.any():
        return 0.0

    beam_dir = cost_fn.simulator.beam_direction
    base_rot = cost_fn.sample.sample_to_lab_matrix[:3, :3]

    obs_g = g_lab_batch[obs_mask]
    obs_mag = torch.norm(obs_g, dim=1, keepdim=True)
    obs_normals = obs_g / obs_mag

    obs_omega1 = bragg_result.omega1[obs_mask]
    obs_omega2 = bragg_result.omega2[obs_mask]
    all_omegas = torch.cat([obs_omega1, obs_omega2])
    all_normals = torch.cat([obs_normals, obs_normals])
    N = all_omegas.shape[0]

    cos_w = torch.cos(all_omegas)
    sin_w = torch.sin(all_omegas)
    Rz = torch.zeros(N, 3, 3)
    Rz[:, 0, 0] = cos_w; Rz[:, 0, 1] = -sin_w
    Rz[:, 1, 0] = sin_w; Rz[:, 1, 1] = cos_w
    Rz[:, 2, 2] = 1.0

    full_rot = Rz @ base_rot
    lab_normals = torch.bmm(full_rot, all_normals.unsqueeze(-1)).squeeze(-1)
    dot = (beam_dir.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)
    reflected = beam_dir.unsqueeze(0) - 2.0 * dot * lab_normals

    # Eta filter
    rd_norms = torch.norm(reflected, dim=1)
    safe_norms = torch.where(rd_norms > 0, rd_norms, torch.ones_like(rd_norms))
    ry = torch.abs(reflected[:, 1]) / safe_norms
    rz_val = torch.abs(reflected[:, 2]) / safe_norms
    eta = torch.atan2(ry, rz_val)
    valid = (eta < cost_fn.eta_limit) & (rd_norms > 0)

    peak_omegas_t = all_omegas[valid]
    peak_normals_t = all_normals[valid]
    M_peaks = peak_omegas_t.shape[0]

    if M_peaks == 0:
        return 0.0

    # Stage A: omega-to-wedge mapping
    low = cost_fn.range_map.low
    width = cost_fn.range_map.width
    bin_indices = ((peak_omegas_t.numpy() - low) / width).astype(int)
    index_list = cost_fn.range_map.index_list
    n_bins = len(index_list)

    wedge_indices = np.full(M_peaks, -1, dtype=np.int64)
    for i in range(M_peaks):
        n = bin_indices[i]
        if 0 <= n < n_bins and index_list[n] is not None:
            wedge_indices[i] = index_list[n]

    valid_mask = wedge_indices >= 0
    if not np.any(valid_mask):
        return 0.0

    vi = np.where(valid_mask)[0]
    v_omegas = peak_omegas_t[vi]
    v_normals = peak_normals_t[vi]
    v_wedge = wedge_indices[vi]
    M = len(vi)

    # Stage B: batch geometry
    orig_matrix = cost_fn.sample.sample_to_lab_matrix
    cos_w2 = torch.cos(v_omegas); sin_w2 = torch.sin(v_omegas)
    Rz2 = torch.zeros(M, 3, 3)
    Rz2[:, 0, 0] = cos_w2; Rz2[:, 0, 1] = -sin_w2
    Rz2[:, 1, 0] = sin_w2; Rz2[:, 1, 1] = cos_w2
    Rz2[:, 2, 2] = 1.0

    base_rot2 = orig_matrix[:3, :3]
    full_rot2 = Rz2 @ base_rot2
    lab_normals2 = torch.bmm(full_rot2, v_normals.unsqueeze(-1)).squeeze(-1)
    dot2 = (beam_dir.unsqueeze(0) * lab_normals2).sum(dim=1, keepdim=True)
    reflected2 = beam_dir.unsqueeze(0) - 2.0 * dot2 * lab_normals2

    full_4x4 = torch.zeros(M, 4, 4)
    full_4x4[:, :3, :3] = full_rot2
    full_4x4[:, :3, 3] = orig_matrix[:3, 3]
    full_4x4[:, 3, 3] = 1.0

    verts_4d = torch.cat([vertices, torch.ones(3, 1)], dim=1)
    lab_verts = torch.einsum('mij,vj->mvi', full_4x4, verts_4d)[:, :, :3]

    # Stage C: ray-detector intersection (for all detectors)
    n_det = len(detector_list)
    quality_sum = 0.0
    n_quality = 0

    for det_idx in range(n_det):
        detector = detector_list[det_idx]
        plane = detector._detector_plane
        plane_n = plane.normal
        plane_d = plane.d

        origins = lab_verts.reshape(M * 3, 3)
        dirs = reflected2.unsqueeze(1).expand(M, 3, 3).reshape(M * 3, 3)

        denom = (dirs * plane_n).sum(dim=1)
        numer = -((origins * plane_n).sum(dim=1) + plane_d)
        parallel = torch.abs(denom) < 1e-8
        safe_denom = torch.where(parallel, torch.ones_like(denom), denom)
        t = torch.where(parallel, torch.zeros_like(denom), numer / safe_denom)
        hits = (~parallel) & (t > 0)

        pts = origins + t.unsqueeze(1) * dirs

        relative = pts - detector._position
        pixel_loc = relative - detector._lab_frame_coord_origin
        j = (pixel_loc * detector._lab_frame_basis_j).sum(dim=1)
        k = (pixel_loc * detector._lab_frame_basis_k).sum(dim=1)

        col = (j + detector.pixel_half_width) / detector.pixel_width
        row = (k + detector.pixel_half_height) / detector.pixel_height

        hits_mv = hits.reshape(M, 3)
        all_hit = hits_mv.all(dim=1)
        cols_mv = col.reshape(M, 3)
        rows_mv = row.reshape(M, 3)

        # Centroid of triangle vertices
        centroids_col = cols_mv.mean(dim=1)  # (M,)
        centroids_row = rows_mv.mean(dim=1)  # (M,)

        # Build flat indices for this detector
        flat_idx = torch.from_numpy(v_wedge.astype(np.int64)) * n_det + det_idx

        # Gather images
        batch_imgs = image_stack.images[flat_idx]  # (M, 1, H, W)

        # Normalize to [-1, 1] for grid_sample
        grid_x = (2.0 * centroids_col / image_stack.W - 1.0)
        grid_y = (2.0 * centroids_row / image_stack.H - 1.0)
        grid = torch.stack([grid_x, grid_y], dim=1).unsqueeze(1).unsqueeze(1)  # (M, 1, 1, 2)

        sampled = F.grid_sample(
            batch_imgs, grid, mode='bilinear',
            padding_mode='zeros', align_corners=False
        ).squeeze()  # (M,)

        # Mask by all_hit
        sampled = sampled * all_hit.float()

        quality_sum += sampled.sum().item()
        n_quality += all_hit.sum().item()

    quality = quality_sum / max(n_quality, 1)
    return quality

centroid_quality = centroid_eval()
print(f"  Centroid eval reference: quality={centroid_quality:.4f}")

mean_centroid = benchmark("Centroid grid_sample eval (manual)", 20, centroid_eval)

print()
print("=" * 70)
print(f"  Hard eval:     {mean_hard:>10.1f} us/eval")
print(f"  Centroid eval: {mean_centroid:>10.1f} us/eval")
print(f"  Speedup:       {mean_hard / mean_centroid:.1f}x")
print("=" * 70)
