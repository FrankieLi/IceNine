"""
Cost functions for reconstruction — overlap computation between simulated
voxel projections and experimental detector images.

Ports the C++ OverlapInfo / VoxelCostFunction from Src/OverlapInfo.h,
Src/OverlapInfo.tmpl.cpp, and Src/CostFunctions.h.

The main entry point is `VoxelCostFunction.evaluate()`, which:
1. Computes observable Bragg peaks for a candidate orientation
2. For each peak, rotates sample to the corresponding omega angle
3. Projects voxel triangle onto all detectors
4. Counts pixel overlap with experimental data
5. Aggregates into a quality metric (running mean across peaks)
"""

import math
from dataclasses import dataclass
from typing import List, Tuple, Union

import numpy as np
import torch

from .crystal_structure import CrystalStructure
from .detector import Detector
from .diffraction_core import (
    build_reflected_ray,
    get_illuminated_pixel,
    get_scattering_omegas_torch,
)
from .experimental_data import ExperimentalData
from .sample import Sample
from .simulation import Simulation
from .simulation_range import SimulationRange

try:
    from ._rasterize import stage_d_overlap as _c_stage_d_overlap
    _HAS_C_RASTERIZE = True
except ImportError:
    _HAS_C_RASTERIZE = False


# ---------------------------------------------------------------------------
# OverlapInfo — aggregate overlap metrics
# ---------------------------------------------------------------------------

@dataclass
class OverlapInfo:
    """
    Aggregate overlap metrics between simulation and experiment.

    C++ Reference: Src/OverlapInfo.h SOverlapInfo
    """

    pixel_overlap: int = 0
    pixel_on_detector: int = 0
    peak_overlap: int = 0
    peak_on_detector: int = 0
    detectors_overlap: int = 0
    quality: float = 0.0
    n_quality_points: int = 0

    def update_counts(
        self,
        peak_pixel_overlap: int,
        peak_pixel_on_detector: int,
        peak_overlap_flag: int,
        peak_on_detector_flag: int,
    ) -> None:
        """
        Accumulate counts from a single peak evaluation.

        C++ Reference: OverlapInfo.h SOverlapInfo::UpdateCounts
        """
        self.pixel_overlap += peak_pixel_overlap
        self.pixel_on_detector += peak_pixel_on_detector
        self.peak_overlap += peak_overlap_flag
        self.peak_on_detector += peak_on_detector_flag

    def update_quality(
        self,
        peak_pixel_overlap: int,
        peak_pixel_on_detector: int,
        peak_detectors_overlap: int,
        n_detectors_total: int,
    ) -> None:
        """
        Update quality metric using incremental (Welford-style) running mean.

        Uses per-peak values (not accumulated totals) to compute the per-peak
        quality contribution, matching the C++ UpdateQuality(oRHS, ...) signature
        where oRHS is a separate per-peak OverlapInfo.

        quality_i = (peak_pixel_overlap / peak_pixel_on_detector) * (peak_det_overlap / n_det)
        quality = running_mean(quality_0, quality_1, ..., quality_n)

        C++ Reference: OverlapInfo.h SOverlapInfo::UpdateQuality
        """
        if peak_pixel_on_detector == 0 or n_detectors_total == 0:
            return

        pixel_ratio = peak_pixel_overlap / peak_pixel_on_detector
        if peak_detectors_overlap > 0:
            det_ratio = peak_detectors_overlap / n_detectors_total
            cur_quality = pixel_ratio * det_ratio
        else:
            cur_quality = 0.0

        # Incremental mean: q = q + (cur - q) / (n + 1)
        self.quality += (cur_quality - self.quality) / (self.n_quality_points + 1)
        self.n_quality_points += 1

    @property
    def cost(self) -> float:
        """Cost = 1 - quality. Lower is better."""
        return 1.0 - self.quality

    @property
    def hit_ratio(self) -> float:
        """Fraction of simulated pixels overlapping experimental data."""
        if self.pixel_on_detector == 0:
            return 0.0
        return self.pixel_overlap / self.pixel_on_detector

    @property
    def confidence(self) -> float:
        """
        Confidence metric combining hit ratio and quality.

        C++ Reference: CostFunctions.h GetConfidence
        """
        return self.quality


# ---------------------------------------------------------------------------
# Qualified peak counting — contiguous detector validation
# ---------------------------------------------------------------------------

def count_qualified_peaks(
    detector_lit: List[bool],
    spot_overlap: List[bool],
    n_detectors: int,
) -> Tuple[int, int, int]:
    """
    Check if a peak has valid (contiguous) detector coverage.

    A peak is "qualified" as on-detector if it lights up a contiguous set
    of detectors starting from detector 0 (allowing trailing gap).
    An overlap is "qualified" if overlapping detectors are also contiguous.

    Args:
        detector_lit: Per-detector flag: True if simulated peak hits this detector
        spot_overlap: Per-detector flag: True if simulated peak overlaps experiment
        n_detectors: Total number of detectors

    Returns:
        (peak_on_detector, peak_overlap, n_detectors_overlap):
            peak_on_detector: 1 if peak is validly on detector, 0 otherwise
            peak_overlap: 1 if overlap is valid, 0 otherwise
            n_detectors_overlap: count of detectors with valid overlap

    C++ Reference: OverlapInfo.tmpl.cpp CountQualifiedPeaks
    """
    # Check if peak is valid: contiguous detector lit pattern
    # Valid: [True, True, False, False], [True, True, True, True]
    # Invalid: [False, True, True], [True, False, True]
    valid_sim_peak = False
    valid_overlap = False
    n_det_overlap = 0

    if n_detectors == 0:
        return 0, 0, 0

    # Check contiguous lit pattern (must start at detector 0)
    if detector_lit[0]:
        valid_sim_peak = True
        gap_seen = False
        for i in range(1, n_detectors):
            if not detector_lit[i]:
                gap_seen = True
            elif gap_seen:
                # Detector lit after a gap → invalid
                valid_sim_peak = False
                break

    if not valid_sim_peak:
        return 0, 0, 0

    # Check contiguous overlap pattern
    if spot_overlap[0]:
        valid_overlap = True
        n_det_overlap = 1
        gap_seen = False
        for i in range(1, n_detectors):
            if spot_overlap[i]:
                if gap_seen and detector_lit[i]:
                    # Overlap after gap where detector was lit → invalid
                    valid_overlap = False
                    break
                n_det_overlap += 1
            elif detector_lit[i]:
                gap_seen = True
    elif not detector_lit[0]:
        # detector 0 not lit at all — cannot have valid overlap
        pass
    else:
        # detector 0 lit but no overlap — still might be valid if later has overlap
        # C++ logic: bValidOverlap remains false unless bSpotOverlap[0] is true
        pass

    if not valid_overlap:
        n_det_overlap = 0

    return (
        1 if valid_sim_peak else 0,
        1 if valid_overlap else 0,
        n_det_overlap,
    )


# ---------------------------------------------------------------------------
# Diffraction overlap calculation
# ---------------------------------------------------------------------------

def calculate_diffraction_overlap(
    sample: Sample,
    voxel_vertices: torch.Tensor,
    peak_omegas: List[float],
    peak_normals: List[torch.Tensor],
    detector_list: List[Detector],
    range_map: SimulationRange,
    exp_data: ExperimentalData,
    beam_direction: torch.Tensor,
    mode: str = 'hard',
) -> OverlapInfo:
    """
    Calculate overlap between simulated diffraction peaks and experimental data.

    For each observable peak:
    1. Look up omega index in the simulation range
    2. Rotate sample to the omega angle
    3. Project voxel triangle onto all detectors
    4. Count pixel overlap with experimental images
    5. Validate contiguous detector coverage
    6. Update running quality metric

    Args:
        sample: Sample object (will be temporarily rotated per-peak)
        voxel_vertices: Voxel triangle vertices in sample frame, shape (3, 3)
        peak_omegas: List of omega angles (radians) for observable peaks
        peak_normals: List of scattering direction vectors (unit, sample frame)
        detector_list: List of detectors
        range_map: SimulationRange for omega-to-index mapping
        exp_data: ExperimentalData with experimental images
        beam_direction: Beam direction vector
        mode: 'hard' (discrete counting) or 'soft' (differentiable)

    Returns:
        OverlapInfo with accumulated overlap metrics

    C++ Reference: OverlapInfo.tmpl.cpp CalculateDiffractionOverlap
    """
    n_detectors = len(detector_list)
    overlap_info = OverlapInfo()

    # Save original sample orientation for restoration
    orig_matrix = sample.sample_to_lab_matrix.clone()

    for peak_idx in range(len(peak_omegas)):
        omega = peak_omegas[peak_idx]
        normal = peak_normals[peak_idx]

        # Map omega to experimental data index
        omega_index = range_map.angle_to_wedge_index(omega)
        if omega_index is None:
            continue

        # Rotate sample to this omega angle
        sample.sample_to_lab_matrix = orig_matrix.clone()
        sample.rotate_z(omega)

        # Transform scattering direction to lab frame
        rot_matrix = sample.sample_to_lab_matrix[:3, :3]
        lab_normal = rot_matrix @ normal

        # Compute reflected ray direction
        reflected_dir = beam_direction - 2.0 * torch.dot(beam_direction, lab_normal) * lab_normal

        # Per-detector tracking for qualified peak counting
        detector_lit = [False] * n_detectors
        spot_overlap = [False] * n_detectors
        peak_pixel_overlap = 0
        peak_pixel_on_det = 0

        # Project voxel vertices onto each detector using existing helpers
        # Use build_reflected_ray + get_illuminated_pixel (same as forward sim)
        for det_idx in range(n_detectors):
            detector = detector_list[det_idx]
            exp_image = exp_data.get_image(omega_index, det_idx)

            # Project all 3 vertices
            projected = []
            all_hit = True
            for vi in range(3):
                vertex = voxel_vertices[vi]
                reflected_ray = build_reflected_ray(sample, vertex, reflected_dir)
                hit, pixel_col, pixel_row = get_illuminated_pixel(
                    detector, reflected_ray
                )

                if not hit.item():
                    all_hit = False
                    break

                projected.append(
                    torch.tensor([pixel_col.item(), pixel_row.item()],
                                 dtype=torch.float32)
                )

            if not all_hit or len(projected) < 3:
                continue

            # Count pixel overlap with experimental image
            n_overlap, n_lit = exp_image.get_triangle_overlap_property(
                projected[0], projected[1], projected[2],
                mode=mode,
            )

            n_overlap_int = int(n_overlap.item())
            n_lit_int = int(n_lit.item())

            peak_pixel_overlap += n_overlap_int
            peak_pixel_on_det += n_lit_int

            # Match C++ ShapePixelOverlapCounter<SVoxel> (OverlapInfo.h:498-517):
            # detector_lit is set only when overlap > 0 OR triangle has in-bounds
            # pixels (n_lit > 0, equivalent to C++ IsInBound check on vertices).
            if n_overlap_int > 0:
                detector_lit[det_idx] = True
                spot_overlap[det_idx] = True
            elif n_lit_int > 0:
                detector_lit[det_idx] = True

        # Qualified peak counting
        peak_on_det, peak_ovlp, n_det_ovlp = count_qualified_peaks(
            detector_lit, spot_overlap, n_detectors
        )

        # Update overlap info — accumulate counts, then update quality with per-peak values
        # C++ calls UpdateCounts(oCurOverlapInfo) then UpdateQuality(oCurOverlapInfo, ...)
        overlap_info.detectors_overlap = n_det_ovlp
        overlap_info.update_counts(
            peak_pixel_overlap, peak_pixel_on_det,
            peak_ovlp, peak_on_det,
        )
        overlap_info.update_quality(
            peak_pixel_overlap, peak_pixel_on_det,
            n_det_ovlp, n_detectors,
        )

    # Restore original sample orientation
    sample.sample_to_lab_matrix = orig_matrix

    return overlap_info


def calculate_diffraction_overlap_batched(
    sample: Sample,
    voxel_vertices: torch.Tensor,
    peak_omegas: "Union[List[float], torch.Tensor]",
    peak_normals: "Union[List[torch.Tensor], torch.Tensor]",
    detector_list: List[Detector],
    range_map: SimulationRange,
    exp_data: ExperimentalData,
    beam_direction: torch.Tensor,
    mode: str = 'hard',
    pixel_radius: int = 0,
) -> OverlapInfo:
    """
    Batched version of calculate_diffraction_overlap.

    Replaces the per-peak Python loop for geometric computations (Z-rotation,
    normal transform, reflection, vertex projection, ray-plane intersection)
    with vectorized tensor operations. Only the per-peak overlap counting
    (which accesses different experimental images) remains sequential.

    When pixel_radius > 0, uses pixel-based peak overlap (matching C++
    PixelBasedPeakOverlapCounter) instead of full triangle rasterization.
    This checks a ±pixel_radius square around the projected center point
    for brightness, providing tolerance for small orientation errors during
    discrete search.

    C++ Reference: OverlapInfo.h PixelBasedPeakOverlapCounter (pixel_radius > 0)
                   OverlapInfo.h ShapePixelOverlapCounter (pixel_radius == 0)
    """
    n_detectors = len(detector_list)
    overlap_info = OverlapInfo()
    n_peaks = len(peak_omegas)

    if n_peaks == 0:
        return overlap_info

    # ---- Stage A: Vectorized omega filtering ----
    # Accept both tensor and list inputs for backward compatibility
    if isinstance(peak_omegas, torch.Tensor):
        omegas = peak_omegas
    else:
        omegas = torch.tensor(peak_omegas, dtype=torch.float32)
    if isinstance(peak_normals, torch.Tensor):
        normals = peak_normals
    else:
        normals = torch.stack(peak_normals)  # (N, 3)

    # Map omega to wedge index using range_map internals
    low = range_map.low
    width = range_map.width
    bin_indices = ((omegas.numpy() - low) / width).astype(int)

    # Look up wedge indices, filtering invalid bins
    wedge_indices = np.full(n_peaks, -1, dtype=np.int64)
    index_list = range_map.index_list
    n_bins = len(index_list)
    for i in range(n_peaks):
        n = bin_indices[i]
        if 0 <= n < n_bins and index_list[n] is not None:
            wedge_indices[i] = index_list[n]

    valid_mask = wedge_indices >= 0
    if not np.any(valid_mask):
        return overlap_info

    valid_indices = np.where(valid_mask)[0]
    v_omegas = omegas[valid_indices]          # (M,)
    v_normals = normals[valid_indices]        # (M, 3)
    v_wedge = wedge_indices[valid_indices]    # (M,)
    M = len(valid_indices)

    # ---- Stage B: Batch geometric computations ----
    orig_matrix = sample.sample_to_lab_matrix  # (4, 4) tensor

    # Build M rotation matrices Rz(omega) — (M, 3, 3)
    cos_w = torch.cos(v_omegas)
    sin_w = torch.sin(v_omegas)
    Rz = torch.zeros(M, 3, 3)
    Rz[:, 0, 0] = cos_w
    Rz[:, 0, 1] = -sin_w
    Rz[:, 1, 0] = sin_w
    Rz[:, 1, 1] = cos_w
    Rz[:, 2, 2] = 1.0

    base_rot = orig_matrix[:3, :3]       # (3, 3)
    full_rot = Rz @ base_rot             # (M, 3, 3)

    # Transform normals to lab frame: (M, 3)
    lab_normals = torch.bmm(full_rot, v_normals.unsqueeze(-1)).squeeze(-1)

    # Reflected directions: beam - 2*(beam·n)*n — (M, 3)
    dot = (beam_direction.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)  # (M, 1)
    reflected = beam_direction.unsqueeze(0) - 2.0 * dot * lab_normals  # (M, 3)

    # Build full 4x4 matrices for vertex transform — (M, 4, 4)
    full_4x4 = torch.zeros(M, 4, 4)
    full_4x4[:, :3, :3] = full_rot
    full_4x4[:, :3, 3] = orig_matrix[:3, 3]  # translation unchanged by Rz
    full_4x4[:, 3, 3] = 1.0

    # Transform 3 vertices for all M peaks: (3, 4) homogeneous
    verts_4d = torch.cat([voxel_vertices, torch.ones(3, 1)], dim=1)  # (3, 4)
    # einsum: for each peak m, for each vertex v, compute full_4x4[m] @ verts_4d[v]
    lab_verts = torch.einsum('mij,vj->mvi', full_4x4, verts_4d)[:, :, :3]  # (M, 3, 3)

    # ---- Stage C: Batch ray-detector intersection per detector ----
    # For each detector, compute pixel coordinates for all M peaks × 3 vertices
    # Store results: per_det_hit[det_idx] = (M,) bool, per_det_pixels[det_idx] = (M, 3, 2)
    per_det_all_hit = []
    per_det_pixels = []

    for det_idx in range(n_detectors):
        detector = detector_list[det_idx]
        plane = detector._detector_plane

        plane_n = plane.normal   # (3,)
        plane_d = plane.d        # scalar

        # Ray origins: lab_verts (M, 3, 3) reshaped to (M*3, 3)
        origins = lab_verts.reshape(M * 3, 3)
        # Ray directions: reflected (M, 3) expanded to (M, 3, 3) then (M*3, 3)
        dirs = reflected.unsqueeze(1).expand(M, 3, 3).reshape(M * 3, 3)

        # t = -(n·origin + d) / (n·dir)
        denom = (dirs * plane_n).sum(dim=1)                     # (M*3,)
        numer = -((origins * plane_n).sum(dim=1) + plane_d)     # (M*3,)
        parallel = torch.abs(denom) < 1e-8
        safe_denom = torch.where(parallel, torch.ones_like(denom), denom)
        t = torch.where(parallel, torch.zeros_like(denom), numer / safe_denom)
        hits = (~parallel) & (t > 0)

        # Intersection points → pixel coordinates (inline detector math)
        pts = origins + t.unsqueeze(1) * dirs  # (M*3, 3)

        # lab_to_detector_coordinate inlined:
        relative = pts - detector._position       # (M*3, 3)
        pixel_loc = relative - detector._lab_frame_coord_origin  # (M*3, 3)
        j = (pixel_loc * detector._lab_frame_basis_j).sum(dim=1)  # (M*3,)
        k = (pixel_loc * detector._lab_frame_basis_k).sum(dim=1)  # (M*3,)

        # to_col_pixel / to_row_pixel inlined:
        col = (j + detector.pixel_half_width) / detector.pixel_width
        row = (k + detector.pixel_half_height) / detector.pixel_height

        # Reshape to (M, 3): per peak, per vertex
        # Match serial: detector_lit uses ray-plane hit only (no bounds check)
        hits_mv = hits.reshape(M, 3)
        all_hit = hits_mv.all(dim=1)  # (M,) — peaks where all 3 vertices intersect plane

        cols_mv = col.reshape(M, 3)
        rows_mv = row.reshape(M, 3)

        # Store pixel coords as (M, 3, 2) — [col, row]
        pixels = torch.stack([cols_mv, rows_mv], dim=2)  # (M, 3, 2)
        per_det_all_hit.append(all_hit)
        per_det_pixels.append(pixels)

    # ---- Stage D: Per-peak overlap accumulation ----
    # NOTE: Stage D intentionally breaks autograd. The overlap check is
    # inherently non-differentiable (binary pixel test, integer counting,
    # contiguity validation). Stages A-C preserve the autograd graph;
    # if gradient-based refinement is ever needed, differentiate through
    # those stages and use a differentiable renderer for image comparison.

    if _HAS_C_RASTERIZE and mode == 'hard':
        # ---- Fast path: batch C extension ----
        # Build compact image lookup — only unique wedges needed
        v_wedge_np = np.asarray(v_wedge, dtype=np.int32)
        unique_wedges = np.unique(v_wedge_np)
        wedge_to_local = np.full(int(v_wedge_np.max()) + 1, -1, dtype=np.int32)
        for i, w in enumerate(unique_wedges):
            wedge_to_local[w] = i
        local_wedge = wedge_to_local[v_wedge_np]  # (M,) re-indexed

        # Build image list: [local_wedge * n_det + det] -> binary numpy
        image_list = []
        for w in unique_wedges:
            for d in range(n_detectors):
                image_list.append(exp_data.get_image(int(w), d).get_binary_numpy())

        # Build det_hit_mask: (M, n_det) uint8
        det_hit = np.stack(
            [per_det_all_hit[d].numpy().astype(np.uint8) for d in range(n_detectors)],
            axis=1,
        )  # (M, n_det)

        # Get image dimensions from first image
        first_image = exp_data.get_image(int(unique_wedges[0]), 0)
        img_rows = first_image.num_rows
        img_cols = first_image.num_cols

        if pixel_radius > 0:
            # pixel_centers: (M, n_det, 2) int32 — first vertex as center
            centers = np.stack(
                [per_det_pixels[d][:, 0, :].numpy().astype(np.int32)
                 for d in range(n_detectors)],
                axis=1,
            )  # (M, n_det, 2) — [col, row]
            centers = np.ascontiguousarray(centers)
            result = _c_stage_d_overlap(
                image_list, local_wedge, det_hit,
                centers, None,
                n_detectors, pixel_radius, M, img_rows, img_cols,
            )
        else:
            # pixel_vertices: (M, n_det, 3, 2) float32
            verts = np.stack(
                [per_det_pixels[d].numpy().astype(np.float32)
                 for d in range(n_detectors)],
                axis=1,
            )  # (M, n_det, 3, 2)
            verts = np.ascontiguousarray(verts)
            result = _c_stage_d_overlap(
                image_list, local_wedge, det_hit,
                None, verts,
                n_detectors, pixel_radius, M, img_rows, img_cols,
            )

        # Unpack results
        overlap_info.pixel_overlap = result[0]
        overlap_info.pixel_on_detector = result[1]
        overlap_info.peak_overlap = result[2]
        overlap_info.peak_on_detector = result[3]
        overlap_info.quality = result[4]
        overlap_info.n_quality_points = result[5]
        overlap_info.detectors_overlap = result[6]

    else:
        # ---- Fallback: sequential Python loop ----
        for p in range(M):
            wedge_idx = int(v_wedge[p])

            detector_lit = [False] * n_detectors
            spot_overlap = [False] * n_detectors
            peak_pixel_overlap = 0
            peak_pixel_on_det = 0

            for det_idx in range(n_detectors):
                if not per_det_all_hit[det_idx][p].item():
                    continue

                exp_image = exp_data.get_image(wedge_idx, det_idx)

                if pixel_radius > 0:
                    pixels = per_det_pixels[det_idx][p]
                    cx = int(pixels[0, 0].item())
                    cy = int(pixels[0, 1].item())

                    found_in_bounds = False
                    found_bright = False
                    num_rows = exp_image.num_rows
                    num_cols = exp_image.num_cols
                    binary_img = exp_image.get_binary_numpy()
                    for dx in range(-pixel_radius, pixel_radius + 1):
                        if found_bright:
                            break
                        for dy in range(-pixel_radius, pixel_radius + 1):
                            px, py = cx + dx, cy + dy
                            if 0 <= px < num_cols and 0 <= py < num_rows:
                                found_in_bounds = True
                                if binary_img[py, px]:
                                    found_bright = True
                                    break

                    if found_in_bounds:
                        detector_lit[det_idx] = True
                        peak_pixel_on_det += 1
                    if found_bright:
                        peak_pixel_overlap += 1
                        spot_overlap[det_idx] = True
                else:
                    pixels = per_det_pixels[det_idx][p]
                    v0 = pixels[0]
                    v1 = pixels[1]
                    v2 = pixels[2]

                    n_overlap, n_lit = exp_image.get_triangle_overlap_property(
                        v0, v1, v2, mode=mode,
                    )

                    n_overlap_int = int(n_overlap.item())
                    n_lit_int = int(n_lit.item())

                    peak_pixel_overlap += n_overlap_int
                    peak_pixel_on_det += n_lit_int

                    if n_overlap_int > 0:
                        detector_lit[det_idx] = True
                        spot_overlap[det_idx] = True
                    elif n_lit_int > 0:
                        detector_lit[det_idx] = True

            peak_on_det, peak_ovlp, n_det_ovlp = count_qualified_peaks(
                detector_lit, spot_overlap, n_detectors
            )

            overlap_info.detectors_overlap = n_det_ovlp
            overlap_info.update_counts(
                peak_pixel_overlap, peak_pixel_on_det,
                peak_ovlp, peak_on_det,
            )
            overlap_info.update_quality(
                peak_pixel_overlap, peak_pixel_on_det,
                n_det_ovlp, n_detectors,
            )

    return overlap_info


# ---------------------------------------------------------------------------
# VoxelCostFunction — high-level cost evaluator
# ---------------------------------------------------------------------------

class VoxelCostFunction:
    """
    Evaluate how well a candidate voxel orientation matches experimental data.

    Given a voxel with a candidate crystal orientation:
    1. Get reciprocal lattice vectors for the voxel's crystal phase
    2. Compute observable Bragg peaks (via diffraction_core)
    3. For each peak, project voxel onto detectors and count overlap
    4. Return OverlapInfo with quality metric

    C++ Reference: Src/CostFunctions.h VoxelCostFunction
    """

    def __init__(
        self,
        simulator: Simulation,
        detector_list: List[Detector],
        range_map: SimulationRange,
        exp_data: ExperimentalData,
        sample: Sample,
        structure_list: List[CrystalStructure],
        mode: str = 'hard',
        eta_limit: float = math.pi / 2.0,
        pixel_radius: int = 0,
        max_q: float = 0.0,
    ):
        """
        Args:
            simulator: Initialized Simulation object
            detector_list: List of detectors
            range_map: SimulationRange for omega-to-index mapping
            exp_data: ExperimentalData with experimental images
            sample: Sample object (used for rotations)
            structure_list: List of CrystalStructure (one per phase)
            mode: 'hard' for discrete counting, 'soft' for differentiable
            eta_limit: Maximum eta angle for peak acceptance (radians).
                       Peaks with eta >= eta_limit are skipped, matching C++
                       XDMEtaAcceptFn in Simulation.tmpl.cpp GetProjectedVertices.
            pixel_radius: Pixel radius for overlap check expansion. When > 0,
                          uses pixel-based peak overlap (±pixel_radius square)
                          instead of triangle rasterization. Used for discrete
                          search to widen the cost landscape.
                          C++ Reference: OverlapInfo.h PixelBasedPeakOverlapCounter
            max_q: Maximum Q for reciprocal vector filtering (Å⁻¹). When > 0,
                   only reciprocal vectors with |q| <= max_q are used. When 0,
                   uses all vectors. C++ uses nQMaxDiscrete=5 for discrete search.
                   C++ Reference: DiscreteSearch.h:298-300
        """
        self.simulator = simulator
        self.detector_list = detector_list
        self.range_map = range_map
        self.exp_data = exp_data
        self.sample = sample
        self.structure_list = structure_list
        self.mode = mode
        self.eta_limit = eta_limit
        self.pixel_radius = pixel_radius

        # Pre-compute reciprocal vectors per phase
        self._phase_recip_vecs = {}
        for phase_idx, structure in enumerate(structure_list):
            recp_vecs = structure.get_reflection_vectors()
            if recp_vecs:
                # Filter by max_q if specified (C++ nQMaxDiscrete)
                if max_q > 0:
                    recp_vecs = [rv for rv in recp_vecs if rv.q_mag <= max_q]
                if not recp_vecs:
                    continue
                g_hkl_batch = torch.stack(
                    [torch.from_numpy(rv.q_vec).float() for rv in recp_vecs]
                )
                g_mag_batch = torch.tensor(
                    [rv.q_mag for rv in recp_vecs], dtype=torch.float32
                )
                self._phase_recip_vecs[phase_idx] = (g_hkl_batch, g_mag_batch)

    def evaluate(
        self,
        orientation: np.ndarray,
        voxel_vertices: torch.Tensor,
        phase_index: int = 0,
    ) -> OverlapInfo:
        """
        Evaluate cost for a candidate voxel orientation.

        Args:
            orientation: 3x3 rotation matrix (crystal → sample frame)
            voxel_vertices: Triangle vertices in sample frame, shape (3, 3)
            phase_index: Crystal phase index

        Returns:
            OverlapInfo with quality, cost, hit_ratio, etc.

        C++ Reference: CostFunctions.h VoxelCostFunction::operator()
        """
        if phase_index not in self._phase_recip_vecs:
            return OverlapInfo()

        g_hkl_batch, g_mag_batch = self._phase_recip_vecs[phase_index]

        # NOTE: Do NOT call set_orientation_matrix() here.
        # The orientation is applied only to reciprocal vectors (g_lab = orient @ g_hkl).
        # The sample-to-lab matrix stays as the base (identity) matrix, matching the
        # forward simulation where vertices are in the sample frame and only undergo
        # Rz(omega) rotation, not crystal orientation rotation.

        # Transform reciprocal vectors to sample frame
        orientation_t = torch.from_numpy(orientation).float()
        g_lab_batch = (orientation_t @ g_hkl_batch.T).T

        # Compute Bragg angles for all reflections
        bragg_result = get_scattering_omegas_torch(
            g_lab_batch,
            g_mag_batch,
            self.simulator.beam_energy,
            self.simulator.beam_deflection_chi,
        )

        # Collect observable peaks using vectorized eta filtering
        obs_mask = bragg_result.observable
        if not obs_mask.any():
            return OverlapInfo()

        beam_dir = self.simulator.beam_direction
        base_rot = self.sample.sample_to_lab_matrix[:3, :3]

        # Extract observable peaks as tensors — no Python loop needed
        obs_g = g_lab_batch[obs_mask]  # (K, 3)
        obs_mag = torch.norm(obs_g, dim=1, keepdim=True)  # (K, 1)
        obs_normals = obs_g / obs_mag  # (K, 3) unit vectors

        # Expand: each observable peak has 2 omega solutions → 2K rows
        obs_omega1 = bragg_result.omega1[obs_mask]  # (K,)
        obs_omega2 = bragg_result.omega2[obs_mask]  # (K,)
        all_omegas = torch.cat([obs_omega1, obs_omega2])  # (2K,)
        all_normals = torch.cat([obs_normals, obs_normals])  # (2K, 3)
        N = all_omegas.shape[0]

        # Batch rotation matrices Rz(omega) — (N, 3, 3)
        cos_w = torch.cos(all_omegas)
        sin_w = torch.sin(all_omegas)
        Rz = torch.zeros(N, 3, 3)
        Rz[:, 0, 0] = cos_w
        Rz[:, 0, 1] = -sin_w
        Rz[:, 1, 0] = sin_w
        Rz[:, 1, 1] = cos_w
        Rz[:, 2, 2] = 1.0

        # Batch: full_rot = Rz @ base_rot, lab_normals = full_rot @ normals
        full_rot = Rz @ base_rot  # (N, 3, 3)
        lab_normals = torch.bmm(full_rot, all_normals.unsqueeze(-1)).squeeze(-1)  # (N, 3)

        # Batch reflection: reflected = beam - 2*(beam·n)*n
        dot = (beam_dir.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)  # (N, 1)
        reflected = beam_dir.unsqueeze(0) - 2.0 * dot * lab_normals  # (N, 3)

        # Batch eta filter: eta = atan2(|ry|/|r|, |rz|/|r|) < eta_limit
        rd_norms = torch.norm(reflected, dim=1)  # (N,)
        safe_norms = torch.where(rd_norms > 0, rd_norms, torch.ones_like(rd_norms))
        ry = torch.abs(reflected[:, 1]) / safe_norms
        rz_val = torch.abs(reflected[:, 2]) / safe_norms
        eta = torch.atan2(ry, rz_val)  # (N,)
        valid = (eta < self.eta_limit) & (rd_norms > 0)  # (N,) bool

        peak_omegas_t = all_omegas[valid]  # (M,) tensor
        peak_normals_t = all_normals[valid]  # (M, 3) tensor

        if peak_omegas_t.shape[0] == 0:
            return OverlapInfo()

        # Calculate diffraction overlap (batched for performance)
        return calculate_diffraction_overlap_batched(
            sample=self.sample,
            voxel_vertices=voxel_vertices,
            peak_omegas=peak_omegas_t,
            peak_normals=peak_normals_t,
            detector_list=self.detector_list,
            range_map=self.range_map,
            exp_data=self.exp_data,
            beam_direction=self.simulator.beam_direction,
            mode=self.mode,
            pixel_radius=self.pixel_radius,
        )
