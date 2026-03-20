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
from typing import List, Tuple

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

    def update_quality(self, n_detectors_total: int) -> None:
        """
        Update quality metric using incremental (Welford-style) running mean.

        quality_i = (pixel_overlap / pixel_on_detector) * (detectors_overlap / n_detectors)
        quality = running_mean(quality_0, quality_1, ..., quality_n)

        C++ Reference: OverlapInfo.h SOverlapInfo::UpdateQuality
        """
        if self.pixel_on_detector == 0 or n_detectors_total == 0:
            return

        pixel_ratio = self.pixel_overlap / self.pixel_on_detector
        det_ratio = self.detectors_overlap / n_detectors_total
        cur_quality = pixel_ratio * det_ratio

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

            detector_lit[det_idx] = True

            # Count pixel overlap with experimental image
            n_overlap, n_lit = exp_image.get_triangle_overlap_property(
                projected[0], projected[1], projected[2],
                mode=mode,
            )

            n_overlap_int = int(n_overlap.item())
            n_lit_int = int(n_lit.item())

            peak_pixel_overlap += n_overlap_int
            peak_pixel_on_det += n_lit_int

            if n_overlap_int > 0:
                spot_overlap[det_idx] = True

        # Qualified peak counting
        peak_on_det, peak_ovlp, n_det_ovlp = count_qualified_peaks(
            detector_lit, spot_overlap, n_detectors
        )

        # Update overlap info
        overlap_info.detectors_overlap = n_det_ovlp
        overlap_info.update_counts(
            peak_pixel_overlap, peak_pixel_on_det,
            peak_ovlp, peak_on_det,
        )
        overlap_info.update_quality(n_detectors)

    # Restore original sample orientation
    sample.sample_to_lab_matrix = orig_matrix

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
        """
        self.simulator = simulator
        self.detector_list = detector_list
        self.range_map = range_map
        self.exp_data = exp_data
        self.sample = sample
        self.structure_list = structure_list
        self.mode = mode

        # Pre-compute reciprocal vectors per phase
        self._phase_recip_vecs = {}
        for phase_idx, structure in enumerate(structure_list):
            recp_vecs = structure.get_reflection_vectors()
            if recp_vecs:
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

        # Collect observable peaks
        peak_omegas = []
        peak_normals = []

        obs_mask = bragg_result.observable
        if not obs_mask.any():
            return OverlapInfo()

        obs_indices = torch.where(obs_mask)[0].tolist()
        for idx in obs_indices:
            g_vec = g_lab_batch[idx]
            g_mag = torch.norm(g_vec)
            normal = g_vec / g_mag  # scattering direction (unit)

            # Both omega solutions
            peak_omegas.append(bragg_result.omega1[idx].item())
            peak_normals.append(normal)

            peak_omegas.append(bragg_result.omega2[idx].item())
            peak_normals.append(normal)

        # Calculate diffraction overlap
        return calculate_diffraction_overlap(
            sample=self.sample,
            voxel_vertices=voxel_vertices,
            peak_omegas=peak_omegas,
            peak_normals=peak_normals,
            detector_list=self.detector_list,
            range_map=self.range_map,
            exp_data=self.exp_data,
            beam_direction=self.simulator.beam_direction,
            mode=self.mode,
        )
