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
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch

from .crystal_structure import CrystalStructure
from .detector import Detector
from .diffraction_core import (
    BraggConstants,
    build_reflected_ray,
    get_illuminated_pixel,
    get_scattering_omegas_torch,
    scattering_omegas_fixed_beam,
)
from .experimental_data import ExperimentalData
from .ray_overlap import RayRule, apply_cases_python, ray_flags, ray_flags_c
from .sample import Sample
from .simulation import Simulation
from .simulation_range import SimulationRange

try:
    from ._rasterize import stage_d_overlap as _c_stage_d_overlap
    _HAS_C_RASTERIZE = True
except ImportError:
    _HAS_C_RASTERIZE = False
try:
    from ._rasterize import HAS_CPP_RULE as _C_HAS_CPP_RULE
except ImportError:
    _C_HAS_CPP_RULE = False  # extension missing or built before the cpp_rule argument existed
try:
    from ._rasterize import HAS_RAY_RULE as _C_HAS_RAY_RULE
except ImportError:
    _C_HAS_RAY_RULE = False
try:
    from ._rasterize import HAS_VERTEX_LIT as _C_VERTEX_LIT_LEVEL

    # 2 = the exact C++ vertex-lit rule (1 was the superset rule of the first Task A build)
    _C_HAS_VERTEX_LIT = _C_VERTEX_LIT_LEVEL >= 2
except ImportError:
    _C_HAS_VERTEX_LIT = False  # extension missing or built before the exact cpp_edge_rules
try:
    from ._rasterize import HAS_FUSED as _C_FUSED_LEVEL
    from ._rasterize import fused_context as _c_fused_context
    from ._rasterize import fused_evaluate as _c_fused_evaluate

    # 1 = the fused evaluate of the FusedCost key (S2, _fused_eval.c)
    _C_HAS_FUSED = _C_FUSED_LEVEL >= 1
except ImportError:
    _C_HAS_FUSED = False  # extension missing, or built before (or without) _fused_eval.c
_WARNED_STALE_EXT = False
# Observed-ray flags in C (False = the per-peak Python reference ray_overlap.ray_flags; tests).
USE_C_RAY_FLAGS = True


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
    # Only with a ray rule (ray_overlap.RayRule): counts of the cases a, b, c, d, e, f_dark,
    # f_bright and of detector-0 overlaps, over this evaluation's peaks. None otherwise.
    ray_cases: Optional[Tuple[int, ...]] = None

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
        Fraction of peaks with overlap (peak-level hit ratio).

        C++ Reference: CostFunctions.cpp GetConfidence
        Returns: peak_overlap / peak_on_detector
        """
        if self.peak_on_detector == 0:
            return 0.0
        return self.peak_overlap / self.peak_on_detector


# ---------------------------------------------------------------------------
# Qualified peak counting — contiguous detector validation
# ---------------------------------------------------------------------------

def count_qualified_peaks(
    detector_lit: List[bool],
    spot_overlap: List[bool],
    n_detectors: int,
    cpp_rule: bool = False,
) -> Tuple[int, int, int]:
    """
    Check if a peak has valid (contiguous) detector coverage.

    With ``cpp_rule=True`` the literal C++ rule is applied instead (``count_qualified_peaks_cpp``):
    a later detector that is lit but shows no overlap makes the overlap invalid. The default
    (False) is the original Python rule.

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
    if cpp_rule:
        return count_qualified_peaks_cpp(detector_lit, spot_overlap, n_detectors)

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


def count_qualified_peaks_cpp(
    detector_lit: List[bool],
    spot_overlap: List[bool],
    n_detectors: int,
) -> Tuple[int, int, int]:
    """
    Literal port of C++ ``Utilities::CountQualifiedPeaks`` (OverlapInfo.tmpl.cpp:78-139).

    Peak on detector: detector 0 is lit; the first unlit later detector is tolerated once
    (``bOnPreviousDetector``), and afterwards a further unlit detector invalidates the peak, while
    a lit one does not (lit, unlit, lit is valid; lit, unlit, unlit is not). Overlap: detector 0
    overlaps and every later detector either overlaps or is unlit; one unlit, non-overlapping
    detector is tolerated (``bPreviousOverlap``), and a lit non-overlapping one invalidates the
    overlap. The count of overlapping detectors is reported only when both are valid.

    Returns:
        (peak_on_detector, peak_overlap, n_detectors_overlap)
    """
    if n_detectors == 0:
        return 0, 0, 0
    valid_sim = bool(detector_lit[0])
    on_previous = bool(detector_lit[0])
    i = 1
    while i < n_detectors and valid_sim:
        if detector_lit[i]:
            valid_sim = True
        elif on_previous:
            valid_sim = True
            on_previous = False
        else:
            valid_sim = False
        i += 1

    previous_overlap = bool(spot_overlap[0])
    valid_overlap = bool(spot_overlap[0])
    n_overlap = 1 if spot_overlap[0] else 0
    i = 1
    while i < n_detectors and valid_overlap:
        if spot_overlap[i]:
            n_overlap += 1
            valid_overlap = True
        elif previous_overlap and not detector_lit[i]:
            valid_overlap = True
            previous_overlap = False
        else:
            valid_overlap = False
        i += 1

    return (
        1 if valid_sim else 0,
        1 if valid_overlap else 0,
        n_overlap if (valid_sim and valid_overlap) else 0,
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
    cpp_rule: bool = False,
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
            detector_lit, spot_overlap, n_detectors, cpp_rule
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


def _omega_bins(omegas: np.ndarray, low: float, width: float, cpp_edge_rules: bool) -> np.ndarray:
    """Uniform omega bin of each angle. Default: truncation toward zero (historical), so an omega in
    [low - width, low) falls in bin 0. ``cpp_edge_rules``: floor, so it is bin -1 (no frame), as C++
    ``CSimulationRange::AngleToIndex``."""
    f = (omegas - low) / width
    return (np.floor(f) if cpp_edge_rules else f).astype(int)


def frame_index_array(range_map: SimulationRange) -> np.ndarray:
    """Frame (wedge) index of each omega bin of ``range_map`` as an ``int64`` array of length
    ``n_bins + 1``: ``-1`` where ``index_list`` has ``None`` (a gap), and a final ``-1`` entry that
    ``frame_of_omegas`` uses for every bin outside the range."""
    return np.array(
        [-1 if i is None else int(i) for i in range_map.index_list] + [-1], dtype=np.int64
    )


def frame_of_omegas(
    omegas: np.ndarray, low: float, width: float, index: np.ndarray, cpp_edge_rules: bool
) -> np.ndarray:
    """Frame index of each omega (``-1``: no frame), one vectorised lookup.

    ``index`` is ``frame_index_array(range_map)``; the bin is ``_omega_bins`` (truncation, or floor
    with ``cpp_edge_rules``). A bin below 0, at or above the number of bins, or on a gap has no
    frame. The result equals the loop it replaces (``0 <= n < n_bins and index_list[n] is not
    None``) for every input. The bin is clipped to ``[-1, n_bins]``: both ends read the trailing
    ``-1`` entry of ``index`` (index ``-1`` is its last element)."""
    bins = _omega_bins(omegas, low, width, cpp_edge_rules)
    return index[bins.clip(-1, len(index) - 1)]


_PLACEHOLDER_FRAME = np.zeros((1, 1), dtype=np.uint8)
_PLACEHOLDER_FRAME.flags.writeable = False


class FrameTable:
    """The binary frames of one ``ExperimentalData`` for the C stage: ``images[frame * n_det +
    det]``, each filled on first use.

    The C extension (``stage_d_overlap``, ``ray_flags``) only indexes this list, so passing the
    whole table with the global frame indices gives the results of a per-call list of just the
    frames in use, without building that list on every call. A frame that has not been used holds
    a 1x1 placeholder that the C code never reads (it reads a frame only for a peak that maps to
    it, and ``ensure`` fills those first).

    **Contract.** The table belongs to the ``exp_data`` it was built for (``source``) and keeps
    each frame's ``get_binary_numpy()`` array from the first fill: the images must not change
    afterwards. Swapping the whole ``exp_data`` object of a cost function is handled
    (``VoxelCostFunction.frame_table`` rebuilds the table). The frame size is read from detector 0
    of the first frame filled and every later detector-0 frame must have the same size.
    ``source`` is not touched until the first fill, so a cost function may be built with
    ``exp_data=None`` and have its data attached later.
    """

    __slots__ = ("source", "n_frames", "n_det", "images", "rows", "cols", "_filled")

    def __init__(self, exp_data: ExperimentalData, n_frames: int, n_det: int) -> None:
        self.source = exp_data
        self.n_frames = n_frames
        self.n_det = n_det
        self.images: List[np.ndarray] = [_PLACEHOLDER_FRAME] * (n_frames * n_det)
        self.rows = 0
        self.cols = 0
        self._filled = np.zeros(n_frames, dtype=bool)

    def ensure(self, frames: np.ndarray) -> None:
        """Fill the frames in ``frames`` (integer frame indices, all in ``[0, n_frames)``)."""
        filled = self._filled[frames]
        if filled.all():
            return
        exp = self.source
        n_det = self.n_det
        for w in np.unique(frames[~filled]):
            w = int(w)
            for d in range(n_det):
                im = exp.get_image(w, d)
                self.images[w * n_det + d] = im.get_binary_numpy()
                if d == 0:
                    size = (im.num_rows, im.num_cols)
                    if self.rows == 0:
                        self.rows, self.cols = size
                    elif size != (self.rows, self.cols):
                        raise ValueError(
                            f"frame {w} of detector 0 has size {size}, the table's frames have "
                            f"{(self.rows, self.cols)}"
                        )
            self._filled[w] = True


def _pixel_coordinate(
    d: torch.Tensor, half: float, size: float, cpp_edge_rules: bool
) -> torch.Tensor:
    """Continuous pixel coordinate of a detector coordinate ``d`` (mm). ``cpp_edge_rules``: a
    negative ``d`` is pixel -1, as C++ ``ToColPixel`` / ``ToRowPixel``."""
    pix = (d + half) / size
    if cpp_edge_rules:
        pix = torch.where(d < 0, torch.full_like(pix, -1.0), pix)
    return pix


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
    cpp_rule: bool = False,
    ray_rule: "Optional[RayRule]" = None,
    voxel_size: float = 0.0,
    cpp_edge_rules: bool = False,
    frame_index: Optional[np.ndarray] = None,
    frame_table: "Optional[FrameTable]" = None,
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

    ``cpp_edge_rules`` (opt-in, default False): three C++ conventions at the edges.
    (1) Frame lookup as C++ ``CSimulationRange`` (``AngleToIndex``: an omega below the range's low
    edge has no frame). The default keeps the historical truncation toward zero, which maps an
    omega in [low - width, low) to the first bin (about 0.3% of the peaks of a 180-frame,
    1 degree scan are compared with frame 0 instead of being dropped). C++ Reference:
    SimulationData.h:317-322 and 446-453.
    (2) Detector coordinate below zero is pixel -1 (C++ ``ToRowPixel``/``ToColPixel`` test the
    millimetre coordinate ``d < 0`` before adding the half pixel), where the default floors
    ``(d + half) / width`` and so gives pixel 0 for a coordinate in [0, 0.5). C++ Reference:
    Detector.cpp:528-548.
    (3) Detector-lit rule as C++ ``ShapePixelOverlapCounter`` (OverlapInfo.h:501-514), exactly: a
    detector is lit iff the pixel overlap is > 0 or one of the three vertices is in bounds. There is
    no "pixels produced" clause, so a triangle that covers pixels with all three vertices out of
    bounds stays unlit, and a triangle that produced no pixel (the fill stops one row early at the
    last image row or column) with a vertex in bounds is lit. The default lights on any produced
    pixel (unchanged). Pixel counts (``peak_pixel_on_det``) are the rasterizer's and need not match
    C++.
    Needs the extension built with ``HAS_VERTEX_LIT`` >= 2 (exact rule), else the slow per-peak
    path is used.

    ``ray_rule`` (opt-in, default None): apply the ray-consistency rule of ``ray_overlap`` on top
    of the C++ CountQualifiedPeaks rule (it implies ``cpp_rule``). It is meant for the coarse
    search levels only; the final cost stays the strict C++ rule. See ray_overlap.py.

    ``frame_index`` / ``frame_table`` (optional, S1): ``frame_index_array(range_map)`` and a
    ``FrameTable`` of ``exp_data`` that the caller keeps across calls (``VoxelCostFunction`` does).
    With neither, the index array is built from ``range_map`` and the fast path fills a table
    private to this call, so the results are the same; with them the per-call work of mapping
    omegas to frames and of collecting the binary frames for the C stage is skipped. A
    ``frame_table`` built for another ``exp_data`` or detector count is ignored.

    C++ Reference: OverlapInfo.h PixelBasedPeakOverlapCounter (pixel_radius > 0)
                   OverlapInfo.h ShapePixelOverlapCounter (pixel_radius == 0)
    """
    n_detectors = len(detector_list)
    overlap_info = OverlapInfo()
    n_peaks = len(peak_omegas)
    if ray_rule is not None:
        cpp_rule = True

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

    # Map omega to wedge index (one vectorised lookup; -1: no frame)
    if frame_index is None:
        frame_index = frame_index_array(range_map)
    wedge_indices = frame_of_omegas(
        omegas.numpy(), range_map.low, range_map.width, frame_index, cpp_edge_rules
    )

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
        col = _pixel_coordinate(j, detector.pixel_half_width, detector.pixel_width, cpp_edge_rules)
        row = _pixel_coordinate(
            k, detector.pixel_half_height, detector.pixel_height, cpp_edge_rules
        )

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

    obs_in = obs_bright = None
    if ray_rule is not None:
        if _HAS_C_RASTERIZE and not _C_HAS_RAY_RULE:
            raise RuntimeError(
                "ray_rule needs the _rasterize extension built with HAS_RAY_RULE. Rebuild: cd "
                "icenine_py && uv run --no-project --with setuptools --with numpy python setup.py "
                "build_ext --inplace"
            )
        pix0 = per_det_pixels[0].numpy().astype(np.float64)  # (M, 3, 2) [col, row]
        ray_origins = lab_verts.mean(dim=1).double().numpy()
        ray_omega_step = abs(float(range_map.width))
        if not (_HAS_C_RASTERIZE and mode == 'hard' and USE_C_RAY_FLAGS):
            obs_in, obs_bright = ray_flags(
                ray_rule,
                detector_list,
                exp_data,
                np.asarray(v_wedge),
                ray_origins,
                per_det_all_hit[0].numpy(),
                pix0,
                pixel_radius,
                voxel_size,
                ray_omega_step,
            )

    global _WARNED_STALE_EXT
    if cpp_rule and _HAS_C_RASTERIZE and not _C_HAS_CPP_RULE and not _WARNED_STALE_EXT:
        _WARNED_STALE_EXT = True
        warnings.warn(
            "CppOverlapRule is on but the _rasterize extension was built before cpp_rule; "
            "using the much slower per-peak Python loop. Rebuild: cd icenine_py && uv run "
            "--no-project --with setuptools --with numpy python setup.py build_ext --inplace",
            RuntimeWarning,
            stacklevel=2,
        )

    if cpp_edge_rules and _HAS_C_RASTERIZE and not _C_HAS_VERTEX_LIT and not _WARNED_STALE_EXT:
        _WARNED_STALE_EXT = True
        warnings.warn(
            "CppEdgeRules is on but the _rasterize extension was built before the exact "
            "cpp_edge_rules; using the much slower per-peak Python loop. Rebuild: cd icenine_py "
            "&& uv run --no-project --with setuptools --with numpy python setup.py "
            "build_ext --inplace",
            RuntimeWarning,
            stacklevel=2,
        )

    # bit field for the C extension: 1 = C++ CountQualifiedPeaks rule, 2 = vertex-lit rule
    rule_bits = (1 if cpp_rule else 0) | (2 if cpp_edge_rules else 0)

    if (
        _HAS_C_RASTERIZE
        and mode == 'hard'
        and (_C_HAS_CPP_RULE or not cpp_rule)
        and (_C_HAS_VERTEX_LIT or not cpp_edge_rules)
    ):
        # ---- Fast path: batch C extension ----
        # Binary frames [frame * n_det + det], filled on first use and kept across calls by the
        # caller's FrameTable; the C code indexes it with the global frame indices
        v_wedge_np = np.asarray(v_wedge, dtype=np.int32)
        table = frame_table
        if table is None or table.source is not exp_data or table.n_det != n_detectors:
            table = FrameTable(exp_data, int(frame_index.max()) + 1, n_detectors)
        table.ensure(v_wedge_np)
        image_list = table.images
        local_wedge = v_wedge_np  # (M,) frame indices into the table

        if ray_rule is not None and obs_in is None:
            # observed-ray check for all peaks of the evaluation in one C call (no Python loop)
            obs_in, obs_bright = ray_flags_c(
                ray_rule,
                detector_list,
                image_list,
                local_wedge,
                ray_origins,
                per_det_all_hit[0].numpy(),
                pix0,
                pixel_radius,
                voxel_size,
                ray_omega_step,
            )

        # Build det_hit_mask: (M, n_det) uint8
        det_hit = np.stack(
            [per_det_all_hit[d].numpy().astype(np.uint8) for d in range(n_detectors)],
            axis=1,
        )  # (M, n_det)

        # Image dimensions: detector 0's frame size (the table checks every frame against it)
        img_rows = table.rows
        img_cols = table.cols

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
                *((rule_bits,) if rule_bits else ()),
                *((obs_in, obs_bright) if ray_rule is not None else ()),
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
                *((rule_bits,) if rule_bits else ()),
                *((obs_in, obs_bright) if ray_rule is not None else ()),
            )

        # Unpack results
        overlap_info.pixel_overlap = result[0]
        overlap_info.pixel_on_detector = result[1]
        overlap_info.peak_overlap = result[2]
        overlap_info.peak_on_detector = result[3]
        overlap_info.quality = result[4]
        overlap_info.n_quality_points = result[5]
        overlap_info.detectors_overlap = result[6]
        if ray_rule is not None:
            overlap_info.ray_cases = tuple(int(x) for x in result[7])

    else:
        # ---- Fallback: sequential Python loop ----
        case_counts = np.zeros(8, dtype=np.int64) if ray_rule is not None else None
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
                    elif cpp_edge_rules:
                        # CppEdgeRules, exact C++ ShapePixelOverlapCounter (OverlapInfo.h:501-514):
                        # lit iff overlap > 0 (above) or a vertex is in bounds; no "pixels produced"
                        # clause, so a triangle covering pixels with all vertices out of bounds is
                        # unlit
                        for vtx in (v0, v1, v2):
                            vx = -1.0 if float(vtx[0]) < 0 else math.floor(float(vtx[0]))
                            vy = -1.0 if float(vtx[1]) < 0 else math.floor(float(vtx[1]))
                            if 0 <= vx < exp_image.num_cols and 0 <= vy < exp_image.num_rows:
                                detector_lit[det_idx] = True
                                break
                    elif n_lit_int > 0:
                        detector_lit[det_idx] = True

            if ray_rule is not None:
                d_ov, d_on = apply_cases_python(
                    detector_lit, spot_overlap, obs_in[p], obs_bright[p], case_counts
                )
                peak_pixel_overlap += d_ov
                peak_pixel_on_det += d_on

            peak_on_det, peak_ovlp, n_det_ovlp = count_qualified_peaks(
                detector_lit, spot_overlap, n_detectors, cpp_rule
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
        if case_counts is not None:
            overlap_info.ray_cases = tuple(int(x) for x in case_counts)

    return overlap_info


# ---------------------------------------------------------------------------
# Fused evaluate (S2, opt-in key FusedCost): one C call for VoxelCostFunction.evaluate
# ---------------------------------------------------------------------------

#: How the fused path of ``VoxelCostFunction.evaluate`` (key FusedCost) treats each argument of
#: ``VoxelCostFunction.__init__``: "context" (read into the C context, rebuilt when it changes),
#: "covered" (every value is implemented in C), "fallback" (a value the C code does not implement
#: falls back to the torch path, see ``fused_unsupported_reason``), "key" (the switch itself).
#: tests/test_fused_cost.py fails when ``__init__`` gains an argument that is not listed here, so a
#: new cost mode cannot reach the fused path untested: classify it, and make
#: ``fused_unsupported_reason`` return a reason for its values that C does not implement.
FUSED_ARGUMENTS: Dict[str, str] = {
    "simulator": "context",  # beam direction, energy, chi (one chi for every reflection)
    "detector_list": "fallback",  # geometry; more than 8 detectors fall back
    "range_map": "context",  # the frame lookup
    "exp_data": "context",  # the FrameTable
    "sample": "context",  # sample-to-lab rotation and translation
    "structure_list": "context",  # the reciprocal vectors of each phase
    "mode": "fallback",  # 'hard' only; soft mode falls back
    "eta_limit": "fallback",  # None falls back
    "pixel_radius": "covered",  # 0 (triangle) and > 0 (pixel radius)
    "max_q": "covered",  # filters the reflections before the context is built
    "min_sin_eta": "covered",
    "cpp_overlap_rule": "covered",
    "center_screen": "covered",  # centroid or vertex 0 of the pixel-radius screen
    "ray_rule": "fallback",  # RayOverlapRule falls back
    "cpp_edge_rules": "covered",
    "fused_cost": "key",
}

_FUSED_MAX_DETECTORS = 8  # MAX_DETECTORS of _stage_d.h
_FUSED_BUILD_HINT = (
    "build it: cd icenine_py && uv run --no-project --with setuptools --with numpy python setup.py "
    "build_ext --inplace"
)
_FUSED_VERTICES_REASON = "voxel vertices that are not a (3, 3) float32 tensor on the CPU"
_WARNED_FUSED: set = set()  # reasons already warned about (one RuntimeWarning per reason)


def fused_unsupported_reason(cost_fn: "VoxelCostFunction") -> Optional[str]:
    """Why the fused evaluate (key FusedCost) cannot compute ``cost_fn``'s costs (None: it can).
    Checks every "fallback" argument of ``FUSED_ARGUMENTS`` and the build of the extension; the
    dtypes and the beam deflection are checked when the context is built
    (``_FusedContext.reason``)."""
    if not _HAS_C_RASTERIZE:
        return "the _rasterize C extension is not built; " + _FUSED_BUILD_HINT
    if not _C_HAS_FUSED:
        return (
            "the _rasterize extension was built without the fused evaluate (HAS_FUSED); rebuild: "
            + _FUSED_BUILD_HINT
        )
    if cost_fn.mode != "hard":
        return "soft mode"
    if cost_fn.ray_rule is not None:
        return "the ray rule (RayOverlapRule)"
    if cost_fn.eta_limit is None:
        return "no eta limit"
    n_det = len(cost_fn.detector_list)
    if n_det < 1 or n_det > _FUSED_MAX_DETECTORS:
        return f"{n_det} detectors (1 to {_FUSED_MAX_DETECTORS} supported)"
    return None


def _warn_fused_fallback(reason: str, stacklevel: int = 5) -> None:
    """One ``RuntimeWarning`` per reason; ``stacklevel`` points at the caller of
    ``VoxelCostFunction.evaluate`` (5 from ``_fused_context``, 4 from ``_evaluate_fused``)."""
    if reason not in _WARNED_FUSED:
        _WARNED_FUSED.add(reason)
        warnings.warn(
            f"FusedCost is not supported for this cost function ({reason}): using the torch "
            "evaluate",
            RuntimeWarning,
            stacklevel=stacklevel,
        )


def _numpy_frame_mode(low: object, width: object) -> int:
    """How numpy evaluates ``(omegas - low) / width`` (``_omega_bins``) for float32 omegas, as the
    C frame lookup must repeat it: 0 both steps in float64, 1 the subtraction in float32 and the
    division in float64, 2 both in float32 (NEP 50: a numpy float64 scalar promotes, a Python float
    does not)."""
    a = np.zeros(1, np.float32) - low
    if a.dtype == np.float64:
        return 0
    return 1 if (a / width).dtype == np.float64 else 2


#: g = R g_hkl in the fused C code: 1 = the plain order ``(r0 g0 + r1 g1) + r2 g2`` (the default),
#: 0 = the fmaf chain of the S2 Step 0 prototype. Neither reproduces torch's sgemm bit for bit; on
#: S1's parity set the plain order gave fewer pixel-edge flips against torch (S2 Step 1,
#: D-20261027-I5; scripts/speed_parity/s2/s2_core_parity.py, benchmarks/s2_core/). Read when a
#: context is built; 0 is kept for that comparison and the bit-for-bit check against the prototype.
_FUSED_GLAB_ORDER = 1


class _FusedContext:
    """The C context of the fused evaluate for one (cost function, phase), the arrays it points
    into, and the inputs it was built from. ``is_current`` is False after any of them changed: a
    new FrameTable (``exp_data``, ``range_map`` or the detector count replaced; the checks of
    ``VoxelCostFunction.frame_table``), new reciprocal vectors, beam energy or chi (the checks of
    ``BraggConstants.is_for``), a cost mode changed on the live object, or the geometry: the
    sample-to-lab matrix or the beam direction replaced or modified in place (tensor version), a
    detector's plane, position, frame origin or basis tensors replaced (``Detector.set_position``
    and ``set_orientation`` replace them) or its pixel size changed.

    **Contract (D-20261027-I8).** The detector tensors are compared by identity only: an in-place
    edit of a detector's plane, position, frame origin or basis tensors is not detected (checking
    their versions would add cost to every call). Move a detector with ``Detector.set_position`` or
    ``set_orientation``, which replace the tensors; no library code edits them in place."""

    __slots__ = (
        "capsule",
        "reason",
        "table",
        "g_hkl_t",
        "g_mag_t",
        "energy",
        "chi",
        "modes",
        "s2l",
        "s2l_version",
        "beam",
        "beam_version",
        "detectors",
        "det_state",
        "_keep",
    )

    def __init__(self, cf: "VoxelCostFunction", phase: int) -> None:
        self.capsule: object = None
        self.reason: Optional[str] = None
        g_hkl_t, g_mag_t = cf._phase_recip_vecs[phase]
        sim = cf.simulator
        energy, chi = sim.beam_energy, sim.beam_deflection_chi
        self.table = cf.frame_table()
        self.g_hkl_t, self.g_mag_t = g_hkl_t, g_mag_t
        self.energy, self.chi = energy, chi
        self.modes = cf._fused_modes()
        s2l = cf.sample.sample_to_lab_matrix
        beam = sim.beam_direction
        self.s2l, self.s2l_version = s2l, s2l._version
        self.beam, self.beam_version = beam, beam._version
        self.detectors = cf.detector_list
        self.det_state = [_detector_state(det) for det in cf.detector_list]
        self._keep: Tuple[np.ndarray, ...] = ()

        # the torch path's Bragg constants (the same object, so the same values)
        consts = cf._bragg_consts.get(phase)
        if consts is None or not consts.is_for(g_mag_t, energy, chi):
            consts = cf._bragg_consts[phase] = BraggConstants(g_mag_t, energy, chi)
        sin_chi = consts.sin_chi_laue.reshape(-1)
        cos_chi = consts.cos_chi_laue.reshape(-1)
        tensors = [g_hkl_t, consts.sin_theta, consts.g_mag_eps, s2l, beam]
        for det in cf.detector_list:
            pl = det._detector_plane
            tensors += [
                pl.normal,
                det._position,
                det._lab_frame_coord_origin,
                det._lab_frame_basis_j,
                det._lab_frame_basis_k,
            ]
        if any(t.dtype != torch.float32 or t.device.type != "cpu" for t in tensors):
            self.reason = "a geometry or reflection tensor that is not float32 on the CPU"
            return
        if not bool((sin_chi == sin_chi[0]).all()) or not bool((cos_chi == cos_chi[0]).all()):
            self.reason = "a beam deflection chi that differs between reflections"
            return

        g_hkl = np.ascontiguousarray(g_hkl_t.numpy(), dtype=np.float32)
        sin_theta = np.ascontiguousarray(consts.sin_theta.numpy(), dtype=np.float32)
        g_mag_eps = np.ascontiguousarray(consts.g_mag_eps.numpy(), dtype=np.float32)
        scal = np.array(
            [
                float(sin_chi[0]),
                float(cos_chi[0]),
                np.float32(consts.epsilon),
                np.float32(cf.eta_limit),
                np.float32(cf.min_sin_eta),
            ],
            np.float32,
        )
        rule_bits = (1 if cf.cpp_overlap_rule else 0) | (2 if cf.cpp_edge_rules else 0)
        rmap = cf.range_map
        ints = np.array(
            [
                cf.pixel_radius,
                int(bool(cf.center_screen)),
                rule_bits,
                int(cf.min_sin_eta > 0.0),
                _numpy_frame_mode(rmap.low, rmap.width),
                len(cf.detector_list),
                _FUSED_GLAB_ORDER,
            ],
            np.int32,
        )
        geom = np.stack(
            [
                np.concatenate(
                    [
                        det._detector_plane.normal.numpy(),
                        np.array([float(det._detector_plane.d)], np.float32),
                        det._position.numpy(),
                        det._lab_frame_coord_origin.numpy(),
                        det._lab_frame_basis_j.numpy(),
                        det._lab_frame_basis_k.numpy(),
                        np.array(
                            [
                                det.pixel_half_width,
                                det.pixel_width,
                                det.pixel_half_height,
                                det.pixel_height,
                            ],
                            np.float32,
                        ),
                    ]
                ).astype(np.float32)
                for det in cf.detector_list
            ]
        )
        frame_index = cf._frame_index
        assert frame_index is not None  # set by frame_table()
        self._keep = (g_hkl, sin_theta, g_mag_eps, frame_index)
        self.capsule = _c_fused_context(
            g_hkl,
            sin_theta,
            g_mag_eps,
            scal,
            beam.numpy(),
            np.ascontiguousarray(s2l[:3, :3].numpy()),
            np.ascontiguousarray(s2l[:3, 3].numpy()),
            ints,
            np.array([float(rmap.low), float(rmap.width)], np.float64),
            np.ascontiguousarray(frame_index, dtype=np.int64),
            np.ascontiguousarray(geom, dtype=np.float32),
            self.table,
            self.table.images,
            self.table._filled,
        )

    def is_current(self, cf: "VoxelCostFunction", phase: int) -> bool:
        sim = cf.simulator
        s2l = cf.sample.sample_to_lab_matrix
        beam = sim.beam_direction
        recip = cf._phase_recip_vecs.get(phase)
        if not (
            self.table is cf.frame_table()
            and recip is not None
            and recip[0] is self.g_hkl_t
            and recip[1] is self.g_mag_t
            and self.energy == sim.beam_energy
            and self.chi == sim.beam_deflection_chi
            and s2l is self.s2l
            and s2l._version == self.s2l_version
            and beam is self.beam
            and beam._version == self.beam_version
            and cf.detector_list is self.detectors
            and self.modes == cf._fused_modes()
        ):
            return False
        if len(cf.detector_list) != len(self.det_state):
            return False
        for det, st in zip(cf.detector_list, self.det_state):
            pl = det._detector_plane
            if not (
                pl is st[0]
                and pl.coeffs is st[1]
                and det._position is st[2]
                and det._lab_frame_coord_origin is st[3]
                and det._lab_frame_basis_j is st[4]
                and det._lab_frame_basis_k is st[5]
                and det.pixel_half_width == st[6]
                and det.pixel_width == st[7]
                and det.pixel_half_height == st[8]
                and det.pixel_height == st[9]
            ):
                return False
        return True


def _detector_state(det: Detector) -> Tuple[object, ...]:
    """What the fused context reads from a detector (``_FusedContext.is_current`` compares the
    tensors by identity, the pixel sizes by value)."""
    pl = det._detector_plane
    return (
        pl,
        pl.coeffs,
        det._position,
        det._lab_frame_coord_origin,
        det._lab_frame_basis_j,
        det._lab_frame_basis_k,
        det.pixel_half_width,
        det.pixel_width,
        det.pixel_half_height,
        det.pixel_height,
    )


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
        min_sin_eta: float = 0.0,
        cpp_overlap_rule: bool = False,
        center_screen: bool = True,
        ray_rule: "Optional[RayRule]" = None,
        cpp_edge_rules: bool = False,
        fused_cost: bool = False,
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
            min_sin_eta: Optional lower bound on |sin eta| of the diffracted beam
                   (eta as in the eta filter above): peaks nearer the rotation axis
                   are skipped. 0 (default) = no filter (the C++ behaviour). Used to
                   give the optimizers the same eligible spots as the orientation
                   network (generate_toy_orientation_dataset.build_problem).
            cpp_overlap_rule: Apply the literal C++ CountQualifiedPeaks overlap rule (a lit,
                   non-overlapping later detector invalidates a peak's overlap). Default False
                   (the original Python rule). C++ Reference: OverlapInfo.tmpl.cpp:78-139.
            center_screen: For pixel_radius > 0 project the voxel centroid (C++ VoxelToCenter)
                   instead of the first triangle vertex (vertex 0). Default True for every
                   rule (identical start point to C++); False is an explicit opt-out kept for
                   history.
            ray_rule: Opt-in ray-consistency rule (ray_overlap.RayRule), default None = off. It
                   implies the C++ rule. The reconstructor gives it only to the coarse levels
                   0-2 (screen, re-cost, quick MC); FindOptimal, VarianceMinimizing and the final
                   cost use a cost function without it (two costs by design, see ray_overlap.py).
            cpp_edge_rules: Opt-in, default False: the C++ conventions at the edges (an omega
                   below the range has no frame; a detector coordinate below 0 is pixel -1; a
                   pixel-less triangle with an in-bounds vertex lights the detector; see
                   calculate_diffraction_overlap_batched).
            fused_cost: Opt-in, default False (config key FusedCost, S2): ``evaluate`` computes
                   the cost in one C call (``_rasterize.fused_evaluate``, icenine/_fused_eval.c)
                   instead of the torch operations. Float32 in torch's operation order, but not
                   bit for bit: a pixel coordinate within float32 rounding of an integer can floor
                   to the other pixel (parity rule D-20261027-I2). Covered: hard mode, pixel radius
                   0 and > 0 (centroid and vertex 0), the default rule, ``cpp_overlap_rule``,
                   ``cpp_edge_rules``, ``min_sin_eta``. Everything else (soft mode, ``ray_rule``,
                   no eta limit, more than 8 detectors, an extension built without
                   ``HAS_FUSED``) uses the torch path, with one ``RuntimeWarning`` per reason
                   (``fused_unsupported_reason``, ``FUSED_ARGUMENTS``).
        """
        self.fused_cost = fused_cost
        self._fused_ctx: Dict[int, _FusedContext] = {}
        self._fused_vt: Optional[torch.Tensor] = None  # vertices of the cached array, its version
        self._fused_vv = -1
        self._fused_vn: Optional[np.ndarray] = None
        self.cpp_edge_rules = cpp_edge_rules
        self.ray_rule = ray_rule
        self.cpp_overlap_rule = cpp_overlap_rule or ray_rule is not None
        self.center_screen = center_screen
        self.simulator = simulator
        self.min_sin_eta = min_sin_eta
        self.detector_list = detector_list
        self.range_map = range_map
        self.exp_data = exp_data
        self.sample = sample
        self.structure_list = structure_list
        self.mode = mode
        self.eta_limit = eta_limit
        self.pixel_radius = pixel_radius
        self.eval_count = 0
        # S1: the frame lookup array and the binary-frame table are built on first use
        # (``frame_table``); the Bragg constants per phase on its first ``evaluate``
        self._frame_index: Optional[np.ndarray] = None
        self._frame_range: Optional[SimulationRange] = None
        self._frame_table: Optional[FrameTable] = None
        self._bragg_consts: Dict[int, BraggConstants] = {}

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

    def frame_table(self) -> FrameTable:
        """The binary-frame table of ``self.exp_data`` (see ``FrameTable``), built on first use and
        rebuilt, together with the frame lookup array, when ``exp_data`` (or ``range_map``, or the
        detector count) has been replaced since: swapping the data of a live cost function keeps
        working. Changing the images of the same ``exp_data`` object after the first ``evaluate`` is
        not supported (``FrameTable`` contract)."""
        t = self._frame_table
        if (
            t is not None
            and t.source is self.exp_data
            and self._frame_range is self.range_map
            and t.n_det == len(self.detector_list)
        ):
            return t
        if self._frame_index is None or self._frame_range is not self.range_map:
            self._frame_index = frame_index_array(self.range_map)
            self._frame_range = self.range_map
        t = self._frame_table = FrameTable(
            self.exp_data, int(self._frame_index.max()) + 1, len(self.detector_list)
        )
        return t

    def _fused_modes(self) -> Tuple[object, ...]:
        """The cost modes the fused C context is built for (a change rebuilds it)."""
        return (
            self.mode,
            self.ray_rule,
            self.pixel_radius,
            self.center_screen,
            self.cpp_overlap_rule,
            self.cpp_edge_rules,
            self.eta_limit,
            self.min_sin_eta,
        )

    def _fused_context(self, phase_index: int) -> Optional[_FusedContext]:
        """The current fused context of a phase (built or rebuilt as needed), or None when the
        fused path does not cover this cost function (one RuntimeWarning per reason)."""
        ctx = self._fused_ctx.get(phase_index)
        if ctx is not None and ctx.is_current(self, phase_index):
            return ctx if ctx.reason is None else None
        reason = fused_unsupported_reason(self)
        if reason is None:
            ctx = self._fused_ctx[phase_index] = _FusedContext(self, phase_index)
            reason = ctx.reason
        if reason is not None:
            _warn_fused_fallback(reason)
            return None
        return ctx

    def _fused_vertices(self, voxel_vertices: torch.Tensor) -> Optional[np.ndarray]:
        """``voxel_vertices`` as the (3, 3) float32 array of the C call (cached while the same
        tensor is passed unmodified); None for anything the torch path would treat otherwise."""
        vt = voxel_vertices
        if vt is self._fused_vt and vt._version == self._fused_vv:
            return self._fused_vn
        if (
            not isinstance(vt, torch.Tensor)
            or vt.dtype != torch.float32
            or vt.device.type != "cpu"
            or tuple(vt.shape) != (3, 3)
        ):
            return None
        self._fused_vn = np.ascontiguousarray(vt.detach().numpy())
        self._fused_vt, self._fused_vv = vt, vt._version
        return self._fused_vn

    def _evaluate_fused(
        self, orientation: np.ndarray, voxel_vertices: torch.Tensor, phase_index: int
    ) -> Optional[OverlapInfo]:
        """``evaluate`` through the fused C call (key FusedCost); None: use the torch path."""
        if phase_index not in self._phase_recip_vecs:
            return OverlapInfo()
        ctx = self._fused_context(phase_index)
        if ctx is None:
            return None
        verts = self._fused_vertices(voxel_vertices)
        if verts is None:
            _warn_fused_fallback(_FUSED_VERTICES_REASON, stacklevel=4)
            return None
        r = _c_fused_evaluate(ctx.capsule, orientation, verts)
        if type(r) is list:
            # lazy FrameTable fill: C names the frames it needs, Python fills them, C runs again
            ctx.table.ensure(np.asarray(r, dtype=np.int64))
            r = _c_fused_evaluate(ctx.capsule, orientation, verts)
            if type(r) is list:
                raise RuntimeError(f"fused evaluate: frames {r} still missing after ensure")
        return OverlapInfo(r[0], r[1], r[2], r[3], r[6], r[4], r[5])

    def fused_debug(
        self, orientation: np.ndarray, voxel_vertices: torch.Tensor, phase_index: int = 0
    ) -> Dict[str, np.ndarray]:
        """The intermediates of the fused evaluate for one orientation (parity checks): per
        reflection ``g_lab``, ``obs``, ``omega1``, ``omega2``; per (solution, reflection) row
        ``r = s * K + k`` ``omega``, ``cos``, ``sin``, ``eta``, ``valid``, ``frame``,
        ``lab_normal``, ``reflected``, ``pix`` (n_det, 3, 2), ``hit``, ``t``; and ``result``, the
        7 values ``evaluate`` returns (frames are filled first). Does not count an evaluation.
        Raises RuntimeError where the fused path does not cover the cost function."""
        if phase_index not in self._phase_recip_vecs:
            raise RuntimeError(f"no reflections for phase {phase_index}")
        ctx = self._fused_context(phase_index)
        verts = self._fused_vertices(voxel_vertices)
        if ctx is None or verts is None:
            raise RuntimeError("the fused evaluate does not cover this cost function or input")
        d: Dict[str, np.ndarray] = _c_fused_evaluate(ctx.capsule, orientation, verts, 1)
        if "missing" in d:
            ctx.table.ensure(np.asarray(d["missing"], dtype=np.int64))
            d = _c_fused_evaluate(ctx.capsule, orientation, verts, 1)
        return d

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
        self.eval_count += 1
        if self.fused_cost:
            fused = self._evaluate_fused(orientation, voxel_vertices, phase_index)
            if fused is not None:
                return fused
        voxel_size = 0.0
        if self.ray_rule is not None:
            edges = voxel_vertices - voxel_vertices.roll(1, dims=0)
            voxel_size = float(torch.linalg.norm(edges, dim=1).max())
        if self.center_screen and self.pixel_radius > 0:
            voxel_vertices = voxel_vertices.mean(dim=0, keepdim=True).repeat(3, 1)

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

        # Compute Bragg angles for all reflections: get_scattering_omegas_torch with the
        # call-independent constants hoisted (BraggConstants; the omegas of unobservable rows are
        # not zeroed, they are dropped below). Rebuilt when |g|, the beam energy or chi changed.
        energy = self.simulator.beam_energy
        chi = self.simulator.beam_deflection_chi
        consts = self._bragg_consts.get(phase_index)
        if consts is None or not consts.is_for(g_mag_batch, energy, chi):
            consts = self._bragg_consts[phase_index] = BraggConstants(g_mag_batch, energy, chi)
        bragg_result = scattering_omegas_fixed_beam(g_lab_batch, consts)

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
        if self.min_sin_eta > 0.0:
            valid = valid & (torch.sin(eta) >= self.min_sin_eta)

        peak_omegas_t = all_omegas[valid]  # (M,) tensor
        peak_normals_t = all_normals[valid]  # (M, 3) tensor

        if peak_omegas_t.shape[0] == 0:
            return OverlapInfo()

        # Calculate diffraction overlap (batched for performance)
        table = self.frame_table()
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
            frame_index=self._frame_index,
            frame_table=table,
            **({"cpp_rule": True} if self.cpp_overlap_rule else {}),
            **({"cpp_edge_rules": True} if self.cpp_edge_rules else {}),
            **(
                {"ray_rule": self.ray_rule, "voxel_size": voxel_size}
                if self.ray_rule is not None
                else {}
            ),
        )
