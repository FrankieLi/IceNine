"""
Core simulation engine for X-ray diffraction forward model.

Provides physics calculations for projecting diffraction peaks onto detectors,
including Bragg condition evaluation, ray tracing, and peak generation.

Python port of Src/Simulation.h/cpp

Author: S. F. Li
"""

from dataclasses import dataclass
from typing import List, Optional, Callable, Tuple
import numpy as np
import torch

from .experiment_setup import ExperimentSetup
from .detector import Detector
from .sample import Sample
from .image_data import ImageData
from .crystal_structure import CrystalStructure
from .geometry import Ray
from .diffraction_core import (
    get_scattering_omegas_torch,
    get_reflected_ray,
    get_illuminated_pixel
)
from .constants import KEV_OVER_HBAR_C_IN_ANG


# ============================================================================
# Data Structures
# ============================================================================

@dataclass
class PeakInfo:
    """
    Information about a single observable diffraction peak.

    C++ Reference:
        Simulation.h:~30-50 SPeakInfo struct

    Attributes:
        omega: Sample rotation angle (radians) where peak is observable
        g_vector: Scattering vector in sample frame, shape (3,)
        g_magnitude: Magnitude of scattering vector (Å⁻¹)
        intensity: Peak intensity (arbitrary units)
        detector_index: Which detector observes this peak
    """
    omega: float
    g_vector: torch.Tensor
    g_magnitude: float
    intensity: float
    detector_index: int = 0


# Type alias for peak filter callables
# Returns: (accept: bool, intensity: float)
PeakFilterFn = Callable[[torch.Tensor], Tuple[bool, float]]


# ============================================================================
# Core Simulation Class
# ============================================================================

class Simulation:
    """
    Core diffraction simulation engine.

    Handles fundamental physics calculations for X-ray diffraction including:
    - Solving Bragg condition for observable peaks
    - Ray tracing from sample vertices to detector pixels
    - Projecting voxels onto detector images

    C++ Reference:
        Simulation.h:59-317 template <class DetectorListT> class CSimulation

    Attributes:
        beam_direction: X-ray beam direction vector in lab frame (normalized)
        beam_energy: Beam energy in keV
        beam_deflection_chi: Beam deflection angle (radians)
        initialized: Whether simulation has been initialized

    Example:
        >>> from icenine.experiment_setup import ExperimentSetup
        >>> from icenine.config_file import ConfigFile
        >>>
        >>> config = ConfigFile.from_file("test.config")
        >>> exp_setup = ExperimentSetup(config)
        >>> simulator = Simulation(exp_setup)
    """

    def __init__(self, experiment_setup: Optional[ExperimentSetup] = None):
        """
        Initialize simulation from experiment setup.

        Args:
            experiment_setup: Experiment parameters (optional)

        C++ Reference:
            Simulation.cpp:49-61 CSimulation constructor
        """
        if experiment_setup is not None:
            self.beam_direction = torch.from_numpy(
                experiment_setup.get_beam_direction()
            ).float()
            self.beam_energy = experiment_setup.get_beam_energy()
            self.beam_deflection_chi = experiment_setup.get_beam_deflection_chi_laue()
            self.initialized = True
        else:
            self.beam_direction = torch.tensor([0.0, 0.0, 1.0])
            self.beam_energy = 0.0
            self.beam_deflection_chi = 0.0
            self.initialized = False

    def get_observable_peaks(
        self,
        voxel_orientation: torch.Tensor,
        reciprocal_vectors: List[torch.Tensor],
        detector_list: Optional[List[Detector]] = None
    ) -> List[PeakInfo]:
        """
        Generate list of observable peaks for a voxel orientation.

        For each reciprocal lattice vector (reflection), calculates the omega
        angles where the Bragg condition is satisfied and creates PeakInfo
        structures for both solutions.

        Args:
            voxel_orientation: Voxel orientation matrix, shape (3, 3)
            reciprocal_vectors: List of reciprocal lattice vectors G_hkl, each shape (3,)
            detector_list: List of detectors (optional, for filtering)

        Returns:
            List of PeakInfo structures for all observable peaks

        Algorithm:
            1. For each reciprocal vector G:
                - Transform to lab frame: G' = orientation @ G
                - Calculate omega angles via get_scattering_omegas_torch
                - Create PeakInfo for both ω₁ and ω₂ if observable

        C++ Reference:
            Simulation.cpp:130-152 CSimulation::GetObservablePeaks

        Example:
            >>> # For a voxel with identity orientation
            >>> orientation = torch.eye(3)
            >>> reflections = [torch.tensor([2.668, 0., 0.])]  # (111) for Au
            >>> peaks = simulator.get_observable_peaks(orientation, reflections)
        """
        if not self.initialized:
            raise RuntimeError("Simulation not initialized")

        peaks = []

        # Stack all reciprocal vectors for batched processing
        if len(reciprocal_vectors) == 0:
            return peaks

        g_hkl_batch = torch.stack(reciprocal_vectors)  # (M, 3)

        # Transform to lab frame: G' = orientation @ G
        # (3, 3) @ (M, 3).T -> (3, M) -> (M, 3)
        g_lab_batch = (voxel_orientation @ g_hkl_batch.T).T

        # Calculate magnitudes
        g_magnitudes = torch.norm(g_lab_batch, dim=1)

        # Batch calculate omega angles
        result = get_scattering_omegas_torch(
            g_lab_batch,
            g_magnitudes,
            self.beam_energy,
            self.beam_deflection_chi
        )

        # Extract observable peaks
        for i in range(len(reciprocal_vectors)):
            if not result.observable[i]:
                continue

            # Create PeakInfo for first omega solution
            peaks.append(PeakInfo(
                omega=result.omega1[i].item(),
                g_vector=g_lab_batch[i],
                g_magnitude=g_magnitudes[i].item(),
                intensity=1.0,  # Will be updated by peak filter
                detector_index=0
            ))

            # Create PeakInfo for second omega solution
            peaks.append(PeakInfo(
                omega=result.omega2[i].item(),
                g_vector=g_lab_batch[i],
                g_magnitude=g_magnitudes[i].item(),
                intensity=1.0,
                detector_index=0
            ))

        return peaks

    def project_vertex(
        self,
        detector: Detector,
        sample: Sample,
        vertex: torch.Tensor,
        normal: torch.Tensor
    ) -> Tuple[bool, float, float]:
        """
        Project single vertex onto detector via ray tracing.

        Builds reflected ray from vertex with given normal, intersects with
        detector plane, and converts to pixel coordinates.

        Args:
            detector: Detector to project onto
            sample: Sample containing vertex
            vertex: Vertex position in sample frame, shape (3,)
            normal: Surface normal in sample frame, shape (3,)

        Returns:
            Tuple of (hit, pixel_col, pixel_row):
                - hit: True if ray hits detector
                - pixel_col: Column coordinate (J direction)
                - pixel_row: Row coordinate (K direction)

        C++ Reference:
            Simulation.cpp:187-194 CSimulation::ProjectVertex

        Example:
            >>> vertex = torch.tensor([0., 0., 0.])
            >>> normal = torch.tensor([0., 0., 1.])
            >>> hit, col, row = simulator.project_vertex(detector, sample, vertex, normal)
        """
        # Build reflected ray from vertex
        # C++: GetReflectedRay(oSample, oVertex, oNormal, oBeamDirection)
        reflected_ray = get_reflected_ray(
            sample,
            vertex,
            normal,
            self.beam_direction
        )

        # Project ray onto detector
        # C++: return GetIlluminatedPixel(p, oDetector, oReflectedRay)
        hit, pixel_col, pixel_row = get_illuminated_pixel(detector, reflected_ray)

        return hit.item(), pixel_col.item(), pixel_row.item()

    def project_voxel(
        self,
        image: ImageData,
        detector: Detector,
        sample: Sample,
        voxel_vertices: torch.Tensor,
        normal: torch.Tensor,
        peak_filter: PeakFilterFn
    ) -> bool:
        """
        Project voxel onto detector image by rasterizing triangular peak.

        Projects the voxel's 3 vertices onto the detector. If all hit, rasterizes
        the resulting triangle with intensity determined by the peak filter.

        Args:
            image: ImageData to rasterize into
            detector: Detector geometry
            sample: Sample containing voxel
            voxel_vertices: Voxel vertices in sample frame, shape (3, 3)
                           [vertex0, vertex1, vertex2]
            normal: Surface normal in sample frame, shape (3,)
                   (same for all vertices - constant interpolation)
            peak_filter: Callable that takes scattering direction and returns
                        (accept: bool, intensity: float)

        Returns:
            True if voxel was successfully projected

        Algorithm:
            1. Calculate reflected ray direction once (shared by all vertices)
            2. For each of 3 vertices:
                - Build reflected ray
                - Intersect with detector
                - Get pixel coordinates
            3. Apply peak acceptance filter
            4. If all vertices hit and filter accepts:
                - Rasterize triangle with given intensity

        C++ Reference:
            Simulation.tmpl.cpp:62-104 CSimulation::ProjectVoxel (template)

        Notes:
            - Template parameter FPeakFilter in C++ → callable in Python
            - Optimization: reflected direction calculated once, reused for all vertices

        Example:
            >>> from icenine.peak_filters import TrivialAcceptFn
            >>> vertices = torch.tensor([[0., 0., 0.],
            ...                          [0.1, 0., 0.],
            ...                          [0., 0.1, 0.]])
            >>> normal = torch.tensor([0., 0., 1.])
            >>> filter_fn = TrivialAcceptFn(intensity=1.0)
            >>> success = simulator.project_voxel(
            ...     image, detector, sample, vertices, normal, filter_fn
            ... )
        """
        # Calculate reflected ray direction once (optimization)
        # C++: GetReflectedRayDir(oSample, oNormal, oBeamDirection)
        from .diffraction_core import get_reflected_ray_dir, build_reflected_ray

        reflected_dir = get_reflected_ray_dir(sample, normal, self.beam_direction)

        # Apply peak acceptance filter
        # Normalize reflected direction for filter
        reflected_dir_normalized = reflected_dir / torch.norm(reflected_dir)
        accept, intensity = peak_filter(reflected_dir_normalized)

        if not accept:
            return False

        # Project all 3 vertices
        pixels = []
        for i in range(3):
            vertex = voxel_vertices[i]

            # Build reflected ray from this vertex
            # C++: BuildReflectedRay(oSample, oVertex, oRefDir)
            reflected_ray = build_reflected_ray(sample, vertex, reflected_dir)

            # Get illuminated pixel
            # C++: GetIlluminatedPixel(oPixels[i], oDetector, oReflectedRay)
            hit, pixel_col, pixel_row = get_illuminated_pixel(detector, reflected_ray)

            if not hit.item():
                return False  # Voxel partially off detector

            pixels.append((pixel_col.item(), pixel_row.item()))

        # All 3 vertices hit detector - rasterize triangle
        # C++: oOutImage.AddTriangle(oPixels[0], oPixels[1], oPixels[2], fIntensity)
        v0 = torch.tensor(pixels[0])
        v1 = torch.tensor(pixels[1])
        v2 = torch.tensor(pixels[2])

        image.add_triangle(v0, v1, v2, intensity)

        return True

    def __repr__(self):
        """String representation for debugging."""
        return (
            f"Simulation("
            f"beam_energy={self.beam_energy:.2f} keV, "
            f"beam_dir={self.beam_direction.numpy()}, "
            f"initialized={self.initialized})"
        )
