"""
Diffraction core calculations for X-ray crystallography.

This module provides fundamental diffraction physics calculations including:
- Scattering omega angle calculations (Bragg condition satisfaction)
- Observable peak generation

PyTorch implementation for differentiable physics and neural network integration.

Based on C++ implementation in Src/DiffractionCore.h and Src/Simulation.cpp
"""

import numpy as np
import torch
from typing import Tuple, Optional, Union
from dataclasses import dataclass

from .constants import KEV_OVER_HBAR_C_IN_ANG


@dataclass
class ScatteringResult:
    """
    Result from scattering omega calculation.

    Attributes:
        observable: Whether the peak can be observed
        omega1: First omega angle (radians), None if not observable
        omega2: Second omega angle (radians), None if not observable
        g_vector: Scattering vector (Å⁻¹)
        g_magnitude: Magnitude of scattering vector (Å⁻¹)
    """
    observable: bool
    omega1: Optional[float]
    omega2: Optional[float]
    g_vector: np.ndarray
    g_magnitude: float


@dataclass
class BatchScatteringResult:
    """
    Batched result from scattering omega calculation.

    Attributes:
        observable: Boolean tensor indicating observability, shape (N,)
        omega1: First omega angles (radians), shape (N,)
        omega2: Second omega angles (radians), shape (N,)
        g_vectors: Scattering vectors (Å⁻¹), shape (N, 3)
        g_magnitudes: Magnitudes of scattering vectors (Å⁻¹), shape (N,)
    """
    observable: torch.Tensor  # bool tensor
    omega1: torch.Tensor      # float tensor
    omega2: torch.Tensor      # float tensor
    g_vectors: torch.Tensor   # (N, 3)
    g_magnitudes: torch.Tensor  # (N,)


def get_scattering_omegas_torch(
    g_vectors: torch.Tensor,
    g_magnitudes: torch.Tensor,
    beam_energy: Union[float, torch.Tensor],
    beam_deflection_chi: Union[float, torch.Tensor] = 0.0,
    epsilon: float = 1e-10
) -> BatchScatteringResult:
    """
    Calculate omega angles for Bragg condition (PyTorch batched version).

    This is the core differentiable physics function for neural network integration.
    Supports batched operations and automatic differentiation.

    Physics:
        - Wavenumber k = 2π/λ = E/(ℏc)
        - Bragg angle: sin(θ) = |G| / (2k)
        - Chi angle: tilt of G relative to z-axis
        - Two solutions for omega correspond to the two positions where
          the scattering vector can satisfy Bragg condition

    Args:
        g_vectors: Scattering vectors (G) in sample frame, shape (N, 3) in Å⁻¹
                   Can also be (3,) for single vector
        g_magnitudes: Magnitudes of scattering vectors |G|, shape (N,) in Å⁻¹
                     Can also be scalar for single vector
        beam_energy: X-ray beam energy in keV, scalar or shape (N,)
        beam_deflection_chi: Beam deflection angle in radians, scalar or shape (N,)
        epsilon: Small value for numerical stability

    Returns:
        BatchScatteringResult containing:
            - observable: Boolean tensor, shape (N,)
            - omega1, omega2: Rotation angles in radians, shape (N,), range [-π, π]
            - g_vectors, g_magnitudes: Echo of inputs

    Example:
        >>> import torch
        >>> # Single reflection
        >>> g_vec = torch.tensor([[2.668, 0.0, 0.0]])
        >>> g_mag = torch.tensor([2.668])
        >>> result = get_scattering_omegas_torch(g_vec, g_mag, beam_energy=50.02)
        >>> print(f"Observable: {result.observable[0]}")
        >>> print(f"Omega1: {result.omega1[0]:.4f} rad")

        >>> # Batched reflections
        >>> g_vecs = torch.tensor([[2.668, 0.0, 0.0],
        ...                        [3.081, 0.0, 0.0],
        ...                        [4.357, 0.0, 0.0]])
        >>> g_mags = torch.norm(g_vecs, dim=1)
        >>> result = get_scattering_omegas_torch(g_vecs, g_mags, beam_energy=50.02)
        >>> observable_count = result.observable.sum()

    Notes:
        - Fully differentiable for autograd
        - Non-observable reflections have omega values set to 0
        - Works on CPU or GPU (follows input tensor device)
        - Numerically stable with epsilon for edge cases
    """
    # Handle single vector case - ensure batch dimension
    if g_vectors.dim() == 1:
        g_vectors = g_vectors.unsqueeze(0)
        g_magnitudes = g_magnitudes.unsqueeze(0) if g_magnitudes.dim() == 0 else g_magnitudes.unsqueeze(0)
        single_input = True
    else:
        single_input = False

    device = g_vectors.device
    batch_size = g_vectors.shape[0]

    # Convert scalars to tensors if needed
    if not isinstance(beam_energy, torch.Tensor):
        beam_energy = torch.tensor(beam_energy, dtype=g_vectors.dtype, device=device)
    if not isinstance(beam_deflection_chi, torch.Tensor):
        beam_deflection_chi = torch.tensor(beam_deflection_chi, dtype=g_vectors.dtype, device=device)

    # Broadcast to batch size if needed
    if beam_energy.dim() == 0:
        beam_energy = beam_energy.expand(batch_size)
    if beam_deflection_chi.dim() == 0:
        beam_deflection_chi = beam_deflection_chi.expand(batch_size)

    # Wavenumber k = E/(ℏc) where E is in keV, result in Å⁻¹
    wavenumber = KEV_OVER_HBAR_C_IN_ANG * beam_energy

    # Bragg angle: sin(θ) = |G| / (2k)
    sin_theta = g_magnitudes / (2.0 * wavenumber)

    # Chi angle: tilt of G relative to z-axis
    # cos(χ) = G_z / |G|
    cos_chi = g_vectors[:, 2] / (g_magnitudes + epsilon)
    # Clamp to avoid numerical issues with sqrt
    cos_chi = torch.clamp(cos_chi, -1.0, 1.0)
    sin_chi = torch.sqrt(1.0 - cos_chi * cos_chi + epsilon)

    # Beam deflection angles
    sin_chi_laue = torch.sin(beam_deflection_chi)
    cos_chi_laue = torch.cos(beam_deflection_chi)

    # Check if solution exists
    # |sin(θ) + cos(χ)sin(χ_Laue)| ≤ |sin(χ)cos(χ_Laue)|
    numerator = sin_theta + cos_chi * sin_chi_laue
    denominator = sin_chi * cos_chi_laue + epsilon

    # Observable condition
    observable = torch.abs(numerator) <= torch.abs(denominator)

    # Δω₀: angle to bring G to nominal position along +y-axis
    delta_omega_0 = torch.atan2(g_vectors[:, 0], g_vectors[:, 1])

    # Δω_b: Bragg rotation angle
    # Clamp to avoid arcsin domain errors
    arcsin_arg = torch.clamp(numerator / denominator, -1.0, 1.0)
    delta_omega_b1 = torch.asin(arcsin_arg)
    delta_omega_b2 = torch.pi - delta_omega_b1

    # Final omega angles
    omega1 = delta_omega_b1 + delta_omega_0
    omega2 = delta_omega_b2 + delta_omega_0

    # Normalize to [-π, π] range
    omega1 = torch.where(omega1 > torch.pi, omega1 - 2.0 * torch.pi, omega1)
    omega1 = torch.where(omega1 < -torch.pi, omega1 + 2.0 * torch.pi, omega1)
    omega2 = torch.where(omega2 > torch.pi, omega2 - 2.0 * torch.pi, omega2)
    omega2 = torch.where(omega2 < -torch.pi, omega2 + 2.0 * torch.pi, omega2)

    # Set non-observable peaks to zero
    omega1 = torch.where(observable, omega1, torch.zeros_like(omega1))
    omega2 = torch.where(observable, omega2, torch.zeros_like(omega2))

    return BatchScatteringResult(
        observable=observable,
        omega1=omega1,
        omega2=omega2,
        g_vectors=g_vectors,
        g_magnitudes=g_magnitudes
    )


def get_scattering_omegas(
    g_vector: Union[np.ndarray, torch.Tensor],
    g_magnitude: float,
    beam_energy: float,
    beam_deflection_chi: float = 0.0
) -> ScatteringResult:
    """
    Calculate the omega angles at which the Bragg condition is satisfied.

    NumPy/PyTorch compatible version for backward compatibility with tests.
    For new code and neural network applications, use get_scattering_omegas_torch.

    This function determines the rotation angles (around the z-axis) at which
    a given scattering vector satisfies the Bragg diffraction condition for
    the specified beam energy and geometry.

    Based on C++ implementation in Simulation.cpp lines 65-119.

    Args:
        g_vector: Scattering vector (G) in sample frame, shape (3,) in Å⁻¹
                 Accepts numpy array or torch tensor
        g_magnitude: Magnitude of scattering vector |G| in Å⁻¹
        beam_energy: X-ray beam energy in keV
        beam_deflection_chi: Beam deflection angle in radians (default 0)

    Returns:
        ScatteringResult containing:
            - observable: True if peak can be observed
            - omega1, omega2: Rotation angles in radians, range [-π, π]
            - g_vector, g_magnitude: Echo of inputs (as numpy arrays)

    Example:
        >>> import numpy as np
        >>> g_vec = np.array([2.668, 0.0, 0.0])  # (111) reflection
        >>> g_mag = 2.668
        >>> result = get_scattering_omegas(g_vec, g_mag, beam_energy=50.02)
        >>> if result.observable:
        ...     print(f"Omega angles: {result.omega1:.3f}, {result.omega2:.3f}")
    """
    # Convert to torch tensor if needed
    if isinstance(g_vector, np.ndarray):
        g_vec_torch = torch.from_numpy(g_vector).float()
        return_numpy = True
    else:
        g_vec_torch = g_vector.float()
        return_numpy = True  # Always return numpy for consistency

    g_mag_torch = torch.tensor([g_magnitude], dtype=torch.float32)

    # Use batched PyTorch implementation
    result = get_scattering_omegas_torch(
        g_vec_torch.unsqueeze(0),
        g_mag_torch,
        beam_energy,
        beam_deflection_chi
    )

    # Extract single result
    observable = result.observable[0].item()

    if observable:
        omega1 = float(result.omega1[0].item())
        omega2 = float(result.omega2[0].item())
    else:
        omega1 = None
        omega2 = None

    # Return numpy array for consistency
    if return_numpy:
        g_vec_out = g_vector if isinstance(g_vector, np.ndarray) else g_vector.cpu().numpy()
    else:
        g_vec_out = g_vector.copy()

    return ScatteringResult(
        observable=observable,
        omega1=omega1,
        omega2=omega2,
        g_vector=g_vec_out.copy() if isinstance(g_vec_out, np.ndarray) else g_vec_out,
        g_magnitude=float(g_magnitude)
    )


def calculate_bragg_angle(
    g_magnitude: Union[float, torch.Tensor],
    beam_energy: Union[float, torch.Tensor]
) -> Union[float, torch.Tensor]:
    """
    Calculate Bragg angle for a given scattering vector magnitude.

    Supports both scalar and batched tensor inputs for neural network integration.

    Args:
        g_magnitude: Magnitude of scattering vector |G| in Å⁻¹
                    Can be float or torch.Tensor
        beam_energy: X-ray beam energy in keV
                    Can be float or torch.Tensor

    Returns:
        Bragg angle θ in radians (float or torch.Tensor)

    Example:
        >>> # Scalar version
        >>> theta = calculate_bragg_angle(g_mag=2.668, beam_energy=50.02)
        >>> print(f"Bragg angle: {np.degrees(theta):.2f}°")

        >>> # Batched tensor version
        >>> import torch
        >>> g_mags = torch.tensor([2.668, 3.081, 4.357])
        >>> thetas = calculate_bragg_angle(g_mags, beam_energy=50.02)
    """
    if isinstance(g_magnitude, torch.Tensor) or isinstance(beam_energy, torch.Tensor):
        # PyTorch path
        if not isinstance(g_magnitude, torch.Tensor):
            g_magnitude = torch.tensor(g_magnitude, dtype=torch.float32)
        if not isinstance(beam_energy, torch.Tensor):
            beam_energy = torch.tensor(beam_energy, dtype=torch.float32)

        wavenumber = KEV_OVER_HBAR_C_IN_ANG * beam_energy
        sin_theta = g_magnitude / (2.0 * wavenumber)
        sin_theta = torch.clamp(sin_theta, -1.0, 1.0)
        return torch.asin(sin_theta)
    else:
        # NumPy path
        wavenumber = KEV_OVER_HBAR_C_IN_ANG * beam_energy
        sin_theta = g_magnitude / (2.0 * wavenumber)
        return float(np.arcsin(sin_theta))


def wavelength_to_wavenumber(wavelength: Union[float, torch.Tensor]) -> Union[float, torch.Tensor]:
    """
    Convert wavelength to wavenumber.

    Args:
        wavelength: X-ray wavelength in Å (float or torch.Tensor)

    Returns:
        Wavenumber k = 2π/λ in Å⁻¹ (float or torch.Tensor)
    """
    if isinstance(wavelength, torch.Tensor):
        return 2.0 * torch.pi / wavelength
    else:
        return 2.0 * np.pi / wavelength


def wavenumber_to_wavelength(wavenumber: Union[float, torch.Tensor]) -> Union[float, torch.Tensor]:
    """
    Convert wavenumber to wavelength.

    Args:
        wavenumber: Wavenumber k in Å⁻¹ (float or torch.Tensor)

    Returns:
        Wavelength λ = 2π/k in Å (float or torch.Tensor)
    """
    if isinstance(wavenumber, torch.Tensor):
        return 2.0 * torch.pi / wavenumber
    else:
        return 2.0 * np.pi / wavenumber


# ============================================================================
# Neural Network Integration Utilities
# ============================================================================

def batch_scattering_omegas_from_reflections(
    hkl_indices: torch.Tensor,
    reciprocal_lattice_param: Union[float, torch.Tensor],
    orientation_matrix: torch.Tensor,
    beam_energy: Union[float, torch.Tensor],
    beam_deflection_chi: Union[float, torch.Tensor] = 0.0,
) -> BatchScatteringResult:
    """
    Calculate omega angles for batched reflections with crystal orientation.

    This is a high-level function for neural network forward simulation that
    combines Miller indices, reciprocal lattice parameters, and orientation
    to compute scattering angles.

    Args:
        hkl_indices: Miller indices, shape (N, 3) or (B, N, 3) for batch
        reciprocal_lattice_param: Reciprocal lattice parameter a* in Å⁻¹
        orientation_matrix: Rotation matrix, shape (3, 3) or (B, 3, 3)
        beam_energy: X-ray beam energy in keV
        beam_deflection_chi: Beam deflection angle in radians

    Returns:
        BatchScatteringResult for all reflections

    Example:
        >>> import torch
        >>> # Define reflections
        >>> hkl = torch.tensor([[1, 1, 1], [2, 0, 0], [2, 2, 0]])
        >>> a_recip = 2 * torch.pi / 4.0782  # Au FCC
        >>>
        >>> # Identity orientation
        >>> orientation = torch.eye(3)
        >>>
        >>> result = batch_scattering_omegas_from_reflections(
        ...     hkl, a_recip, orientation, beam_energy=50.02
        ... )
        >>> print(f"Observable: {result.observable.sum()} / {len(hkl)}")
    """
    # Convert hkl to reciprocal vectors
    g_vectors_crystal = hkl_indices.float() * reciprocal_lattice_param

    # Apply orientation matrix
    # If batched: (B, N, 3) @ (B, 3, 3) -> (B, N, 3)
    # If single: (N, 3) @ (3, 3) -> (N, 3)
    if orientation_matrix.dim() == 3:
        # Batched orientation: (B, 3, 3)
        g_vectors_lab = torch.bmm(
            g_vectors_crystal,
            orientation_matrix.transpose(-2, -1)
        )
    else:
        # Single orientation: (3, 3)
        g_vectors_lab = torch.mm(g_vectors_crystal, orientation_matrix.T)

    # Calculate magnitudes
    g_magnitudes = torch.norm(g_vectors_lab, dim=-1)

    # Flatten for batch processing if needed
    if g_vectors_lab.dim() == 3:
        B, N, _ = g_vectors_lab.shape
        g_vectors_flat = g_vectors_lab.reshape(B * N, 3)
        g_magnitudes_flat = g_magnitudes.reshape(B * N)
    else:
        g_vectors_flat = g_vectors_lab
        g_magnitudes_flat = g_magnitudes

    # Calculate omega angles
    result = get_scattering_omegas_torch(
        g_vectors_flat,
        g_magnitudes_flat,
        beam_energy,
        beam_deflection_chi
    )

    # Reshape back if batched
    if g_vectors_lab.dim() == 3:
        result.observable = result.observable.reshape(B, N)
        result.omega1 = result.omega1.reshape(B, N)
        result.omega2 = result.omega2.reshape(B, N)
        result.g_vectors = result.g_vectors.reshape(B, N, 3)
        result.g_magnitudes = result.g_magnitudes.reshape(B, N)

    return result


# ============================================================================
# Ray Projection and Reflection Functions
# ============================================================================

def get_reflection_vector(
    r_in: torch.Tensor,
    normal: torch.Tensor
) -> torch.Tensor:
    """
    Calculate geometric reflection vector.

    Based on geometric optics: R_out = R_in - 2(R_in·n)n

    Args:
        r_in: Incident ray direction, shape (3,) or (N, 3)
        normal: Surface normal vector, shape (3,) or (N, 3)
                Must be normalized

    Returns:
        Reflected ray direction, same shape as inputs

    C++ Reference:
        DiffractionCore.h:89-98 GetReflectionVector

    Example:
        >>> # Reflect downward ray off horizontal surface
        >>> r_in = torch.tensor([0., 0., -1.])
        >>> normal = torch.tensor([0., 0., 1.])
        >>> r_out = get_reflection_vector(r_in, normal)
        >>> r_out
        tensor([0., 0., 1.])
    """
    # R_out = R_in - 2 * (R_in · normal) * normal
    if r_in.dim() == 1:
        # Single vector
        dot_product = torch.dot(r_in, normal)
    else:
        # Batched vectors
        dot_product = torch.sum(r_in * normal, dim=-1, keepdim=True)

    r_out = r_in - 2.0 * dot_product * normal
    return r_out


def get_reflected_ray_dir(
    sample,
    normal: torch.Tensor,
    beam_dir: torch.Tensor
) -> torch.Tensor:
    """
    Calculate reflected ray direction in lab frame.

    Transforms normal from sample frame to lab frame, then applies
    geometric reflection of beam direction.

    Purpose: Avoids redundant calculation when projecting multiple vertices
    with the same normal (constant interpolation across voxel).

    Args:
        sample: Sample object with to_lab_frame() method
        normal: Surface normal in sample frame, shape (3,)
        beam_dir: Incident beam direction in lab frame, shape (3,)

    Returns:
        Reflected ray direction in lab frame, shape (3,)

    C++ Reference:
        DiffractionCore.h:114-120 GetReflectedRayDir

    Example:
        >>> from icenine.sample import Sample
        >>> sample = Sample()
        >>> sample.set_location(torch.tensor([0., 0., 0.]))
        >>> normal = torch.tensor([0., 0., 1.])  # Sample frame
        >>> beam_dir = torch.tensor([0., 0., 1.])  # Lab frame
        >>> reflected_dir = get_reflected_ray_dir(sample, normal, beam_dir)
    """
    # Transform normal to lab frame
    lab_normal = sample.to_lab_frame(normal)

    # Apply geometric reflection
    reflected_dir = get_reflection_vector(beam_dir, lab_normal)

    return reflected_dir


def build_reflected_ray(
    sample,
    vertex: torch.Tensor,
    ref_dir: torch.Tensor
):
    """
    Build reflected ray from vertex in lab frame.

    Creates a Ray object with origin at vertex (transformed to lab frame)
    and direction along ref_dir.

    Args:
        sample: Sample object with to_lab_frame() method
        vertex: Vertex position in sample frame, shape (3,)
        ref_dir: Reflected ray direction in lab frame, shape (3,)

    Returns:
        Ray object with origin and direction in lab frame

    C++ Reference:
        DiffractionCore.h:132-140 BuildReflectedRay

    Example:
        >>> from icenine.sample import Sample
        >>> from icenine.geometry import Ray
        >>> sample = Sample()
        >>> vertex = torch.tensor([1., 0., 0.])  # Sample frame
        >>> ref_dir = torch.tensor([0., 0., 1.])  # Lab frame
        >>> ray = build_reflected_ray(sample, vertex, ref_dir)
        >>> ray.direction
        tensor([0., 0., 1.])
    """
    from .geometry import Ray

    # Transform vertex to lab frame
    lab_vertex = sample.to_lab_frame(vertex)

    # Create ray
    reflected_ray = Ray(origin=lab_vertex, direction=ref_dir)

    return reflected_ray


def get_reflected_ray(
    sample,
    vertex: torch.Tensor,
    normal: torch.Tensor,
    beam_dir: torch.Tensor
):
    """
    Calculate reflected ray from vertex with given normal (convenience wrapper).

    Combines get_reflected_ray_dir() and build_reflected_ray() into single call.

    Args:
        sample: Sample object with to_lab_frame() method
        vertex: Vertex position in sample frame, shape (3,)
        normal: Surface normal in sample frame, shape (3,)
        beam_dir: Incident beam direction in lab frame, shape (3,)

    Returns:
        Ray object representing reflected ray in lab frame

    C++ Reference:
        DiffractionCore.h:150-156 GetReflectedRay

    Example:
        >>> from icenine.sample import Sample
        >>> sample = Sample()
        >>> vertex = torch.tensor([1., 0., 0.])
        >>> normal = torch.tensor([0., 0., 1.])
        >>> beam_dir = torch.tensor([0., 0., 1.])
        >>> ray = get_reflected_ray(sample, vertex, normal, beam_dir)
    """
    # Calculate reflected direction
    ref_dir = get_reflected_ray_dir(sample, normal, beam_dir)

    # Build ray from vertex
    reflected_ray = build_reflected_ray(sample, vertex, ref_dir)

    return reflected_ray


def get_illuminated_pixel(
    detector,
    incident_ray
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Calculate pixel coordinates where ray intersects detector.

    Performs ray-plane intersection with detector, then converts intersection
    point from lab coordinates to pixel coordinates (row, col).

    Args:
        detector: Detector object with intersect_ray() and lab_to_pixel() methods
        incident_ray: Ray object (incident upon the detector)

    Returns:
        Tuple of (hit, pixel_col, pixel_row):
            - hit: Boolean tensor indicating if ray intersects detector
            - pixel_col: Column coordinate (J direction), valid if hit=True
            - pixel_row: Row coordinate (K direction), valid if hit=True

    C++ Reference:
        DiffractionCore.h:165-179 GetIlluminatedPixel

    Notes:
        - C++ returns Point with (x=col, y=row) convention
        - Python returns separate (col, row) tensors
        - Non-intersecting rays return (False, 0, 0)

    Example:
        >>> from icenine.detector import Detector
        >>> from icenine.geometry import Ray
        >>> import torch
        >>>
        >>> # Create detector
        >>> detector = Detector(num_rows=1024, num_cols=1024,
        ...                     beam_center_j=512, beam_center_k=512)
        >>>
        >>> # Ray hitting detector center
        >>> ray = Ray(
        ...     origin=torch.tensor([0., 0., -1.]),
        ...     direction=torch.tensor([0., 0., 1.])
        ... )
        >>> hit, col, row = get_illuminated_pixel(detector, ray)
        >>> hit
        tensor(True)
    """
    # Perform ray-detector intersection
    intersects, t = detector.intersect_ray(incident_ray)

    if not intersects.item():
        # No intersection
        return (
            torch.tensor(False),
            torch.tensor(0.0),
            torch.tensor(0.0)
        )

    # Calculate intersection point in lab frame
    # C++: oIntersectLoc = oIncidentRay.Evaluate(fT)
    intersect_loc = incident_ray.at(t)

    # Convert to pixel coordinates
    # C++: oDetector.LabToPixel(ImageRow, ImageCol, oIntersectLoc)
    pixel_row, pixel_col = detector.lab_to_pixel(intersect_loc)

    # C++ convention: p.x = ImageCol, p.y = ImageRow
    return intersects, pixel_col, pixel_row
