"""
Geometric primitives and transformations for X-ray diffraction.

This module provides:
- Euler angle ↔ rotation matrix conversions (Bunge/ZXZ convention)
- Geometric primitives (Plane, Ray) for ray-detector intersection
- PyTorch-based implementations for differentiability and GPU support

All implementations match the C++ IceNine code for compatibility.

Coordinate Conventions:
- Euler angles: Bunge convention (ZXZ intrinsic rotations)
- Rotation matrices: Active rotation convention
- C++ reference: XDM++/libXDM/3dMath.cpp

Author: S. F. Li
"""

import numpy as np
import torch
from scipy.spatial.transform import Rotation
from typing import Tuple
from dataclasses import dataclass


# =============================================================================
# Euler Angle Conversions (Bunge Convention: ZXZ Intrinsic)
# =============================================================================


def euler_to_matrix(phi1_deg: float, Phi_deg: float, phi2_deg: float) -> np.ndarray:
    """
    Convert Bunge Euler angles (degrees) to rotation matrix.

    Uses ZXZ (intrinsic) convention matching C++ BuildActiveEulerMatrix.
    This is the "active" rotation convention used in the IceNine C++ code.

    C++ reference: XDM++/libXDM/3dMath.cpp BuildActiveEulerMatrix (lines 152-174)

    The Bunge convention defines rotations as:
    R = Rz(φ₁) @ Rx(Φ) @ Rz(φ₂)

    Args:
        phi1_deg: First Euler angle φ₁ in degrees [0, 360]
        Phi_deg: Second Euler angle Φ in degrees [0, 180]
        phi2_deg: Third Euler angle φ₂ in degrees [0, 360]

    Returns:
        3x3 rotation matrix (numpy array, float32)

    Example:
        >>> # Identity rotation
        >>> R = euler_to_matrix(0, 0, 0)
        >>> np.allclose(R, np.eye(3))
        True

        >>> # 90 degree rotation around Z
        >>> R = euler_to_matrix(90, 0, 0)
        >>> np.allclose(R @ np.array([1, 0, 0]), np.array([0, 1, 0]), atol=1e-6)
        True

    Note:
        This is a NumPy/scipy implementation for I/O and compatibility.
        For differentiable operations, use euler_to_matrix_torch().
    """
    # Convert degrees to radians
    phi1 = np.deg2rad(phi1_deg)
    Phi = np.deg2rad(Phi_deg)
    phi2 = np.deg2rad(phi2_deg)

    # Use scipy Rotation with ZXZ extrinsic = ZXZ intrinsic
    # scipy's 'ZXZ' is extrinsic, which equals intrinsic in reverse order
    rot = Rotation.from_euler("ZXZ", [phi1, Phi, phi2], degrees=False)
    return rot.as_matrix().astype(np.float32)


def matrix_to_euler(rotation_matrix: np.ndarray) -> Tuple[float, float, float]:
    """
    Convert rotation matrix to Bunge Euler angles (degrees).

    Uses ZXZ (intrinsic) convention matching C++ GetEulerAngles.

    C++ reference: XDM++/libXDM/3dMath.cpp GetEulerAngles (lines 180-209)

    Args:
        rotation_matrix: 3x3 rotation matrix (numpy array)

    Returns:
        Tuple of (φ₁, Φ, φ₂) in degrees
        - φ₁ ∈ [0, 360]
        - Φ ∈ [0, 180]
        - φ₂ ∈ [0, 360]

    Example:
        >>> R = np.eye(3)
        >>> phi1, Phi, phi2 = matrix_to_euler(R)
        >>> phi1, Phi, phi2
        (0.0, 0.0, 0.0)

        >>> # Round-trip test
        >>> R = euler_to_matrix(45, 30, 60)
        >>> phi1, Phi, phi2 = matrix_to_euler(R)
        >>> R2 = euler_to_matrix(phi1, Phi, phi2)
        >>> np.allclose(R, R2)
        True

    Note:
        This is a NumPy/scipy implementation for I/O and compatibility.
        For differentiable operations, use matrix_to_euler_torch().
    """
    rot = Rotation.from_matrix(rotation_matrix)
    euler_angles_rad = rot.as_euler("ZXZ", degrees=False)

    # Convert to degrees
    phi1_deg = float(np.rad2deg(euler_angles_rad[0]))
    Phi_deg = float(np.rad2deg(euler_angles_rad[1]))
    phi2_deg = float(np.rad2deg(euler_angles_rad[2]))

    # Ensure positive angles [0, 360] for phi1 and phi2, [0, 180] for Phi
    phi1_deg = phi1_deg % 360.0
    phi2_deg = phi2_deg % 360.0
    Phi_deg = Phi_deg % 180.0

    return phi1_deg, Phi_deg, phi2_deg


def euler_to_matrix_torch(
    phi1_deg: torch.Tensor, Phi_deg: torch.Tensor, phi2_deg: torch.Tensor
) -> torch.Tensor:
    """
    Convert batched Euler angles (degrees) to rotation matrices (PyTorch).

    Differentiable version for neural network integration and gradient-based
    optimization. Supports both scalar and batched inputs.

    Uses ZXZ (intrinsic) Bunge convention matching C++ BuildActiveEulerMatrix.

    C++ reference: XDM++/libXDM/3dMath.cpp BuildActiveEulerMatrix (lines 152-174)

    The implementation directly evaluates the matrix formula:
    R = Rz(φ₁) @ Rx(Φ) @ Rz(φ₂)

    Where the matrix elements are:
    m[0][0] = cos(φ₁)cos(φ₂) - sin(φ₁)cos(Φ)sin(φ₂)
    m[1][0] = sin(φ₁)cos(φ₂) + cos(φ₁)cos(Φ)sin(φ₂)
    m[2][0] = sin(Φ)sin(φ₂)
    ... (9 elements total)

    Args:
        phi1_deg: First Euler angles in degrees, shape () or (N,)
        Phi_deg: Second Euler angles in degrees, shape () or (N,)
        phi2_deg: Third Euler angles in degrees, shape () or (N,)

    Returns:
        Rotation matrices:
        - Shape (3, 3) if scalar input
        - Shape (N, 3, 3) if batched input

    Example:
        >>> # Scalar input
        >>> phi1 = torch.tensor(0.)
        >>> Phi = torch.tensor(0.)
        >>> phi2 = torch.tensor(0.)
        >>> R = euler_to_matrix_torch(phi1, Phi, phi2)
        >>> R.shape
        torch.Size([3, 3])
        >>> torch.allclose(R, torch.eye(3))
        True

        >>> # Batched input
        >>> phi1 = torch.tensor([0., 90.])
        >>> Phi = torch.tensor([0., 0.])
        >>> phi2 = torch.tensor([0., 0.])
        >>> R = euler_to_matrix_torch(phi1, Phi, phi2)
        >>> R.shape
        torch.Size([2, 3, 3])

        >>> # Differentiable
        >>> phi1 = torch.tensor(45., requires_grad=True)
        >>> Phi = torch.tensor(0.)
        >>> phi2 = torch.tensor(0.)
        >>> R = euler_to_matrix_torch(phi1, Phi, phi2)
        >>> loss = R.sum()
        >>> loss.backward()
        >>> phi1.grad is not None
        True

    Note:
        All inputs must have the same batch dimension or be scalars.
        The function preserves gradients for automatic differentiation.
    """
    # Convert to radians
    phi1 = torch.deg2rad(phi1_deg)
    Phi = torch.deg2rad(Phi_deg)
    phi2 = torch.deg2rad(phi2_deg)

    # Build rotation matrix using ZXZ convention
    # R = Rz(phi1) @ Rx(Phi) @ Rz(phi2)
    cos_phi1 = torch.cos(phi1)
    sin_phi1 = torch.sin(phi1)
    cos_Phi = torch.cos(Phi)
    sin_Phi = torch.sin(Phi)
    cos_phi2 = torch.cos(phi2)
    sin_phi2 = torch.sin(phi2)

    # Build matrix (matching C++ BuildActiveEulerMatrix line-by-line)
    # C++ indices: m[row][col], Python: m[row, col]
    m00 = cos_phi1 * cos_phi2 - sin_phi1 * cos_Phi * sin_phi2
    m10 = sin_phi1 * cos_phi2 + cos_phi1 * cos_Phi * sin_phi2
    m20 = sin_Phi * sin_phi2

    m01 = -cos_phi1 * sin_phi2 - sin_phi1 * cos_Phi * cos_phi2
    m11 = -sin_phi1 * sin_phi2 + cos_phi1 * cos_Phi * cos_phi2
    m21 = sin_Phi * cos_phi2

    m02 = sin_phi1 * sin_Phi
    m12 = -cos_phi1 * sin_Phi
    m22 = cos_Phi

    # Check if input is scalar or batched
    if phi1.dim() == 0:
        # Scalar input - stack tensors to preserve gradients
        return torch.stack(
            [
                torch.stack([m00, m01, m02]),
                torch.stack([m10, m11, m12]),
                torch.stack([m20, m21, m22]),
            ]
        )
    else:
        # Batched input - stack along batch dimension
        return torch.stack(
            [
                torch.stack([m00, m01, m02], dim=-1),
                torch.stack([m10, m11, m12], dim=-1),
                torch.stack([m20, m21, m22], dim=-1),
            ],
            dim=-2,
        )


def matrix_to_euler_torch(
    rotation_matrix: torch.Tensor,
    epsilon: float = 1e-6
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Convert rotation matrix to Bunge Euler angles (degrees) - PyTorch version.

    Differentiable inverse of euler_to_matrix_torch(). Supports both scalar
    and batched inputs.

    C++ reference: XDM++/libXDM/3dMath.cpp GetEulerAngles (lines 180-209)

    Args:
        rotation_matrix: Rotation matrix
            - Shape (3, 3) for single matrix
            - Shape (N, 3, 3) for batched matrices
        epsilon: Threshold for gimbal lock detection (default: 1e-6)

    Returns:
        Tuple of (phi1_deg, Phi_deg, phi2_deg) in degrees
        - phi1_deg ∈ [0, 360], shape () or (N,)
        - Phi_deg ∈ [0, 180], shape () or (N,)
        - phi2_deg ∈ [0, 360], shape () or (N,)

    Example:
        >>> # Scalar input
        >>> R = torch.eye(3)
        >>> phi1, Phi, phi2 = matrix_to_euler_torch(R)
        >>> phi1.item(), Phi.item(), phi2.item()
        (0.0, 0.0, 0.0)

        >>> # Round-trip test
        >>> R = euler_to_matrix_torch(
        ...     torch.tensor(45.), torch.tensor(30.), torch.tensor(60.)
        ... )
        >>> phi1, Phi, phi2 = matrix_to_euler_torch(R)
        >>> R2 = euler_to_matrix_torch(phi1, Phi, phi2)
        >>> torch.allclose(R, R2, atol=1e-5)
        True

        >>> # Batched input
        >>> R = torch.eye(3).unsqueeze(0).repeat(5, 1, 1)
        >>> phi1, Phi, phi2 = matrix_to_euler_torch(R)
        >>> phi1.shape
        torch.Size([5])

    Note:
        Handles gimbal lock singularities when Φ ≈ 0 or Φ ≈ π using
        differentiable torch.where() operations to maintain gradient flow.

    Raises:
        NotImplementedError: This function needs to be implemented
            (see DETECTOR_IMPLEMENTATION_PLAN.md Phase 1.2)

    Note:
        This function is not currently needed for Detector implementation
        since we don't differentiate through matrix → Euler conversions.
        It can be implemented later if required.
    """
    raise NotImplementedError(
        "matrix_to_euler_torch() not yet implemented. "
        "Not needed for current Detector implementation. "
        "See DETECTOR_IMPLEMENTATION_PLAN.md Phase 1.2 if implementation is required."
    )


# =============================================================================
# Geometric Primitives for Ray-Detector Intersection
# =============================================================================


@dataclass
class Plane:
    """
    Plane in 3D space: A*x + B*y + C*z + D = 0

    Represents a plane using the implicit equation where [A, B, C] is the
    normal vector and D is the offset from the origin.

    Attributes:
        coeffs: Plane coefficients [A, B, C, D]
                Shape (4,) for single plane
                Shape (N, 4) for batched planes

    Example:
        >>> # XY plane (z = 0)
        >>> plane = Plane(torch.tensor([0., 0., 1., 0.]))
        >>> plane.normal
        tensor([0., 0., 1.])

        >>> # Normalized plane
        >>> plane = Plane(torch.tensor([1., 1., 1., 5.]))
        >>> normalized = plane.normalize()
        >>> torch.allclose(torch.norm(normalized.normal), torch.tensor(1.))
        True
    """
    coeffs: torch.Tensor

    @property
    def normal(self) -> torch.Tensor:
        """
        Return plane normal vector [A, B, C].

        Returns:
            Normal vector, shape (3,) or (N, 3)
        """
        return self.coeffs[..., :3]

    @property
    def d(self) -> torch.Tensor:
        """
        Return plane offset D.

        Returns:
            Offset value, shape () or (N,)
        """
        return self.coeffs[..., 3]

    def normalize(self) -> 'Plane':
        """
        Normalize plane equation so |normal| = 1.

        Returns:
            New Plane with normalized coefficients

        Example:
            >>> plane = Plane(torch.tensor([2., 0., 0., 10.]))
            >>> norm_plane = plane.normalize()
            >>> norm_plane.coeffs
            tensor([1., 0., 0., 5.])
        """
        # Compute norm of normal vector
        norm = torch.norm(self.normal, dim=-1, keepdim=True)

        # Avoid division by zero
        norm = torch.where(norm > 1e-10, norm, torch.ones_like(norm))

        # Normalize coefficients
        if self.coeffs.dim() == 1:
            # Single plane
            normalized_coeffs = self.coeffs / norm.squeeze()
        else:
            # Batched planes
            normalized_coeffs = self.coeffs / norm

        return Plane(normalized_coeffs)

    def distance_to_point(self, point: torch.Tensor) -> torch.Tensor:
        """
        Calculate signed distance from point(s) to plane.

        Positive distance means point is on the side of the plane that
        the normal vector points to.

        Args:
            point: 3D point(s), shape (3,) or (N, 3)

        Returns:
            Signed distance, shape () or (N,)

        Example:
            >>> # Point above XY plane
            >>> plane = Plane(torch.tensor([0., 0., 1., 0.]))
            >>> point = torch.tensor([0., 0., 5.])
            >>> plane.distance_to_point(point)
            tensor(5.)
        """
        # Normalize plane first for accurate distance
        normalized = self.normalize()

        # Distance = N·P + D
        if point.dim() == 1:
            # Single point
            distance = torch.dot(normalized.normal, point) + normalized.d
        else:
            # Batched points
            distance = torch.sum(normalized.normal * point, dim=-1) + normalized.d

        return distance


@dataclass
class Ray:
    """
    Ray in 3D space: P(t) = origin + t * direction

    Represents a ray with an origin point and a direction vector.

    Attributes:
        origin: Ray origin point, shape (3,) or (N, 3)
        direction: Ray direction vector, shape (3,) or (N, 3)
                  (does not need to be normalized)

    Example:
        >>> # Ray along X-axis from origin
        >>> ray = Ray(
        ...     origin=torch.tensor([0., 0., 0.]),
        ...     direction=torch.tensor([1., 0., 0.])
        ... )
        >>> ray.origin
        tensor([0., 0., 0.])
    """
    origin: torch.Tensor
    direction: torch.Tensor

    def intersect_plane(
        self,
        plane: Plane,
        epsilon: float = 1e-8
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute ray-plane intersection.

        Solves for t in the ray equation P(t) = O + t*D where P lies on
        the plane: N·P + d = 0

        Substituting: N·(O + t*D) + d = 0
        Solving for t: t = -(N·O + d) / (N·D)

        Args:
            plane: Plane to intersect with
            epsilon: Threshold for parallel detection (default: 1e-8)

        Returns:
            Tuple of (intersects, t):
            - intersects: Boolean tensor indicating if intersection exists
                         Shape () or (N,)
            - t: Intersection parameter (valid only where intersects=True)
                 Shape () or (N,)

        Example:
            >>> # Ray along Z-axis intersecting XY plane
            >>> ray = Ray(
            ...     origin=torch.tensor([0., 0., -5.]),
            ...     direction=torch.tensor([0., 0., 1.])
            ... )
            >>> plane = Plane(torch.tensor([0., 0., 1., 0.]))  # z = 0
            >>> intersects, t = ray.intersect_plane(plane)
            >>> intersects
            tensor(True)
            >>> t
            tensor(5.)

            >>> # Parallel ray (no intersection)
            >>> ray = Ray(
            ...     origin=torch.tensor([0., 0., 1.]),
            ...     direction=torch.tensor([1., 0., 0.])
            ... )
            >>> plane = Plane(torch.tensor([0., 0., 1., 0.]))
            >>> intersects, t = ray.intersect_plane(plane)
            >>> intersects
            tensor(False)
        """
        # Compute denominator: N·D
        if self.direction.dim() == 1:
            # Single ray
            denominator = torch.dot(plane.normal, self.direction)
        else:
            # Batched rays
            denominator = torch.sum(plane.normal * self.direction, dim=-1)

        # Check if ray is parallel to plane (denominator ≈ 0)
        intersects = torch.abs(denominator) > epsilon

        # Compute numerator: -(N·O + d)
        if self.origin.dim() == 1:
            # Single ray
            numerator = -(torch.dot(plane.normal, self.origin) + plane.d)
        else:
            # Batched rays
            numerator = -(torch.sum(plane.normal * self.origin, dim=-1) + plane.d)

        # Compute t (set to 0 where no intersection)
        t = torch.where(
            intersects,
            numerator / (denominator + epsilon * (~intersects).float()),
            torch.zeros_like(denominator)
        )

        # Also check that intersection is in front of ray (t > 0)
        intersects = intersects & (t > 0)

        return intersects, t

    def at(self, t: torch.Tensor) -> torch.Tensor:
        """
        Evaluate ray position at parameter t: P(t) = O + t*D

        Args:
            t: Ray parameter, shape () or (N,)

        Returns:
            Position on ray, shape (3,) or (N, 3)

        Example:
            >>> ray = Ray(
            ...     origin=torch.tensor([1., 0., 0.]),
            ...     direction=torch.tensor([0., 1., 0.])
            ... )
            >>> ray.at(torch.tensor(5.))
            tensor([1., 5., 0.])
        """
        if t.dim() == 0:
            # Scalar t
            return self.origin + t * self.direction
        else:
            # Batched t
            return self.origin + t.unsqueeze(-1) * self.direction


def passive_euler_matrix(phi: float, theta: float, psi: float) -> torch.Tensor:
    """
    Build passive Euler rotation matrix (sample → lab frame transformation).

    Matches C++ SMatrix4x4::SetPassiveEulerMatrix from XDM++/libXDM/3dMath.cpp:679-698.

    This is the EXACT convention used by CSample::SetOrientation() for transforming
    vectors from the sample frame to the lab frame.

    Args:
        phi: First Euler angle (radians)
        theta: Second Euler angle (radians)
        psi: Third Euler angle (radians)

    Returns:
        3x3 rotation matrix as torch.Tensor

    Note:
        Unlike the active Euler matrix (ZXZ Bunge convention used for voxel orientations),
        this uses a passive rotation convention for global sample orientation.
    """
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)

    # Matrix elements from C++ lines 687-697
    m00 = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi
    m10 = -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi
    m20 = sin_theta * sin_phi

    m01 = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi
    m11 = -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi
    m21 = -sin_theta * cos_phi

    m02 = sin_psi * sin_theta
    m12 = cos_psi * sin_theta
    m22 = cos_theta

    return torch.tensor([
        [m00, m01, m02],
        [m10, m11, m12],
        [m20, m21, m22]
    ], dtype=torch.float32)
