"""
Sample class for IceNine microstructure and coordinate transformations.

Python port of CSample from Src/Sample.h and Src/Sample.cpp.

The Sample class manages:
- Global sample position and orientation in lab frame
- Microstructure voxel grid (.mic file I/O)
- Crystal structure(s) for diffraction calculations
- Coordinate transformations (sample → lab frame)
- Multi-phase material support

Author: S. F. Li
"""

from typing import List, Optional
import torch
import numpy as np

from .mic_file import MicFile
from .crystal_structure import CrystalStructure
from .geometry import euler_to_matrix_torch, matrix_to_euler, passive_euler_matrix
from .symmetry import CrystalSymmetry


class Sample:
    """
    Sample with voxel grid and coordinate transformations.

    This class represents a polycrystalline sample with:
    - A global orientation and position in the lab frame
    - A microstructure described by a voxel grid (each voxel has local orientation)
    - One or more crystal structures (multi-phase materials supported)

    The key operation is transforming vectors from the sample frame to the lab frame,
    which is needed for every diffraction calculation.

    Python port of C++ CSample class (Src/Sample.h, Src/Sample.cpp).

    Attributes:
        sample_to_lab_matrix: 4x4 homogeneous transformation matrix (sample → lab)
        location: Sample center position in lab frame [x, y, z] (meters)
        orientation_euler: Euler angles [phi, theta, psi] (radians)
        mic_file: Microstructure voxel grid (MicFile object)
        crystal_structures: List of crystal structures (multi-phase support)
        symmetry: Crystal symmetry for the sample

    Example:
        >>> sample = Sample()
        >>> sample.set_orientation(90, 0, 0)  # Rotate 90° around Z
        >>> sample.set_location(np.array([0, 0, 0.1]))  # 10 cm above origin
        >>> sample.load_sample("data/sample.mic")
        >>> sample.add_crystal_structure(gold_structure)
    """

    def __init__(self):
        """Initialize sample with identity transformation.

        C++ Reference: Sample.cpp:33-40 (CSample constructor)
        """
        # Transformation matrix (4x4 homogeneous coordinates)
        # [R | t]  where R is 3x3 rotation, t is 3x1 translation
        # [0 | 1]
        self.sample_to_lab_matrix: torch.Tensor = torch.eye(4, dtype=torch.float32)

        # Sample position and orientation
        self.location: torch.Tensor = torch.zeros(3, dtype=torch.float32)
        self.orientation_euler: torch.Tensor = torch.zeros(3, dtype=torch.float32)

        # Microstructure voxel grid
        self.mic_file: Optional[MicFile] = None

        # Crystal structure(s) - support multi-phase materials
        self.crystal_structures: List[CrystalStructure] = []

        # Crystal symmetry
        self.symmetry: Optional[CrystalSymmetry] = None

    # =========================================================================
    # Coordinate Transformations
    # =========================================================================

    def to_lab_frame(self, vectors: torch.Tensor) -> torch.Tensor:
        """
        Transform vector(s) from sample frame to lab frame.

        This is the MOST FREQUENTLY CALLED method in the entire IceNine codebase.
        It is called in the innermost loop of diffraction calculations for every
        voxel, every vertex, and every reflection.

        C++ Reference: Sample.cpp:181-187 (ToLabFrame)

        Performance notes:
        - Supports batched operations (N vectors at once)
        - Uses PyTorch for GPU acceleration
        - Preserves gradients for optimization

        Args:
            vectors: Shape (3,) or (N, 3) - vectors in sample frame

        Returns:
            Transformed vectors in lab frame, same shape as input

        Example:
            >>> sample = Sample()
            >>> sample.set_orientation(90, 0, 0)  # 90° rotation around Z
            >>> v = torch.tensor([1.0, 0.0, 0.0])
            >>> v_lab = sample.to_lab_frame(v)
            >>> print(v_lab)  # Should be approximately [0, 1, 0]
        """
        # Extract 3x3 rotation part from 4x4 transformation matrix
        rotation = self.sample_to_lab_matrix[:3, :3]

        if vectors.dim() == 1:
            # Single vector (3,): result = R @ v
            return torch.matmul(rotation, vectors)
        else:
            # Batched vectors (N, 3): result = v @ R^T
            # This is more efficient than looping and gives same result:
            # (v @ R^T)[i] = v[i] @ R^T = (R @ v[i])^T^T = R @ v[i]
            return torch.matmul(vectors, rotation.T)

    def to_lab_frame_4d(self, vectors: torch.Tensor) -> torch.Tensor:
        """
        Transform 4D homogeneous vector(s) from sample frame to lab frame.

        Includes both rotation and translation.

        C++ Reference: Sample.cpp:171-179 (ToLabFrame 4D version)

        Args:
            vectors: Shape (4,) or (N, 4) - homogeneous vectors in sample frame

        Returns:
            Transformed vectors in lab frame, same shape as input
        """
        if vectors.dim() == 1:
            return torch.matmul(self.sample_to_lab_matrix, vectors)
        else:
            return torch.matmul(vectors, self.sample_to_lab_matrix.T)

    def set_location(self, location: np.ndarray) -> None:
        """
        Set sample location in lab frame.

        C++ Reference: Sample.cpp:66-72 (SetLocation)

        Args:
            location: Position [x, y, z] in meters
        """
        self.location = torch.from_numpy(location.astype(np.float32))
        # Update translation part of 4x4 matrix
        self.sample_to_lab_matrix[:3, 3] = self.location

    def translate(self, displacement: np.ndarray) -> None:
        """
        Apply translation to current location.

        C++ Reference: Sample.cpp:157-163 (Translate)

        Args:
            displacement: Translation vector [dx, dy, dz] in meters
        """
        displacement_tensor = torch.from_numpy(displacement.astype(np.float32))
        self.location += displacement_tensor
        self.sample_to_lab_matrix[:3, 3] = self.location

    def set_orientation(self, phi: float, theta: float, psi: float) -> None:
        """
        Set global sample orientation using Euler angles.

        Uses passive Euler angle convention matching C++ SetPassiveEulerMatrix.
        This is the convention used for global sample orientation (different from
        the active ZXZ Bunge convention used for individual voxel orientations).

        C++ Reference: Sample.cpp:80-94 (SetOrientation - Euler angles version)

        Args:
            phi: First Euler angle (degrees)
            theta: Second Euler angle (degrees)
            psi: Third Euler angle (degrees)

        Example:
            >>> sample = Sample()
            >>> sample.set_orientation(90, 0, 0)  # Rotate 90° around Z
        """
        # Convert to radians
        phi_rad = np.deg2rad(phi)
        theta_rad = np.deg2rad(theta)
        psi_rad = np.deg2rad(psi)

        # Build rotation matrix using passive Euler angle convention
        # This matches C++ SetPassiveEulerMatrix exactly
        rotation = passive_euler_matrix(phi_rad, theta_rad, psi_rad)

        # Update 3x3 rotation part of 4x4 matrix
        self.sample_to_lab_matrix[:3, :3] = rotation

        # Store Euler angles
        self.orientation_euler = torch.tensor(
            [phi_rad, theta_rad, psi_rad], dtype=torch.float32
        )

    def set_orientation_matrix(self, rotation_matrix: np.ndarray) -> None:
        """
        Set global sample orientation using rotation matrix.

        C++ Reference: Sample.cpp:102-115 (SetOrientation - matrix version)

        Args:
            rotation_matrix: 3x3 rotation matrix (sample → lab)
        """
        rotation = torch.from_numpy(rotation_matrix.astype(np.float32))
        self.sample_to_lab_matrix[:3, :3] = rotation

        # Extract Euler angles from rotation matrix
        phi, theta, psi = matrix_to_euler(rotation.numpy())
        self.orientation_euler = torch.tensor(
            [phi, theta, psi], dtype=torch.float32
        )

    def rotate(self, phi: float, theta: float, psi: float) -> None:
        """
        Apply rotation to existing orientation (composition).

        Composes the new rotation with the existing transformation:
        M_new = R(phi, theta, psi) @ M_old

        C++ Reference: Sample.cpp:123-136 (Rotate - Euler version)

        Args:
            phi: Rotation around Z-axis (degrees)
            theta: Rotation around X-axis (degrees)
            psi: Rotation around Z-axis (degrees)
        """
        # Convert to radians
        phi_rad = np.deg2rad(phi)
        theta_rad = np.deg2rad(theta)
        psi_rad = np.deg2rad(psi)

        # Build rotation matrix using passive convention
        R_new = passive_euler_matrix(phi_rad, theta_rad, psi_rad)

        # Compose with existing rotation
        rotation_old = self.sample_to_lab_matrix[:3, :3]
        self.sample_to_lab_matrix[:3, :3] = torch.matmul(R_new, rotation_old)

        # Update Euler angles (Note: extraction may not be unique due to gimbal lock)
        # For now, we don't update orientation_euler to avoid gimbal lock issues
        # The authoritative source is the rotation matrix itself

    def rotate_axis_angle(self, axis: np.ndarray, angle_deg: float) -> None:
        """
        Apply axis-angle rotation to existing orientation.

        C++ Reference: Sample.cpp:144-149 (Rotate - axis-angle version)

        Args:
            axis: Rotation axis (will be normalized)
            angle_deg: Rotation angle in degrees
        """
        # Normalize axis
        axis = axis / np.linalg.norm(axis)
        angle_rad = np.deg2rad(angle_deg)

        # Rodrigues' rotation formula for PASSIVE rotation (transpose of active)
        # Active: R = I + sin(θ)K + (1-cos(θ))K²
        # Passive: R = I - sin(θ)K + (1-cos(θ))K²  (note the minus sign)
        # where K is the skew-symmetric cross-product matrix
        K = torch.tensor([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ], dtype=torch.float32)

        I = torch.eye(3, dtype=torch.float32)
        # Use negative angle for passive rotation (equivalent to transpose)
        R = I - np.sin(angle_rad) * K + (1 - np.cos(angle_rad)) * torch.matmul(K, K)

        # Compose with existing rotation
        rotation_old = self.sample_to_lab_matrix[:3, :3]
        self.sample_to_lab_matrix[:3, :3] = torch.matmul(R, rotation_old)

    def rotate_z(self, omega_deg: float) -> None:
        """
        Apply rotation around Z-axis (optimized version).

        This is a performance-optimized version for the common case of
        rotating around the Z-axis, which is frequently used during
        sample rotation in experiments.

        C++ Reference: Sample.cpp:165-169 (RotateZ - hand-optimized)

        Args:
            omega_deg: Rotation angle around Z-axis in degrees
        """
        # For now, use general rotate() method
        # Can optimize later with hand-coded matrix if needed
        self.rotate(omega_deg, 0, 0)

    # =========================================================================
    # Accessors
    # =========================================================================

    def get_location(self) -> np.ndarray:
        """
        Get sample location in lab frame.

        C++ Reference: Sample.h:96 (GetLocation)

        Returns:
            Position [x, y, z] in meters
        """
        return self.location.numpy()

    def get_orientation(self) -> np.ndarray:
        """
        Get sample orientation as Euler angles.

        C++ Reference: Sample.h:97 (GetOrientation)

        Returns:
            Euler angles [phi, theta, psi] in radians
        """
        return self.orientation_euler.numpy()

    def get_orientation_matrix(self) -> np.ndarray:
        """
        Get sample orientation as rotation matrix.

        C++ Reference: Sample.h:98 (GetOrientationMatrix)

        Returns:
            3x3 rotation matrix (sample → lab)
        """
        return self.sample_to_lab_matrix[:3, :3].numpy()

    def get_transformation_matrix(self) -> np.ndarray:
        """
        Get full 4x4 transformation matrix (sample → lab).

        Returns:
            4x4 homogeneous transformation matrix
        """
        return self.sample_to_lab_matrix.numpy()

    # =========================================================================
    # File I/O
    # =========================================================================

    def load_sample(self, filename: str) -> bool:
        """
        Load microstructure from .mic file.

        C++ Reference: Sample.cpp:189-195 (LoadSample)

        Args:
            filename: Path to .mic file

        Returns:
            True if successful, False otherwise

        Example:
            >>> sample = Sample()
            >>> success = sample.load_sample("data/Au1007_small.mic")
            >>> if success:
            ...     print(f"Loaded {len(sample.mic_file.voxels)} voxels")
        """
        try:
            self.mic_file = MicFile.read(filename)
            return True
        except Exception as e:
            print(f"Failed to load sample from {filename}: {e}")
            return False

    def save_sample(self, filename: str) -> bool:
        """
        Save microstructure to .mic file.

        Args:
            filename: Path to output .mic file

        Returns:
            True if successful, False otherwise
        """
        if self.mic_file is None:
            print("No mic file to save")
            return False

        try:
            self.mic_file.save(filename)
            return True
        except Exception as e:
            print(f"Failed to save sample to {filename}: {e}")
            return False

    # =========================================================================
    # Crystal Structure Management
    # =========================================================================

    def add_crystal_structure(self, structure: CrystalStructure) -> None:
        """
        Add crystal structure for a phase.

        Supports multi-phase materials. Each voxel in the MicFile has a
        phase ID that indexes into this list.

        C++ Reference: Sample.cpp:207-211 (AddCrystalStructure)

        Args:
            structure: CrystalStructure object for this phase

        Example:
            >>> sample = Sample()
            >>> gold = CrystalStructure.fcc(4.0782)
            >>> sample.add_crystal_structure(gold)
        """
        self.crystal_structures.append(structure)

    def get_structure_list(self) -> List[CrystalStructure]:
        """
        Get list of all crystal structures.

        C++ Reference: Sample.h:117 (GetStructureList)

        Returns:
            List of CrystalStructure objects
        """
        return self.crystal_structures

    def set_sample_symmetry(self, symmetry: CrystalSymmetry) -> None:
        """
        Set crystal symmetry for the sample.

        C++ Reference: Sample.cpp:213-217 (SetSampleSymmetry)

        Args:
            symmetry: CrystalSymmetry object
        """
        self.symmetry = symmetry

    def get_sample_symmetry(self) -> Optional[CrystalSymmetry]:
        """
        Get crystal symmetry for the sample.

        C++ Reference: Sample.h:119-122 (GetSampleSymmetry)

        Returns:
            CrystalSymmetry object or None
        """
        return self.symmetry

    # =========================================================================
    # Microstructure Access
    # =========================================================================

    def get_mic(self) -> Optional[MicFile]:
        """
        Get microstructure voxel grid.

        C++ Reference: Sample.h:124-125 (GetMic)

        Returns:
            MicFile object or None if not loaded
        """
        return self.mic_file

    def num_voxels(self) -> int:
        """
        Get number of voxels in microstructure.

        Returns:
            Number of voxels, or 0 if no mic file loaded
        """
        if self.mic_file is None:
            return 0
        return len(self.mic_file.voxels)

    def __repr__(self) -> str:
        """String representation of Sample."""
        location = self.get_location()
        orientation = self.get_orientation()
        num_voxels = self.num_voxels()
        num_phases = len(self.crystal_structures)

        return (
            f"Sample(\n"
            f"  location={location},\n"
            f"  orientation_euler={orientation} rad,\n"
            f"  num_voxels={num_voxels},\n"
            f"  num_phases={num_phases}\n"
            f")"
        )
