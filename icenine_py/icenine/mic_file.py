"""
MIC file I/O for microstructure voxel data.

This module provides read/write functionality for .mic files, the standard
format for storing 3D microstructure data with crystal orientations.

File format: Text-based format storing voxel positions, orientations (Euler angles),
and quality metrics. Maintains backward compatibility with C++ implementation in
XDM++/libXDM/MicIO.h.

Key features:
- Read/write C++ compatible .mic files
- Automatic Euler angle ↔ rotation matrix conversion
- PyTorch tensor storage for differentiability
- Batch operations support

Based on C++ implementation: XDM++/libXDM/MicIO.h
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Set, Dict
import numpy as np
import torch
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree  # C-optimized KDTree for 10-100x speedup
from pathlib import Path

# Import geometry utilities (moved to geometry.py for reusability)
from .geometry import euler_to_matrix, matrix_to_euler, euler_to_matrix_torch


class ReconstructionState:
    """
    Voxel state machine for BFS reconstruction.

    C++ Reference: ReconstructionStrategies.h:257-263 MultiStagedDetails
    """

    NOT_VISITED = -1
    VISITED = 0
    FITTED = 1
    REFIT = 2


@dataclass
class Voxel:
    """
    Single voxel in a microstructure.

    Represents a spatial element with crystal orientation and quality metrics.

    Attributes:
        position: 3D coordinates (x, y, z) in sample frame [m]
        orientation: 3x3 rotation matrix (crystal → sample frame)
        side_length: Voxel edge length [m]
        generation: Refinement level (0 = coarsest)
        phase: Material phase ID (0 = unfitted, 1+ = fitted)
        confidence: Fitting quality [0, 1]
        cost: Fitting cost function value
        overlap_ratio: Fraction of simulated peaks matching experimental
        points_up: Triangle orientation (triangular mesh only)
        id: Unique voxel identifier
        deformation: 3x3 deformation tensor (optional, for strain)
        reconstruction_id: BFS state (-1=NOT_VISITED, 0=VISITED, 1=FITTED, 2=REFIT)
    """

    position: np.ndarray  # (3,) - x, y, z
    orientation: np.ndarray  # (3, 3) - rotation matrix
    side_length: float = 0.0
    generation: int = 0
    phase: int = 1
    confidence: float = 0.0
    cost: float = 0.0
    overlap_ratio: float = 0.0
    points_up: bool = True
    id: int = -1
    deformation: Optional[np.ndarray] = None  # (3, 3) optional
    reconstruction_id: int = -1  # BFS state, see ReconstructionState

    def __post_init__(self):
        """Validate voxel data after initialization."""
        self.position = np.asarray(self.position, dtype=np.float32)
        self.orientation = np.asarray(self.orientation, dtype=np.float32)
        assert self.position.shape == (3,), f"Position must be (3,), got {self.position.shape}"
        assert self.orientation.shape == (
            3,
            3,
        ), f"Orientation must be (3,3), got {self.orientation.shape}"
        if self.deformation is not None:
            self.deformation = np.asarray(self.deformation, dtype=np.float32)
            assert self.deformation.shape == (
                3,
                3,
            ), f"Deformation must be (3,3), got {self.deformation.shape}"


class MicFile:
    """
    MIC file reader/writer with PyTorch tensor storage.

    Stores microstructure voxel data in batched PyTorch tensors for efficient
    computation and automatic differentiation. Maintains backward compatibility
    with C++ .mic file format.

    Attributes:
        voxels: List of Voxel objects
        initial_side_length: Initial voxel side length (generation 0) [m]
        positions: Batched positions, shape (N, 3)
        orientations: Batched orientation matrices, shape (N, 3, 3)
        confidence: Batched confidence values, shape (N,)
        phase: Batched phase IDs, shape (N,)
    """

    def __init__(
        self,
        voxels: Optional[List[Voxel]] = None,
        initial_side_length: float = 1.2e-2,
    ):
        """
        Initialize MIC file data structure.

        Args:
            voxels: List of Voxel objects (optional)
            initial_side_length: Initial voxel side length in meters
        """
        self.initial_side_length = initial_side_length
        self.voxels = voxels if voxels is not None else []

        # Create batched PyTorch tensors from voxels
        self._create_tensors()

    def _create_tensors(self):
        """Convert voxel list to batched PyTorch tensors."""
        if len(self.voxels) == 0:
            # Empty tensors
            self.positions = torch.empty(0, 3, dtype=torch.float32)
            self.orientations = torch.empty(0, 3, 3, dtype=torch.float32)
            self.confidence = torch.empty(0, dtype=torch.float32)
            self.cost = torch.empty(0, dtype=torch.float32)
            self.overlap_ratio = torch.empty(0, dtype=torch.float32)
            self.phase = torch.empty(0, dtype=torch.int32)
            self.generation = torch.empty(0, dtype=torch.int32)
        else:
            # Stack voxel data into tensors
            self.positions = torch.from_numpy(
                np.stack([v.position for v in self.voxels])
            ).float()
            self.orientations = torch.from_numpy(
                np.stack([v.orientation for v in self.voxels])
            ).float()
            self.confidence = torch.tensor([v.confidence for v in self.voxels]).float()
            self.cost = torch.tensor([v.cost for v in self.voxels]).float()
            self.overlap_ratio = torch.tensor(
                [v.overlap_ratio for v in self.voxels]
            ).float()
            self.phase = torch.tensor([v.phase for v in self.voxels], dtype=torch.int32)
            self.generation = torch.tensor(
                [v.generation for v in self.voxels], dtype=torch.int32
            )

    @staticmethod
    def _parse_header(header_line: str) -> float:
        """
        Parse header line to extract initial side length.

        Args:
            header_line: First line of .mic file

        Returns:
            Initial side length in meters

        Raises:
            ValueError: If header format is invalid
        """
        tokens = header_line.strip().split()
        if len(tokens) != 1:
            raise ValueError(
                f"Invalid MIC file format (line 1): expected 1 token, got {len(tokens)}"
            )
        return float(tokens[0])

    @staticmethod
    def _parse_deformation_tensor(tokens: List[str]) -> np.ndarray:
        """
        Parse symmetric 3x3 deformation tensor from 6 values.

        The deformation tensor is stored as [m00, m11, m22, m01, m12, m02]
        and reconstructed as a symmetric matrix.

        Args:
            tokens: List of 6 string values representing tensor components

        Returns:
            3x3 symmetric deformation tensor
        """
        # Convert to floats
        vals = [float(t) for t in tokens]

        # Build symmetric matrix: D[i,j] = D[j,i]
        D = np.array(
            [[vals[0], vals[3], vals[5]], [vals[3], vals[1], vals[4]], [vals[5], vals[4], vals[2]]],
            dtype=np.float32,
        )

        return D

    @classmethod
    def _parse_voxel_from_tokens(
        cls, tokens: List[str], line_num: int, initial_side_length: float
    ) -> Voxel:
        """
        Parse a single voxel from line tokens.

        Args:
            tokens: List of string tokens from one line
            line_num: Line number (for error reporting)
            initial_side_length: Initial voxel side length from header

        Returns:
            Parsed Voxel object

        Raises:
            ValueError: If token count is invalid
        """
        if len(tokens) < 9:
            raise ValueError(
                f"Invalid MIC file format (line {line_num}): "
                f"expected >= 9 columns, got {len(tokens)}"
            )

        # Parse required fields using tuple unpacking (more Pythonic)
        try:
            x, y, z = float(tokens[0]), float(tokens[1]), float(tokens[2])
            direction = int(tokens[3])  # 1 = UP, 2 = DOWN
            generation = int(tokens[4])
            phase = int(tokens[5])
            phi1_deg, Phi_deg, phi2_deg = float(tokens[6]), float(tokens[7]), float(tokens[8])

            # Parse optional fields with defaults
            confidence = float(tokens[9]) if len(tokens) > 9 else 0.0
            cost = float(tokens[10]) if len(tokens) > 10 else 0.0
            overlap_ratio = float(tokens[11]) if len(tokens) > 11 else 0.0
        except ValueError as e:
            # Show the offending line when parsing fails
            line_text = ' '.join(tokens)
            raise ValueError(
                f"Failed to parse MIC file (line {line_num}):\n"
                f"  Error: {e}\n"
                f"  Line: {line_text}\n"
                f"  Tokens: {tokens}"
            )

        # Parse deformation tensor if present (19 total columns)
        deformation = None
        if len(tokens) == 19:
            deformation = cls._parse_deformation_tensor(tokens[13:19])

        # Convert Euler angles (degrees) to rotation matrix
        # C++ uses BuildActiveEulerMatrix(phi1, Phi, phi2) with ZXZ convention
        orientation_matrix = euler_to_matrix(phi1_deg, Phi_deg, phi2_deg)

        # Calculate side length from generation: side = initial / 2^gen
        side_length = initial_side_length / (2**generation)

        return Voxel(
            position=np.array([x, y, z], dtype=np.float32),
            orientation=orientation_matrix,
            side_length=side_length,
            generation=generation,
            phase=phase,
            confidence=confidence,
            cost=cost,
            overlap_ratio=overlap_ratio,
            points_up=(direction == 1),
            deformation=deformation,
        )

    @classmethod
    def read(cls, filename: str) -> "MicFile":
        """
        Read .mic file (triangular mesh format).

        C++ compatible reader for triangular mesh .mic files. Automatically
        converts Euler angles (degrees) to rotation matrices.

        File format (from C++ MicIO.h lines 263-356):
            Line 1: <SideLength>
            Line 2+: <x> <y> <z> <dir> <gen> <phase> <φ1°> <Φ°> <φ2°> <conf> [<cost> <overlap> <time> <deformation...>]

        Args:
            filename: Path to .mic file

        Returns:
            MicFile object with voxel data

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid

        Example:
            >>> mic = MicFile.read("sample.mic")
            >>> print(f"Loaded {len(mic.voxels)} voxels")
            >>> print(f"Side length: {mic.initial_side_length}")
        """
        path = Path(filename)
        if not path.exists():
            raise FileNotFoundError(f"MIC file not found: {filename}")

        with open(filename, "r") as f:
            lines = f.readlines()

        if len(lines) == 0:
            raise ValueError(f"Empty MIC file: {filename}")

        # Parse header (line 1)
        initial_side_length = cls._parse_header(lines[0])

        # Parse voxels (lines 2+)
        voxels = []
        for line_num, line in enumerate(lines[1:], start=2):
            tokens = line.strip().split()

            # Skip empty lines
            if len(tokens) == 0:
                continue

            voxel = cls._parse_voxel_from_tokens(tokens, line_num, initial_side_length)
            voxels.append(voxel)

        return cls(voxels=voxels, initial_side_length=initial_side_length)

    def write(self, filename: str):
        """
        Write .mic file (triangular mesh format).

        C++ compatible writer for triangular mesh .mic files. Automatically
        converts rotation matrices to Euler angles (degrees).

        Writes format matching C++ MicIO.h lines 363-429.

        Args:
            filename: Path to output .mic file

        Example:
            >>> mic = MicFile.read("input.mic")
            >>> mic.write("output.mic")
        """
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)

        with open(filename, "w") as f:
            # Line 1: initial side length
            f.write(f"{self.initial_side_length:13.7E}\n")

            # Lines 2+: voxel data
            for voxel in self.voxels:
                # Get left-most vertex position for triangular mesh
                # (In triangular mesh, position is the left vertex)
                x, y, z = voxel.position

                # Direction: 1 = UP, 2 = DOWN
                direction = 1 if voxel.points_up else 2

                # Convert rotation matrix to Euler angles (degrees)
                phi1_deg, Phi_deg, phi2_deg = matrix_to_euler(voxel.orientation)

                # Write core fields (columns 1-10)
                f.write(
                    f" {x:13.7E} {y:13.7E} {z:13.7E} "
                    f"{direction:11d} {voxel.generation:11d} {voxel.phase:11d} "
                    f"{phi1_deg:11.4f} {Phi_deg:14.6f} {phi2_deg:14.6f} "
                    f"{voxel.confidence:14.7f}"
                )

                # Write extended fields (cost, overlap, fitting time, deformation)
                f.write(f"    {voxel.cost:13.6E}")
                f.write(f" {voxel.overlap_ratio:13.6E}")
                f.write(f" {0.0:13.6E}")  # fitting_time placeholder

                # Write deformation tensor (symmetric 3x3)
                if voxel.deformation is not None:
                    D = voxel.deformation
                    f.write(f" {D[0,0]:13.6E}")
                    f.write(f" {D[1,1]:13.6E}")
                    f.write(f" {D[2,2]:13.6E}")
                    f.write(f" {D[1,0]:13.6E}")  # off-diagonal
                    f.write(f" {D[1,2]:13.6E}")
                    f.write(f" {D[2,0]:13.6E}")
                else:
                    # Identity matrix
                    f.write(f" {1.0:13.6E}")
                    f.write(f" {1.0:13.6E}")
                    f.write(f" {1.0:13.6E}")
                    f.write(f" {0.0:13.6E}")
                    f.write(f" {0.0:13.6E}")
                    f.write(f" {0.0:13.6E}")

                f.write("\n")

    def save_torch(self, filename: str):
        """
        Save to PyTorch native format (fast, preserves gradients).

        Args:
            filename: Path to output .pt file
        """
        torch.save(
            {
                "positions": self.positions,
                "orientations": self.orientations,
                "confidence": self.confidence,
                "cost": self.cost,
                "overlap_ratio": self.overlap_ratio,
                "phase": self.phase,
                "generation": self.generation,
                "initial_side_length": self.initial_side_length,
            },
            filename,
        )

    @classmethod
    def load_torch(cls, filename: str) -> "MicFile":
        """
        Load from PyTorch native format.

        Args:
            filename: Path to .pt file

        Returns:
            MicFile object
        """
        data = torch.load(filename)

        # Reconstruct voxels from tensors
        n_voxels = data["positions"].shape[0]
        voxels = []

        for i in range(n_voxels):
            voxel = Voxel(
                position=data["positions"][i].numpy(),
                orientation=data["orientations"][i].numpy(),
                side_length=data["initial_side_length"] / (2 ** data["generation"][i].item()),
                generation=data["generation"][i].item(),
                phase=data["phase"][i].item(),
                confidence=data["confidence"][i].item(),
                cost=data["cost"][i].item(),
                overlap_ratio=data["overlap_ratio"][i].item(),
            )
            voxels.append(voxel)

        return cls(voxels=voxels, initial_side_length=data["initial_side_length"])

    # ===========================================================================
    # Spatial Indexing and Neighbor Queries
    # ===========================================================================

    def build_spatial_index(self) -> None:
        """
        Build KDTree for fast spatial queries.

        This method constructs a KD-tree from voxel positions for efficient
        neighbor finding. The tree is cached and automatically rebuilt if
        positions change.

        Example:
            >>> mic = MicFile.read("sample.mic")
            >>> mic.build_spatial_index()
            >>> neighbors = mic.get_neighbors(0, radius=0.02)
        """
        if len(self.voxels) == 0:
            self._kdtree = None
            self._built_index = False
            return

        # Build cKDTree from positions (C-optimized for speed)
        positions_np = self.positions.numpy() if isinstance(self.positions, torch.Tensor) else self.positions
        self._kdtree = cKDTree(positions_np)
        self._built_index = True

    def _ensure_spatial_index(self) -> None:
        """Ensure spatial index is built before queries."""
        if not hasattr(self, "_built_index") or not self._built_index:
            self.build_spatial_index()

    def get_neighbors(
        self, voxel_idx: int, radius: float, max_neighbors: Optional[int] = None
    ) -> List[int]:
        """
        Find neighboring voxels within radius.

        Compatible with C++ MicGrid::GetNeighbors() but using radius-based query.

        Args:
            voxel_idx: Index of query voxel
            radius: Search radius in meters
            max_neighbors: Maximum number of neighbors to return (optional)

        Returns:
            List of voxel indices within radius (excluding query voxel itself)

        Example:
            >>> mic = MicFile.read("sample.mic")
            >>> # Find all neighbors within 2cm
            >>> neighbors = mic.get_neighbors(0, radius=0.02)
            >>> print(f"Found {len(neighbors)} neighbors")
        """
        self._ensure_spatial_index()

        if voxel_idx < 0 or voxel_idx >= len(self.voxels):
            raise IndexError(f"Voxel index {voxel_idx} out of range [0, {len(self.voxels)})")

        query_pos = self.positions[voxel_idx].numpy()

        # Query KDTree for neighbors within radius
        indices = self._kdtree.query_ball_point(query_pos, r=radius)

        # Remove query voxel itself
        neighbors = [i for i in indices if i != voxel_idx]

        # Limit number of neighbors if requested
        if max_neighbors is not None and len(neighbors) > max_neighbors:
            # Keep closest neighbors
            distances = [np.linalg.norm(self.positions[i].numpy() - query_pos) for i in neighbors]
            sorted_indices = np.argsort(distances)
            neighbors = [neighbors[i] for i in sorted_indices[:max_neighbors]]

        return neighbors

    def get_k_nearest_neighbors(self, voxel_idx: int, k: int) -> Tuple[List[int], List[float]]:
        """
        Find k nearest neighboring voxels.

        Args:
            voxel_idx: Index of query voxel
            k: Number of neighbors to find

        Returns:
            Tuple of (neighbor_indices, distances) where distances are in meters

        Example:
            >>> mic = MicFile.read("sample.mic")
            >>> # Find 6 nearest neighbors
            >>> neighbors, distances = mic.get_k_nearest_neighbors(0, k=6)
            >>> print(f"Nearest neighbor at {distances[0]:.6f} m")
        """
        self._ensure_spatial_index()

        if voxel_idx < 0 or voxel_idx >= len(self.voxels):
            raise IndexError(f"Voxel index {voxel_idx} out of range [0, {len(self.voxels)})")

        query_pos = self.positions[voxel_idx].numpy()

        # Query for k+1 neighbors (including the query point itself)
        distances, indices = self._kdtree.query(query_pos, k=k + 1)

        # Handle scalar return for k=1
        if k == 1:
            distances = np.array([distances])
            indices = np.array([indices])

        # Filter out query voxel explicitly by index value
        neighbor_indices = []
        neighbor_distances = []
        for i, d in zip(indices, distances):
            if i != voxel_idx:
                neighbor_indices.append(int(i))
                neighbor_distances.append(float(d))

        return neighbor_indices, neighbor_distances

    def query_region(self, center: np.ndarray, radius: float) -> List[int]:
        """
        Find all voxels within radius of a point.

        Args:
            center: 3D point (x, y, z) in meters
            radius: Search radius in meters

        Returns:
            List of voxel indices in region

        Example:
            >>> mic = MicFile.read("sample.mic")
            >>> # Find all voxels near origin
            >>> indices = mic.query_region(np.array([0, 0, 0]), radius=0.05)
        """
        self._ensure_spatial_index()

        center = np.asarray(center, dtype=np.float32)
        if center.shape != (3,):
            raise ValueError(f"Center must be (3,) array, got {center.shape}")

        indices = self._kdtree.query_ball_point(center, r=radius)
        return indices

    def get_boundary_voxels(
        self, fitted_phase: Optional[int] = None, unfitted_phase: int = 0
    ) -> List[int]:
        """
        Find boundary voxels (fitted voxels with unfitted neighbors).

        This is essential for reconstruction algorithms that propagate
        orientations from fitted to unfitted voxels.

        Args:
            fitted_phase: Phase ID for fitted voxels (None = any phase > 0)
            unfitted_phase: Phase ID for unfitted voxels (default: 0)

        Returns:
            List of voxel indices on boundary

        Example:
            >>> mic = MicFile.read("partial_reconstruction.mic")
            >>> boundary = mic.get_boundary_voxels()
            >>> print(f"Boundary has {len(boundary)} voxels ready to propagate")
        """
        self._ensure_spatial_index()

        boundary_indices = []

        for idx, voxel in enumerate(self.voxels):
            # Skip unfitted voxels
            if voxel.phase == unfitted_phase:
                continue

            # Check if fitted_phase filter applies
            if fitted_phase is not None and voxel.phase != fitted_phase:
                continue

            # Find neighbors within 2 * side_length
            search_radius = 2.0 * voxel.side_length
            neighbors = self.get_neighbors(idx, radius=search_radius)

            # Check if any neighbor is unfitted
            has_unfitted_neighbor = any(
                self.voxels[n].phase == unfitted_phase for n in neighbors
            )

            if has_unfitted_neighbor:
                boundary_indices.append(idx)

        return boundary_indices

    def is_boundary_voxel(
        self, voxel_idx: int, radius: Optional[float] = None
    ) -> bool:
        """
        Check if voxel is on boundary (has unfitted neighbors).

        Args:
            voxel_idx: Index of voxel to check
            radius: Search radius (default: 2 * side_length)

        Returns:
            True if voxel is on boundary

        Example:
            >>> mic = MicFile.read("sample.mic")
            >>> if mic.is_boundary_voxel(10):
            ...     print("Voxel 10 is on the boundary")
        """
        if voxel_idx < 0 or voxel_idx >= len(self.voxels):
            raise IndexError(f"Voxel index {voxel_idx} out of range")

        voxel = self.voxels[voxel_idx]

        # Unfitted voxels are not on boundary
        if voxel.phase == 0:
            return False

        # Determine search radius
        if radius is None:
            radius = 2.0 * voxel.side_length

        # Find neighbors
        neighbors = self.get_neighbors(voxel_idx, radius=radius)

        # Check if any neighbor is unfitted (phase == 0)
        return any(self.voxels[n].phase == 0 for n in neighbors)

    def __len__(self) -> int:
        """Return number of voxels."""
        return len(self.voxels)

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"MicFile(n_voxels={len(self.voxels)}, "
            f"side_length={self.initial_side_length:.6e})"
        )


# =============================================================================
# Euler Angle Conversion Utilities
# =============================================================================
# NOTE: These functions have been moved to geometry.py for reusability.
# They are imported at the top of this file and re-exported for backward
# compatibility with existing code that imports from mic_file.
#
# See icenine/geometry.py for implementations and documentation.
# =============================================================================
