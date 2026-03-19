"""
Crystal symmetry operations using pymatgen.

This module replaces the C++ Symmetry.h/cpp implementation with
a modern pymatgen-based approach that automatically handles all
230 space groups instead of hardcoded matrices.

Key improvements over C++:
- Automatic symmetry generation from space group
- No manual matrix entry (error-prone)
- Support for all 230 space groups (not just cubic/hex/tetragonal)
- Industry-standard crystallographic database
"""

from typing import List, Tuple, Optional
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp


class CrystalSymmetry:
    """
    Crystal symmetry operations wrapper around pymatgen.

    Replaces C++ CCubicSymmetry, CHexagonalSymmetry, CTetragonalSymmetry
    with a single unified class.

    Attributes:
        structure: pymatgen Structure object
        space_group_number: International space group number (1-230)
        point_group: Point group symbol (e.g., 'm-3m' for cubic)
        symmetry_ops: List of symmetry operations (SymmOp objects)
    """

    def __init__(self, structure: Structure):
        """
        Initialize symmetry operations for a crystal structure.

        Args:
            structure: pymatgen Structure object defining the crystal
        """
        self.structure = structure
        self.analyzer = SpacegroupAnalyzer(structure)

        # Get symmetry operations (point group for diffraction)
        self.symmetry_ops = self.analyzer.get_point_group_operations()

        # Cache symmetry information
        self.space_group_number = self.analyzer.get_space_group_number()
        self.space_group_symbol = self.analyzer.get_space_group_symbol()
        self.point_group = self.analyzer.get_point_group_symbol()

    def get_rotation_matrices(self) -> List[np.ndarray]:
        """
        Get all symmetry rotation matrices.

        Equivalent to C++: CCubicSymmetry::Get().GetOperatorList()

        Returns:
            List of 3x3 rotation matrices (numpy arrays)

        Example:
            >>> sym = CrystalSymmetry(fcc_structure)
            >>> matrices = sym.get_rotation_matrices()
            >>> len(matrices)  # 48 for cubic (24 rotations + inversion)
            48
        """
        return [op.rotation_matrix for op in self.symmetry_ops]

    def get_symmetry_ops(self) -> List[SymmOp]:
        """
        Get all symmetry operations (rotation + translation).

        Returns:
            List of pymatgen SymmOp objects
        """
        return self.symmetry_ops

    def vectors_equivalent(
        self, v1: np.ndarray, v2: np.ndarray, tolerance: float = 1e-5
    ) -> bool:
        """
        Check if two vectors are symmetry-equivalent.

        Equivalent to C++: Equivilent(oSymOps, v1, v2)

        Args:
            v1: First vector (h1, k1, l1) or (x, y, z)
            v2: Second vector (h2, k2, l2) or (x, y, z)
            tolerance: Numerical tolerance for equivalence

        Returns:
            True if v1 and v2 are related by a symmetry operation

        Example:
            >>> sym = CrystalSymmetry(fcc_structure)
            >>> sym.vectors_equivalent([1, 1, 1], [1, -1, -1])
            True
            >>> sym.vectors_equivalent([1, 0, 0], [0, 1, 0])
            True
        """
        v1_array = np.asarray(v1, dtype=float)
        v2_array = np.asarray(v2, dtype=float)

        # Try all symmetry operations
        for sym_op in self.symmetry_ops:
            rotated = sym_op.rotation_matrix @ v1_array
            if np.allclose(rotated, v2_array, atol=tolerance):
                return True

        return False

    def get_unique_vectors(
        self, vectors: List[np.ndarray], tolerance: float = 1e-5
    ) -> List[np.ndarray]:
        """
        Remove symmetry-equivalent vectors from a list.

        Equivalent to C++: GetUniqueVectors(oSymOps, oVectorList)

        Args:
            vectors: List of vectors (can be Miller indices or Cartesian)
            tolerance: Numerical tolerance for equivalence

        Returns:
            List of unique vectors (one representative per equivalence class)

        Example:
            >>> vectors = [[1,1,1], [1,-1,-1], [2,0,0], [0,2,0]]
            >>> unique = sym.get_unique_vectors(vectors)
            >>> len(unique)  # Fewer than input due to symmetry
            2
        """
        unique_vectors: List[np.ndarray] = []

        for vec in vectors:
            vec_array = np.asarray(vec, dtype=float)
            is_unique = True

            # Check against all unique vectors found so far
            for unique_vec in unique_vectors:
                if self.vectors_equivalent(vec_array, unique_vec, tolerance):
                    is_unique = False
                    break

            if is_unique:
                unique_vectors.append(vec_array)

        return unique_vectors

    def get_equivalent_hkls(self, hkl: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        """
        Generate all symmetry-equivalent Miller indices.

        Uses pymatgen's built-in function for generating equivalent reflections.

        Args:
            hkl: Miller index tuple (h, k, l)

        Returns:
            List of all equivalent (h, k, l) tuples

        Example:
            >>> sym = CrystalSymmetry(fcc_structure)
            >>> equivalents = sym.get_equivalent_hkls((1, 1, 1))
            >>> len(equivalents)  # 8 for cubic {111}
            8
        """
        from pymatgen.core.surface import get_symmetrically_equivalent_miller_indices

        equivalent_hkls = get_symmetrically_equivalent_miller_indices(
            self.structure, hkl
        )
        return list(equivalent_hkls)

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"CrystalSymmetry("
            f"space_group={self.space_group_number}, "
            f"symbol='{self.space_group_symbol}', "
            f"point_group='{self.point_group}', "
            f"n_ops={len(self.symmetry_ops)})"
        )


def create_cubic_symmetry(lattice_a: float, element: str = "Au") -> CrystalSymmetry:
    """
    Create cubic symmetry for FCC/BCC structures.

    Convenience function equivalent to C++: CCubicSymmetry::Get()

    Args:
        lattice_a: Cubic lattice parameter in Angstroms
        element: Element symbol (default: Au)

    Returns:
        CrystalSymmetry object for cubic crystal

    Example:
        >>> sym = create_cubic_symmetry(4.0782, "Au")
        >>> len(sym.get_rotation_matrices())
        48
    """
    lattice = Lattice.cubic(lattice_a)
    structure = Structure(lattice, [element], [[0, 0, 0]])
    return CrystalSymmetry(structure)


def create_fcc_symmetry(lattice_a: float, element: str = "Au") -> CrystalSymmetry:
    """
    Create FCC (face-centered cubic) symmetry.

    FCC has space group 225 (Fm-3m).

    Args:
        lattice_a: Cubic lattice parameter in Angstroms
        element: Element symbol (default: Au for gold)

    Returns:
        CrystalSymmetry object for FCC crystal

    Example:
        >>> sym = create_fcc_symmetry(4.0782)  # Au lattice parameter
        >>> sym.space_group_number
        225
    """
    from pymatgen.core.structure import Structure
    from pymatgen.core.lattice import Lattice

    # FCC lattice with face-centered atoms
    # Conventional FCC cell has 4 atoms: (0,0,0), (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5)
    lattice = Lattice.cubic(lattice_a)
    structure = Structure(
        lattice,
        [element] * 4,
        [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    )

    return CrystalSymmetry(structure)


# Convenience functions matching C++ API
def vectors_equivalent(
    symmetry: CrystalSymmetry, v1: np.ndarray, v2: np.ndarray, tolerance: float = 1e-5
) -> bool:
    """
    Standalone function matching C++ API: Equivilent(oSym, v1, v2)
    """
    return symmetry.vectors_equivalent(v1, v2, tolerance)


def get_unique_vectors(
    symmetry: CrystalSymmetry, vectors: List[np.ndarray], tolerance: float = 1e-5
) -> List[np.ndarray]:
    """
    Standalone function matching C++ API: GetUniqueVectors(oSym, vectors)
    """
    return symmetry.get_unique_vectors(vectors, tolerance)
