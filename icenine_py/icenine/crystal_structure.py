"""
Crystal structure and reciprocal lattice calculations.

This module provides Miller index generation, reciprocal lattice vectors,
and symmetry reduction for diffraction calculations.

Replaces parts of CrystalStructure.h/cpp from the C++ implementation.
"""

from typing import List, Tuple, Optional
import numpy as np
from dataclasses import dataclass
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from .symmetry import CrystalSymmetry


@dataclass
class Reflection:
    """
    A reciprocal lattice reflection (Miller indices with Q magnitude).

    Attributes:
        h, k, l: Miller indices
        q_mag: Magnitude of scattering vector |Q| in Å⁻¹
        q_vec: Cartesian scattering vector (optional)
        intensity: Structure factor intensity |F|² (optional)
    """
    h: int
    k: int
    l: int
    q_mag: float
    q_vec: Optional[np.ndarray] = None
    intensity: float = 1.0

    def __repr__(self) -> str:
        return f"({self.h} {self.k} {self.l}) Q={self.q_mag:.4f}"

    def as_tuple(self) -> Tuple[int, int, int]:
        """Return Miller indices as tuple."""
        return (self.h, self.k, self.l)


class CrystalStructure:
    """
    Crystal structure with reciprocal lattice calculations.

    This class handles:
    - Reciprocal lattice generation
    - Miller index enumeration
    - Q-vector calculations
    - Symmetry reduction

    Example:
        >>> structure = CrystalStructure.create_fcc("Au", 4.0782)
        >>> reflections = structure.generate_reflections(max_q=8.0)
        >>> unique = structure.get_unique_reflections(reflections)
    """

    def __init__(self, structure: Structure):
        """
        Initialize crystal structure.

        Args:
            structure: pymatgen Structure object
        """
        self.structure = structure
        self.lattice = structure.lattice
        self.reciprocal_lattice = structure.lattice.reciprocal_lattice

        # Create symmetry object
        self.symmetry = CrystalSymmetry(structure)

        # Cache for reflection vectors
        self._reflection_vectors: Optional[List[Reflection]] = None
        self._max_q: Optional[float] = None

    @classmethod
    def create_fcc(cls, element: str, a: float) -> "CrystalStructure":
        """
        Create FCC (face-centered cubic) crystal structure.

        Args:
            element: Element symbol (e.g., "Au")
            a: Cubic lattice parameter in Angstroms

        Returns:
            CrystalStructure object

        Example:
            >>> au = CrystalStructure.create_fcc("Au", 4.0782)
            >>> au.lattice.a
            4.0782
        """
        lattice = Lattice.cubic(a)
        # FCC conventional cell has 4 atoms at face-centered positions
        structure = Structure(
            lattice,
            [element] * 4,
            [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        )
        return cls(structure)

    @classmethod
    def create_cubic(cls, element: str, a: float) -> "CrystalStructure":
        """Alias for create_fcc."""
        return cls.create_fcc(element, a)

    def calculate_q_magnitude(self, h: int, k: int, l: int) -> float:
        """
        Calculate magnitude of reciprocal lattice vector |Q|.

        For cubic crystals: |Q| = (2π/a) * sqrt(h² + k² + l²)

        Args:
            h, k, l: Miller indices

        Returns:
            |Q| in Å⁻¹

        Example:
            >>> structure = CrystalStructure.create_fcc("Au", 4.0782)
            >>> q = structure.calculate_q_magnitude(1, 1, 1)
            >>> round(q, 6)
            2.668530
        """
        # Get reciprocal lattice vector in Cartesian coordinates
        q_vec = self.reciprocal_lattice.get_cartesian_coords([h, k, l])
        return float(np.linalg.norm(q_vec))

    def calculate_q_vector(self, h: int, k: int, l: int) -> np.ndarray:
        """
        Calculate reciprocal lattice vector Q in Cartesian coordinates.

        Args:
            h, k, l: Miller indices

        Returns:
            Q vector [Qx, Qy, Qz] in Å⁻¹
        """
        return self.reciprocal_lattice.get_cartesian_coords([h, k, l])

    def passes_systematic_absences(self, h: int, k: int, l: int) -> bool:
        """
        Check if reflection passes systematic absence rules.

        For FCC: h,k,l all even OR all odd
        For BCC: h+k+l even
        For simple cubic: no restrictions

        Args:
            h, k, l: Miller indices

        Returns:
            True if reflection is allowed (not systematically absent)

        Example:
            >>> structure = CrystalStructure.create_fcc("Au", 4.0782)
            >>> structure.passes_systematic_absences(1, 1, 1)  # All odd - allowed
            True
            >>> structure.passes_systematic_absences(1, 0, 0)  # Mixed - forbidden
            False
        """
        # Check if FCC by examining space group
        if self.symmetry.space_group_number == 225:  # Fm-3m (FCC)
            # FCC rule: h,k,l all even OR all odd
            all_even = (h % 2 == 0) and (k % 2 == 0) and (l % 2 == 0)
            all_odd = (h % 2 != 0) and (k % 2 != 0) and (l % 2 != 0)
            return all_even or all_odd

        # For simple cubic or other structures, no restrictions by default
        return True

    def generate_reflections(
        self,
        max_q: float,
        include_systematic_absences: bool = False
    ) -> List[Reflection]:
        """
        Generate all possible reflections up to max_q.

        This generates all (h,k,l) combinations within the Q sphere,
        optionally filtering systematic absences.

        Args:
            max_q: Maximum |Q| in Å⁻¹
            include_systematic_absences: If False, filter forbidden reflections

        Returns:
            List of Reflection objects

        Example:
            >>> structure = CrystalStructure.create_fcc("Au", 4.0782)
            >>> reflections = structure.generate_reflections(max_q=8.0)
            >>> len(reflections)  # Should match C++: 136 for Au FCC
            136
        """
        reflections = []

        # Estimate maximum h,k,l to check
        # For cubic: a* = 2π/a, so h_max ~ max_q * a / (2π)
        a_recip = self.reciprocal_lattice.a  # pymatgen already includes 2π
        max_index = int(np.ceil(max_q / a_recip)) + 1

        # Generate all combinations
        for h in range(-max_index, max_index + 1):
            for k in range(-max_index, max_index + 1):
                for l in range(-max_index, max_index + 1):
                    # Skip (0,0,0)
                    if h == 0 and k == 0 and l == 0:
                        continue

                    # Check systematic absences
                    if not include_systematic_absences:
                        if not self.passes_systematic_absences(h, k, l):
                            continue

                    # Calculate |Q|
                    q_mag = self.calculate_q_magnitude(h, k, l)

                    # Keep if within max_q
                    if q_mag <= max_q:
                        q_vec = self.calculate_q_vector(h, k, l)

                        # Calculate structure factor intensity
                        # C++: oRefRecipVector.fIntensity = CalculateIntensity( oRefRecipVector.v )
                        intensity = self.calculate_intensity(q_vec)

                        reflection = Reflection(h, k, l, q_mag, q_vec, intensity)
                        reflections.append(reflection)

        return reflections

    def get_unique_reflections(
        self,
        reflections: List[Reflection],
        tolerance: float = 1e-5
    ) -> List[Reflection]:
        """
        Remove symmetry-equivalent reflections.

        Uses crystal symmetry to reduce the reflection list to unique
        representatives (one per symmetry-equivalent family).

        Args:
            reflections: List of all reflections
            tolerance: Numerical tolerance for equivalence

        Returns:
            List of unique reflections (one per family)

        Example:
            >>> structure = CrystalStructure.create_fcc("Au", 4.0782)
            >>> all_refs = structure.generate_reflections(max_q=8.0)
            >>> unique_refs = structure.get_unique_reflections(all_refs)
            >>> len(all_refs)  # C++ ground truth
            136
            >>> len(unique_refs)  # C++ ground truth
            9
        """
        unique_reflections: List[Reflection] = []

        for reflection in reflections:
            # Check if this reflection is equivalent to any unique one
            hkl = np.array([reflection.h, reflection.k, reflection.l], dtype=float)
            is_unique = True

            for unique_ref in unique_reflections:
                unique_hkl = np.array(
                    [unique_ref.h, unique_ref.k, unique_ref.l], dtype=float
                )

                # Use symmetry to check equivalence
                if self.symmetry.vectors_equivalent(hkl, unique_hkl, tolerance):
                    is_unique = False
                    break

            if is_unique:
                unique_reflections.append(reflection)

        return unique_reflections

    def generate_unique_reflections(
        self,
        max_q: float,
        include_systematic_absences: bool = False,
        tolerance: float = 1e-5
    ) -> List[Reflection]:
        """
        Generate unique reflections in one step.

        Convenience method that combines generate_reflections() and
        get_unique_reflections().

        Args:
            max_q: Maximum |Q| in Å⁻¹
            include_systematic_absences: Include forbidden reflections
            tolerance: Numerical tolerance for symmetry equivalence

        Returns:
            List of unique Reflection objects

        Example:
            >>> structure = CrystalStructure.create_fcc("Au", 4.0782)
            >>> unique = structure.generate_unique_reflections(max_q=8.0)
            >>> len(unique)
            9
        """
        all_reflections = self.generate_reflections(
            max_q, include_systematic_absences
        )
        return self.get_unique_reflections(all_reflections, tolerance)

    def calculate_intensity(self, q_vec: np.ndarray) -> float:
        """
        Calculate structure factor intensity |F|² for a reciprocal vector.

        This implements the kinematical diffraction structure factor:
            F = Σ_j f_j * exp(i * Q · r_j)
            I = |F|² = Real(F)² + Imag(F)²

        where:
            f_j = atomic scattering factor (approximated as atomic number)
            Q = reciprocal lattice vector
            r_j = atom position in Cartesian coordinates

        C++ Reference: CrystalStructure.cpp:404-425 CalculateIntensity()

        Args:
            q_vec: Reciprocal lattice vector Q in Cartesian coordinates (Å⁻¹)

        Returns:
            Structure factor intensity |F|²

        Example:
            >>> structure = CrystalStructure.create_fcc("Au", 4.0782)
            >>> q_vec = structure.calculate_q_vector(1, 1, 1)
            >>> intensity = structure.calculate_intensity(q_vec)
        """
        s_real = 0.0
        s_imaginary = 0.0

        # C++: for(Size_Type i = 0; i < oTranslationVector.size(); i ++)
        for site in self.structure.sites:
            # Get atom position in Cartesian coordinates
            # C++: SVector3 oAtomPosition = oTranslationVector[i].v.m_fX * oPrimitiveVector[0] + ...
            atom_position = site.coords  # pymatgen already provides Cartesian coords

            # Calculate Q · r
            # C++: Float fRDotK = Dot( oReciprocalVector, oAtomPosition )
            r_dot_k = np.dot(q_vec, atom_position)

            # Atomic scattering factor (approximated as atomic number Z)
            # C++: oTranslationVector[i].fEffectiveZ
            # For more accuracy, could use proper scattering factors from tables
            effective_z = float(site.specie.Z)

            # Structure factor: F = Σ f_j * e^(i*k·r) = Σ f_j * (cos(k·r) + i*sin(k·r))
            # C++: fSReal += oTranslationVector[i].fEffectiveZ * cos( fRDotK )
            # C++: fSImaginary += oTranslationVector[i].fEffectiveZ * sin( fRDotK )
            s_real += effective_z * np.cos(r_dot_k)
            s_imaginary += effective_z * np.sin(r_dot_k)

        # Intensity = |F|² = Real² + Imag²
        # C++: return ( fSReal * fSReal + fSImaginary * fSImaginary )
        intensity = s_real * s_real + s_imaginary * s_imaginary

        return intensity

    def get_d_spacing(self, h: int, k: int, l: int) -> float:
        """
        Calculate d-spacing for a reflection.

        d = 2π / |Q|  (note: pymatgen uses 2π convention)

        Args:
            h, k, l: Miller indices

        Returns:
            d-spacing in Angstroms
        """
        q_mag = self.calculate_q_magnitude(h, k, l)
        return 2 * np.pi / q_mag if q_mag > 0 else np.inf

    def get_reciprocal_lattice_parameters(self) -> Tuple[float, float, float]:
        """
        Get reciprocal lattice parameters (a*, b*, c*).

        Returns:
            (a_star, b_star, c_star) in Å⁻¹
        """
        return (
            self.reciprocal_lattice.a,
            self.reciprocal_lattice.b,
            self.reciprocal_lattice.c,
        )

    def set_reflection_limits(
        self,
        max_h: int,
        max_k: int,
        max_l: int,
        max_q: float,
        min_intensity_fraction: float = 0.0
    ) -> None:
        """
        Set reflection generation limits and cache reflections.

        Compatible with C++ CUnitCell::SetReflectionVectorLimits().

        Args:
            max_h, max_k, max_l: Maximum Miller indices (currently unused, using max_q instead)
            max_q: Maximum scattering vector magnitude in Å⁻¹
            min_intensity_fraction: Minimum intensity fraction to include (currently unused)
        """
        self._max_q = max_q
        # Generate and cache reflections
        self._reflection_vectors = self.generate_reflections(max_q)

    def get_reflection_vectors(self) -> List[Reflection]:
        """
        Get cached reflection vectors.

        Compatible with C++ CUnitCell::GetReflectionVectorList().
        Must call set_reflection_limits() first.

        Returns:
            List of Reflection objects with q_vec populated

        Raises:
            RuntimeError: If set_reflection_limits() not called
        """
        if self._reflection_vectors is None:
            raise RuntimeError(
                "Reflection vectors not initialized. "
                "Call set_reflection_limits() first."
            )
        return self._reflection_vectors

    def __repr__(self) -> str:
        """String representation."""
        formula = self.structure.composition.reduced_formula
        a = self.lattice.a
        sg = self.symmetry.space_group_number
        return f"CrystalStructure({formula}, a={a:.4f}Å, SG={sg})"


# Convenience functions for common crystal types

def create_fcc_structure(element: str, a: float) -> CrystalStructure:
    """
    Create FCC crystal structure.

    Args:
        element: Element symbol
        a: Lattice parameter in Angstroms

    Returns:
        CrystalStructure object
    """
    return CrystalStructure.create_fcc(element, a)


def create_gold_fcc(a: float = 4.0782) -> CrystalStructure:
    """
    Create Gold FCC structure with default lattice parameter.

    Default value from ConfigFiles/ReconstructTest.config

    Args:
        a: Lattice parameter (default: 4.0782 Å)

    Returns:
        CrystalStructure for Au FCC
    """
    return CrystalStructure.create_fcc("Au", a)
