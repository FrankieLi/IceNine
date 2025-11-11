"""
Test suite for Miller indices generation and symmetry reduction.

These tests validate the Python crystal structure implementation
against C++ ground truth data for Au FCC with Q_max = 8.0 Å⁻¹.

Test Data Sources:
    - cpp_outputs/au_fcc_all_hkl_qmax8.json (136 reflections)
    - cpp_outputs/au_fcc_unique_hkl_qmax8.json (9 unique reflections)
    - cpp_outputs/q_magnitudes.json (Q-values for test cases)

Critical Tests:
    - Test 3: Generated (h,k,l) match C++ count exactly (136)
    - Test 4: Unique reflections match C++ count exactly (9)
    - Test 5: Q-magnitudes match C++ (tolerance < 1e-10)
"""

import json
import numpy as np
import pytest
from pathlib import Path
from typing import List, Dict, Tuple

# Import modules to test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from icenine.crystal_structure import (
    CrystalStructure,
    Reflection,
    create_gold_fcc,
)
from icenine.constants import TEST_PARAMS


# Path to C++ ground truth data
CPP_OUTPUTS_DIR = Path(__file__).parent.parent / "cpp_outputs"


def load_json(filename: str) -> dict:
    """Load JSON ground truth data from C++."""
    filepath = CPP_OUTPUTS_DIR / filename
    with open(filepath, "r") as f:
        return json.load(f)


@pytest.fixture
def cpp_all_hkl() -> dict:
    """Load all Miller indices from C++ (before symmetry reduction)."""
    return load_json("au_fcc_all_hkl_qmax8.json")


@pytest.fixture
def cpp_unique_hkl() -> dict:
    """Load unique Miller indices from C++ (after symmetry reduction)."""
    return load_json("au_fcc_unique_hkl_qmax8.json")


@pytest.fixture
def cpp_q_magnitudes() -> dict:
    """Load test Q-magnitudes from C++."""
    return load_json("q_magnitudes.json")


@pytest.fixture
def au_fcc() -> CrystalStructure:
    """Create Au FCC crystal structure for testing."""
    return create_gold_fcc(a=TEST_PARAMS["a"])


class TestCrystalStructureBasics:
    """Basic tests for crystal structure properties."""

    def test_lattice_parameter(self, au_fcc):
        """Verify lattice parameter is set correctly."""
        assert np.isclose(au_fcc.lattice.a, TEST_PARAMS["a"])

    def test_space_group(self, au_fcc):
        """Verify Au FCC has correct space group (Fm-3m = 225)."""
        assert au_fcc.symmetry.space_group_number == 225

    def test_reciprocal_lattice_params(self, au_fcc):
        """Verify reciprocal lattice parameter calculation."""
        a_star, b_star, c_star = au_fcc.get_reciprocal_lattice_parameters()

        # For cubic: a* = b* = c* = 2π/a
        expected_a_star = 2 * np.pi / TEST_PARAMS["a"]

        np.testing.assert_allclose(a_star, expected_a_star, rtol=1e-10)
        np.testing.assert_allclose(b_star, expected_a_star, rtol=1e-10)
        np.testing.assert_allclose(c_star, expected_a_star, rtol=1e-10)


class TestSystematicAbsences:
    """Test FCC systematic absence rules."""

    def test_fcc_allowed_reflections(self, au_fcc):
        """Test that FCC allows h,k,l all even or all odd."""
        # All even - allowed
        assert au_fcc.passes_systematic_absences(2, 2, 2) is True
        assert au_fcc.passes_systematic_absences(4, 0, 0) is True

        # All odd - allowed
        assert au_fcc.passes_systematic_absences(1, 1, 1) is True
        assert au_fcc.passes_systematic_absences(3, 1, 1) is True

    def test_fcc_forbidden_reflections(self, au_fcc):
        """Test that FCC forbids mixed even/odd."""
        # Mixed - forbidden
        assert au_fcc.passes_systematic_absences(1, 0, 0) is False
        assert au_fcc.passes_systematic_absences(2, 1, 0) is False
        assert au_fcc.passes_systematic_absences(1, 2, 3) is False


class TestQMagnitudeCalculations:
    """
    Test 5: Q-magnitude calculations match C++ ground truth.

    Validates that |Q| = (2π/a) * sqrt(h² + k² + l²) matches C++
    with numerical tolerance < 1e-10.
    """

    def test_q_magnitudes_match_cpp(self, au_fcc, cpp_q_magnitudes):
        """
        CRITICAL TEST: Q-magnitudes must match C++ exactly.

        Validates 5 test reflections with high precision.
        """
        test_reflections = cpp_q_magnitudes["reflections"]
        failures = []

        for ref in test_reflections:
            h, k, l = ref["h"], ref["k"], ref["l"]
            cpp_q = ref["q_mag"]

            # Calculate Python Q
            py_q = au_fcc.calculate_q_magnitude(h, k, l)

            # Check match within tight tolerance
            if not np.isclose(py_q, cpp_q, atol=1e-10):
                failures.append({
                    "hkl": (h, k, l),
                    "cpp_q": cpp_q,
                    "py_q": py_q,
                    "diff": abs(cpp_q - py_q),
                })

        assert len(failures) == 0, (
            f"Q-magnitude mismatch in {len(failures)} reflections:\n"
            + "\n".join(str(f) for f in failures)
        )

    def test_q_magnitude_for_111(self, au_fcc):
        """Specific test for (111) reflection."""
        q = au_fcc.calculate_q_magnitude(1, 1, 1)

        # For Au FCC: a = 4.0782, so a* = 2π/4.0782 = 1.54043...
        # Q(111) = a* * sqrt(3) = 1.54043 * 1.732 = 2.6685...
        expected_q = (2 * np.pi / 4.0782) * np.sqrt(3)

        np.testing.assert_allclose(q, expected_q, rtol=1e-10)

    def test_q_magnitude_for_200(self, au_fcc):
        """Specific test for (200) reflection."""
        q = au_fcc.calculate_q_magnitude(2, 0, 0)

        # Q(200) = a* * 2
        expected_q = (2 * np.pi / 4.0782) * 2

        np.testing.assert_allclose(q, expected_q, rtol=1e-10)

    def test_q_vector_calculation(self, au_fcc):
        """Test Q-vector has correct magnitude."""
        h, k, l = 1, 1, 1
        q_vec = au_fcc.calculate_q_vector(h, k, l)

        # Check magnitude
        q_mag_from_vec = np.linalg.norm(q_vec)
        q_mag_direct = au_fcc.calculate_q_magnitude(h, k, l)

        np.testing.assert_allclose(q_mag_from_vec, q_mag_direct, rtol=1e-10)


class TestMillerIndexGeneration:
    """
    Test 3: Generated Miller indices match C++ exactly.

    Validates that generate_reflections() produces the same 136
    reflections as C++ for Au FCC with Q_max = 8.0 Å⁻¹.
    """

    def test_reflection_count_matches_cpp(self, au_fcc, cpp_all_hkl):
        """
        CRITICAL TEST: Must generate exactly 136 reflections like C++.
        """
        max_q = cpp_all_hkl["max_q"]
        cpp_count = cpp_all_hkl["count"]

        # Generate reflections with Python
        py_reflections = au_fcc.generate_reflections(max_q=max_q)

        assert len(py_reflections) == cpp_count, (
            f"Python generated {len(py_reflections)} reflections, "
            f"C++ generated {cpp_count}"
        )

    def test_all_cpp_reflections_present(self, au_fcc, cpp_all_hkl):
        """
        CRITICAL TEST: Every C++ reflection must be present in Python.

        Checks that all 136 (h,k,l) from C++ are generated by Python.
        """
        max_q = cpp_all_hkl["max_q"]
        cpp_reflections = cpp_all_hkl["reflections"]

        # Generate Python reflections
        py_reflections = au_fcc.generate_reflections(max_q=max_q)

        # Convert to sets of (h,k,l) tuples
        cpp_hkl_set = {(r["h"], r["k"], r["l"]) for r in cpp_reflections}
        py_hkl_set = {(r.h, r.k, r.l) for r in py_reflections}

        # Check all C++ reflections are in Python
        missing = cpp_hkl_set - py_hkl_set
        extra = py_hkl_set - cpp_hkl_set

        assert len(missing) == 0, f"Missing {len(missing)} C++ reflections: {missing}"
        assert len(extra) == 0, f"Extra {len(extra)} Python reflections: {extra}"

    def test_q_values_match_cpp(self, au_fcc, cpp_all_hkl):
        """
        Verify Q-magnitudes match for all 136 reflections.

        For each (h,k,l), check that Python Q matches C++ Q.
        """
        max_q = cpp_all_hkl["max_q"]
        cpp_reflections = cpp_all_hkl["reflections"]

        # Generate Python reflections
        py_reflections = au_fcc.generate_reflections(max_q=max_q)

        # Create lookup dict for Python: (h,k,l) -> q_mag
        py_q_dict = {(r.h, r.k, r.l): r.q_mag for r in py_reflections}

        failures = []
        for cpp_ref in cpp_reflections:
            h, k, l = cpp_ref["h"], cpp_ref["k"], cpp_ref["l"]
            cpp_q = cpp_ref["q_mag"]

            # Get Python Q
            py_q = py_q_dict.get((h, k, l))

            if py_q is None:
                failures.append(f"Missing: ({h},{k},{l})")
            elif not np.isclose(py_q, cpp_q, atol=1e-10):
                failures.append(
                    f"({h},{k},{l}): C++={cpp_q:.15f}, Py={py_q:.15f}, "
                    f"diff={abs(cpp_q - py_q):.2e}"
                )

        assert len(failures) == 0, (
            f"Q-magnitude mismatch in {len(failures)} reflections:\n"
            + "\n".join(failures[:10])  # Show first 10
        )

    def test_no_systematic_absences_in_output(self, au_fcc):
        """Verify that forbidden reflections are not generated."""
        reflections = au_fcc.generate_reflections(max_q=8.0)

        # Check that no forbidden reflections are present
        for ref in reflections:
            # All must pass systematic absence rules
            assert au_fcc.passes_systematic_absences(ref.h, ref.k, ref.l), (
                f"Forbidden reflection found: ({ref.h}, {ref.k}, {ref.l})"
            )


class TestUniqueReflections:
    """
    Test 4: Unique reflections after symmetry reduction match C++.

    Validates that get_unique_reflections() reduces 136 → 9
    reflections exactly as C++ does.
    """

    def test_unique_count_matches_cpp(self, au_fcc, cpp_all_hkl, cpp_unique_hkl):
        """
        CRITICAL TEST: Must reduce to exactly 9 unique reflections like C++.

        This is the most important test - validates the complete symmetry
        reduction pipeline matches C++ exactly.
        """
        max_q = cpp_all_hkl["max_q"]
        cpp_unique_count = cpp_unique_hkl["count"]

        # Generate all reflections
        all_reflections = au_fcc.generate_reflections(max_q=max_q)

        # Apply symmetry reduction
        unique_reflections = au_fcc.get_unique_reflections(all_reflections)

        assert len(unique_reflections) == cpp_unique_count, (
            f"Python found {len(unique_reflections)} unique reflections, "
            f"C++ found {cpp_unique_count}"
        )

    def test_unique_reflections_match_cpp(self, au_fcc, cpp_unique_hkl):
        """
        Verify that each unique reflection has a C++ match.

        Allows for different representative choices (e.g., (1,1,1) vs (-1,-1,-1)),
        but every Python unique reflection must be equivalent to some C++ one.
        """
        max_q = cpp_unique_hkl["max_q"]
        cpp_unique = cpp_unique_hkl["reflections"]

        # Generate Python unique reflections
        py_unique = au_fcc.generate_unique_reflections(max_q=max_q)

        # For each Python unique, check it matches some C++ unique
        failures = []
        for py_ref in py_unique:
            py_hkl = np.array([py_ref.h, py_ref.k, py_ref.l])

            # Check if equivalent to any C++ unique
            found_match = False
            for cpp_ref in cpp_unique:
                cpp_hkl = np.array([cpp_ref["h"], cpp_ref["k"], cpp_ref["l"]])

                if au_fcc.symmetry.vectors_equivalent(py_hkl, cpp_hkl):
                    found_match = True
                    break

            if not found_match:
                failures.append(
                    f"({py_ref.h},{py_ref.k},{py_ref.l}) has no C++ equivalent"
                )

        assert len(failures) == 0, (
            f"Unmatched Python unique reflections:\n" + "\n".join(failures)
        )

    def test_no_duplicates_in_unique_set(self, au_fcc):
        """Verify unique reflections contain no symmetry-equivalent pairs."""
        unique_reflections = au_fcc.generate_unique_reflections(max_q=8.0)

        # Check every pair is NOT equivalent
        for i, ref_i in enumerate(unique_reflections):
            for j, ref_j in enumerate(unique_reflections):
                if i >= j:
                    continue  # Skip self and already-checked pairs

                hkl_i = np.array([ref_i.h, ref_i.k, ref_i.l])
                hkl_j = np.array([ref_j.h, ref_j.k, ref_j.l])

                # Should NOT be equivalent
                assert not au_fcc.symmetry.vectors_equivalent(hkl_i, hkl_j), (
                    f"Duplicate found in unique set: "
                    f"({ref_i.h},{ref_i.k},{ref_i.l}) ≡ "
                    f"({ref_j.h},{ref_j.k},{ref_j.l})"
                )

    def test_reduction_percentage(self, au_fcc):
        """Verify reduction ratio matches expected value."""
        all_refs = au_fcc.generate_reflections(max_q=8.0)
        unique_refs = au_fcc.get_unique_reflections(all_refs)

        # Should reduce 136 → 9 (93.4% reduction)
        reduction_ratio = 1 - (len(unique_refs) / len(all_refs))

        assert len(all_refs) == 136
        assert len(unique_refs) == 9
        assert np.isclose(reduction_ratio, 0.9338, atol=0.01)


class TestReflectionDataclass:
    """Test the Reflection dataclass."""

    def test_reflection_creation(self):
        """Test Reflection object creation."""
        q_vec = np.array([1.0, 2.0, 3.0])
        ref = Reflection(h=1, k=1, l=1, q_mag=3.74, q_vec=q_vec)

        assert ref.h == 1
        assert ref.k == 1
        assert ref.l == 1
        assert ref.q_mag == 3.74
        assert ref.as_tuple() == (1, 1, 1)

    def test_reflection_repr(self):
        """Test Reflection string representation."""
        ref = Reflection(h=2, k=0, l=0, q_mag=3.0816)
        repr_str = repr(ref)

        assert "(2 0 0)" in repr_str
        assert "3.0816" in repr_str


class TestDSpacing:
    """Test d-spacing calculations."""

    def test_d_spacing_calculation(self, au_fcc):
        """Verify d = 2π/|Q|."""
        h, k, l = 1, 1, 1
        q = au_fcc.calculate_q_magnitude(h, k, l)
        d = au_fcc.get_d_spacing(h, k, l)

        expected_d = 2 * np.pi / q
        np.testing.assert_allclose(d, expected_d, rtol=1e-10)


# ==================== Integration Tests ====================

class TestIntegration:
    """End-to-end integration tests."""

    def test_complete_workflow(self, au_fcc):
        """Test complete Miller index generation → symmetry reduction workflow."""
        # Step 1: Generate all reflections
        all_refs = au_fcc.generate_reflections(max_q=8.0)
        assert len(all_refs) == 136

        # Step 2: Apply symmetry reduction
        unique_refs = au_fcc.get_unique_reflections(all_refs)
        assert len(unique_refs) == 9

        # Step 3: Verify Q-values are sorted or reasonable
        for ref in unique_refs:
            assert 0 < ref.q_mag <= 8.0

    def test_convenience_function(self):
        """Test generate_unique_reflections convenience method."""
        au = create_gold_fcc()
        unique = au.generate_unique_reflections(max_q=8.0)

        assert len(unique) == 9


# ==================== Performance Benchmarks ====================

@pytest.mark.benchmark
class TestPerformance:
    """Performance benchmarks (optional)."""

    def test_generation_performance(self, au_fcc, benchmark):
        """Benchmark reflection generation."""
        result = benchmark(au_fcc.generate_reflections, max_q=8.0)
        assert len(result) == 136

    def test_reduction_performance(self, au_fcc, benchmark):
        """Benchmark symmetry reduction."""
        all_refs = au_fcc.generate_reflections(max_q=8.0)
        result = benchmark(au_fcc.get_unique_reflections, all_refs)
        assert len(result) == 9


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
