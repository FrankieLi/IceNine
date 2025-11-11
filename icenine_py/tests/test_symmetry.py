"""
Test suite for crystal symmetry operations.

These tests validate the Python/pymatgen symmetry implementation
against C++ ground truth data generated from the original IceNine code.

Test Data Source: cpp_outputs/cubic_symmetry_24ops.json
                  cpp_outputs/vector_equivalence_tests.json
"""

import json
import numpy as np
import pytest
from pathlib import Path

# Import modules to test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from icenine.symmetry import (
    CrystalSymmetry,
    create_fcc_symmetry,
    vectors_equivalent,
    get_unique_vectors,
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
def cpp_cubic_symmetry() -> dict:
    """Load C++ cubic symmetry matrices."""
    return load_json("cubic_symmetry_24ops.json")


@pytest.fixture
def cpp_vector_tests() -> dict:
    """Load C++ vector equivalence test cases."""
    return load_json("vector_equivalence_tests.json")


@pytest.fixture
def au_fcc_symmetry() -> CrystalSymmetry:
    """Create Au FCC symmetry object for testing."""
    return create_fcc_symmetry(
        lattice_a=TEST_PARAMS["a"], element=TEST_PARAMS["element"]
    )


class TestCubicSymmetryMatrices:
    """
    Test 1: Verify cubic symmetry operators match C++ implementation.

    C++ generates 24 rotation matrices for cubic point group m-3m.
    pymatgen may generate 48 (including inversion), so we need to
    verify that all 24 C++ matrices are present in the Python set.
    """

    def test_symmetry_operator_count(self, au_fcc_symmetry, cpp_cubic_symmetry):
        """Verify we have at least 24 symmetry operations."""
        py_matrices = au_fcc_symmetry.get_rotation_matrices()
        cpp_count = cpp_cubic_symmetry["count"]

        assert len(py_matrices) >= cpp_count, (
            f"Python has {len(py_matrices)} ops, "
            f"expected at least {cpp_count} from C++"
        )

    def test_all_cpp_matrices_present(self, au_fcc_symmetry, cpp_cubic_symmetry):
        """
        CRITICAL TEST: All 24 C++ matrices must be present in Python set.

        This test allows for different orderings and additional operations
        (like inversion in pymatgen), but ensures every C++ matrix has a match.
        """
        py_matrices = au_fcc_symmetry.get_rotation_matrices()
        cpp_matrices_data = cpp_cubic_symmetry["matrices"]

        # Convert C++ matrices to numpy arrays
        cpp_matrices = [np.array(m, dtype=float) for m in cpp_matrices_data]

        # Check each C++ matrix has a match in Python
        missing_matrices = []
        for i, cpp_mat in enumerate(cpp_matrices):
            found_match = False
            for py_mat in py_matrices:
                if np.allclose(cpp_mat, py_mat, atol=1e-10):
                    found_match = True
                    break

            if not found_match:
                missing_matrices.append((i, cpp_mat))

        assert len(missing_matrices) == 0, (
            f"Missing {len(missing_matrices)} C++ matrices in Python output:\n"
            f"Indices: {[i for i, _ in missing_matrices]}"
        )

    def test_matrix_properties(self, au_fcc_symmetry):
        """Verify symmetry matrices have correct mathematical properties."""
        matrices = au_fcc_symmetry.get_rotation_matrices()

        for i, mat in enumerate(matrices):
            # Rotation matrices must be orthogonal: R^T @ R = I
            product = mat.T @ mat
            identity = np.eye(3)
            assert np.allclose(product, identity, atol=1e-10), (
                f"Matrix {i} is not orthogonal:\n{mat}"
            )

            # Determinant must be ±1
            det = np.linalg.det(mat)
            assert np.isclose(abs(det), 1.0, atol=1e-10), (
                f"Matrix {i} has invalid determinant: {det}"
            )

    def test_identity_operator_present(self, au_fcc_symmetry):
        """Identity matrix must be present in symmetry operators."""
        matrices = au_fcc_symmetry.get_rotation_matrices()
        identity = np.eye(3)

        found_identity = any(np.allclose(m, identity, atol=1e-10) for m in matrices)
        assert found_identity, "Identity operator not found in symmetry operations"


class TestVectorEquivalence:
    """
    Test 2: Verify vector equivalence checking matches C++.

    Uses test cases from C++ where pairs of vectors are checked
    for symmetry equivalence under cubic symmetry.
    """

    def test_all_cpp_equivalence_cases(self, au_fcc_symmetry, cpp_vector_tests):
        """
        CRITICAL TEST: All C++ vector equivalence tests must match.

        Compares Python results against C++ ground truth for:
        - Equivalent vectors (should return True)
        - Non-equivalent vectors (should return False)
        """
        test_cases = cpp_vector_tests["test_cases"]
        failures = []

        for i, case in enumerate(test_cases):
            v1 = np.array(case["v1"], dtype=float)
            v2 = np.array(case["v2"], dtype=float)
            cpp_expected = case["expected"]
            cpp_computed = case["computed"]

            # Verify C++ computed matches expected (sanity check)
            assert cpp_expected == cpp_computed, (
                f"C++ test case {i} has inconsistent expected vs computed"
            )

            # Compare Python result to C++ ground truth
            py_result = au_fcc_symmetry.vectors_equivalent(v1, v2)

            if py_result != cpp_expected:
                failures.append(
                    {
                        "case_index": i,
                        "v1": v1.tolist(),
                        "v2": v2.tolist(),
                        "cpp_result": cpp_expected,
                        "py_result": py_result,
                    }
                )

        assert len(failures) == 0, (
            f"Vector equivalence mismatch in {len(failures)} cases:\n"
            + "\n".join(str(f) for f in failures)
        )

    def test_specific_cubic_equivalences(self, au_fcc_symmetry):
        """Test known cubic symmetry equivalences."""
        # {111} family - all permutations and sign changes should be equivalent
        assert au_fcc_symmetry.vectors_equivalent([1, 1, 1], [1, -1, -1])
        assert au_fcc_symmetry.vectors_equivalent([1, 1, 1], [-1, 1, -1])
        assert au_fcc_symmetry.vectors_equivalent([1, 1, 1], [-1, -1, 1])

        # {100} family - all axes are equivalent
        assert au_fcc_symmetry.vectors_equivalent([1, 0, 0], [0, 1, 0])
        assert au_fcc_symmetry.vectors_equivalent([1, 0, 0], [0, 0, 1])
        assert au_fcc_symmetry.vectors_equivalent([2, 0, 0], [0, 2, 0])

        # Non-equivalent vectors
        assert not au_fcc_symmetry.vectors_equivalent([1, 0, 0], [1, 1, 0])
        assert not au_fcc_symmetry.vectors_equivalent([1, 1, 1], [2, 0, 0])

    def test_equivalence_tolerance(self, au_fcc_symmetry):
        """Test numerical tolerance for equivalence checking."""
        v1 = np.array([1.0, 1.0, 1.0])
        v2_exact = np.array([1.0, -1.0, -1.0])
        v2_approx = np.array([1.0, -1.0 + 1e-8, -1.0])

        # Exact match
        assert au_fcc_symmetry.vectors_equivalent(v1, v2_exact)

        # Within default tolerance (1e-5)
        assert au_fcc_symmetry.vectors_equivalent(v1, v2_approx)

        # Outside tight tolerance
        assert not au_fcc_symmetry.vectors_equivalent(v1, v2_approx, tolerance=1e-10)


class TestUniqueVectors:
    """Test symmetry reduction of vector lists."""

    def test_remove_cubic_duplicates(self, au_fcc_symmetry):
        """Test removal of symmetry-equivalent vectors."""
        # {111} family: all 8 should reduce to 1
        vectors_111 = [
            [1, 1, 1],
            [1, -1, -1],
            [-1, 1, -1],
            [-1, -1, 1],
            [1, 1, -1],
            [1, -1, 1],
            [-1, 1, 1],
            [-1, -1, -1],
        ]

        unique = au_fcc_symmetry.get_unique_vectors(vectors_111)
        assert len(unique) == 1, f"Expected 1 unique, got {len(unique)}"

    def test_mixed_families(self, au_fcc_symmetry):
        """Test reduction of mixed vector families."""
        vectors = [
            [1, 1, 1],  # {111} family
            [1, -1, -1],  # {111} - should be removed
            [2, 0, 0],  # {200} family
            [0, 2, 0],  # {200} - should be removed
            [2, 2, 0],  # {220} family
        ]

        unique = au_fcc_symmetry.get_unique_vectors(vectors)

        # Should have 3 unique families: {111}, {200}, {220}
        assert len(unique) == 3, f"Expected 3 unique families, got {len(unique)}"


class TestSymmetryInfo:
    """Test symmetry metadata and string representations."""

    def test_space_group_info(self, au_fcc_symmetry):
        """Verify space group information for FCC."""
        assert au_fcc_symmetry.space_group_number == 225  # Fm-3m
        assert "m-3m" in au_fcc_symmetry.point_group  # Cubic point group

    def test_string_representation(self, au_fcc_symmetry):
        """Test __repr__ output."""
        repr_str = repr(au_fcc_symmetry)
        assert "space_group=225" in repr_str
        assert "point_group" in repr_str


# ==================== Performance Benchmarks ====================

@pytest.mark.benchmark
class TestPerformance:
    """Performance tests (optional, requires pytest-benchmark)."""

    def test_equivalence_check_performance(self, au_fcc_symmetry, benchmark):
        """Benchmark vector equivalence checking."""

        def check_equivalence():
            return au_fcc_symmetry.vectors_equivalent([1, 1, 1], [1, -1, -1])

        result = benchmark(check_equivalence)
        assert result is True

    def test_unique_vectors_performance(self, au_fcc_symmetry, benchmark):
        """Benchmark symmetry reduction of 100 vectors."""
        # Generate 100 test vectors
        vectors = [[i % 5, (i + 1) % 5, (i + 2) % 5] for i in range(100)]

        result = benchmark(au_fcc_symmetry.get_unique_vectors, vectors)
        assert len(result) > 0


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
