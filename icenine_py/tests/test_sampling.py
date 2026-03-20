"""
Tests for SO(3) sampling module — Sukharev grids and quaternion operations.

Tests verify:
1. Sukharev grid point counts match expected values
2. Local grid generation produces correct number of orientations
3. Quaternion arithmetic is correct
4. SLERP matches known values
5. FZ reduction works for cubic symmetry
"""

import math

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine.sampling import (
    QuaternionGrid,
    _quat_inverse,
    _quat_multiply,
    generate_local_grid,
    generate_local_grid_multi_level,
    get_layered_sukarev_grid_point,
    get_misorientation,
    get_sukarev_grid_point,
    is_in_fundamental_zone,
    make_sukarev_grid_points,
    matrix_to_quaternion,
    quaternion_to_matrix,
    reduce_to_fundamental_zone,
    slerp,
)


class TestSukharevGrid:
    """Tests for the Sukharev grid sequence."""

    def test_sukarev_grid_point_range(self):
        """Grid points should be in [0, 1)^3."""
        for state in range(100):
            pt = get_sukarev_grid_point(state)
            assert pt.shape == (3,)
            assert np.all(pt >= 0.0)
            assert np.all(pt < 1.0)

    def test_layered_sukarev_grid_point_range(self):
        """Layered grid points should be in [0, 1)^3."""
        for state in range(100):
            pt = get_layered_sukarev_grid_point(state)
            assert pt.shape == (3,)
            assert np.all(pt >= 0.0)
            assert np.all(pt <= 1.0)

    def test_make_sukarev_grid_point_count(self):
        """Grid at level L should have (2^L)^3 points."""
        for level in range(4):
            pts = make_sukarev_grid_points(level)
            expected = (2 ** level) ** 3
            assert pts.shape == (expected, 3), f"Level {level}: expected {expected}, got {pts.shape[0]}"

    def test_make_sukarev_grid_level0(self):
        """Level 0: single point at center of unit cube."""
        pts = make_sukarev_grid_points(0, side_width=1.0)
        assert pts.shape == (1, 3)
        np.testing.assert_allclose(pts[0], [0.5, 0.5, 0.5])

    def test_make_sukarev_grid_level1(self):
        """Level 1: 8 points in 2x2x2 grid."""
        pts = make_sukarev_grid_points(1, side_width=1.0)
        assert pts.shape == (8, 3)
        # Points should be at centers of 8 sub-cubes
        expected_coords = {0.25, 0.75}
        for pt in pts:
            for val in pt:
                assert val in expected_coords or np.isclose(val, 0.25) or np.isclose(val, 0.75)

    def test_make_sukarev_grid_scaled(self):
        """Grid with custom side_width should be scaled."""
        pts = make_sukarev_grid_points(0, side_width=0.5)
        np.testing.assert_allclose(pts[0], [0.25, 0.25, 0.25])


class TestSLERP:
    """Tests for spherical linear interpolation."""

    def test_slerp_endpoints(self):
        """SLERP at t=0 gives q1, t=1 gives q2."""
        q1 = np.array([1.0, 0.0, 0.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0, 0.0])

        result_0 = slerp(q1, q2, 0.0)
        np.testing.assert_allclose(result_0, q1, atol=1e-10)

        result_1 = slerp(q1, q2, 1.0)
        np.testing.assert_allclose(np.abs(result_1), np.abs(q2), atol=1e-10)

    def test_slerp_midpoint(self):
        """SLERP at t=0.5 gives unit quaternion."""
        q1 = np.array([1.0, 0.0, 0.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0, 0.0])
        mid = slerp(q1, q2, 0.5)
        assert abs(np.linalg.norm(mid) - 1.0) < 1e-10

    def test_slerp_near_parallel(self):
        """SLERP handles near-parallel quaternions (linear fallback)."""
        q1 = np.array([1.0, 0.0, 0.0, 0.0])
        q2 = np.array([0.9999, 0.01, 0.0, 0.0])
        q2 = q2 / np.linalg.norm(q2)
        mid = slerp(q1, q2, 0.5)
        assert abs(np.linalg.norm(mid) - 1.0) < 0.01


class TestQuaternionArithmetic:
    """Tests for quaternion multiply, inverse, conversions."""

    def test_identity_multiply(self):
        """q * identity = q."""
        q = np.array([0.5, 0.5, 0.5, 0.5])
        identity = np.array([1.0, 0.0, 0.0, 0.0])
        result = _quat_multiply(q, identity)
        np.testing.assert_allclose(result, q, atol=1e-14)

    def test_inverse(self):
        """q * q_inv = identity for unit quaternion."""
        q = np.array([0.5, 0.5, 0.5, 0.5])
        q_inv = _quat_inverse(q)
        result = _quat_multiply(q, q_inv)
        np.testing.assert_allclose(result, [1, 0, 0, 0], atol=1e-14)

    def test_quaternion_to_matrix_identity(self):
        """Identity quaternion → identity matrix."""
        q = np.array([1.0, 0.0, 0.0, 0.0])
        m = quaternion_to_matrix(q)
        np.testing.assert_allclose(m, np.eye(3), atol=1e-14)

    def test_roundtrip_matrix_quaternion(self):
        """Matrix → quaternion → matrix roundtrip."""
        r = Rotation.random(random_state=42)
        m = r.as_matrix()
        q = matrix_to_quaternion(m)
        m2 = quaternion_to_matrix(q)
        np.testing.assert_allclose(m, m2, atol=1e-10)


class TestQuaternionGrid:
    """Tests for the CQuaternionGrid port."""

    def test_grid_by_count(self):
        """Grid should return exactly n_points quaternions."""
        grid = QuaternionGrid()
        for n in [4, 16, 100]:
            quats = grid.get_grid_by_count(n)
            assert quats.shape == (n, 4)
            # All should be unit quaternions
            norms = np.linalg.norm(quats, axis=1)
            np.testing.assert_allclose(norms, 1.0, atol=1e-6)

    def test_grid_positive_hemisphere(self):
        """All quaternions should be in positive-w hemisphere."""
        grid = QuaternionGrid()
        quats = grid.get_grid_by_count(100)
        assert np.all(quats[:, 0] >= 0), "All w components should be >= 0"

    def test_grid_by_dispersion(self):
        """Grid by dispersion should return reasonable number of points."""
        grid = QuaternionGrid()
        quats = grid.get_grid_by_dispersion(0.5)
        assert quats.shape[1] == 4
        assert quats.shape[0] > 10  # Should have a decent number of points

    def test_structured_local_grid_count(self):
        """Local grid at level L should have (2^L)^3 points."""
        grid = QuaternionGrid()
        for level in range(4):
            quats = grid.get_structured_local_grid(
                math.radians(5.0), level
            )
            expected = (2 ** level) ** 3
            assert quats.shape == (expected, 4), (
                f"Level {level}: expected {expected}, got {quats.shape[0]}"
            )

    def test_structured_local_grid_near_identity(self):
        """Local grid points should be small rotations near identity."""
        grid = QuaternionGrid()
        quats = grid.get_structured_local_grid(math.radians(5.0), 1)
        # All quaternions should have large w component (near identity)
        assert np.all(quats[:, 0] > 0.9), (
            f"Expected w > 0.9, got min w = {quats[:, 0].min()}"
        )


class TestGenerateLocalGrid:
    """Tests for the high-level GenerateLocalGrid wrapper."""

    def test_generates_rotation_matrices(self):
        """Output should be valid rotation matrices."""
        matrices = generate_local_grid(math.radians(5.0), 1)
        assert matrices.shape == (8, 3, 3)
        for m in matrices:
            # Check orthogonality: M @ M.T = I
            np.testing.assert_allclose(m @ m.T, np.eye(3), atol=1e-10)
            # Check det = 1
            np.testing.assert_allclose(np.linalg.det(m), 1.0, atol=1e-10)

    def test_multi_level(self):
        """Multi-level grid concatenates all levels."""
        matrices = generate_local_grid_multi_level(math.radians(5.0), 0, 2)
        # Level 0: 1, Level 1: 8, Level 2: 64
        expected = 1 + 8 + 64
        assert matrices.shape == (expected, 3, 3)


class TestFundamentalZone:
    """Tests for FZ reduction."""

    @staticmethod
    def _get_cubic_symmetry_quats():
        """Get 24 proper cubic symmetry operators as quaternions [w, x, y, z]."""
        from icenine.symmetry import create_cubic_symmetry
        sym = create_cubic_symmetry(4.0)
        matrices = sym.get_rotation_matrices()
        # Filter to proper rotations only (det = +1)
        proper = [m for m in matrices if np.linalg.det(m) > 0]
        quats = np.array([matrix_to_quaternion(np.array(m)) for m in proper])
        return quats

    def test_identity_in_fz(self):
        """Identity quaternion should be in FZ."""
        sym_quats = self._get_cubic_symmetry_quats()
        q_id = np.array([1.0, 0.0, 0.0, 0.0])
        assert is_in_fundamental_zone(q_id, sym_quats)

    def test_reduce_identity(self):
        """Reducing identity should return identity."""
        sym_quats = self._get_cubic_symmetry_quats()
        q_id = np.array([1.0, 0.0, 0.0, 0.0])
        q_fz = reduce_to_fundamental_zone(q_id, sym_quats)
        np.testing.assert_allclose(q_fz, q_id, atol=1e-10)

    def test_reduce_maximizes_w(self):
        """FZ reduction should pick the symmetry-equivalent with largest |w|."""
        sym_quats = self._get_cubic_symmetry_quats()
        # Create a random quaternion
        r = Rotation.random(random_state=42)
        q_scipy = r.as_quat()  # [x, y, z, w]
        q = np.array([q_scipy[3], q_scipy[0], q_scipy[1], q_scipy[2]])

        q_fz = reduce_to_fundamental_zone(q, sym_quats)

        # Verify that no symmetry equivalent has larger |w|
        for sym_q in sym_quats:
            product = _quat_multiply(q, sym_q)
            assert abs(product[0]) <= abs(q_fz[0]) + 1e-10

    def test_misorientation_identity(self):
        """Misorientation between identical orientations should be 0."""
        sym_quats = self._get_cubic_symmetry_quats()
        q = np.array([1.0, 0.0, 0.0, 0.0])
        angle = get_misorientation(q, q, sym_quats)
        assert abs(angle) < 1e-10

    def test_misorientation_symmetric(self):
        """Misorientation should be symmetric."""
        sym_quats = self._get_cubic_symmetry_quats()
        r1 = Rotation.random(random_state=42)
        r2 = Rotation.random(random_state=43)
        q1 = matrix_to_quaternion(r1.as_matrix())
        q2 = matrix_to_quaternion(r2.as_matrix())
        assert abs(
            get_misorientation(q1, q2, sym_quats)
            - get_misorientation(q2, q1, sym_quats)
        ) < 1e-10
