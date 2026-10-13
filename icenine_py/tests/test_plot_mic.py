"""Tests for scripts/plot_mic.py: symmetry handling, colours and voxel matching."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import plot_mic as P  # noqa: E402
from icenine.orientation_search import get_symmetry_quaternions  # noqa: E402
from icenine.symmetry import create_cubic_symmetry  # noqa: E402


@pytest.fixture(scope="module")
def sym() -> np.ndarray:
    return np.asarray(get_symmetry_quaternions(create_cubic_symmetry(1.0)), dtype=np.float64)


def _random_quats(n: int, seed: int = 0) -> np.ndarray:
    q = np.random.default_rng(seed).normal(size=(n, 4))
    return q / np.linalg.norm(q, axis=1, keepdims=True)


def test_symmetry_equivalents_have_zero_misorientation_and_same_colour(sym: np.ndarray) -> None:
    q = _random_quats(50)
    q_eq = P.qmul(q, sym[np.arange(50) % len(sym)])
    assert np.allclose(P.misorientation_deg(q, q_eq, sym), 0.0, atol=1e-5)
    assert np.allclose(P.orientation_colors(q, sym), P.orientation_colors(q_eq, sym), atol=1e-6)
    assert np.allclose(P.orientation_colors(q, sym), P.orientation_colors(-q, sym), atol=1e-6)


def test_misorientation_bounded_by_cubic_maximum(sym: np.ndarray) -> None:
    a, b = _random_quats(2000, 1), _random_quats(2000, 2)
    m = P.misorientation_deg(a, b, sym)
    assert m.max() <= P.CUBIC_MAX_MISORIENTATION_DEG + 1e-3
    assert m.max() > 55.0


def test_known_rotation_angle(sym: np.ndarray) -> None:
    ident = np.array([[1.0, 0.0, 0.0, 0.0]])
    rot = np.array([[np.cos(np.radians(5.0)), 0.0, 0.0, np.sin(np.radians(5.0))]])  # 10 deg about z
    assert P.misorientation_deg(ident, rot, sym)[0] == pytest.approx(10.0, abs=1e-6)


def test_fz_rodrigues_inside_cubic_zone(sym: np.ndarray) -> None:
    rf = P.fz_rodrigues(_random_quats(500, 3), sym)
    assert np.abs(rf).max() <= P.RF_MAX + 1e-9


def _grid(n: int) -> np.ndarray:
    tri = []
    for i in range(n):
        x = float(i)
        up = i % 2 == 0
        h = np.sqrt(3.0) / 2.0 * (1.0 if up else -1.0)
        tri.append([[x, 0.0], [x + 1.0, 0.0], [x + 0.5, h]])
    return np.array(tri)


def test_match_same_order_is_identity() -> None:
    t = _grid(10)
    assert (P.match_to_reference(t, t.copy(), 0.3) == np.arange(10)).all()


def test_match_shuffled_subset() -> None:
    t = _grid(10)
    perm = np.random.default_rng(0).permutation(10)[:7]
    m = P.match_to_reference(t, t[perm], 0.3)
    matched = np.where(m >= 0)[0]
    assert sorted(matched.tolist()) == sorted(perm.tolist())
    assert (perm[m[matched]] == matched).all()


def test_match_rejects_far_voxels() -> None:
    t = _grid(5)
    assert (P.match_to_reference(t, t + 5.0, 0.3) == -1).all()
