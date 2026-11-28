"""Tests of icenine.gradient_sample (AR3 step 3): the spread, level 0, the centre, the seed."""

from typing import Dict, Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine.gradient_sample import LEVELS_DEG, gradient_orientations, spread_deg
from icenine.ray_indexing import misorientation_deg

SYM = np.eye(3)[None]  # identity only: enough for small rotations


def _layer(
    n_grains: int = 6, per: int = 30, seed: int = 3
) -> Tuple[np.ndarray, np.ndarray, Dict[int, np.ndarray]]:
    rg = np.random.default_rng(seed)
    xy, ids = [], []
    for g in range(n_grains):
        ctr = rg.uniform(-1, 1, 2)
        xy.append(ctr + rg.normal(0, 0.05, (per, 2)))
        ids += [g] * per
    pos = np.concatenate(xy)
    pos = np.hstack([pos, np.zeros((len(pos), 1))])
    Rg = {g: Rotation.random(random_state=g).as_matrix() for g in range(n_grains)}
    return pos, np.array(ids), Rg


def _levels(n: int = 6) -> Dict[int, float]:
    return {g: LEVELS_DEG[g % 5] for g in range(n)}


def test_measured_spread_equals_level() -> None:
    xy, ids, Rg = _layer()
    lv = _levels()
    out = gradient_orientations(xy, ids, Rg, lv, np.random.default_rng(0))
    sp = spread_deg(out, ids, Rg, SYM)
    for g, s in sp.items():
        assert s == pytest.approx(lv[g], abs=1e-6)


def test_level_zero_bit_for_bit_and_one_voxel_grain() -> None:
    xy, ids, Rg = _layer()
    ids = ids.copy()
    ids[0] = 99  # a grain of one voxel
    Rg[99] = Rotation.random(random_state=99).as_matrix()
    lv = _levels()
    lv[99] = 1.0
    out = gradient_orientations(xy, ids, Rg, lv, np.random.default_rng(0))
    assert np.array_equal(out[ids == 0], np.broadcast_to(Rg[0], out[ids == 0].shape))
    assert np.array_equal(out[ids == 99][0], Rg[99])


def test_centre_orientation_kept() -> None:
    """The orientation at the grain centre c is the grain's orientation.

    Checked through the mean rotation vector: the model's rotation vector is M (x - c), linear in
    position, and c is the mean of the voxel centres, so the rotation vector at c is zero exactly
    when the mean of the voxels' rotation vectors is zero. For this model the two tests are
    equivalent (no voxel need sit at c)."""
    xy, ids, Rg = _layer()
    out = gradient_orientations(xy, ids, Rg, _levels(), np.random.default_rng(0))
    for g in range(1, 6):
        sel = ids == g
        # the rotation vector is linear in (x - c): its mean over the grain is zero
        rv = Rotation.from_matrix(out[sel] @ Rg[g].T).as_rotvec()
        assert np.linalg.norm(rv.mean(axis=0)) < 1e-9


def test_same_seed_agrees_and_positions_untouched() -> None:
    xy, ids, Rg = _layer()
    xy0, ids0 = xy.copy(), ids.copy()
    a = gradient_orientations(xy, ids, Rg, _levels(), np.random.default_rng(0))
    b = gradient_orientations(xy, ids, Rg, _levels(), np.random.default_rng(0))
    c = gradient_orientations(xy, ids, Rg, _levels(), np.random.default_rng(1))
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)
    assert np.array_equal(xy, xy0) and np.array_equal(ids, ids0)


def test_orientations_are_rotations_and_bounded() -> None:
    xy, ids, Rg = _layer()
    out = gradient_orientations(xy, ids, Rg, _levels(), np.random.default_rng(0))
    assert np.abs(out @ out.transpose(0, 2, 1) - np.eye(3)).max() < 1e-12
    assert np.allclose(np.linalg.det(out), 1.0)
    for g in range(6):
        sel = ids == g
        d = misorientation_deg(out[sel], np.broadcast_to(Rg[g], out[sel].shape), SYM)
        assert d.max() <= _levels()[g] + 1e-6
