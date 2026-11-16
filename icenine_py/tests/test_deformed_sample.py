"""Random voxel-to-voxel misorientation samples (icenine/deformed_sample.py)."""

import numpy as np
from scipy.spatial.transform import Rotation

from icenine.deformed_sample import random_field_orientations, random_rotvec_field


def test_field_has_zero_mean_and_the_requested_rms() -> None:
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 10, (60, 2))
    for corr in (0.0, 3.0):
        w = random_rotvec_field(xy, 2.0, corr, np.random.default_rng(1))
        assert np.abs(w.mean(0)).max() < 1e-12
        rms = np.degrees(np.sqrt((w * w).sum(1).mean()))
        assert abs(rms - 2.0) < 1e-9


def test_correlated_field_is_smoother_than_white_noise() -> None:
    g = np.stack(np.meshgrid(np.arange(10.0), np.arange(10.0)), -1).reshape(-1, 2)
    white = random_rotvec_field(g, 1.0, 0.0, np.random.default_rng(2))
    smooth = random_rotvec_field(g, 1.0, 3.0, np.random.default_rng(2))
    right = [(i, i + 1) for i in range(len(g)) if (i + 1) % 10]

    def step(w: np.ndarray) -> float:
        return float(np.mean([np.linalg.norm(w[i] - w[j]) for i, j in right]))

    assert step(smooth) < 0.3 * step(white)


def test_orientations_left_multiply_and_keep_level_zero_grains() -> None:
    rng = np.random.default_rng(3)
    xy = rng.uniform(0, 5, (30, 2))
    ids = np.array([0] * 15 + [1] * 15)
    Rg = {0: Rotation.random(random_state=4).as_matrix(), 1: Rotation.random(random_state=5).as_matrix()}
    R = random_field_orientations(xy, ids, Rg, {0: 1.0, 1: 0.0}, 0.0, np.random.default_rng(6))
    assert np.array_equal(R[15:], np.broadcast_to(Rg[1], (15, 3, 3)))
    d = R[:15] @ Rg[0].T  # exp([r]x)
    ang = np.degrees(Rotation.from_matrix(d).magnitude())
    assert abs(np.sqrt((ang**2).mean()) - 1.0) < 1e-9
