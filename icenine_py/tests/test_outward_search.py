"""Task F: the outward searches of icenine/outward_search.py.

Shell geometry and ``ShellSearch`` are checked on a synthetic cost (the angle to a hidden truth, no
images); ``RayGaussNewton`` on the ThreeVoxels toy problem (real two-detector images): the exact
Jacobian against finite differences of the residual, recovery from starts 1 and 3 degrees out, the
strict cost evaluated exactly twice, and the callback hook."""

import math
import sys
from pathlib import Path
from typing import Any, Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine import cost_functions as CF
from icenine.outward_search import (
    RayGaussNewton,
    ShellSearch,
    angle_deg,
    rotate_left,
    shell_directions,
    shell_rotations,
    shell_size,
)
from icenine.reconstructor import _get_voxel_vertices

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))


# ------------------------------------------------------------------------------ shells
def test_shell_size_grows_with_the_square_of_the_radius() -> None:
    assert shell_size(0.01, 0.3) == 6  # a floor
    n1, n2 = shell_size(1.5, 0.3), shell_size(3.0, 0.3)
    assert 3.7 < n2 / n1 < 4.3
    # hexagonal packing: area / (sqrt(3)/2 s^2)
    assert n1 == math.ceil(4 * math.pi * 1.5**2 / (math.sqrt(3) / 2 * 0.3**2))


def test_shell_points_lie_on_the_sphere_and_are_deterministic() -> None:
    R0 = Rotation.random(random_state=3).as_matrix()
    Rs = shell_rotations(R0, 2.0, 0.3, 4)
    assert Rs.shape == (shell_size(2.0, 0.3), 3, 3)
    np.testing.assert_allclose([angle_deg(R, R0) for R in Rs], 2.0, atol=1e-9)
    np.testing.assert_array_equal(Rs, shell_rotations(R0, 2.0, 0.3, 4))
    assert not np.allclose(shell_directions(50, 1), shell_directions(50, 2))  # turned per shell
    d = shell_directions(200, 0)
    np.testing.assert_allclose(np.linalg.norm(d, axis=1), 1.0, atol=1e-12)
    assert abs(d[:, 2].mean()) < 1e-9  # equal-area: symmetric in z


def test_shell_covers_its_sphere_within_the_spacing() -> None:
    """Every random direction on a shell of radius r has a shell point within the spacing."""
    R0 = np.eye(3)
    r, s = 2.4, 0.3
    Rs = shell_rotations(R0, r, s, 1)
    probe = np.random.default_rng(0).normal(size=(300, 3))
    probe /= np.linalg.norm(probe, axis=1, keepdims=True)
    P = rotate_left(R0, np.radians(r) * probe)
    worst = max(min(angle_deg(p, R) for R in Rs) for p in P)
    assert worst < s


# ------------------------------------------------------------------------------ ShellSearch
class _Synthetic:
    """Strict cost 0.05 + angle/1 deg (so it passes 0.8 within 0.75 deg); the screen scores the
    negative angle; the smooth cost is the angle. All counted."""

    def __init__(self, R_true: np.ndarray) -> None:
        self.R_true = R_true
        self.n_strict = self.n_screen = self.n_smooth = 0

    def strict(self, R: np.ndarray) -> float:
        self.n_strict += 1
        return min(1.0, 0.05 + angle_deg(R, self.R_true))

    def screen(self, Rs: np.ndarray) -> np.ndarray:
        self.n_screen += len(Rs)
        return np.array([-angle_deg(R, self.R_true) for R in Rs])

    def smooth(self, Rs: np.ndarray) -> np.ndarray:
        self.n_smooth += len(Rs)
        return np.array([angle_deg(R, self.R_true) for R in Rs])


def _search(offset: float, r_max: float = 6.0, with_smooth: bool = True) -> Tuple[Any, Any, Any]:
    R_true = Rotation.from_rotvec(np.radians([3.0, -2.0, 1.0])).as_matrix()
    v = np.array([0.3, -0.5, 0.8])
    R0 = rotate_left(R_true, np.radians(offset) * v / np.linalg.norm(v))
    syn = _Synthetic(R_true)
    ss = ShellSearch(
        syn.screen,
        syn.strict,
        syn.smooth if with_smooth else None,
        spacing_deg=0.3,
        r_max_deg=r_max,
        accept_cost=0.8,
    )
    return ss.search(R0), syn, R_true


def test_shell_search_finds_the_truth_at_the_shell_of_its_offset() -> None:
    res, syn, R_true = _search(2.0)
    assert res.found and angle_deg(res.R, R_true) < 0.75
    assert 1.5 <= res.radius_deg <= 2.5
    assert res.n_strict == syn.n_strict and res.n_screen == syn.n_screen
    assert (
        res.n_smooth == syn.n_smooth and res.n_evals == syn.n_strict + syn.n_screen + syn.n_smooth
    )
    assert res.n_screen == sum(shell_size(k * 0.3, 0.3) for k in range(1, res.n_shells + 1))
    assert res.trace[-1][2] == res.n_evals


def test_shell_search_cost_grows_with_the_offset_and_stops_early_when_close() -> None:
    near, _, _ = _search(0.3)
    mid, _, _ = _search(1.5)
    far, _, _ = _search(3.0)
    assert near.found and near.n_screen == 0 and near.radius_deg == 0.0  # the start passes
    assert near.n_evals == 1
    assert near.n_evals < mid.n_evals < far.n_evals


def test_shell_search_returns_the_best_point_when_nothing_passes() -> None:
    res, syn, R_true = _search(4.0, r_max=1.5)
    assert not res.found and math.isnan(res.radius_deg)
    assert res.n_shells == 5
    assert (
        res.strict_cost > 0.8 and angle_deg(res.R, R_true) < 4.0
    )  # the best it saw, not the start


def test_shell_search_without_the_smooth_stage_still_works() -> None:
    res, syn, R_true = _search(1.2, with_smooth=False)
    assert res.found and res.n_smooth == 0 and syn.n_smooth == 0


def test_shell_search_validates_its_arguments() -> None:
    f = lambda R: 0.0  # noqa: E731
    with pytest.raises(ValueError):
        ShellSearch(f, f, spacing_deg=0.0)
    with pytest.raises(ValueError):
        ShellSearch(f, f, n_smooth=1, n_strict=2)


# ------------------------------------------------------------------------------ RayGaussNewton
def _toy_gn(
    tau_schedule: Tuple[float, ...] = (0.25,)
) -> Tuple[RayGaussNewton, Any, np.ndarray, Any]:
    from test_smooth_cost import _sc, _strict, _toy

    rec, voxel, R_true = _toy()
    sc, cf = _sc(rec), _strict(rec)
    verts = _get_voxel_vertices(voxel)
    calls = {"n": 0}

    def strict(R: np.ndarray) -> float:
        calls["n"] += 1
        return float(cf.evaluate(R.astype(np.float32), verts, voxel.phase).cost)

    gn = RayGaussNewton(sc, verts, strict, voxel.phase, tau_schedule_deg=tau_schedule)
    return gn, calls, R_true, verts


@needs_c
def test_ray_gn_jacobian_matches_finite_differences_of_the_residual() -> None:
    gn, _, R_true, _ = _toy_gn()
    R = rotate_left(R_true, np.radians([0.15, -0.1, 0.2]))
    tau = math.radians(0.25)
    lin = gn.linearize(R, tau)
    assert lin["m"] > 20 and lin["J"].shape == (3 * lin["m"], 3)
    h = 1e-6
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        up = gn.linearize(rotate_left(R, e), tau)
        dn = gn.linearize(rotate_left(R, -e), tau)
        assert up["m"] == dn["m"] == lin["m"]
        fd = ((up["res"] - dn["res"]) / (2 * h)).numpy()
        np.testing.assert_allclose(lin["J"][:, i].numpy(), fd, atol=2e-5)
    # the cost reported by the linearisation is the cost of the value function
    assert lin["cost"] == pytest.approx(gn.e4_value(R, tau), abs=1e-12)


@needs_c
def test_ray_gn_recovers_the_truth_from_one_and_three_degrees() -> None:
    gn, calls, R_true, _ = _toy_gn()
    d = np.array([0.6, -0.3, 0.74])
    d /= np.linalg.norm(d)
    for off in (1.0, 3.0):
        calls["n"] = 0
        seen = []
        res = gn.solve(
            rotate_left(R_true, np.radians(off) * d),
            callback=lambda R, ns, nt: seen.append((ns, nt)),
        )
        assert angle_deg(res.R, R_true) < 0.1, off
        assert res.accepted and res.strict_cost <= 0.8 and res.converged
        assert res.e4_end < res.e4_start
        assert calls["n"] == 2 == res.n_strict  # start and end only
        assert res.n_evals == res.n_jac + res.n_val + 2 and res.n_jac >= 3
        assert seen[0][0] == 1 and seen[-1][0] <= res.n_smooth and len(seen) >= 2
        assert [s[0] for s in seen] == sorted(s[0] for s in seen)


@needs_c
def test_ray_gn_at_the_optimum_does_nothing_and_a_bad_start_is_not_accepted() -> None:
    gn, _, R_true, _ = _toy_gn()
    res = gn.solve(R_true)
    assert angle_deg(res.R, R_true) < 0.05 and res.accepted
    assert res.n_jac <= 4
    with pytest.raises(ValueError):
        RayGaussNewton(gn.sc, gn.verts, gn.strict_fn, 0, tau_schedule_deg=(0.0,))
    # the end result is never worse by the strict cost than the start
    far = rotate_left(R_true, np.radians(40.0) * np.array([1.0, 0.0, 0.0]))
    out = gn.solve(far)
    assert out.strict_cost <= out.start_strict_cost
