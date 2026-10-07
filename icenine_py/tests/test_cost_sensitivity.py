"""scripts/cost_sensitivity: resolution (Fisher / Cramer-Rao) formula, plateau statistics."""

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
for sub in ("scripts/cost_sensitivity", "scripts/common"):
    sys.path.insert(0, str(ROOT / sub))

import resolution as res  # noqa: E402
import summary as sm  # noqa: E402


def test_fisher_single_reflection():
    dw = 0.01
    gw = np.array([[2.0, 0.0, 0.0]])
    gu = np.zeros((1, 2, 3))
    gu[0, 0, 0] = 1.0
    gu[0, 1, 1] = 3.0
    jf, jp = res.fisher_quantisation(gw, gu, dw)
    assert np.allclose(jf, np.diag([12 * 4 / dw**2, 0, 0]))
    assert np.allclose(jp, np.diag([12.0, 12.0 * 9, 0.0]))
    # bin sizes scale the variance as bin^2
    jf2, jp2 = res.fisher_quantisation(gw, gu, dw, omega_bin=2.0, pixel_bin=2.0)
    assert np.allclose(jf2, jf / 4) and np.allclose(jp2, jp / 4)


def test_fisher_adds_over_reflections():
    rng = np.random.default_rng(0)
    gw, gu = rng.normal(size=(5, 3)), rng.normal(size=(5, 2, 3))
    jf, jp = res.fisher_quantisation(gw, gu, 0.002)
    parts = [res.fisher_quantisation(gw[i : i + 1], gu[i : i + 1], 0.002) for i in range(5)]
    assert np.allclose(jf, sum(p[0] for p in parts)) and np.allclose(jp, sum(p[1] for p in parts))


def test_crb_diagonal_and_singular():
    J = np.diag([100.0, 400.0, 900.0])
    out = res.crb(J)
    assert np.allclose(out["sigma_axes_deg"], np.degrees([0.1, 0.05, 1 / 30]))
    assert math.isclose(out["rms3_deg"], math.degrees(math.sqrt(0.01 + 1 / 400 + 1 / 900)))
    assert np.allclose(out["principal_deg"], np.degrees([1 / 30, 0.05, 0.1]))  # ascending sigma
    assert res.crb(np.diag([1.0, 1.0, 0.0]))["rms3_deg"] == math.inf


def test_more_peaks_tighten_bound():
    rng = np.random.default_rng(1)
    gw, gu = rng.normal(size=(40, 3)), rng.normal(size=(40, 2, 3))
    r_small = res.crb(sum(res.fisher_quantisation(gw[:10], gu[:10], 0.002)))
    r_big = res.crb(sum(res.fisher_quantisation(gw, gu, 0.002)))
    assert r_big["rms3_deg"] < r_small["rms3_deg"]


def _plateau(scale: list, radii: np.ndarray, eps: float) -> dict:
    import landscape as ls

    dirs = ls.unit_directions(400, 3)
    cost = np.sqrt(((radii[:, None, None] * dirs[None]) ** 2 * np.array(scale) ** 2).sum(-1))
    return sm.plateau_stats(cost, dirs, radii, 0.0, eps)


def test_plateau_isotropic_and_anisotropic():
    radii = np.linspace(0.01, 0.1, 10)
    iso = _plateau([1.0, 1.0, 1.0], radii, 0.05)  # cost = r, plateau r <= 0.05
    assert iso["r_any"] == radii[4] and iso["ext_med"] == radii[4]
    assert iso["aniso"] < 1.3
    assert iso["r50"] == radii[5]
    # cost = r * (1, 1/3, 1/3 scaled): plateau elongated 3x along two axes
    an = _plateau([1.0, 1 / 3, 1 / 3], radii, 0.05)
    assert 1.5 < an["aniso"] < 4.0
    assert an["r_any"] > iso["r_any"]


def test_min_set_stats_tie_safe():
    # one case: samples (cost, radius, offset). Three samples tie at the minimum 0.1 (< truth 0.5):
    # permuting the samples must not change any statistic.
    cost = np.array([[0.3, 0.1, 0.1, 0.4, 0.1, 0.6]])
    rad = np.array([[0.001, 0.002, 0.002, 0.003, 0.0005, 0.004]])
    pts = np.array(
        [[[0, 0, 0.001], [0.002, 0, 0], [-0.002, 0, 0], [0, 0, 0], [0, 0.0005, 0], [0, 0, 1]]]
    )
    ct = np.array([0.5])
    out = sm.min_set_stats(cost, rad, pts, ct)
    assert out["n_lower"] == 1
    assert math.isclose(out["min_drop_among_lower_quartiles"][1], 0.4)
    assert out["n_tied_at_min_quartiles"][1] == 3
    assert out["min_radius_smallest_tied_quartiles"][1] == 0.0005
    # centroid of the three tied offsets: (0, 0.0005/3, 0)
    assert np.allclose(out["centroid_offset"], 0.0005 / 3)
    perm = np.random.default_rng(0).permutation(6)
    out2 = sm.min_set_stats(cost[:, perm], rad[:, perm], pts[:, perm], ct)
    assert np.allclose(out2["centroid_offset"], out["centroid_offset"])
    assert out2["min_radius_smallest_tied_quartiles"] == out["min_radius_smallest_tied_quartiles"]
    # truth already at the minimum: no case is "lower", the plateau counts the tie
    out3 = sm.min_set_stats(cost, rad, pts, np.array([0.1]))
    assert out3["n_lower"] == 0 and out3["plateau_n_samples_quartiles"][1] == 3


def test_monotone_rays():
    grid = np.array([0.001, 0.005, 0.01])
    ct = np.array([1.0, 1.0])
    # case 0: rays (0 up, 1 dips at the last radius, 2 dips at the first radius)
    c0 = np.array([[1.1, 1.2, 1.3], [1.1, 1.2, 1.15], [0.9, 1.2, 1.3]]).T
    c1 = np.array([[1.1, 1.2, 1.3]] * 3).T
    out = sm.monotone_rays(ct, np.stack([c0, c1]), grid)
    assert math.isclose(out["all_radii"], 4 / 6)  # rays 0 of case 0 and the three of case 1
    # from 0.005 outward the dip at the first radius no longer counts, only the last-radius dip does
    assert math.isclose(out["from_0p005"], 5 / 6)
    assert out["cases_all_rays_all_radii"]["k"] == 1 and out["cases_all_rays_from_0p005"]["k"] == 1


def test_lower_than_result():
    flat = np.array([[0.1, 0.2, 0.3, 0.9], [0.5, 0.6, 0.7, 0.8]])
    rad = np.array([[0.001, 0.02, 0.05, 0.1]] * 2)
    # case 0: result cost 0.5 at 0.01 deg: points 0.1 (r 0.001, closer), 0.2 (0.02), 0.3 (0.05)
    # are lower, two of them farther. Case 1: no lower point. Case 2 unusable (error > 0.1)
    flat = np.vstack([flat, flat[:1]])
    rad = np.vstack([rad, rad[:1]])
    out = sm.lower_than_result(flat, rad, np.array([0.5, 0.4, 0.5]), np.array([0.01, 0.01, 0.5]))
    assert out["n_usable"] == 2 and out["cases_with_lower"] == 1
    assert out["points"] == 3 and out["points_farther"] == 2 and out["cases_with_farther"] == 1
    assert math.isclose(out["per_case_fraction_farther_quartiles"][1], 2 / 3)


def test_small_radius_formatting():
    # radii below 0.001 must not round to 0.001 in the tables
    assert sm.fq([0.0005, 0.0005, 0.002], 4) == "0.0005/0.0005/0.0020"
    assert f"{0.0005:g}" == "0.0005"
