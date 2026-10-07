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


def _plateau(scale, radii, eps):
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
