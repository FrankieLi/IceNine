"""Opt-in ray-consistency overlap rule (icenine/ray_overlap.py): cases a-f on hand-built input,
C stage equals the Python twin, observed-ray geometry on synthetic detectors, the window arms,
config keys, and the default (off) leaving the cost unchanged."""

import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from icenine import cost_functions as CF
from icenine.config_file import ConfigFile
from icenine.cost_functions import count_qualified_peaks
from icenine.orientation_search import SearchParameters
from icenine.ray_overlap import RayRule, apply_cases_python, ray_flags, ray_rule_from_config

needs_ext = pytest.mark.skipif(
    not CF._HAS_C_RASTERIZE or not getattr(CF, "_C_HAS_RAY_RULE", False),
    reason="C ext missing or built before the ray rule (rebuild _rasterize)",
)


def _stage_d(lit_pred, bright_pred, obs_in, obs_bright):
    """Run the C stage on one peak with 2 detectors whose predicted-spot state is set by images.

    lit_pred / bright_pred: predicted spot in bounds / bright on detector 1 (detector 0 is bright).
    """
    rows = cols = 16
    img0 = np.zeros((rows, cols), np.uint8)
    img0[8, 8] = 1
    img1 = np.zeros((rows, cols), np.uint8)
    if bright_pred:
        img1[8, 8] = 1
    wedge = np.zeros(1, np.int32)
    hit = np.ones((1, 2), np.uint8)
    centers = np.array([[[8, 8], [8, 8] if lit_pred else [-50, -50]]], np.int32)
    oi = np.array([[0, int(obs_in)]], np.uint8)
    ob = np.array([[0, int(obs_bright)]], np.uint8)
    return CF._c_stage_d_overlap(
        [img0, img1], wedge, hit, centers, None, 2, 1, 1, rows, cols, 1, oi, ob
    )


@needs_ext
def test_cases_a_to_f_c_stage():
    # (predicted in bounds, predicted bright, observed in bounds, observed bright) -> valid
    table = {
        "a": ((1, 1, 0, 0), 1, 0),
        "b": ((1, 0, 1, 1), 1, 1),
        "c": ((1, 0, 1, 0), 0, 2),
        "d": ((0, 0, 0, 0), 1, 3),
        "e": ((1, 0, 0, 0), 1, 4),
        "f": ((0, 0, 1, 0), 1, 5),  # always forgiven (f1 dropped 2026-10-12)
    }
    for name, (args, valid, case) in table.items():
        res = _stage_d(*args)
        assert res[2] == valid, name  # peak_overlap
        cases = res[7]
        assert cases[7] == 1 and cases[case] == 1, (name, cases)
        assert sum(cases[:7]) == 1, (name, cases)


@needs_ext
def test_c_stage_equals_python_twin():
    rng = np.random.default_rng(1)
    rows = cols = 32
    n_det, n_peaks, rad = 2, 300, 1
    images = [(rng.random((rows, cols)) < 0.3).astype(np.uint8) for _ in range(n_det)]
    wedge = np.zeros(n_peaks, np.int32)
    hit = np.ones((n_peaks, n_det), np.uint8)
    centers = np.ascontiguousarray(rng.integers(-3, 35, (n_peaks, n_det, 2)).astype(np.int32))
    oi = (rng.random((n_peaks, n_det)) < 0.6).astype(np.uint8)
    ob = (oi & (rng.random((n_peaks, n_det)) < 0.5)).astype(np.uint8)
    res = CF._c_stage_d_overlap(
        images, wedge, hit, centers, None, n_det, rad, n_peaks, rows, cols, 1, oi, ob
    )
    po = pon = 0
    cases = np.zeros(8, np.int64)
    for p in range(n_peaks):
        lit, ov = [False] * n_det, [False] * n_det
        for d in range(n_det):
            c, r = centers[p, d]
            sub = images[d][max(r - rad, 0) : r + rad + 1, max(c - rad, 0) : c + rad + 1]
            inb = bool(c + rad >= 0 and r + rad >= 0 and c - rad < cols and r - rad < rows)
            lit[d], ov[d] = inb, bool(inb and sub.any())
        apply_cases_python(lit, ov, oi[p], ob[p], cases)
        on, o, _ = count_qualified_peaks(lit, ov, n_det, cpp_rule=True)
        pon += on
        po += o
    assert (res[2], res[3]) == (po, pon)
    assert tuple(res[7]) == tuple(int(x) for x in cases)


def _fake_detector(x: float, n: int = 40):
    """Detector in the plane x = const: col = y + n/2, row = z + n/2, pixel 1 (arbitrary units)."""
    t = lambda *v: torch.tensor(v, dtype=torch.float64)  # noqa: E731
    return SimpleNamespace(
        _detector_plane=SimpleNamespace(normal=t(1.0, 0.0, 0.0), d=torch.tensor(-x)),
        _position=t(x, 0.0, 0.0),
        _lab_frame_coord_origin=t(0.0, 0.0, 0.0),
        _lab_frame_basis_j=t(0.0, 1.0, 0.0),
        _lab_frame_basis_k=t(0.0, 0.0, 1.0),
        pixel_width=1.0,
        pixel_height=1.0,
        pixel_half_width=n / 2.0,
        pixel_half_height=n / 2.0,
    )


class _Exp:
    def __init__(self, imgs):
        self.imgs = imgs

    def get_image(self, w, d):
        return SimpleNamespace(get_binary_numpy=lambda: self.imgs[d])


def _geometry_case(offset_col: int, rule: RayRule):
    n = 40
    img0 = np.zeros((n, n), np.uint8)
    img0[20, 22] = 1  # detector 0 pixel: col 22, row 20 -> y = 2.5, z = 0.5 at x = 10
    img1 = np.zeros((n, n), np.uint8)
    img1[22, 25 + offset_col] = 1  # the ray reaches (x=20): y = 5, z = 1 -> col 25, row 21
    px = np.zeros((1, 3, 2))
    px[:, :, 0], px[:, :, 1] = 22.5, 20.5
    return ray_flags(
        rule,
        [_fake_detector(10.0), _fake_detector(20.0)],
        _Exp([img0, img1]),
        np.zeros(1, np.int64),
        np.zeros((1, 3)),
        np.ones(1, bool),
        px,
        1,
        0.0,
        0.0,
    )


def test_ray_geometry_and_radius():
    # the bright detector-1 pixel is at row 22, one row below the ray hit (row 21)
    inb, br = _geometry_case(0, RayRule(arm="R", radius_px=0))
    assert (inb[0, 1], br[0, 1]) == (1, 0)
    inb, br = _geometry_case(0, RayRule(arm="R", radius_px=1))
    assert br[0, 1] == 1
    inb, br = _geometry_case(3, RayRule(arm="R", radius_px=2))
    assert br[0, 1] == 0
    inb, br = _geometry_case(3, RayRule(arm="R", radius_px=4))
    assert br[0, 1] == 1
    assert br[0, 0] == 0 and inb[0, 0] == 0  # detector 0 itself gets no flags


def test_anisotropic_window_grows_with_terms():
    # no voxel/arc/registration terms: the detector-0 pixel term only (L1/L0 = 2 -> 2 px half-width)
    _, br = _geometry_case(0, RayRule(arm="A"))
    assert br[0, 1] == 1  # row offset 1 is inside the half-width
    _, br = _geometry_case(3, RayRule(arm="A"))
    assert br[0, 1] == 0
    _, br = _geometry_case(3, RayRule(arm="A", registration_mm=3.0))
    assert br[0, 1] == 1


def test_off_ray_out_of_bounds():
    n = 40
    img = np.zeros((n, n), np.uint8)
    img[20, 39] = (
        1  # detector-0 bright pixel at the far edge: y = 19.5 -> detector 1 at y = 39 (col 59)
    )
    px = np.zeros((1, 3, 2))
    px[:, :, 0], px[:, :, 1] = 39.5, 20.5
    inb, _ = ray_flags(
        RayRule(arm="R", radius_px=0),
        [_fake_detector(10.0), _fake_detector(20.0)],
        _Exp([img, np.zeros((n, n), np.uint8)]),
        np.zeros(1, np.int64),
        np.zeros((1, 3)),
        np.ones(1, bool),
        px,
        1,
        0.0,
        0.0,
    )
    assert inb[0, 1] == 0  # col 59: beyond the 40-pixel detector


def test_off_by_default():
    assert SearchParameters().ray_rule is None
    assert RayRule().levels == 3
    with pytest.raises(ValueError):
        RayRule(arm="X")


def test_config_keys(tmp_path):
    base = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"
    text = (base / "ReconstructPhaseD_mc_clean.config").read_text()
    plain = tmp_path / "p.config"
    plain.write_text(text)
    cfg = ConfigFile.from_file(str(plain))
    assert cfg.ray_overlap_rule == "off" and ray_rule_from_config(cfg) is None
    assert SearchParameters.from_config(cfg).ray_rule is None
    for spec, arm, rad in (("R4", "R", 4), ("A", "A", 2)):
        f = tmp_path / f"{spec}.config"
        f.write_text(text + f"\nRayOverlapRule {spec}\n")
        rule = SearchParameters.from_config(ConfigFile.from_file(str(f))).ray_rule
        assert rule is not None and rule.arm == arm
        if arm == "R":
            assert rule.radius_px == rad
        assert rule.registration_mm == pytest.approx(0.010)
    # ray rule alone warns (unstudied combination); with CppOverlapRule 1 it does not
    f = tmp_path / "R2.config"
    f.write_text(text + "\nRayOverlapRule R2\n")
    with pytest.warns(UserWarning, match="without CppOverlapRule"):
        SearchParameters.from_config(ConfigFile.from_file(str(f)))
    g = tmp_path / "R2cpp.config"
    g.write_text(text + "\nRayOverlapRule R2\nCppOverlapRule 1\n")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sp = SearchParameters.from_config(ConfigFile.from_file(str(g)))
    assert sp.cpp_overlap_rule and sp.ray_rule is not None
    bad = tmp_path / "bad.config"
    bad.write_text(text + "\nRayOverlapRule banana\n")
    with pytest.raises(ValueError):
        ConfigFile.from_file(str(bad))


def test_rule_only_in_coarse_levels():
    """Two costs by design: levels 0-2 get the ray rule, level 3 (and FindOptimal) do not."""
    from icenine.reconstructor import AdaptiveVoxelReconstructor

    rec = object.__new__(AdaptiveVoxelReconstructor)
    rec.setup = SimpleNamespace(search_params=SearchParameters(ray_rule=RayRule(arm="R")))
    assert [rec._ray_rule_for_level(k) is not None for k in range(5)] == [
        True,
        True,
        True,
        False,
        False,
    ]
    rec.setup = SimpleNamespace(search_params=SearchParameters())
    assert all(rec._ray_rule_for_level(k) is None for k in range(5))
    assert rec._cost_flags() == {}
    assert rec._cost_flags(RayRule()) == dict(cpp_overlap_rule=True, center_screen=True)


@needs_ext
@pytest.mark.parametrize(
    "rule",
    [RayRule(arm="R", radius_px=r) for r in (0, 1, 3)] + [RayRule(arm="A", registration_mm=0.5)],
)
@pytest.mark.parametrize("pixel_radius", [0, 2])
def test_c_ray_flags_equal_python_reference(rule, pixel_radius):
    """The C ray check gives the same flags as the Python reference (ray_overlap.ray_flags)."""
    from icenine.ray_overlap import ray_flags_c

    rng = np.random.default_rng(5)
    n, m, n_det = 40, 400, 3
    dets = [_fake_detector(10.0 * (1 + d)) for d in range(n_det)]
    imgs = [(rng.random((n, n)) < 0.25).astype(np.uint8) for _ in range(n_det)]
    px = np.zeros((m, 3, 2))
    px[:, :, 0] = rng.uniform(-3, 42, (m, 1)) + rng.uniform(0, 1.5, (m, 3))
    px[:, :, 1] = rng.uniform(-3, 42, (m, 1)) + rng.uniform(0, 1.5, (m, 3))
    origins = np.concatenate([rng.uniform(-0.5, 0.5, (m, 2)), np.zeros((m, 1))], axis=1)
    hit = rng.random(m) < 0.9
    ref = ray_flags(
        rule, dets, _Exp(imgs), np.zeros(m, np.int64), origins, hit, px, pixel_radius, 0.3, 0.035
    )
    got = ray_flags_c(
        rule, dets, imgs, np.zeros(m, np.int32), origins, hit, px, pixel_radius, 0.3, 0.035
    )
    assert ref[0].sum() > 0 and ref[1].sum() > 0
    assert np.array_equal(ref[0], got[0]) and np.array_equal(ref[1], got[1])
