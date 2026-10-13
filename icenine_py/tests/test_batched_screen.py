"""Opt-in batched / interpolated level screen (BatchedScreen, CppOmegaInterp): the batched exact
screen equals the per-orientation screen, the interpolation follows the C++ formula, the default
discrete search is unchanged and the config keys reach the search."""

from pathlib import Path
from typing import Any, Dict

import numpy as np
import pytest
import torch
from scipy.spatial.transform import Rotation

from icenine import cost_functions as CF
from icenine.batched_screen import ScreenEvaluator, make_screen_evaluator
from icenine.cost_functions import VoxelCostFunction

from test_center_screen import three_voxels  # noqa: F401  (fixture)

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")
CONFIGS = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"


def _cost(c: Dict[str, Any], **kw: Any) -> VoxelCostFunction:
    return VoxelCostFunction(pixel_radius=3, max_q=5.0, **kw, **c["kwargs"])


def _orients(truth: np.ndarray, n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = [truth]
    for deg in (0.2, 1.0, 3.0, 8.0):
        for _ in range(n):
            ax = rng.normal(size=3)
            ax /= np.linalg.norm(ax)
            out.append(Rotation.from_rotvec(ax * np.radians(deg)).as_matrix() @ truth)
    out += [Rotation.random(random_state=int(s)).as_matrix() for s in rng.integers(1 << 30, size=n)]
    return np.array(out)


@needs_c
@pytest.mark.parametrize("flags", [{}, {"cpp_overlap_rule": True, "cpp_edge_rules": True}])
def test_batched_exact_equals_per_orientation(
    three_voxels: Dict[str, Any], flags: Dict[str, Any]
) -> None:  # noqa: F811
    c = three_voxels
    cf = _cost(c, **flags)
    ev = ScreenEvaluator(cf)
    orients = _orients(np.asarray(c["voxel"].orientation, dtype=np.float64), 10)
    po, pon = ev.overlaps_exact(orients, c["vertices"], c["voxel"].phase)
    for i, R in enumerate(orients):
        info = cf.evaluate(R, c["vertices"], c["voxel"].phase)
        assert (po[i], pon[i]) == (info.peak_overlap, info.peak_on_detector), i
    assert po.max() > 0  # the test is not vacuous: the truth is bright


@needs_c
def test_interp_with_identity_local_equals_exact(
    three_voxels: Dict[str, Any],
) -> None:  # noqa: F811
    c = three_voxels
    ev = ScreenEvaluator(_cost(c))
    orients = _orients(np.asarray(c["voxel"].orientation, dtype=np.float64), 6)
    po_i, pon_i = ev.overlaps_interp(orients, np.eye(3)[None], c["vertices"], c["voxel"].phase)
    po_e, pon_e = ev.overlaps_exact(orients, c["vertices"], c["voxel"].phase)
    assert np.array_equal(po_i[:, 0], po_e) and np.array_equal(pon_i[:, 0], pon_e)


def test_interp_omegas_follow_the_cpp_formula(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    """Independent float64 transcription of PrecomputeScatteringVector and
    InterpolateObservablePeaks."""
    c = three_voxels
    cf = _cost(c)
    ev = ScreenEvaluator(cf)
    phase = c["voxel"].phase
    centre = np.asarray(c["voxel"].orientation, dtype=np.float64)
    local = Rotation.from_rotvec(np.radians(2.0) * np.array([0.3, -0.5, 0.8]) / 0.99).as_matrix()
    om, nrm, seg = ev.interp_peaks(centre[None], local[None], phase)
    assert len(om) > 20 and (seg == 0).all()
    g_hkl, g_mag = cf._phase_recip_vecs[phase]
    g = centre @ g_hkl.numpy().astype(np.float64).T  # (3, n)
    E = cf.simulator.beam_energy
    k = 0.506773182 * E
    ref = []
    for r in range(g.shape[1]):
        gv, mag = g[:, r], float(g_mag[r])
        sin_t, cos_c = mag / (2 * k), gv[2] / mag
        sin_c = np.sqrt(1 - cos_c**2)
        if abs(sin_t) > abs(sin_c):
            continue
        b1 = np.arcsin(sin_t / sin_c)
        d0 = np.arctan2(gv[0], gv[1])
        d = gv / np.linalg.norm(gv)
        for om0 in (b1 + d0, np.pi - b1 + d0):
            om0 = (om0 + np.pi) % (2 * np.pi) - np.pi
            rz = np.array(
                [[np.cos(om0), -np.sin(om0), 0], [np.sin(om0), np.cos(om0), 0], [0, 0, 1]]
            )
            x_lab = (rz @ d)[0]
            dn = local @ d
            d_om = ((rz[0] @ dn) - x_lab) / (rz[1] @ dn)
            if abs(d_om) <= 1:
                ref.append((om0 + d_om, dn))
    assert len(ref) == len(om)
    got = sorted(zip(om.numpy().tolist(), map(tuple, nrm.numpy().tolist())))
    want = sorted((o, tuple(n)) for o, n in ref)
    # tolerance: float32 arithmetic against a float64 reference
    assert np.allclose([a[0] for a in got], [a[0] for a in want], atol=2e-5)


@needs_c
def test_discrete_search_default_unchanged_and_batched_equal(
    three_voxels: Dict[str, Any],
) -> None:  # noqa: F811
    from icenine.orientation_search import run_discrete_search_spaced
    from icenine.sampling import generate_local_grid_multi_level

    c = three_voxels
    g_cost = _cost(c)
    l_cost = VoxelCostFunction(pixel_radius=0, **c["kwargs"])
    fz = _orients(np.asarray(c["voxel"].orientation, dtype=np.float64), 8)[:24]
    grid = generate_local_grid_multi_level(np.radians(5.0), 0, 1)
    common = dict(
        local_cost_fn=l_cost,
        fz_orientations=fz,
        local_grid=grid,
        voxel_vertices=c["vertices"],
        angular_radius=np.radians(5.0),
        symmetry_quats=np.zeros((0, 4)),
        phase_index=c["voxel"].phase,
    )
    base = run_discrete_search_spaced(global_cost_fn=_cost(c), **common)
    batched = run_discrete_search_spaced(
        global_cost_fn=g_cost, screen=ScreenEvaluator(g_cost), screen_chunk=5, **common
    )
    assert len(base) == len(batched) > 0
    assert all(np.array_equal(a.orientation, b.orientation) for a, b in zip(base, batched))
    assert [a.cost for a in base] == [b.cost for b in batched]
    assert g_cost.eval_count == len(fz) * len(grid)  # the batched screen counts its orientations


def test_config_keys_and_default_off(tmp_path: Path) -> None:
    from icenine.config_file import ConfigFile
    from icenine.orientation_search import SearchParameters

    text = (CONFIGS / "ReconstructPhaseD_mc_clean.config").read_text()
    off, on = tmp_path / "off.config", tmp_path / "on.config"
    off.write_text(text)
    on.write_text(text + "\nCppOmegaInterp 1\nBatchedScreen 1\n")
    for f, want in ((off, False), (on, True)):
        cfg = ConfigFile.from_file(str(f))
        sp = SearchParameters.from_config(cfg)
        assert cfg.cpp_omega_interp is want and cfg.batched_screen is want
        assert sp.cpp_omega_interp is want and sp.batched_screen is want
    sp0 = SearchParameters()
    assert sp0.cpp_omega_interp is False and sp0.batched_screen is False


@needs_c
def test_screen_kwargs_empty_when_off_and_falls_back_with_ray_rule(
    three_voxels: Dict[str, Any],
) -> None:  # noqa: F811
    from icenine.ray_overlap import RayRule
    from icenine.reconstructor import AdaptiveVoxelReconstructor

    c = three_voxels
    cf = _cost(c)

    class _Setup:
        search_params = type(
            "SP", (), {"cpp_omega_interp": False, "batched_screen": False, "screen_chunk": 7}
        )()

    rec = AdaptiveVoxelReconstructor.__new__(AdaptiveVoxelReconstructor)
    rec.setup = _Setup()  # type: ignore[assignment]
    assert rec._screen_kwargs(cf) == {}
    rec.setup.search_params.batched_screen = True  # type: ignore[attr-defined]
    kw = rec._screen_kwargs(cf)
    assert kw["screen_chunk"] == 7 and kw["screen_interp"] is False and kw["screen"] is not None
    rec.setup.search_params.cpp_omega_interp = True  # type: ignore[attr-defined]
    rec.setup.search_params.batched_screen = False  # type: ignore[attr-defined]
    kw = rec._screen_kwargs(cf)
    assert kw["screen_chunk"] == 1 and kw["screen_interp"] is True
    ray_cf = VoxelCostFunction(pixel_radius=3, max_q=5.0, ray_rule=RayRule(), **c["kwargs"])
    assert make_screen_evaluator(ray_cf) is None
    with pytest.warns(RuntimeWarning, match="not supported"):
        assert rec._screen_kwargs(ray_cf) == {}
