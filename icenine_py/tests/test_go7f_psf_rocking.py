"""Opt-in detector PSF and rocking curve in the forward model (GO7f).

Defaults bit-identical (a fixed small render and fixed predictions against
``tests/data/go7f_baseline.npz``, made before the options existed); mass conservation under PSF
and rocking in ``SpotModel`` and in the stock renderer; the rocking frame weights against the
analytic window width (``docs/omega_peak_width_derivation.md``); the strict-cost rocking window."""

import math
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest
import torch

from go7f_fixture import EX, fingerprint

import icenine.grain_fit as GF
from icenine import cost_functions as CF
from icenine.cost_functions import VoxelCostFunction
from scipy.spatial.transform import Rotation

from test_center_screen import three_voxels  # noqa: F401  (fixture)

DT = torch.float64
BASE = Path(__file__).resolve().parent / "data" / "go7f_baseline.npz"
needs_ex = pytest.mark.skipif(
    not (EX / "SimInput").exists(), reason="Example2.ManyGrains not present"
)
needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")


# ------------------------------------------------------------------ bit-identical defaults
@needs_ex
def test_defaults_bit_identical_to_before() -> None:
    ref = np.load(BASE)
    got = fingerprint()
    assert set(ref.files) == set(got)
    assert len(ref["tri_int"]) > 100 and ref["rows_fpos"].size > 50
    for k in ref.files:
        assert np.array_equal(ref[k], got[k]), k


# ------------------------------------------------------------------ SpotModel
def _model_rows(**kw: Any) -> Tuple[Any, Any, Any]:
    from test_strain_render import _setup

    from icenine.reconstructor import _get_voxel_vertices
    from icenine.smooth_cost import SmoothCost

    cwd = os.getcwd()
    os.chdir(EX)
    try:
        sim, sample, dets, rmap, _n = _setup(3)
        sc = SmoothCost(
            sim.simulator,
            dets,
            rmap,
            None,
            sample,
            sample.get_structure_list(),
            eta_limit=sim.exp_setup.get_eta_limit(),
        )
        mic = sample.get_mic()
        model = GF.SpotModel(sc, phase=int(mic.voxels[0].phase), **kw)
        Rs = torch.as_tensor(np.stack([v.orientation for v in mic.voxels]), dtype=torch.float64)
        Vs = torch.as_tensor(
            np.stack([_get_voxel_vertices(v).numpy() for v in mic.voxels]), dtype=torch.float64
        )
        model.test_RV = (Rs.numpy(), Vs.numpy())  # type: ignore[attr-defined]
        return model, model.predict(Rs, Vs), sc
    finally:
        os.chdir(cwd)


@needs_ex
def test_spot_model_default_flags_and_rows_unchanged() -> None:
    model, rows, _sc = _model_rows()
    assert model.psf_sigma_px == 0.0 and model.rocking_deg == 0.0 and not model.spreading
    assert rows.spread is None
    bary = torch.as_tensor(GF.subtriangle_bary(3))
    pc, pr, fp, m, ok = GF._rows_cloud(rows, bary, 0)
    pc0, pr0, m0, ok0 = GF._rows_points(rows, bary, 0)
    assert torch.equal(pc, pc0) and torch.equal(pr, pr0) and torch.equal(m, m0)
    assert torch.equal(fp, rows.fpos[:, None].expand_as(pc0))


@needs_ex
def test_psf_conserves_mass_and_has_gaussian_moments() -> None:
    sigma = 0.8
    _m, rows0, _sc = _model_rows()
    model, rows, _sc = _model_rows(psf_sigma_px=sigma)
    assert rows.spread is model
    bary = torch.as_tensor(GF.subtriangle_bary(3))
    pc0, pr0, m0, ok0 = GF._rows_points(rows0, bary, 0)
    pc, pr, fp, m, ok = GF._rows_cloud(rows, bary, 0)
    assert pc.shape[1] == 9 * pc0.shape[1] and torch.equal(ok, ok0)
    assert torch.allclose(m.sum(1), m0.sum(1), rtol=1e-12, atol=0)  # mass per row
    k = pc0.shape[1]
    pc3, pr3, m3 = pc.reshape(-1, k, 9), pr.reshape(-1, k, 9), m.reshape(-1, k, 9)
    w = m3 / m3.sum(2, keepdim=True).clamp_min(1e-300)
    dx, dy = pc3 - pc0[:, :, None], pr3 - pr0[:, :, None]
    assert torch.allclose((w * dx).sum(2), torch.zeros((), dtype=DT), atol=1e-12)  # centred
    assert torch.allclose(
        (w * dx**2).sum(2), torch.full((), sigma**2, dtype=DT), atol=1e-12
    )  # variance
    assert torch.allclose((w * dy**2).sum(2), torch.full((), sigma**2, dtype=DT), atol=1e-12)
    assert torch.allclose((w * dx * dy).sum(2), torch.zeros((), dtype=DT), atol=1e-12)  # separable
    assert torch.allclose((w * dx**4).sum(2), torch.full((), 3 * sigma**4, dtype=DT), atol=1e-12)
    assert torch.equal(fp, rows.fpos[:, None, None].expand(-1, k, 9).reshape(len(fp), -1))


def _independent_overlap(lo: float, hi: float, b: int) -> float:
    """Overlap of [lo, hi] with frame bin [b, b + 1] by a fine midpoint sum."""
    xs = lo + (np.arange(200000) + 0.5) * (hi - lo) / 200000
    return float(((xs >= b) & (xs < b + 1)).mean())


@needs_ex
def test_rocking_weights_against_analytic_width() -> None:
    rocking = 1.7
    model, rows, sc = _model_rows(rocking_deg=rocking)
    bins, w = model.rocking_weights(rows.fpos, rows.sin_eta)
    width_deg = math.degrees(float(sc.width))
    live = (rows.frame >= 0) & (bins >= 0).all(1) & (bins < sc.n_bins).all(1)
    assert int(live.sum()) > 20
    n_multi = 0
    for i in torch.nonzero(live)[:, 0].tolist():
        sin_eta = float(rows.sin_eta[i].abs())
        # analytic full window width in frames: rocking / |sin eta| (degrees) / frame width
        full = rocking / max(sin_eta, 0.01) / width_deg
        f = float(rows.fpos[i])
        lo, hi = f - full / 2, f + full / 2
        assert float(w[i].sum()) == pytest.approx(1.0, abs=1e-12)
        for j in range(bins.shape[1]):
            ref = _independent_overlap(lo, hi, int(bins[i, j])) if j % 3 == 0 else None
            if ref is not None:
                assert float(w[i, j]) == pytest.approx(ref, abs=2e-5)
        # the bins span the window: no overlapped frame is missing
        assert int(bins[i].min()) == math.floor(lo) and int(bins[i].max()) >= math.floor(hi)
        n_multi += int((w[i] > 0).sum() > 1)
    assert n_multi > 5  # peaks near the rotation axis do spread over several frames


@needs_ex
def test_rocking_spread_conserves_mass_and_moves_to_frame_centres() -> None:
    _m, rows0, _sc = _model_rows()
    model, rows, sc = _model_rows(psf_sigma_px=0.5, rocking_deg=1.0)
    bary = torch.as_tensor(GF.subtriangle_bary(3))
    pc0, pr0, m0, ok0 = GF._rows_points(rows0, bary, 0)
    pc, pr, fp, m, ok = GF._rows_cloud(rows, bary, 0)
    inside = (rows.frame >= 0) & (rows.fpos >= 2) & (rows.fpos < sc.n_bins - 6)
    assert int(inside.sum()) > 20
    assert torch.allclose(m[inside].sum(1), m0[inside].sum(1), rtol=1e-10)  # PSF x rocking
    live = m > 0
    assert torch.allclose((fp[live] % 1.0), torch.full((), 0.5, dtype=fp.dtype))
    # rocking only: each point's mass is split, positions repeated
    mdl2, rows2, _ = _model_rows(rocking_deg=1.0)
    pc2, pr2, fp2, m2, _ok = GF._rows_cloud(rows2, bary, 0)
    k = pc2.shape[1] // pc0.shape[1]
    assert k >= 2 and torch.allclose(m2[inside].sum(1), m0[inside].sum(1), rtol=1e-10)
    assert torch.equal(pc2.reshape(len(pc2), -1, k)[:, :, 0], pc0)


@needs_ex
def test_pack_predicted_uses_the_spread_points() -> None:
    model0, rows0, _ = _model_rows()
    model, rows, _ = _model_rows(psf_sigma_px=0.7, rocking_deg=0.5)
    bary = torch.as_tensor(GF.subtriangle_bary(3))
    key = (int(rows0.refl[0]), int(rows0.branch[0]), 0)
    z = np.zeros((0, 3))
    region = GF.GrainRegion(
        [0], [], [GF.Spot(key, 0, 2048, 0, 2048, 0, 360, z, np.zeros(0))], {key: 0}
    )
    region.shell_pts, region.shell_mass = [z], [np.zeros(0)]
    x0, a0 = GF.pack_predicted(region, rows0, bary, model0.n_det)
    x1, a1 = GF.pack_predicted(region, rows, bary, model.n_det)
    assert x1.shape[1] > x0.shape[1]
    assert float(a1.sum()) == pytest.approx(float(a0.sum()), rel=1e-6)  # mass conserved


# ------------------------------------------------------------------ stock renderer
@needs_ex
def test_render_rocking_conserves_intensity_per_peak() -> None:
    from test_strain_render import _render, _setup

    cwd = os.getcwd()
    os.chdir(EX)
    try:
        sim, sample, dets, rmap, n_om = _setup(3)
        base = _render(sim, sample, dets, rmap, n_om)
        sim.rocking_deg = 1.5
        on = _render(sim, sample, dets, rmap, n_om)
    finally:
        os.chdir(cwd)

    # sum of intensity per (voxel, reflection, omega, detector) over frames
    def tot(log: List[Any]) -> Dict[Any, Tuple[float, int]]:
        out: Dict[Any, Tuple[float, int]] = {}
        for (v, r, _oi, om), det, _tri, inten in log:
            s, n = out.get((v, r, round(om, 9), det), (0.0, 0))
            out[(v, r, round(om, 9), det)] = (s + inten, n + 1)
        return out

    b, o = tot(base), tot(on)
    assert set(b) == set(o)
    n_split = 0
    for k, (s, n) in b.items():
        s2, n2 = o[k]
        assert n == 1 and s2 <= s * (1 + 1e-9)
        n_split += n2 > 1
    full = [k for k in b if o[k][0] == pytest.approx(b[k][0], rel=1e-9)]
    assert len(full) > 0.8 * len(b)  # the rest lost a share to frames outside the scan
    assert n_split > 5
    # triangles drawn are those of the default render
    keys_b = {(x[0][0], x[0][1], round(x[0][3], 9), x[1]) for x in base}
    keys_o = {(x[0][0], x[0][1], round(x[0][3], 9), x[1]) for x in on}
    assert keys_b == keys_o


def test_psf_blur_conserves_mass_and_thresholds() -> None:
    from icenine.forward_simulation import ForwardSimulation

    img = np.zeros((64, 64), np.float32)
    img[30:34, 28:36] = 5.0
    out = ForwardSimulation.psf_blur(img, 1.0, threshold=0.0)
    assert out.dtype == img.dtype and float(out.sum()) == pytest.approx(float(img.sum()), rel=1e-4)
    assert (out > 0).sum() > (img > 0).sum()  # the footprint grows
    cut = ForwardSimulation.psf_blur(img, 1.0, threshold=0.05)
    assert (cut > 0).sum() < (out > 0).sum() and cut.min() >= 0
    assert cut[cut > 0].min() >= 0.05 * cut.max() * (1 - 1e-6)
    # a delta blurs to a Gaussian of the requested sigma
    d = np.zeros((65, 65), np.float64)
    d[32, 32] = 1.0
    g = ForwardSimulation.psf_blur(d, 2.0, threshold=0.0)
    xs = np.arange(65) - 32
    var = float((g.sum(0) * xs**2).sum() / g.sum())
    assert var == pytest.approx(4.0, rel=2e-2)


def test_apply_psf_blurs_dense_frames_and_skips_sinks() -> None:
    from icenine.forward_simulation import ForwardSimulation
    from icenine.image_data import ImageData

    sim = ForwardSimulation.__new__(ForwardSimulation)
    sim.psf_sigma_px, sim.psf_threshold = 1.0, 0.0
    a, b = ImageData(32, 32), ImageData(32, 32)
    a._pixels_dense[10, 10] = 4.0
    images: List[List[Any]] = [[a, object()], [b, None]]
    sim._apply_psf(images)  # type: ignore[arg-type]
    assert float(a._pixels_dense.sum()) == pytest.approx(4.0, rel=1e-4)
    assert int((a._pixels_dense > 0).sum()) > 1 and float(b._pixels_dense.sum()) == 0.0


def test_render_options_default_off_and_serial_path_refuses() -> None:
    from icenine.forward_simulation import ForwardSimulation

    sim = ForwardSimulation.__new__(ForwardSimulation)
    sim.rocking_deg, sim.psf_sigma_px = 0.5, 0.0
    with pytest.raises(NotImplementedError):
        sim._simulate_peaks([], [], None, None)  # type: ignore[arg-type]


# ------------------------------------------------------------------ strict cost
def test_rocking_window_frames_analytic() -> None:
    from icenine.cost_functions import rocking_window_frames

    n = 20
    frame_index = np.concatenate([np.arange(n), [-1]])
    width = math.radians(1.0)
    low = 0.0
    # eta = 90 deg (ray along y, sin eta = 1): window +- rocking / 2
    refl = np.array([[0.0, 1.0, 0.0]] * 4)
    om = np.radians([5.5, 5.0, 5.9, 0.2])
    # 1.0 deg rocking: +- 0.5 frames: 5.5 -> [5.0, 6.0] -> frames 5, 6; 5.0 -> [4.5, 5.5] -> 4, 5;
    # 5.9 -> [5.4, 6.4] -> 5, 6; 0.2 -> [-0.3, 0.7] -> 0 (frame -1 is outside the scan)
    got = rocking_window_frames(om, refl, 1.0, low, width, frame_index)
    assert got[1] == (4, 5) and got[2] == (5, 6) and got[3] == (0,)
    assert got[0] in ((5, 6), (5,))  # window edge exactly on the bin boundary: both are exact
    # near the rotation axis the window widens by 1 / |sin eta|: sin eta = 0.25 -> +- 2 frames
    refl2 = np.array([[0.0, 0.25, math.sqrt(1 - 0.0625)]])
    got2 = rocking_window_frames(np.radians([10.5]), refl2, 1.0, low, width, frame_index)
    assert got2[0] == (8, 9, 10, 11, 12)
    # beyond the scan edge frames are dropped
    got3 = rocking_window_frames(np.radians([19.9]), refl, 1.0, low, width, frame_index)
    assert got3[0] == (19,)


@needs_c
def test_strict_cost_rocking_window(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    from scipy.spatial.transform import Rotation

    c = three_voxels
    phase, verts = c["voxel"].phase, c["vertices"]
    R0 = np.asarray(c["voxel"].orientation, dtype=np.float64)
    width = math.degrees(float(c["kwargs"]["range_map"].width))
    cf0 = VoxelCostFunction(pixel_radius=0, mode="hard", **c["kwargs"])
    assert cf0.rocking_deg == 0.0
    tiny = VoxelCostFunction(pixel_radius=0, mode="hard", rocking_deg=1e-9, **c["kwargs"])
    wide = VoxelCostFunction(pixel_radius=0, mode="hard", rocking_deg=1.5 * width, **c["kwargs"])
    # a window of one frame is the default cost; a wider window only adds lit frames
    rng = np.random.default_rng(0)
    rots = [R0]
    for deg in (0.05, 0.2, 0.6):
        for _ in range(3):
            ax = rng.normal(size=3)
            ax /= np.linalg.norm(ax)
            rots.append(Rotation.from_rotvec(ax * np.radians(deg)).as_matrix() @ R0)
    for dz in (0.7, -0.7):  # a turn about the rotation axis shifts every spot in omega
        rots.append(Rotation.from_rotvec([0, 0, np.radians(dz * width)]).as_matrix() @ R0)
    gain = 0
    for R in rots:
        a = cf0.evaluate(R, verts, phase)
        t = tiny.evaluate(R, verts, phase)
        w = wide.evaluate(R, verts, phase)
        assert (t.pixel_overlap, t.peak_overlap, t.quality, t.cost) == (
            a.pixel_overlap,
            a.peak_overlap,
            a.quality,
            a.cost,
        )
        assert w.peak_on_detector == a.peak_on_detector
        assert w.pixel_overlap >= a.pixel_overlap and w.peak_overlap >= a.peak_overlap
        gain += int(w.pixel_overlap > a.pixel_overlap)
    assert gain >= 2  # the shifted orientations score in the window, not in the strict frame
    shifted = rots[-1], rots[-2]
    for R in shifted:
        a, w = cf0.evaluate(R, verts, phase), wide.evaluate(R, verts, phase)
        assert w.pixel_overlap > a.pixel_overlap


def test_strict_cost_rocking_rules() -> None:
    from icenine import batched_cost as BC

    with pytest.raises(ValueError):
        VoxelCostFunction(
            None, [], None, None, None, [], rocking_deg=-1.0  # type: ignore[arg-type]
        )
    cf = VoxelCostFunction(
        None,  # type: ignore[arg-type]
        [],
        None,
        None,
        None,
        [],
        rocking_deg=0.3,
        eta_limit=math.pi / 2,
    )
    assert cf.rocking_deg == 0.3
    # the batched evaluator and the fused path do not implement the window: they refuse it
    assert "rocking" in str(BC.unsupported_reason(cf)) or not CF._HAS_C_RASTERIZE
    if CF._C_HAS_FUSED:
        assert "rocking" in str(CF.fused_unsupported_reason(cf))
    # soft mode has no window
    with pytest.raises(ValueError):
        CF.calculate_diffraction_overlap_batched(
            None,
            None,
            [0.1],
            [torch.zeros(3)],
            [],
            None,
            None,
            None,  # type: ignore[arg-type]
            mode="soft",
            rocking_deg=0.3,
        )


@needs_c
def test_rocking_set_after_fused_evaluate_is_honoured(
    three_voxels: Dict[str, Any],  # noqa: F811
) -> None:
    """Setting ``rocking_deg`` after the fused context exists must rebuild it (the fused path has
    no window), so the cost equals that of a cost function built with the window."""
    if not CF._C_HAS_FUSED:
        pytest.skip("fused evaluate not built")
    c = three_voxels
    phase, verts = c["voxel"].phase, c["vertices"]
    R0 = np.asarray(c["voxel"].orientation, dtype=np.float64)
    width = math.degrees(float(c["kwargs"]["range_map"].width))
    R = Rotation.from_rotvec([0, 0, np.radians(0.7 * width)]).as_matrix() @ R0
    kw = dict(pixel_radius=0, mode="hard", eta_limit=math.pi / 2, **c["kwargs"])
    ref = VoxelCostFunction(rocking_deg=1.5 * width, **kw)
    base = VoxelCostFunction(**kw).evaluate(R, verts, phase)
    cf = VoxelCostFunction(fused_cost=True, **kw)
    assert cf.evaluate(R, verts, phase).pixel_overlap == base.pixel_overlap
    cf.rocking_deg = 1.5 * width
    with pytest.warns(RuntimeWarning, match="rocking"):
        got = cf.evaluate(R, verts, phase)
    want = ref.evaluate(R, verts, phase)
    assert got.pixel_overlap == want.pixel_overlap > base.pixel_overlap


def _tiny_table(n: int) -> "CF.FrameTable":
    table = CF.FrameTable.__new__(CF.FrameTable)
    table.unions = CF.OrderedDict()
    table._lit = {}
    table.n_det = 1
    table.n_frames = n
    table.rows = table.cols = 4
    rng = np.random.default_rng(3)
    table.images = [(rng.random((4, 4)) < 0.2).astype(np.uint8) for _ in range(n)]
    table._filled = np.ones(n, dtype=bool)
    return table


def test_frame_table_unions_are_compact_exact_and_hold_a_large_working_set() -> None:
    table = _tiny_table(300)
    keys = [(i, i + 1, i + 2) for i in range(297)]  # more than the old 64-key cap
    for k in keys:
        got = table.union_images(k)[0]
        want = np.maximum(np.maximum(table.images[k[0]], table.images[k[1]]), table.images[k[2]])
        assert got.dtype == np.uint8 and np.array_equal(got, want)
    assert len(table.unions) == len(keys)  # nothing evicted: the working set is held
    stored = table.unions[keys[5]][0]
    assert stored.ndim == 1 and stored.dtype == np.int32  # lit-pixel indices, not a dense image
    again = table.union_images(keys[5])[0]  # a hit: same cached indices, a fresh dense image
    assert table.unions[keys[5]][0] is stored
    again[:] = 7  # the caller owns the dense image; the cache is unaffected
    assert np.array_equal(table.union_images(keys[5])[0], _dense(stored))


def _dense(idx: np.ndarray) -> np.ndarray:
    out = np.zeros(16, dtype=np.uint8)
    out[idx] = 1
    return out.reshape(4, 4)


def test_frame_table_union_cache_lru_keeps_recent_key() -> None:
    table = _tiny_table(40)
    old = CF.FrameTable.MAX_UNIONS
    CF.FrameTable.MAX_UNIONS = 5
    try:
        for i in range(5):
            table.union_images((i, i + 1))
        table.union_images((0, 1))  # touch the oldest
        table.union_images((10, 11))  # evicts the least recently used, (1, 2)
        assert (0, 1) in table.unions and (1, 2) not in table.unions
        assert len(table.unions) == 5
    finally:
        CF.FrameTable.MAX_UNIONS = old


def test_touching_components_frame_dilation() -> None:
    pc = torch.full((1, 1), 5.5)
    pr = torch.full((1, 1), 5.5)
    fr = torch.tensor([10])
    obs = np.array([[5.5, 5.5, 12.5]])  # (col, row, frame + 0.5): two frames from home
    box = (6, 15, 0, 12, 0, 12)
    assert not GF._touching_components(obs, box, pc, pr, fr, 3).any()  # default: one frame
    assert GF._touching_components(obs, box, pc, pr, fr, 3, 3).all()
    obs1 = np.array([[5.5, 5.5, 11.5]])
    assert GF._touching_components(obs1, box, pc, pr, fr, 3).all()


@needs_ex
def test_build_region_keeps_pixels_two_frames_from_home_with_rocking() -> None:
    model0, rows0, _sc = _model_rows()
    width_deg = math.degrees(model0.sc.width)
    model, rows, _sc2 = _model_rows(rocking_deg=4.0 * width_deg)
    R, V = model.test_RV  # type: ignore[attr-defined]
    bary = torch.as_tensor(GF.subtriangle_bary(3))
    pc, pr, _m, ok = GF._rows_points(rows0, bary, 0)
    cand = torch.nonzero(ok & (rows0.frame >= 4) & (rows0.frame < model0.sc.n_bins - 8))[:, 0]
    i = int(cand[0])
    col, row, home = float(pc[i].mean()), float(pr[i].mean()), int(rows0.frame[i])
    H = W = int(max(col, row)) + 40
    n_om = model0.sc.n_bins
    stack = np.zeros((model0.n_det * n_om, H, W), dtype=np.uint8)
    stack[home + 2, int(row), int(col)] = 1  # detector 0, two frames from home
    src = GF.ImageSource(stack, n_om)
    members, shell = [int(rows0.vox[i])], []
    want = (int(rows0.refl[i]), int(rows0.branch[i]), 0)

    def run(m: Any) -> int:
        reg = GF.build_region(m, src, R, V, members, shell, bary, ratio=None)
        return sum(len(s.obs) for s in reg.spots if s.key == want)

    assert run(model0) == 0  # no spread: the pixel is outside the one-frame footprint
    assert run(model) >= 1  # rocking half-window >= 2 frames: kept
