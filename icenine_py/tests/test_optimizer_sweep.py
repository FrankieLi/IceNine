"""scripts/optimizer_sweep.py: pixel-set images, the realism edit, the coarse image stack,
case alignment with the network sweep, and (physics, skipped without the ManyGrains example) the
per-voxel images: a target-only image at the truth gives perfect overlap, the realistic image
contains the distractor pixels."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).parent.parent / "benchmarks"))
import optimizer_sweep as osw  # noqa: E402
import perturbation_sweep as ps  # noqa: E402

PROJECT = Path(__file__).parent.parent.parent
SWEEP_RAW = (
    Path(__file__).parent.parent / "benchmarks/toy_orientation_sweep/perturbation_sweep_raw.npz"
)
GEO = osw.Geometry(n_det=2, n_omega=6, H=32, W=48)


def test_pixel_keys_round_trip_and_grouping():
    rng = np.random.default_rng(0)
    f = rng.integers(0, GEO.n_omega, 200)
    d = rng.integers(0, GEO.n_det, 200)
    r = rng.integers(0, GEO.H, 200)
    c = rng.integers(0, GEO.W, 200)
    k = osw.encode_pixels(f, d, r, c, GEO)
    for a, b in zip(osw.decode_pixels(k, GEO), (f, d, r, c)):
        assert (a == b).all()
    keys = np.unique(k)
    groups = osw.group_pixels(keys, GEO)
    assert sum(len(v) for v in groups.values()) == len(keys)
    fi, di, ri, ci = osw.decode_pixels(keys, GEO)
    for img, pos in groups.items():
        sel = fi * GEO.n_det + di == img
        assert sorted(pos.tolist()) == sorted((ri[sel] * GEO.W + ci[sel]).tolist())


def _spec(col0, row0, frame0):
    return SimpleNamespace(col0=np.array(col0), row0=np.array(row0), frame0=np.array(frame0))


def test_window_pixels_land_at_detector_position_and_frame():
    w = np.zeros((2, 8, 8), dtype=np.uint8)
    w[0, 2, 3] = 3  # code 3 with K = 2 -> frame0 + 0
    w[1, 0, 0] = 1  # frame0 - 2
    spec = _spec([10, 20], [4, 5], [3, 4])
    keys = osw.window_pixel_keys(w, spec, np.array([0, 1]), 2, GEO)
    f, d, r, c = osw.decode_pixels(keys, GEO)
    got = sorted(zip(f.tolist(), d.tolist(), r.tolist(), c.tolist()))
    assert got == sorted([(3, 0, 6, 13), (2, 1, 5, 20)])


def test_realism_edit_reproduces_the_realistic_windows():
    rng = np.random.default_rng(1)
    clean = (rng.random((3, 8, 8)) < 0.15).astype(np.uint8) * 3
    dis = (rng.random((3, 8, 8)) < 0.15).astype(np.uint8) * 2
    real = np.where(rng.random((3, 8, 8)) < 0.2, 0, np.where(clean > 0, clean, dis))
    real[2, 7, 7] = 4  # an added (hot) pixel
    real = real.astype(np.uint8)
    spec = _spec([0, 8, 16], [0, 8, 16], [3, 3, 3])
    det = np.array([0, 0, 1])
    rem, add = osw.realism_edit(clean, dis, real, spec, det, 2, GEO)
    before = osw.window_pixel_keys(np.where(clean > 0, clean, dis), spec, det, 2, GEO)
    after = osw.window_pixel_keys(real, spec, det, 2, GEO)
    assert (np.union1d(np.setdiff1d(before, rem), add) == after).all()
    assert len(np.intersect1d(rem, after)) == 0 and len(add) >= 1


def test_coarse_stack_equals_the_multiscale_stack():
    from icenine.differentiable_cost import MultiScaleImageStack, SparseImageStack

    geo = osw.Geometry(2, 6, 32, 32)
    rng = np.random.default_rng(2)
    keys = np.unique(
        osw.encode_pixels(
            rng.integers(0, 6, 150), rng.integers(0, 2, 150), rng.integers(0, 32, 150),
            rng.integers(0, 32, 150), geo,
        )
    )  # fmt: skip
    groups = osw.group_pixels(keys, geo)
    sp = SparseImageStack.__new__(SparseImageStack)
    sp.n_omega, sp.n_det, sp.H, sp.W = 6, 2, 32, 32
    sp.binary, sp.dtype = True, torch.float32
    sp._pixel_coords, sp._pixel_values = [], []
    for i in range(12):
        pos = groups.get(i)
        coords = (
            torch.empty(0, 2, dtype=torch.int16)
            if pos is None
            else torch.from_numpy(np.stack([pos // 32, pos % 32], 1)).to(torch.int16)
        )
        sp._pixel_coords.append(coords)
        sp._pixel_values.append(None)
    ref = MultiScaleImageStack(sp, [1, 8], omega_window=1).get_at_scale(1)
    mine = osw.CoarseStack(groups, geo, factor=8, omega_window=1, scale_index=1)
    idx = torch.arange(12)
    assert torch.equal(ref.get_images_batch(idx), mine.get_images_batch(idx))


def test_pixel_set_data_gives_binary_images():
    zeros = np.zeros((GEO.H, GEO.W), dtype=np.uint8)
    data = osw.PixelSetData({2 * 1 + 1: np.array([5, 70])}, GEO, zeros)
    img = data.get_image(1, 1)
    assert img.num_rows == GEO.H and img.get_binary_numpy().sum() == 2
    assert data.get_image(0, 0).get_binary_numpy().sum() == 0


def test_case_draws_are_seeded_and_have_the_exact_radius():
    a = ps.case_draws(0, 7336, 3, 0.5, 20, 3, 0.17, 0.5)
    b = ps.case_draws(0, 7336, 3, 0.5, 20, 3, 0.17, 0.5)
    assert (a[0] == b[0]).all() and (a[1][4][0] == b[1][4][0]).all()
    assert np.allclose(np.linalg.norm(a[0], axis=1), 0.5)
    assert a[1][0][0].shape == (3, 3) and a[1][0][1].shape == (3,)


# --- physics (ManyGrains example) ------------------------------------------------------------


@pytest.fixture(scope="module")
def physics(tmp_path_factory):
    example = PROJECT / "Examples" / "Example2.ManyGrains"
    if not (example / "ConfigFiles").exists() or not SWEEP_RAW.exists():
        pytest.skip("ManyGrains example or the sweep's raw results not found")
    cwd = Path.cwd()
    cfg, sweep = osw.load_sweep_config(SWEEP_RAW)
    cfg.update(methods=osw.METHODS, opt_dirs=1, huber_c=1.0)
    try:
        osw.init_worker(cfg)
        yield osw._CTX, sweep
    finally:
        import os

        os.chdir(cwd)


def test_target_only_image_at_truth_has_perfect_overlap(physics):
    ctx, sweep = physics
    vctx = osw.voxel_context(ctx, int(sweep["voxel_indices"][0]))
    groups = osw.group_pixels(vctx.target_keys, ctx.geo)
    ctx.hard_fn.exp_data = osw.PixelSetData(groups, ctx.geo, ctx.zeros)
    R = vctx.R_true.astype(np.float32)
    info = ctx.hard_fn.evaluate(R, vctx.vertices, vctx.voxel.phase)
    assert 1.0 - info.cost == pytest.approx(1.0) and info.peak_on_detector > 20
    # and the cost drops fast away from the truth
    from scipy.spatial.transform import Rotation

    Rp = (Rotation.from_rotvec([0.0, 0.0, np.radians(0.5)]).as_matrix() @ vctx.R_true).astype(
        np.float32
    )
    assert 1.0 - ctx.hard_fn.evaluate(Rp, vctx.vertices, vctx.voxel.phase).cost < 0.5


def test_realistic_image_contains_the_distractors_and_edits(physics):
    ctx, sweep = physics
    a = ctx.args
    v = int(sweep["voxel_indices"][0])
    vctx = osw.voxel_context(ctx, v)
    D, ri = 20, 5
    delta0, draws = ps.case_draws(
        a.sweep_seed, v, ri, a.radii[ri], D, len(vctx.sources), a.neighbor_sigma_deg / np.sqrt(3.0),
        a.neighbor_p,
    )  # fmt: skip
    R_nom = ps.perturbed_nominal(vctx.R_true, delta0)
    preps = [ps.prepare_nominal(ctx, v, R_nom[j])[0] for j in range(D)]
    # alignment with the network sweep: same ROI counts per case
    assert (np.array([p.n for p in preps]) == sweep["n_roi"][0, ri]).all()
    layers = {}
    seed = a.realism_seed + 1009 * ri
    b = ps.render_batch(preps, delta0, draws, vctx.sources, "all", seed, a, layers=layers)
    j = next(j for j in range(D) if draws[j][1].any())
    p = preps[j]
    det = p.obs.det_idx.numpy()
    edit = osw.realism_edit(
        layers["clean"][j, : p.n].numpy(), layers["dis"][j, : p.n].numpy(),
        b["windows"][j, : p.n].numpy(), p.spec, det, a.frame_half_width, ctx.geo,
    )  # fmt: skip
    img = osw.case_image_keys(ctx, vctx, "all", draws[j], edit)
    clean_img = osw.case_image_keys(ctx, vctx, "clean", draws[j], None)
    unedited = osw.case_image_keys(ctx, vctx, "all", draws[j], (edit[0][:0], edit[1][:0]))
    # distractor spots are on the detector beyond the target's own pixels
    assert len(np.setdiff1d(unedited, clean_img)) > 0
    # every pixel the net's realistic windows hold is in the image, and removed ones are gone
    real = osw.window_pixel_keys(
        b["windows"][j, : p.n].numpy(), p.spec, det, a.frame_half_width, ctx.geo
    )
    assert len(np.setdiff1d(real, img)) == 0
    assert len(np.intersect1d(edit[0], img)) == 0 and len(edit[0]) + len(edit[1]) > 0


def test_min_sin_eta_filter_only_removes_peaks(physics):
    from icenine.cost_functions import VoxelCostFunction

    ctx, sweep = physics
    vctx = osw.voxel_context(ctx, int(sweep["voxel_indices"][0]))
    h = ctx.hard_fn
    kw = dict(
        simulator=h.simulator, detector_list=h.detector_list, range_map=h.range_map,
        exp_data=h.exp_data, sample=h.sample, structure_list=h.structure_list,
        eta_limit=h.eta_limit, max_q=8.0,
    )  # fmt: skip
    R = vctx.R_true.astype(np.float32)
    n0 = VoxelCostFunction(**kw).evaluate(R, vctx.vertices, vctx.voxel.phase).n_quality_points
    n1 = VoxelCostFunction(min_sin_eta=0.3, **kw).evaluate(R, vctx.vertices, vctx.voxel.phase)
    assert 0 < n1.n_quality_points <= n0
