"""
Tests for the Stage 0 evaluation tools (orientation_eval.py) and the offset-head
loss/model helpers in orientation_nn.py / toy_orientation_model.py.

Physics tests use Example2.ThreeVoxels and are skipped if its files are absent.
"""

import os
from pathlib import Path

import numpy as np
import pytest
import torch
from scipy.spatial.transform import Rotation

from icenine.config_file import ConfigFile
from icenine.experiment_setup import XDMExperimentSetup
from icenine.mic_file import MicFile
from icenine.orientation_eval import (
    DEG,
    BatchedObserver,
    ExactBayes,
    error_summary,
    lit_pixel_set,
    offsets_to_matrices,
    offsets_to_quaternions,
    quaternions_to_offsets_deg,
    rotvec_to_matrix,
    sample_fixed_magnitude_offsets,
    sample_prior_offsets,
)
from icenine.orientation_nn import (
    _project_all_detectors,
    _project_peak_on_detector,
    _restore_and_rotate,
    cholesky_from_raw,
    decoupled_nll_loss,
    define_roi_set,
    gaussian_nll_loss,
    mse_deg_loss,
    spot_overlaps_grid,
)
from icenine.image_data import ImageData
from icenine.peak_filters import XDMEtaAcceptFn
from icenine.reconstructor import _get_voxel_vertices
from icenine.sample import Sample
from icenine.simulation import Simulation
from icenine.toy_orientation_model import ToyOffsetNet


# ============================================================================
# Pure-math tests (no example data)
# ============================================================================


class TestRotations:
    def test_rodrigues_matches_scipy(self):
        rng = np.random.default_rng(0)
        v = rng.normal(size=(20, 3)) * 0.3
        got = rotvec_to_matrix(torch.from_numpy(v)).numpy()
        assert np.allclose(got, Rotation.from_rotvec(v).as_matrix(), atol=1e-12)

    def test_rodrigues_zero_vector_is_identity(self):
        got = rotvec_to_matrix(torch.zeros(2, 3, dtype=torch.float64)).numpy()
        assert np.allclose(got, np.eye(3))

    def test_offset_quaternion_roundtrip(self):
        rng = np.random.default_rng(1)
        R_nom = Rotation.random(random_state=3).as_matrix()
        offsets = rng.normal(size=(15, 3))
        q = offsets_to_quaternions(offsets, R_nom)
        assert (q[:, 0] >= 0).all()
        assert np.allclose(np.linalg.norm(q, axis=1), 1.0)
        assert np.allclose(quaternions_to_offsets_deg(q, R_nom), offsets, atol=1e-9)

    def test_offsets_apply_in_sample_frame(self):
        """R(delta) = exp([delta]x) R_nom: a pure-z offset rotates about lab z."""
        R_nom = Rotation.random(random_state=4).as_matrix()
        R = offsets_to_matrices(np.array([[0.0, 0.0, 30.0]]), R_nom)[0]
        assert np.allclose(R @ R_nom.T, Rotation.from_euler("z", 30, degrees=True).as_matrix())

    def test_prior_sampling_stays_in_ball(self):
        d = sample_prior_offsets(2000, 2.5, np.random.default_rng(0))
        assert np.linalg.norm(d, axis=1).max() <= 2.5 + 1e-12
        # uniform in a ball: the median radius is R * 0.5^(1/3)
        assert np.median(np.linalg.norm(d, axis=1)) == pytest.approx(2.5 * 0.5 ** (1 / 3), rel=0.05)

    def test_fixed_magnitude_sampling(self):
        d = sample_fixed_magnitude_offsets(50, 0.75, np.random.default_rng(0))
        assert np.allclose(np.linalg.norm(d, axis=1), 0.75)


class TestErrorSummary:
    def test_zero_error(self):
        d = np.random.default_rng(0).normal(size=(10, 3))
        s = error_summary(d, d)
        assert (
            s["rms_z"] == 0
            and s["rms_perp"] == 0
            and s["median_angle"] == pytest.approx(0, abs=1e-9)
        )

    def test_success_rates(self):
        truth = np.zeros((4, 3))
        hat = np.array([[0.05, 0, 0], [0.3, 0, 0], [0.7, 0, 0], [0, 0, 0.02]])
        s = error_summary(hat, truth)
        assert s["success_0p5"] == pytest.approx(0.75)  # 0.05, 0.3, 0.02 pass; 0.7 fails
        assert s["success_0p1"] == pytest.approx(0.5)  # 0.05 and 0.02 pass

    def test_axis_separation(self):
        truth = np.zeros((4, 3))
        z_err = truth + np.array([0.0, 0.0, 0.3])
        s = error_summary(z_err, truth)
        assert s["rms_z"] == pytest.approx(0.3) and s["rms_perp"] == 0
        xy_err = truth + np.array([0.3, 0.3, 0.0])
        s = error_summary(xy_err, truth)
        assert s["rms_z"] == 0 and s["rms_perp"] == pytest.approx(0.3)


class TestPixelGrid:
    def test_spot_overlaps_grid_rule(self):
        assert spot_overlaps_grid([(10.2, 10.7), (11.0, 10.1), (10.5, 11.9)], 2048, 2048)
        assert not spot_overlaps_grid(
            [(-5.2, 10.0), (-3.0, 11.0), (-4.0, 12.0)], 2048, 2048
        )  # left of the grid
        assert not spot_overlaps_grid(
            [(10.0, 2100.0), (11.0, 2101.0), (10.5, 2102.0)], 2048, 2048
        )  # below it
        assert spot_overlaps_grid(
            [(-2.0, 10.0), (3.0, 11.0), (1.0, 12.0)], 2048, 2048
        )  # straddles the edge

    def test_lit_pixel_set_matches_the_rasteriser(self):
        """lit_pixel_set must reproduce ImageData.add_triangle_scanline exactly,
        including spots clipped by, or entirely outside, the grid."""
        rng = np.random.default_rng(0)
        ncols, nrows = 40, 30
        n_partial = n_empty = 0
        for _ in range(400):
            centre = rng.uniform([-6, -6], [ncols + 6, nrows + 6])
            v = centre + rng.uniform(-3.5, 3.5, size=(3, 2))
            img = ImageData(nrows, ncols)
            img.add_triangle_scanline(
                torch.tensor(v[0]), torch.tensor(v[1]), torch.tensor(v[2]), 1.0
            )
            drawn = {(int(c), int(r)) for r, c in zip(*np.nonzero(img.to_numpy()))}
            key = tuple(int(-1 if x < 0 else int(x)) for xy in v for x in xy)
            assert lit_pixel_set(key, ncols, nrows) == drawn
            n_empty += not drawn
            n_partial += bool(drawn) and any(x < 0 or x > ncols - 1 for x in v[:, 0])
        assert n_empty > 20 and n_partial > 10  # the test does exercise clipping

    def test_roi_peaks_are_recorded_peaks(self, problem):
        """Every ROI peak must overlap the pixel grid at the nominal orientation."""
        obs = _observer(problem)
        keys = obs.vertex_keys(obs.observe(torch.zeros(1, 3, dtype=torch.float64)))[0]
        cols, rows = keys[:, 0::2], keys[:, 1::2]
        ncols, nrows = obs.d_ncols[:, None], obs.d_nrows[:, None]
        on_grid = (
            (cols.amax(-1) >= 0)
            & (cols.amin(-1) <= ncols[:, 0] - 1)
            & (rows.amax(-1) >= 0)
            & (rows.amin(-1) <= nrows[:, 0] - 1)
        )
        assert bool(on_grid.all())


class TestOffsetHead:
    def test_cholesky_is_lower_triangular_positive_diagonal(self):
        L = cholesky_from_raw(torch.randn(7, 6))
        assert torch.allclose(L, torch.tril(L))
        assert (torch.diagonal(L, dim1=-2, dim2=-1) > 0).all()

    def test_nll_matches_torch_distribution(self):
        torch.manual_seed(0)
        L = cholesky_from_raw(torch.randn(6, 6))
        mean, tgt = torch.randn(6, 3), torch.randn(6, 3)
        ref = -torch.distributions.MultivariateNormal(mean, scale_tril=L).log_prob(
            tgt
        ) - 1.5 * np.log(2 * np.pi)
        assert torch.allclose(gaussian_nll_loss(mean, L, tgt), ref.mean(), atol=1e-5)

    def test_beta_nll_is_finite_and_differs(self):
        torch.manual_seed(0)
        raw = torch.randn(6, 6, requires_grad=True)
        L = cholesky_from_raw(raw)
        mean, tgt = torch.randn(6, 3), torch.randn(6, 3)
        plain = gaussian_nll_loss(mean, L, tgt)
        weighted = gaussian_nll_loss(mean, L, tgt, beta=0.5)
        weighted.backward()
        assert torch.isfinite(weighted) and torch.isfinite(raw.grad).all()
        assert not torch.isclose(plain, weighted)

    def test_decoupled_loss_mean_gradient_is_mse_gradient(self):
        torch.manual_seed(0)
        raw = torch.randn(6, 6)
        mean = torch.randn(6, 3, requires_grad=True)
        tgt = torch.randn(6, 3)
        decoupled_nll_loss(mean, cholesky_from_raw(raw), tgt, 0.3).backward()
        g_dec = mean.grad.clone()
        mean.grad = None
        mse_deg_loss(mean, tgt, 0.3).backward()
        assert torch.allclose(g_dec, mean.grad, atol=1e-6)
        # independent of sigma: a very different covariance gives the same mean gradient
        mean.grad = None
        decoupled_nll_loss(mean, cholesky_from_raw(raw * 3.0 - 2.0), tgt, 0.3).backward()
        assert torch.allclose(g_dec, mean.grad, atol=1e-6)

    def test_decoupled_loss_covariance_gradient_is_nll_gradient(self):
        torch.manual_seed(1)
        raw = torch.randn(5, 6, requires_grad=True)
        mean, tgt = torch.randn(5, 3, requires_grad=True), torch.randn(5, 3)
        decoupled_nll_loss(mean, cholesky_from_raw(raw), tgt).backward()
        g_dec = raw.grad.clone()
        raw.grad = None
        gaussian_nll_loss(mean.detach(), cholesky_from_raw(raw), tgt).backward()
        assert torch.allclose(g_dec, raw.grad, atol=1e-6)

    def test_nll_mean_gradient_scales_with_inverse_variance(self):
        # the pathology under test: the plain NLL gradient on the mean is Sigma^-1 (mean - y)
        L = torch.diag(torch.tensor([0.02, 0.02, 0.5]))[None]
        mean, tgt = torch.zeros(1, 3, requires_grad=True), torch.ones(1, 3)
        gaussian_nll_loss(mean, L, tgt).backward()
        g = mean.grad[0].abs()
        assert torch.isclose(g[0] / g[2], torch.tensor((0.5 / 0.02) ** 2), rtol=1e-4)

    def test_model_shapes_and_gradients(self):
        net = ToyOffsetNet(n_peaks=3, window_size=8, hidden=(16, 8, 8))
        mean, chol = net(torch.rand(2, 3, 8, 8))
        assert mean.shape == (2, 3) and chol.shape == (2, 3, 3)
        gaussian_nll_loss(mean, chol, torch.zeros(2, 3)).backward()
        assert all(torch.isfinite(p.grad).all() for p in net.parameters())


# ============================================================================
# Physics tests (Example2.ThreeVoxels)
# ============================================================================


@pytest.fixture
def project_root():
    return Path(__file__).parent.parent.parent


@pytest.fixture(autouse=True)
def chdir_to_project_root(project_root, monkeypatch):
    monkeypatch.chdir(project_root)


@pytest.fixture
def problem(project_root):
    example_dir = project_root / "Examples" / "Example2.ThreeVoxels"
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    mic_path = example_dir / "SimInput" / "three_voxels.mic"
    if not config_path.exists() or not mic_path.exists():
        pytest.skip("Example2.ThreeVoxels files not found")

    os.chdir(example_dir)
    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = "3Grains.sim"
    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()
    detector_list = exp_setup.get_detector_list()
    range_map = exp_setup.get_range_to_index_map()
    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])
    simulator = Simulation(exp_setup)
    structure_list = sample.get_structure_list()
    voxel = MicFile.read(str(mic_path)).voxels[0]
    vertices = _get_voxel_vertices(voxel)
    R_nom = voxel.orientation.astype(np.float64)
    roi = define_roi_set(
        torch.from_numpy(R_nom).float(),
        vertices,
        sample,
        detector_list,
        range_map,
        exp_setup,
        structure_list,
        simulator,
        phase_index=voxel.phase,
    )
    if not roi:
        pytest.skip("empty ROI set")
    roi_all = define_roi_set(
        torch.from_numpy(R_nom).float(),
        vertices,
        sample,
        detector_list,
        range_map,
        exp_setup,
        structure_list,
        simulator,
        phase_index=voxel.phase,
        detectors="all",
    )
    return dict(
        R_nom=R_nom,
        vertices=vertices,
        roi=roi,
        roi_all=roi_all,
        structure_list=structure_list,
        voxel=voxel,
        config_path=config_path,
        sample=sample,
        detector_list=detector_list,
        range_map=range_map,
        exp_setup=exp_setup,
        simulator=simulator,
    )


def _observer(problem, roi=None):
    return BatchedObserver(
        problem["R_nom"],
        problem["vertices"],
        problem["sample"],
        problem["detector_list"],
        problem["range_map"],
        problem["exp_setup"],
        roi if roi is not None else problem["roi"],
    )


class TestBatchedObserver:
    @pytest.mark.parametrize("mode", ["first", "all"])
    def test_matches_simulator(self, problem, mode):
        """Presence, frame index and spot centroids agree with the serial simulator
        semantics: drop the peak if any vertex misses any detector plane, then
        require the spot to overlap its own detector's pixel grid."""
        from icenine.diffraction_core import get_scattering_omegas_torch

        roi = problem["roi" if mode == "first" else "roi_all"][
            ::12
        ]  # keeps the serial reference fast
        obs = _observer(problem, roi)
        offsets = np.array([[0.0, 0.0, 0.0], [0.3, -0.2, 0.5], [-1.0, 0.7, 1.2], [1.5, 1.5, -1.5]])
        out = obs.observe(torch.from_numpy(offsets))
        mats = offsets_to_matrices(offsets, problem["R_nom"])
        sample, es = problem["sample"], problem["exp_setup"]
        base = sample.sample_to_lab_matrix[:3, :3].clone()

        n_checked = 0
        for b in range(len(offsets)):
            for m, p in enumerate(roi):
                g_s = mats[b] @ p.g_hkl.double().numpy()
                g = torch.from_numpy(g_s)[None]
                res = get_scattering_omegas_torch(
                    g,
                    torch.norm(g, dim=1),
                    float(es.beam_energy),
                    es.get_beam_deflection_chi_laue(),
                )
                spot, frame = None, None
                if bool(res.observable[0]):
                    w = float((res.omega1 if p.omega_branch == 1 else res.omega2)[0])
                    frame = problem["range_map"].angle_to_wedge_index(w)
                    if frame is not None:
                        _restore_and_rotate(sample, base, w)
                        results = _project_all_detectors(
                            problem["simulator"],
                            sample,
                            problem["detector_list"],
                            problem["vertices"],
                            torch.from_numpy(g_s / np.linalg.norm(g_s)).float(),
                            XDMEtaAcceptFn(0.0, es.get_eta_limit(), p.form_intensity, p.sin_2theta),
                        )
                        _restore_and_rotate(sample, base, 0.0)
                        det = problem["detector_list"][p.detector_index]
                        spot = None if results is None else results[p.detector_index]
                        if spot is not None and not spot_overlaps_grid(
                            spot[3], det.num_cols, det.num_rows
                        ):
                            spot = None
                assert bool(out.present[b, m]) == (spot is not None)
                if spot is not None:
                    assert int(out.frame[b, m]) == frame
                    centroid = out.verts[b, m].mean(dim=0).numpy()
                    assert np.abs(centroid - np.array([spot[1], spot[0]])).max() < 1e-3
                    n_checked += 1
        assert n_checked > 100

    def test_all_detectors_mode(self, problem):
        """'all' keeps the 'first' entries on detector 0, adds entries on other
        detectors, and every entry is present at the nominal orientation."""
        key = lambda p: (p.reflection_index, p.omega_branch, p.detector_index)  # noqa: E731
        first = {key(p) for p in problem["roi"]}
        all_ = {key(p) for p in problem["roi_all"]}
        assert {k for k in all_ if k[2] == 0} == first
        assert any(k[2] == 1 for k in all_)
        out = _observer(problem, problem["roi_all"]).observe(torch.zeros(1, 3, dtype=torch.float64))
        assert bool(out.present.all())

    def test_invalid_detectors_mode(self, problem):
        with pytest.raises(ValueError):
            define_roi_set(
                torch.from_numpy(problem["R_nom"]).float(),
                problem["vertices"],
                problem["sample"],
                problem["detector_list"],
                problem["range_map"],
                problem["exp_setup"],
                problem["structure_list"],
                problem["simulator"],
                phase_index=problem["voxel"].phase,
                detectors="second",
            )

    def test_missing_any_detector_plane_drops_peak_everywhere(self, problem):
        """Adding a detector plane behind the sample (every ray has t < 0 for it)
        must remove every peak on every detector, as in _simulate_peaks."""
        obs = _observer(problem, problem["roi_all"][::20])
        assert bool(obs.observe(torch.zeros(1, 3, dtype=torch.float64)).present.all())
        behind = torch.tensor([[1.0, 0.0, 0.0]], dtype=torch.float64)  # plane x = -1 mm
        obs.all_normals = torch.cat([obs.all_normals, behind])
        obs.all_plane_d = torch.cat([obs.all_plane_d, torch.tensor([1.0], dtype=torch.float64)])
        assert not bool(obs.observe(torch.zeros(1, 3, dtype=torch.float64)).present.any())

    def test_max_q_override_limits_reflections(self, problem):
        """Setting config.max_q before initialize_sample limits the reflection list
        (how the Stage 1 scripts apply Q_max = 8)."""
        config = ConfigFile.from_file(str(problem["config_path"]))
        config.out_file_basename = "3Grains.sim"
        config.max_q = 8.0
        exp_setup = XDMExperimentSetup(config)
        exp_setup.initialize_experiment()
        sample = Sample()
        exp_setup.initialize_sample(sample, exp_setup.get_detector_list()[0])
        phase = problem["voxel"].phase
        limited = sample.get_structure_list()[phase].get_reflection_vectors()
        default = problem["structure_list"][phase].get_reflection_vectors()
        assert 0 < len(limited) < len(default)
        assert max(r.q_mag for r in limited) <= 8.0 + 1e-9

    def test_nominal_offset_reproduces_roi_set(self, problem):
        obs = _observer(problem)
        out = obs.observe(torch.zeros(1, 3, dtype=torch.float64))
        assert bool(out.present.all())  # the ROI set is defined as the peaks present at nominal

    def test_sample_state_untouched(self, problem):
        before = problem["sample"].sample_to_lab_matrix.clone()
        _observer(problem).observe(torch.zeros(2, 3, dtype=torch.float64))
        assert torch.equal(problem["sample"].sample_to_lab_matrix, before)

    def test_rotation_about_stage_axis_shifts_every_omega_exactly(self, problem):
        """Docs 3.1.5: rotating the orientation by beta about z shifts every omega* by -beta."""
        obs = _observer(problem, problem["roi"][::10])
        beta = 0.8
        out0 = obs.observe(torch.zeros(1, 3, dtype=torch.float64))
        outz = obs.observe(torch.tensor([[0.0, 0.0, beta]], dtype=torch.float64))
        both = out0.present[0] & outz.present[0]
        d = outz.omega[0, both] - out0.omega[0, both]
        d = (d + np.pi) % (2 * np.pi) - np.pi
        assert torch.allclose(d, torch.full_like(d, -beta * DEG), atol=1e-9)


class TestExactBayes:
    def test_posterior_contains_truth_and_is_tight(self, problem):
        obs = _observer(problem)
        bayes = ExactBayes(obs, prior_radius_deg=2.5, use_pixels=True)
        truth = np.array([0.6, -0.4, 0.9])
        r = bayes.posterior(truth, np.random.default_rng(0), n_per_round=8000)
        assert r["ess"] > 100
        sd = np.sqrt(np.diag(r["cov"]))
        # the noise-free cell is far smaller than the independent-quantisation estimate (~0.01 deg)
        assert (sd < 0.005).all()
        # truth lies inside the consistent set, so the posterior mean is near it
        assert np.linalg.norm(r["mean"] - truth) < 0.01
        # and the posterior beats predicting nominal by a wide margin
        assert np.sqrt(np.trace(r["cov"])) < 0.01 * np.linalg.norm(truth)

    def test_pixels_tighten_the_posterior(self, problem):
        obs = _observer(problem)
        truth = np.array([-0.5, 0.3, 0.7])
        frames = ExactBayes(obs, 2.5, use_pixels=False).posterior(
            truth, np.random.default_rng(1), n_per_round=8000
        )
        pixels = ExactBayes(obs, 2.5, use_pixels=True).posterior(
            truth, np.random.default_rng(1), n_per_round=8000
        )
        assert np.trace(pixels["cov"]) < np.trace(frames["cov"])


# ============================================================================
# Stage 1: frame-coded windows from the observer
# ============================================================================


class TestRenderWindows:
    def test_decode_windows(self):
        from icenine.orientation_eval import decode_windows

        w = torch.tensor([[0, 1, 5, 9]], dtype=torch.uint8)  # K = 4: codes 1..9 -> offsets -4..4
        out = decode_windows(w, 4)
        assert out.shape == (2, 1, 4)
        assert torch.equal(out[0], torch.tensor([[0.0, 1.0, 1.0, 1.0]]))
        assert torch.allclose(out[1], torch.tensor([[0.0, -1.0, 0.0, 1.0]]))

    def test_windows_match_lit_pixel_sets(self, problem):
        """Window content = the rasteriser's lit pixels shifted by the window origin, with
        code 1 + (frame - frame0 + K); status flags absent / out-of-range spots."""
        from icenine.orientation_eval import WindowSpec, render_windows

        obs = _observer(problem, problem["roi_all"][::5])
        spec = WindowSpec.from_nominal(obs, 32, 4)
        deltas = np.array([[0.0, 0.0, 0.0], [0.2, -0.3, 0.6], [-0.5, 0.4, -1.0]])
        windows, status = render_windows(obs, spec, deltas)
        out = obs.observe(torch.from_numpy(deltas))
        keys = obs.vertex_keys(out)
        n_checked = 0
        for n in range(len(deltas)):
            for m in range(obs.M):
                if not bool(out.present[n, m]):
                    assert int(status[n, m]) == 1 and int(windows[n, m].sum()) == 0
                    continue
                offset = int(out.frame[n, m]) - int(spec.frame0[m])
                if abs(offset) > 4:
                    assert int(status[n, m]) == 2
                    continue
                pix = lit_pixel_set(
                    tuple(keys[n, m].tolist()), int(obs.d_ncols[m]), int(obs.d_nrows[m])
                )
                expected = {
                    (c - spec.col0[m], r - spec.row0[m])
                    for c, r in pix
                    if 0 <= c - spec.col0[m] < 32 and 0 <= r - spec.row0[m] < 32
                }
                rows, cols = np.nonzero(windows[n, m].numpy())
                assert set(zip(cols.tolist(), rows.tolist())) == expected
                assert set(windows[n, m][windows[n, m] > 0].tolist()) <= {1 + offset + 4}
                n_checked += 1
        assert n_checked > 50

    def test_sin_eta_near_axis(self, problem):
        obs = _observer(problem, problem["roi_all"])
        se = obs.sin_eta(torch.zeros(1, 3, dtype=torch.float64))[0]
        assert ((se >= 0) & (se <= 1)).all()
        assert (se < 0.3).any() and (se > 0.9).any()  # both near-axis and far-from-axis spots exist


class TestEndToEndVsForwardSimulation:
    def test_windows_reproduce_simulated_images(self, problem, project_root):
        """At the nominal orientation the 'all' ROI set is every recorded peak, so the
        union of the frame-coded windows must equal, pixel for pixel, the thresholded
        images ForwardSimulation._simulate_peaks produces for that one voxel. At a small
        offset, every rendered pixel must be lit in the simulated images."""
        import copy

        from icenine.forward_simulation import ForwardSimulation
        from icenine.image_data import ImageData
        from icenine.orientation_eval import WindowSpec, render_windows

        config = ConfigFile.from_file(str(problem["config_path"]))
        config.out_file_basename = "3Grains.sim"
        fs = ForwardSimulation(config)
        fs.exp_setup.initialize_experiment()
        detector_list = fs.exp_setup.get_detector_list()
        range_map = fs.exp_setup.get_range_to_index_map()
        fs.simulator = Simulation(fs.exp_setup)
        sample = Sample()
        fs.exp_setup.initialize_sample(sample, detector_list[0])

        obs = _observer(problem, problem["roi_all"])
        spec = WindowSpec.from_nominal(obs, 32, 4)
        n_omega = len(fs.exp_setup.get_omega_range_list())

        def simulated_lit(delta_deg):
            voxel = copy.deepcopy(problem["voxel"])
            voxel.orientation = offsets_to_matrices(np.array([delta_deg]), problem["R_nom"])[0]
            sample.get_mic().voxels = [voxel]
            images = [
                [ImageData(d.num_rows, d.num_cols, mode="sparse") for d in detector_list]
                for _ in range(n_omega)
            ]
            fs._simulate_peaks(images, detector_list, sample, range_map)
            lit = set()
            for w in range(n_omega):
                for di, img in enumerate(images[w]):
                    sp = img._pixels_sparse.coalesce()
                    idx, val = sp.indices(), sp.values()
                    for r, c in idx[:, val > 0].T.tolist():
                        lit.add((w, di, c, r))
            return lit

        def rendered_lit(delta_deg):
            windows, status = render_windows(obs, spec, np.array([delta_deg]))
            lit = set()
            for m in range(obs.M):
                w = windows[0, m]
                rows, cols = np.nonzero(w.numpy())
                for r, c in zip(rows.tolist(), cols.tolist()):
                    frame = int(spec.frame0[m]) + int(w[r, c]) - 1 - 4
                    lit.add(
                        (frame, int(obs.det_idx[m]), c + int(spec.col0[m]), r + int(spec.row0[m]))
                    )
            return lit, status

        sim0 = simulated_lit(np.zeros(3))
        ren0, status0 = rendered_lit(np.zeros(3))
        assert bool((status0 == 0).all())
        assert ren0 == sim0

        delta = np.array([0.15, -0.1, 0.3])
        sim1 = simulated_lit(delta)
        ren1, _ = rendered_lit(delta)
        assert ren1 <= sim1  # new peaks may appear at the offset, but nothing rendered is spurious
        assert len(sim1 - ren1) < 0.02 * len(sim1)
