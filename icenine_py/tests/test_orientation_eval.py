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
    _project_peak_on_detector,
    _restore_and_rotate,
    cholesky_from_raw,
    define_roi_set,
    gaussian_nll_loss,
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
        assert s["rms_z"] == 0 and s["rms_perp"] == 0 and s["median_angle"] == pytest.approx(0, abs=1e-9)

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
        assert not spot_overlaps_grid([(-5.2, 10.0), (-3.0, 11.0), (-4.0, 12.0)], 2048, 2048)  # left of the grid
        assert not spot_overlaps_grid([(10.0, 2100.0), (11.0, 2101.0), (10.5, 2102.0)], 2048, 2048)  # below it
        assert spot_overlaps_grid([(-2.0, 10.0), (3.0, 11.0), (1.0, 12.0)], 2048, 2048)  # straddles the edge

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
            img.add_triangle_scanline(torch.tensor(v[0]), torch.tensor(v[1]), torch.tensor(v[2]), 1.0)
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
        on_grid = (cols.amax(-1) >= 0) & (cols.amin(-1) <= ncols[:, 0] - 1) & (rows.amax(-1) >= 0) & (rows.amin(-1) <= nrows[:, 0] - 1)
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
        ref = -torch.distributions.MultivariateNormal(mean, scale_tril=L).log_prob(tgt) - 1.5 * np.log(2 * np.pi)
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
        torch.from_numpy(R_nom).float(), vertices, sample, detector_list, range_map, exp_setup,
        structure_list, simulator, phase_index=voxel.phase,
    )
    if not roi:
        pytest.skip("empty ROI set")
    return dict(
        R_nom=R_nom, vertices=vertices, roi=roi, sample=sample, detector_list=detector_list,
        range_map=range_map, exp_setup=exp_setup, simulator=simulator,
    )


def _observer(problem, roi=None):
    return BatchedObserver(
        problem["R_nom"], problem["vertices"], problem["sample"], problem["detector_list"],
        problem["range_map"], problem["exp_setup"], roi if roi is not None else problem["roi"],
    )


class TestBatchedObserver:
    def test_matches_simulator(self, problem):
        """Presence, frame index and spot centroids agree with the serial simulator path."""
        from icenine.diffraction_core import get_scattering_omegas_torch

        roi = problem["roi"][::12]  # ~65 peaks keeps the serial reference fast
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
                res = get_scattering_omegas_torch(g, torch.norm(g, dim=1), float(es.beam_energy), es.get_beam_deflection_chi_laue())
                spot, frame = None, None
                if bool(res.observable[0]):
                    w = float((res.omega1 if p.omega_branch == 1 else res.omega2)[0])
                    frame = problem["range_map"].angle_to_wedge_index(w)
                    if frame is not None:
                        _restore_and_rotate(sample, base, w)
                        spot = _project_peak_on_detector(
                            problem["simulator"], sample, problem["detector_list"][p.detector_index], problem["vertices"],
                            torch.from_numpy(g_s / np.linalg.norm(g_s)).float(),
                            XDMEtaAcceptFn(0.0, es.get_eta_limit(), p.form_intensity, p.sin_2theta),
                        )
                        _restore_and_rotate(sample, base, 0.0)
                assert bool(out.present[b, m]) == (spot is not None)
                if spot is not None:
                    assert int(out.frame[b, m]) == frame
                    centroid = out.verts[b, m].mean(dim=0).numpy()
                    assert np.abs(centroid - np.array([spot[1], spot[0]])).max() < 5e-3
                    n_checked += 1
        assert n_checked > 100

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
        frames = ExactBayes(obs, 2.5, use_pixels=False).posterior(truth, np.random.default_rng(1), n_per_round=8000)
        pixels = ExactBayes(obs, 2.5, use_pixels=True).posterior(truth, np.random.default_rng(1), n_per_round=8000)
        assert np.trace(pixels["cov"]) < np.trace(frames["cov"])
