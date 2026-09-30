"""
Tests for the Stage 2 Gauss-Newton baseline (orientation_baselines.py), the per-peak
context features, and the Stage 3 shared per-peak network (PeakSetNet).

Physics tests use Example2.ThreeVoxels with Q_max = 8 and are skipped if it is absent.
"""

import os
from pathlib import Path

import numpy as np
import pytest
import torch

from icenine.config_file import ConfigFile
from icenine.experiment_setup import XDMExperimentSetup
from icenine.mic_file import MicFile
from icenine.orientation_baselines import (
    CentroidGaussNewton,
    extract_measurements,
    frame_center_omega,
)
from icenine.orientation_eval import (
    BatchedObserver,
    WindowSpec,
    error_summary,
    render_windows,
)
from icenine.orientation_nn import define_roi_set
from icenine.reconstructor import _get_voxel_vertices
from icenine.sample import Sample
from icenine.simulation import Simulation
from icenine.orientation_eval import decode_windows
from icenine.toy_orientation_model import PeakSetNet, measurement_features

# ============================================================================
# PeakSetNet (no example data needed)
# ============================================================================


class TestPeakSetNet:
    def _batch(self, B=3, M=7, present=5):
        x = torch.zeros(B, M, 2, 32, 32)
        x[:, :present, 0, 10:12, 12:14] = 1.0
        x[:, :present, 1, 10:12, 12:14] = 0.5
        return x, torch.randn(M, 16)

    def test_shapes_and_positive_definite_covariance(self):
        net = PeakSetNet()
        x, ctx = self._batch()
        mean, chol = net(x, ctx)
        assert mean.shape == (3, 3) and chol.shape == (3, 3, 3)
        assert (torch.diagonal(chol, dim1=-2, dim2=-1) > 0).all()

    def test_parameter_count_independent_of_peak_count(self):
        net = PeakSetNet()
        n = sum(p.numel() for p in net.parameters())
        for M in (5, 40):
            mean, _ = net(torch.zeros(2, M, 2, 32, 32), torch.randn(M, 16))
            assert mean.shape == (2, 3)
        assert n == sum(p.numel() for p in net.parameters())

    def test_permutation_invariance(self):
        torch.manual_seed(0)
        net = PeakSetNet()
        x, ctx = self._batch()
        perm = torch.randperm(x.shape[1])
        a, La = net(x, ctx)
        b, Lb = net(x[:, perm], ctx[perm])
        assert torch.allclose(a, b, atol=1e-5) and torch.allclose(La, Lb, atol=1e-5)

    def test_absent_peaks_are_ignored(self):
        """A peak with an all-zero window must not influence the output, whatever its context."""
        torch.manual_seed(0)
        net = PeakSetNet()
        x, ctx = self._batch(present=5)
        a, _ = net(x, ctx)
        ctx2 = ctx.clone()
        ctx2[5:] = torch.randn(2, 16) * 10
        b, _ = net(x, ctx2)
        assert torch.allclose(a, b, atol=1e-5)

    def test_per_sample_context_gather_matches_per_voxel_calls(self):
        """Multi-voxel batches: context table[voxel_id] (B, M, D) == one call per voxel."""
        torch.manual_seed(0)
        net = PeakSetNet(pool="all")
        x, _ = self._batch(B=4, M=6, present=4)
        table = torch.randn(2, 6, 16)
        vid = torch.tensor([0, 1, 1, 0])
        batched, _ = net(x, table[vid])
        for i in range(4):
            single, _ = net(x[i : i + 1], table[vid[i]])
            assert torch.allclose(batched[i], single[0], atol=1e-5)

    def test_no_peaks_present_is_finite(self):
        net = PeakSetNet()
        mean, chol = net(torch.zeros(2, 4, 2, 32, 32), torch.randn(4, 16))
        assert torch.isfinite(mean).all() and torch.isfinite(chol).all()

    def test_gradients_flow(self):
        from icenine.orientation_nn import gaussian_nll_loss

        net = PeakSetNet()
        x, ctx = self._batch()
        mean, chol = net(x, ctx)
        gaussian_nll_loss(mean, chol, torch.zeros(3, 3)).backward()
        assert all(torch.isfinite(p.grad).all() for p in net.parameters())


class TestMeasurementFeatures:
    def test_synthetic_window(self):
        x = torch.zeros(1, 2, 2, 32, 32)
        x[0, 0, 0, 10:12, 12:15] = 1.0  # rows 10-11, cols 12-14
        x[0, 0, 1, 10:12, 12:15] = 0.5  # frame offset 0.5 * K
        f = measurement_features(x, 4)[0]
        assert torch.allclose(f[0], torch.tensor([1.0, 6 / 20.0, 2.0, 13.5 - 16, 11.0 - 16]))
        assert torch.equal(f[1], torch.zeros(5))  # absent peak

    def test_no_frame_channel_gives_zero_frame_feature(self):
        x = torch.zeros(1, 1, 1, 32, 32)
        x[0, 0, 0, 4:6, 4:6] = 1.0
        assert measurement_features(x, 4)[0, 0, 2] == 0

    @pytest.mark.parametrize("pool", ["meanmax", "meansum", "all"])
    def test_padding_peaks_do_not_change_output(self, pool):
        """Extra all-zero peaks (any context) must not change the prediction."""
        torch.manual_seed(0)
        net = PeakSetNet(pool=pool)
        x = torch.zeros(2, 6, 2, 32, 32)
        x[:, :4, 0, 8:10, 8:10] = 1.0
        x[:, :4, 1, 8:10, 8:10] = -0.25
        ctx = torch.randn(6, 16)
        a, _ = net(x[:, :4], ctx[:4])
        b, _ = net(x, ctx)
        assert torch.allclose(a, b, atol=1e-5)


# ============================================================================
# Physics tests
# ============================================================================


@pytest.fixture
def project_root():
    return Path(__file__).parent.parent.parent


@pytest.fixture(autouse=True)
def chdir_to_project_root(project_root, monkeypatch):
    monkeypatch.chdir(project_root)


@pytest.fixture
def stage1(project_root):
    """Voxel 0, Q_max = 8, both detectors, near-axis spots dropped (the Stage 1 setup)."""
    example_dir = project_root / "Examples" / "Example2.ThreeVoxels"
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    mic_path = example_dir / "SimInput" / "three_voxels.mic"
    if not config_path.exists() or not mic_path.exists():
        pytest.skip("Example2.ThreeVoxels files not found")
    os.chdir(example_dir)
    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = "3Grains.sim"
    config.max_q = 8.0
    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()
    detector_list = exp_setup.get_detector_list()
    range_map = exp_setup.get_range_to_index_map()
    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])
    simulator = Simulation(exp_setup)
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
        sample.get_structure_list(),
        simulator,
        phase_index=voxel.phase,
        detectors="all",
    )
    obs = BatchedObserver(R_nom, vertices, sample, detector_list, range_map, exp_setup, roi)
    keep = obs.sin_eta(torch.zeros(1, 3, dtype=torch.float64))[0].numpy() >= 0.3
    roi = [p for p, k in zip(roi, keep) if k]
    obs = BatchedObserver(R_nom, vertices, sample, detector_list, range_map, exp_setup, roi)
    return dict(obs=obs, spec=WindowSpec.from_nominal(obs, 32, 4))


class TestMeasurements:
    def test_frame_center_omega(self, stage1):
        obs = stage1["obs"]
        w = frame_center_omega(obs, np.array([0, 90, 179]))
        assert np.allclose(np.degrees(w), [-89.5, 0.5, 89.5])

    def test_extract_matches_observer(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        delta = np.array([[0.2, -0.3, 0.5]])
        windows, status = render_windows(obs, spec, delta)
        meas = extract_measurements(windows[0], spec, obs)
        assert np.array_equal(meas.used, (status[0] == 0).numpy())
        out = obs.observe(torch.from_numpy(delta))
        u = meas.used
        # the measured frame centre is within half a frame of the true crossing
        frame_width = abs(obs.range_width)
        d = np.abs((out.omega[0].numpy()[u] - meas.omega[u] + np.pi) % (2 * np.pi) - np.pi)
        assert (d <= 0.5 * frame_width + 1e-9).all()
        # the lit-pixel centroid is within about a pixel of the exact spot centroid
        truth = out.verts[0].mean(dim=1).numpy()[u]
        assert np.abs(meas.centroid[u] - truth).max() < 2.0


class TestMeasurementFeaturesMatchExtraction:
    def test_features_equal_extract_measurements(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        delta = np.array([[0.2, -0.3, 0.5]])
        windows, _ = render_windows(obs, spec, delta)
        meas = extract_measurements(windows[0], spec, obs)
        x = decode_windows(windows, 4)  # (1, M, 2, W, W)
        f = measurement_features(x.double(), 4)[0].numpy()
        u = meas.used
        assert np.array_equal(f[:, 0] > 0, u)
        frame = np.array([spec.frame0[m] for m in range(obs.M)]) + f[:, 2]
        expect = np.full(obs.M, np.nan)
        bin_of = {int(w): b for b, w in enumerate(obs.range_index.tolist()) if w >= 0}
        for m in np.nonzero(u)[0]:
            b = bin_of[int(round(frame[m]))]
            expect[m] = obs.range_low + (b + 0.5) * obs.range_width
        assert np.allclose(expect[u], meas.omega[u])
        cx = spec.col0 + 16 + f[:, 3]
        cy = spec.row0 + 16 + f[:, 4]
        assert np.allclose(cx[u], meas.centroid[u, 0]) and np.allclose(cy[u], meas.centroid[u, 1])


class TestGaussNewton:
    def test_recovers_offsets_far_better_than_nominal(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        rng = np.random.default_rng(0)
        d = rng.normal(size=(12, 3))
        d = d / np.linalg.norm(d, axis=1, keepdims=True) * 0.6
        windows, _ = render_windows(obs, spec, d)
        gn = CentroidGaussNewton(obs)
        est = np.array(
            [gn.solve(extract_measurements(windows[i], spec, obs))["delta"] for i in range(len(d))]
        )
        s = error_summary(est, d)
        assert s["rms_perp"] < 0.02  # pixels fix the perpendicular components
        assert s["rms_z"] < 0.12  # frames fix the stage-axis component to ~a tenth of a degree
        assert s["median_angle"] < 0.1 < 0.6  # vs 0.6 for predicting nominal

    def test_nominal_data_gives_small_offset(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        windows, _ = render_windows(obs, spec, np.zeros((1, 3)))
        r = CentroidGaussNewton(obs).solve(extract_measurements(windows[0], spec, obs))
        assert np.linalg.norm(r["delta"]) < 0.15
        assert r["n_used"] > 50 and np.all(np.isfinite(r["cov"]))

    def test_covariance_reflects_the_axis_anisotropy(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        windows, _ = render_windows(obs, spec, np.zeros((1, 3)))
        r = CentroidGaussNewton(obs).solve(extract_measurements(windows[0], spec, obs))
        sd = np.sqrt(np.diag(r["cov"]))
        assert sd[2] > 3 * max(
            sd[0], sd[1]
        )  # the stage-axis component is the poorly determined one


class TestPeakContext:
    def test_context_features(self, stage1):
        obs = stage1["obs"]
        ctx = obs.peak_context()
        assert ctx.shape == (obs.M, 16) and torch.isfinite(ctx).all()
        # d omega*/d delta_z is exactly -1 for every peak (docs 3.1.5); the z column of the
        # spot Jacobian is ~0 for a voxel this close to the rotation axis.
        assert torch.allclose(ctx[:, 8], torch.full((obs.M,), -1.0), atol=1e-4)
        assert ctx[:, 2].abs().max() < 0.05 and ctx[:, 5].abs().max() < 0.05
        assert ((ctx[:, 9] >= 0.3) & (ctx[:, 9] <= 1.0)).all()  # |sin eta| after the D3 cut
        assert torch.allclose(ctx[:, 11] + ctx[:, 12], torch.ones(obs.M))  # detector one-hot
