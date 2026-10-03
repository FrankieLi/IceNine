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
from icenine.toy_orientation_model import FrameProbeNet, PeakSetNet, measurement_features

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
        assert s["median_angle"] < 0.1  # vs 0.6 for predicting nominal

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


class TestFrameProbeNet:
    def test_uses_only_frame_and_dOmega_and_ignores_pixels(self):
        net = FrameProbeNet()
        x = torch.zeros(2, 6, 2, 32, 32)
        x[:, :4, 0, 10:12, 12:14] = 1.0
        x[:, :4, 1, 10:12, 12:14] = 0.5
        ctx = torch.randn(6, 16)
        mean, chol = net(x, ctx)
        assert mean.shape == (2, 3) and torch.allclose(chol[0], torch.eye(3))
        # moving the lit pixels (same frame) does not change the output
        x2 = torch.roll(x, shifts=(5, 3), dims=(3, 4))
        assert torch.allclose(net(x2, ctx)[0], mean, atol=1e-6)
        # changing the frame offset does
        x3 = x.clone()
        x3[:, :4, 1] *= -1
        assert not torch.allclose(net(x3, ctx)[0], mean)


class TestGaussNewtonStatus:
    def test_nominal_solve_reports_valid_status(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        windows, _ = render_windows(obs, spec, np.zeros((1, 3)))
        r = CentroidGaussNewton(obs).solve(extract_measurements(windows[0], spec, obs))
        assert r["status"] in ("step_below_tol", "no_descent", "too_few_spots", "max_iter")
        assert r["converged"] == (r["status"] in ("step_below_tol", "no_descent"))
        assert r["converged"]


class TestVoxelSelection:
    @pytest.fixture(scope="class")
    def gen(self):
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
        import generate_toy_orientation_dataset as g

        return g

    @staticmethod
    def _fake_mic(n=40):
        from types import SimpleNamespace

        rng = np.random.default_rng(0)
        vox = [
            SimpleNamespace(position=(x * 1e-3, 0.0, 0.0), orientation=rng.normal(size=(3, 3)))
            for x in np.linspace(0, 0.5, n)
        ]
        return SimpleNamespace(voxels=vox)

    def test_select_voxels_single_radius(self, gen):
        targets, cands = gen.select_voxels(self._fake_mic(), 1, 500.0, seed=0)
        assert len(targets) == 1 and len(cands) == 1 and len(cands[0]) > 0

    def test_rejected_then_reused_grain_never_duplicated(self, gen):
        # voxels 1 and 3 share grain B; voxel 0 is unusable, so the first list accepts voxel 1,
        # and the second list (candidates 3 -> grain B again, 2) must skip 3 and take 2.
        grains = {0: "A", 1: "B", 2: "C", 3: "B", 4: "D"}
        orient = {k: np.full((3, 3), ord(g)) for k, g in grains.items()}
        unusable = {0}
        out = gen.accept_voxels(
            [[0, 1], [3, 2], [4]],
            lambda i: orient[i],
            lambda i: None if i in unusable else i,
        )
        assert out == [1, 2, 4]
        assert len({grains[i] for i in out}) == len(out)

    def test_rejected_voxel_grain_remains_available(self, gen):
        # voxel 0 (grain A) is rejected as unusable; a later voxel of grain A is still allowed.
        orient = {0: np.zeros((3, 3)), 1: np.ones((3, 3)), 2: np.zeros((3, 3))}
        out = gen.accept_voxels([[0, 1], [2]], lambda i: orient[i], lambda i: None if i == 0 else i)
        assert out == [1, 2]


# ============================================================================
# Learned Gauss-Newton layer and detector pairing (parallax architecture, Steps 2-3)
# ============================================================================


class TestGNLayer:
    def test_normal_equations_match_weighted_lstsq(self):
        from icenine.toy_orientation_model import gn_normal_equations

        rng = np.random.default_rng(0)
        J = rng.normal(size=(2, 9, 3, 3))
        d_true = rng.normal(size=(2, 3))
        y = np.einsum("bmrc,bc->bmr", J, d_true) + 0.01 * rng.normal(size=(2, 9, 3))
        w = rng.uniform(0.2, 2.0, size=(2, 9, 3))
        w[:, 5] = 0.0  # zero weight removes a peak
        d, Ainv = gn_normal_equations(*(torch.from_numpy(a) for a in (J, y, w)))
        for b in range(2):
            sw = np.sqrt(w[b]).reshape(-1)
            Jb = J[b].reshape(-1, 3)
            ref = np.linalg.lstsq(Jb * sw[:, None], y[b].reshape(-1) * sw, rcond=None)[0]
            assert np.allclose(d[b].numpy(), ref, atol=1e-10)
            assert np.allclose(
                Ainv[b].numpy(), np.linalg.inv((Jb * sw[:, None]).T @ (Jb * sw[:, None]))
            )

    def test_unit_weights_reproduce_linear_gauss_newton_step(self, stage1):
        """W = 1, dy = 0 (the initial state): the layer is CentroidGaussNewton's undamped step."""
        from icenine.orientation_eval import nominal_offsets
        from icenine.toy_orientation_model import GNLayerNet

        obs, spec = stage1["obs"], stage1["spec"]
        delta = np.array([[0.2, -0.3, 0.5], [-0.4, 0.1, -0.6]])
        windows, _ = render_windows(obs, spec, delta)
        net = GNLayerNet(frame_width_rad=obs.frame_width_rad, ridge=0.0, delta_scale=0.05)
        x = decode_windows(windows, 4).float()
        aux = dict(nom_off=torch.from_numpy(nominal_offsets(obs)).float())
        mean, chol = net(x, obs.peak_context(), aux)
        gn = CentroidGaussNewton(obs)
        for i in range(len(delta)):
            ref = gn.solve_linear(extract_measurements(windows[i], spec, obs))
            assert np.allclose(mean[i].detach().numpy(), ref, atol=2e-3), (mean[i], ref)
            # the covariance is (J^T W J)^-1 of the same model
            cov = (chol[i] @ chol[i].T).detach().numpy()
            meas = extract_measurements(windows[i], spec, obs)
            ref_cov = np.linalg.inv(gn.information(meas))
            assert np.allclose(np.diag(cov), np.diag(ref_cov), rtol=0.03)

    @pytest.mark.parametrize("pairing", [False, True])
    def test_permutation_and_padding_invariance(self, pairing):
        from icenine.toy_orientation_model import GNLayerNet

        torch.manual_seed(0)
        net = GNLayerNet(pairing=pairing, n_iter=2)
        for p in net.parameters():  # leave the zero-initialised heads
            if p.abs().sum() == 0:
                torch.nn.init.normal_(p, std=0.1)
        M = 8
        x = torch.zeros(2, M, 2, 32, 32)
        x[:, :6, 0, 10:12, 12:14] = 1.0
        x[:, :6, 1, 10:12, 12:14] = torch.rand(2, 6, 1, 1) - 0.5
        x[1, 2, 0, 12:14, 10:12] = 1.0
        ctx = torch.randn(M, 16)
        nom = torch.rand(M, 3)
        pidx = torch.tensor([3, 2, 1, 0, 5, 4, -1, -1])
        a, La = net(x, ctx, dict(nom_off=nom, pair_index=pidx))
        perm = torch.randperm(M)
        inv = torch.argsort(perm)
        pidx_p = torch.where(pidx[perm] >= 0, inv[pidx[perm].clamp(min=0)], pidx[perm])
        b, Lb = net(x[:, perm], ctx[perm], dict(nom_off=nom[perm], pair_index=pidx_p))
        assert torch.allclose(a, b, atol=1e-5) and torch.allclose(La, Lb, atol=1e-5)
        # padding: extra absent peaks with arbitrary context/offsets do not change the output
        xp = torch.cat([x, torch.zeros(2, 3, 2, 32, 32)], dim=1)
        ctxp = torch.cat([ctx, torch.randn(3, 16) * 10])
        nomp = torch.cat([nom, torch.rand(3, 3)])
        pp = torch.cat([pidx, torch.tensor([-1, -1, -1])])
        c, Lc = net(xp, ctxp, dict(nom_off=nomp, pair_index=pp))
        assert torch.allclose(a, c, atol=1e-5) and torch.allclose(La, Lc, atol=1e-5)

    @staticmethod
    def _gn_batch(M=8):
        torch.manual_seed(1)
        x = torch.rand(3, M, 2, 32, 32) * (torch.rand(3, M, 2, 32, 32) > 0.9).float()
        x[:, :6, 0, 10:12, 12:14] = 1.0
        x[:, :6, 1, 10:12, 12:14] = torch.rand(3, 6, 1, 1) - 0.5
        x[:, 6:] = 0.0  # two absent entries
        ctx = torch.randn(M, 16)
        aux = dict(
            nom_off=torch.rand(M, 3) - 0.5,
            pair_index=torch.tensor([1, 0, 3, 2, 5, 4, -1, -1]),
        )
        return x, ctx, aux

    @staticmethod
    def _backward(net, x, ctx, aux):
        from icenine.orientation_nn import decoupled_nll_loss

        net.zero_grad()
        mean, chol = net(x, ctx, aux) if aux is not None else net(x, ctx)
        decoupled_nll_loss(mean, chol, torch.zeros(3, 3)).backward()

    @pytest.mark.parametrize("kind", ["gn", "gn_paired", "set"])
    def test_gradients_flow_and_are_finite(self, kind):
        """After a few tiny steps off the zero-initialised layers (gradient reaches pair1 only
        after head2/cov2 and then pair2 have left zero), every parameter gets a non-zero, finite
        gradient (regression: the pair MLP once never trained)."""
        from icenine.toy_orientation_model import GNLayerNet

        x, ctx, aux = self._gn_batch()
        if kind == "set":
            net, aux = PeakSetNet(), None
        else:
            net = GNLayerNet(pairing=(kind == "gn_paired"), n_iter=3)
        self._backward(net, x, ctx, aux)
        assert all(torch.isfinite(p.grad).all() for p in net.parameters())
        for _ in range(3):  # zero-initialised layers pass no gradient to their inputs yet
            with torch.no_grad():
                for p in net.parameters():
                    p -= 1e-2 * torch.sign(p.grad)
            self._backward(net, x, ctx, aux)
        dead = [n for n, p in net.named_parameters() if not (p.grad.abs().sum() > 0)]
        assert not dead, dead
        assert all(torch.isfinite(p.grad).all() for p in net.parameters())

    def test_pair_mlp_trains_from_the_zero_init(self):
        """pair2 is zero-initialised and linear: it leaves zero under SGD, and then pair1 trains."""
        from icenine.toy_orientation_model import GNLayerNet

        x, ctx, aux = self._gn_batch()
        net = GNLayerNet(pairing=True, n_iter=3)
        assert net.pair2.weight.abs().sum() == 0
        opt = torch.optim.SGD(net.parameters(), lr=1e-2)
        for _ in range(60):
            self._backward(net, x, ctx, aux)
            opt.step()
        assert net.pair2.weight.abs().sum() > 0
        self._backward(net, x, ctx, aux)
        assert net.pair2.weight.grad.abs().sum() > 0 and net.pair1.weight.grad.abs().sum() > 0

    def test_paired_net_equals_unpaired_net_at_initialisation(self):
        """pair2 = 0: the mixing is off. With the head's extra partner-residual inputs zeroed
        and the shared weights copied, the paired net is the unpaired one; turning pair2 on
        changes the output of paired entries only."""
        import copy

        from icenine.toy_orientation_model import GNLayerNet

        x, ctx, aux = self._gn_batch()
        torch.manual_seed(0)
        plain, paired = GNLayerNet(n_iter=3), GNLayerNet(pairing=True, n_iter=3)
        for q in (plain, paired):  # leave the zero-initialised heads
            for p in (q.head2.weight, q.cov2.weight, q.cov2.bias):
                torch.nn.init.normal_(p, std=0.1)
        sd = plain.state_dict()
        paired.load_state_dict(
            {k: v for k, v in sd.items() if "head1.weight" not in k}, strict=False
        )
        with torch.no_grad():
            paired.head1.weight.zero_()
            paired.head1.weight[:, : sd["head1.weight"].shape[1]] = sd["head1.weight"]
        a, La = plain(x, ctx, aux)
        b, Lb = paired(x, ctx, aux)
        assert torch.allclose(a, b, atol=1e-6) and torch.allclose(La, Lb, atol=1e-6)
        # pair1 is irrelevant while pair2 = 0
        other = copy.deepcopy(paired)
        torch.nn.init.normal_(other.pair1.weight, std=1.0)
        c, _ = other(x, ctx, aux)
        assert torch.allclose(b, c, atol=1e-6)
        # the path is live once pair2 is non-zero
        torch.nn.init.normal_(other.pair2.weight, std=0.3)
        d, _ = other(x, ctx, aux)
        assert not torch.allclose(b, d, atol=1e-4)


class TestPairingAndNominalOffsets:
    def test_pair_index_links_the_two_detector_entries_of_a_ray(self, stage1):
        from icenine.orientation_eval import pair_index

        obs = stage1["obs"]
        roi = obs.roi_list
        pidx = pair_index(roi)
        assert (pidx >= 0).sum() > 40  # most entries are paired
        for i, j in enumerate(pidx):
            if j < 0:
                continue
            assert pidx[j] == i  # symmetric
            assert roi[i].reflection_index == roi[j].reflection_index
            assert roi[i].omega_branch == roi[j].omega_branch
            assert roi[i].detector_index != roi[j].detector_index
            assert abs(roi[i].nominal_omega - roi[j].nominal_omega) < 1e-9  # same ray, same omega

    def test_pair_index_unpaired_and_three_way_cases(self):
        from types import SimpleNamespace as NS

        from icenine.orientation_eval import pair_index

        def p(r, b, d):
            return NS(reflection_index=r, omega_branch=b, detector_index=d)

        roi = [p(0, 1, 0), p(0, 1, 1), p(1, 1, 0), p(0, 2, 1), p(2, 1, 1), p(2, 1, 0)]
        assert pair_index(roi).tolist() == [1, 0, -1, -1, 5, 4]

    def test_features_minus_nominal_offsets_are_measurement_minus_exact_nominal(self, stage1):
        """At the nominal orientation the corrected features are the pure quantisation error."""
        from icenine.orientation_eval import nominal_offsets

        obs, spec = stage1["obs"], stage1["spec"]
        windows, _ = render_windows(obs, spec, np.zeros((1, 3)))
        off = torch.from_numpy(nominal_offsets(obs))
        x = decode_windows(windows, 4).double()
        f = measurement_features(x, 4, off)[0].numpy()
        meas = extract_measurements(windows[0], spec, obs)
        nom = obs.observe(torch.zeros(1, 3, dtype=torch.float64))
        exact = nom.verts[0].mean(dim=1).numpy()
        u = meas.used
        cx = meas.centroid[u, 0] - exact[u, 0]
        cy = meas.centroid[u, 1] - exact[u, 1]
        assert np.allclose(f[u, 3], cx) and np.allclose(f[u, 4], cy)
        assert np.abs(cx).max() < 1.0 and np.abs(cy).max() < 1.0
        # the uncorrected features carry the 0-1 px sub-pixel offset (a large, fixed bias)
        g = measurement_features(x, 4)[0].numpy()
        assert np.abs(g[u, 3]).mean() > np.abs(f[u, 3]).mean()
        # frame: (omega_centre - omega_nom)/width == -(f)... the corrected feature equals that
        wid = obs.range_width
        expect = (meas.omega[u] - nom.omega[0].numpy()[u]) / wid
        assert np.allclose(f[u, 2], expect, atol=1e-9)

    def test_gn_detector_mask_drops_the_other_detector(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        windows, _ = render_windows(obs, spec, np.zeros((1, 3)))
        all_ = extract_measurements(windows[0], spec, obs)
        d0 = extract_measurements(windows[0], spec, obs, detectors=[0])
        d1 = extract_measurements(windows[0], spec, obs, detectors=[1])
        assert d0.used.sum() + d1.used.sum() == all_.used.sum()
        assert not (d0.used & d1.used).any()


class TestDistractors:
    def _shifted_source(self, stage1, shift_mm):
        obs = stage1["obs"]
        verts = obs.vertices.numpy() + np.array([shift_mm, 0.0, 0.0])
        # same orientation/ROI identity, displaced voxel: a same-grain neighbour
        src = BatchedObserver.__new__(BatchedObserver)
        src.__dict__.update(obs.__dict__)
        src.vertices = torch.as_tensor(verts, dtype=obs.dtype)
        return src

    def test_identical_source_reproduces_the_target_windows(self, stage1):
        from icenine.orientation_eval import render_distractor_windows

        obs, spec = stage1["obs"], stage1["spec"]
        delta = np.array([[0.2, -0.3, 0.5], [0.0, 0.1, -0.2]])
        windows, _ = render_windows(obs, spec, delta)
        layer = render_distractor_windows(obs, spec, [obs], [delta])
        assert torch.equal(layer, windows)

    def test_displaced_neighbour_adds_pixels_and_never_alters_the_target(self, stage1):
        from icenine.orientation_eval import combine_windows, render_distractor_windows

        obs, spec = stage1["obs"], stage1["spec"]
        delta = np.array([[0.2, -0.3, 0.5], [0.0, 0.1, -0.2]])
        windows, _ = render_windows(obs, spec, delta)
        src = self._shifted_source(stage1, 0.012)
        layer = render_distractor_windows(obs, spec, [src], [delta + 0.05])
        assert (layer > 0).sum() > 100  # the neighbour lands in the target's windows
        comb = combine_windows(windows, layer)
        lit = windows > 0
        assert torch.equal(comb[lit], windows[lit])  # target pixels (and frame codes) unchanged
        assert ((comb > 0) & ~lit).sum() > 0  # the distractor adds lit pixels elsewhere
        assert torch.equal(comb[~lit & (layer == 0)], windows[~lit & (layer == 0)])
        # the target's own window rendering is independent of the distractor layer
        windows2, _ = render_windows(obs, spec, delta)
        assert torch.equal(windows, windows2)

    def test_realism_variants_are_deterministic_and_bounded(self):
        from icenine.orientation_eval import (
            RealismConfig,
            make_realistic_dataset,
            make_realistic_windows,
        )

        rng = torch.Generator().manual_seed(0)
        w = torch.zeros(6, 10, 32, 32, dtype=torch.uint8)
        w[:, :8, 10:13, 12:15] = 5
        dis = torch.zeros_like(w)
        dis[:, :, 20:22, 20:22] = 3
        assert torch.equal(make_realistic_windows(w, dis, None, 4), w)
        only_nb = make_realistic_windows(w, dis, RealismConfig.named("neighbours"), 4, rng)
        assert torch.equal(only_nb[w > 0], w[w > 0]) and (only_nb[:, :, 20:22, 20:22] == 3).all()
        a = make_realistic_dataset(w, dis, "all", 4, seed=3)
        b = make_realistic_dataset(w, dis, "all", 4, seed=3)
        assert torch.equal(a, b) and int(a.max()) <= 9 and a.dtype == torch.uint8
        dropped = make_realistic_windows(w, None, RealismConfig(False, 1.0, 0.0, 0.0, 0.0), 4)
        assert int(dropped.sum()) == 0
        noisy = make_realistic_dataset(w, None, "noise", 4, seed=1)
        assert (noisy != w).any()


class TestRobustGaussNewton:
    def test_huge_threshold_is_plain_gauss_newton_and_huber_resists_gross_outliers(self, stage1):
        obs, spec = stage1["obs"], stage1["spec"]
        delta = np.array([[0.2, -0.3, 0.5]])
        windows, _ = render_windows(obs, spec, delta)
        meas = extract_measurements(windows[0], spec, obs)
        plain = CentroidGaussNewton(obs).solve(meas)["delta"]
        same = CentroidGaussNewton(obs, huber_c=1e6).solve(meas)["delta"]
        assert np.allclose(plain, same, atol=1e-5)
        # displace the lit spot of 20 entries by 6 px (gross outliers): the robust fit moves
        # less from its own clean-data fit than the plain fit does
        bad = windows.clone()
        present = torch.nonzero((bad[0] > 0).flatten(1).any(1)).flatten()[:20]
        bad[0, present] = torch.roll(bad[0, present], shifts=6, dims=-1)
        mb = extract_measurements(bad[0], spec, obs)
        shifts = []
        for c, clean_fit in (
            (None, plain),
            (1.0, CentroidGaussNewton(obs, huber_c=1.0).solve(meas)["delta"]),
        ):
            shifts.append(
                np.linalg.norm(CentroidGaussNewton(obs, huber_c=c).solve(mb)["delta"] - clean_fit)
            )
        assert shifts[1] < shifts[0]


class TestReviewFixes:
    """Guards and helpers added after the branch review."""

    def test_ema_update_and_checkpoint_selection(self):
        from icenine.orientation_nn import ema_update, select_checkpoint_state

        torch.manual_seed(0)
        model = torch.nn.Sequential(torch.nn.Linear(2, 2), torch.nn.BatchNorm1d(2))
        ema = torch.nn.Sequential(torch.nn.Linear(2, 2), torch.nn.BatchNorm1d(2))
        ema.load_state_dict(model.state_dict())
        w0 = model[0].weight.detach().clone()
        with torch.no_grad():
            model[0].weight.add_(1.0)
            model[1].running_mean.fill_(3.0)
        ema_update(ema, model, 0.9)
        assert torch.allclose(ema[0].weight, 0.9 * w0 + 0.1 * (w0 + 1.0), atol=1e-6)
        assert torch.equal(ema[1].running_mean, model[1].running_mean)  # buffers are copied
        assert not any(p.requires_grad is False for p in model.parameters())
        # decay 1 freezes the EMA
        frozen = ema[0].weight.detach().clone()
        ema_update(ema, model, 1.0)
        assert torch.equal(ema[0].weight, frozen)

        best = {k: v + 5 for k, v in model.state_dict().items() if v.is_floating_point()}
        assert select_checkpoint_state("best", model, best, ema) is best
        s = select_checkpoint_state("ema", model, best, ema)
        assert torch.equal(s["0.weight"], ema[0].weight)
        s = select_checkpoint_state("last", model, best, ema)
        assert torch.equal(s["0.weight"], model[0].weight)
        with pytest.raises(ValueError):
            select_checkpoint_state("ema", model, best, None)
        with pytest.raises(ValueError):
            select_checkpoint_state("best", model, None, ema)
        with pytest.raises(ValueError):
            select_checkpoint_state("other", model, best, ema)

    def test_distractor_frame_filter_and_coding(self, stage1):
        """A source spot on frame f lands in a window iff |f - frame0| <= K, with code
        1 + (f - frame0) + K."""
        import dataclasses

        from icenine.orientation_eval import render_distractor_windows

        obs, spec = stage1["obs"], stage1["spec"]
        K = spec.frame_half_width
        delta = np.zeros((1, 3))  # source spots sit exactly on their nominal frame frame0
        base = render_distractor_windows(obs, spec, [obs], [delta])[0]
        lit = base > 0
        assert lit.any() and (base[lit] == K + 1).all()  # same frame: the centre code
        for s in (1, K):  # target window centred s frames away from the source spot: inside
            sp = dataclasses.replace(spec, frame0=spec.frame0 + s)
            out = render_distractor_windows(obs, sp, [obs], [delta])[0]
            assert (out > 0).sum() == lit.sum()
            assert (out[out > 0] == K + 1 - s).all()
        sp = dataclasses.replace(spec, frame0=spec.frame0 + K + 1)  # outside +-K: filtered out
        assert int(render_distractor_windows(obs, sp, [obs], [delta]).sum()) == 0

    def test_default_ridge_is_close_to_the_undamped_gn_step(self, stage1):
        from icenine.orientation_eval import nominal_offsets
        from icenine.toy_orientation_model import GNLayerNet

        obs, spec = stage1["obs"], stage1["spec"]
        windows, _ = render_windows(obs, spec, np.array([[0.2, -0.3, 0.5], [-0.4, 0.1, -0.6]]))
        x = decode_windows(windows, 4).float()
        aux = dict(nom_off=torch.from_numpy(nominal_offsets(obs)).float())
        out = []
        for ridge in (0.0, 1e-3):
            torch.manual_seed(0)
            net = GNLayerNet(frame_width_rad=obs.frame_width_rad, ridge=ridge, delta_scale=0.05)
            out.append(net(x, obs.peak_context(), aux)[0].detach())
        assert torch.allclose(out[0], out[1], atol=5e-3), (out[0], out[1])
        assert GNLayerNet(frame_width_rad=obs.frame_width_rad).ridge == 1e-3

    def test_realism_of_padded_entries_and_the_valid_mask(self):
        from icenine.orientation_eval import (
            make_realistic_dataset,
            make_realistic_windows,
            RealismConfig,
        )

        w = torch.zeros(8, 6, 32, 32, dtype=torch.uint8)
        w[:, :4, 10:13, 12:15] = 5  # entries 4, 5 are zero padding
        valid = torch.zeros(8, 6, dtype=torch.bool)
        valid[:, :4] = True
        cfg = RealismConfig(False, 0.0, 1.0, 1.0, 0.0)  # a hot pixel and a blob everywhere
        default = make_realistic_windows(w, None, cfg, 4, torch.Generator().manual_seed(0))
        assert (default[:, 4:] > 0).any()  # default: padded entries can be made realistic
        masked = make_realistic_windows(
            w, None, cfg, 4, torch.Generator().manual_seed(0), valid=valid
        )
        assert int(masked[:, 4:].sum()) == 0
        # the draws are unchanged: valid entries are identical with and without the mask
        assert torch.equal(masked[:, :4], default[:, :4])
        a = make_realistic_dataset(w, None, "all", 4, seed=2, chunk=3)
        b = make_realistic_dataset(w, None, "all", 4, seed=2, chunk=3, valid=valid)
        assert (a[:, 4:] > 0).any() and int(b[:, 4:].sum()) == 0

    def test_historical_corruption_names_still_resolve(self):
        import subprocess
        import sys
        from pathlib import Path

        import icenine.orientation_eval as oe

        assert oe.CorruptionConfig is oe.RealismConfig
        assert oe.corrupt_windows is oe.make_realistic_windows
        assert oe.corrupt_dataset is oe.make_realistic_dataset
        scripts = Path(__file__).resolve().parents[1] / "scripts"
        for script, new, old in [
            ("gauss_newton_baseline.py", "--realistic", "--corrupt"),
            ("train_toy_orientation_nn.py", "--realistic-train", "--corrupt-train"),
        ]:
            out = subprocess.run(
                [sys.executable, str(scripts / script), "--help"], capture_output=True, text=True
            ).stdout
            assert new in out and old in out  # both option strings share one argument

    def test_inv3_matches_linalg_inv_and_survives_singular_input(self):
        from icenine.toy_orientation_model import inv3

        torch.manual_seed(0)
        A = torch.randn(50, 3, 3, dtype=torch.float64)
        A = A @ A.transpose(-1, -2) + torch.eye(3, dtype=torch.float64)
        assert torch.allclose(inv3(A), torch.linalg.inv(A), atol=1e-9)
        assert torch.isfinite(inv3(torch.zeros(1, 3, 3, dtype=torch.float64))).all()
