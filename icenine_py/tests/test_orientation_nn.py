"""
Tests for the toy-NN windowed local-refinement pipeline (orientation_nn.py).

Uses Example2.ThreeVoxels' ground-truth .mic as a real voxel/geometry, but
never touches ScatteringData_Python -- these tests only exercise the
synthetic windowed renderer, not experimental-data loading.

Not a full training test: see scripts/generate_toy_orientation_dataset.py and
scripts/train_toy_orientation_nn.py for the smoke-test / scale-up workflow.
"""

import os
from pathlib import Path

import numpy as np
import pytest
import torch

from icenine.config_file import ConfigFile
from icenine.experiment_setup import XDMExperimentSetup
from icenine.mic_file import MicFile
from icenine.orientation_nn import (
    ROIPeak,
    define_roi_set,
    quat_misorientation_deg_batch,
    quaternion_regression_loss,
    render_local_windows,
    sample_local_perturbations,
)
from icenine.reconstructor import _get_voxel_vertices
from icenine.sample import Sample
from icenine.simulation import Simulation


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def project_root():
    return Path(__file__).parent.parent.parent


@pytest.fixture(autouse=True)
def chdir_to_project_root(project_root, monkeypatch):
    monkeypatch.chdir(project_root)


@pytest.fixture
def example_dir(project_root):
    return project_root / "Examples" / "Example2.ThreeVoxels"


@pytest.fixture
def physics_setup(example_dir):
    """Load config/experiment/sample physics -- no experimental data needed."""
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    if not config_path.exists():
        pytest.skip(f"Config not found: {config_path}")

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

    return sample, detector_list, range_map, exp_setup, simulator, structure_list


@pytest.fixture
def ground_truth_mic(example_dir):
    mic_path = example_dir / "SimInput" / "three_voxels.mic"
    if not mic_path.exists():
        pytest.skip(f"Ground truth .mic not found: {mic_path}")
    return MicFile.read(str(mic_path))


@pytest.fixture
def roi_fixture(physics_setup, ground_truth_mic):
    """First voxel with a non-empty ROI set, plus its fixed ROI peaks."""
    sample, detector_list, range_map, exp_setup, simulator, structure_list = physics_setup

    for voxel in ground_truth_mic.voxels:
        vertices = _get_voxel_vertices(voxel)
        orientation = torch.from_numpy(voxel.orientation).float()
        roi_list = define_roi_set(
            orientation,
            vertices,
            sample,
            detector_list,
            range_map,
            exp_setup,
            structure_list,
            simulator,
            phase_index=voxel.phase,
        )
        if roi_list:
            return (
                voxel,
                vertices,
                orientation,
                roi_list,
                (sample, detector_list, range_map, exp_setup, simulator),
            )

    pytest.skip("No voxel in three_voxels.mic produced a non-empty ROI set")


# ============================================================================
# ROI definition
# ============================================================================


class TestDefineROISet:
    def test_roi_set_nonempty(self, roi_fixture):
        _voxel, _vertices, _orientation, roi_list, _ctx = roi_fixture
        assert len(roi_list) > 0
        assert all(isinstance(roi, ROIPeak) for roi in roi_list)

    def test_roi_peaks_have_valid_branch_and_detector(self, roi_fixture):
        _voxel, _vertices, _orientation, roi_list, ctx = roi_fixture
        _sample, detector_list, _range_map, _exp_setup, _simulator = ctx
        for roi in roi_list:
            assert roi.omega_branch in (1, 2)
            assert 0 <= roi.detector_index < len(detector_list)

    def test_roi_set_deterministic(self, roi_fixture, physics_setup):
        """Re-running define_roi_set at the same orientation reproduces the same
        peak identities and doesn't leave the sample's rotation state mutated."""
        voxel, vertices, orientation, roi_list, ctx = roi_fixture
        sample, detector_list, range_map, exp_setup, simulator = ctx

        _sample2, _detector_list2, _range_map2, _exp_setup2, _simulator2, structure_list = (
            physics_setup
        )

        rotation_before = sample.sample_to_lab_matrix.clone()
        roi_list_2 = define_roi_set(
            orientation,
            vertices,
            sample,
            detector_list,
            range_map,
            exp_setup,
            structure_list,
            simulator,
            phase_index=voxel.phase,
        )
        assert torch.allclose(sample.sample_to_lab_matrix, rotation_before)

        ids_1 = [(r.reflection_index, r.omega_branch, r.detector_index) for r in roi_list]
        ids_2 = [(r.reflection_index, r.omega_branch, r.detector_index) for r in roi_list_2]
        assert ids_1 == ids_2


# ============================================================================
# Windowed rendering
# ============================================================================


class TestRenderLocalWindows:
    def test_window_shape(self, roi_fixture):
        voxel, vertices, orientation, roi_list, ctx = roi_fixture
        sample, detector_list, range_map, exp_setup, simulator = ctx
        window_size = 32

        windows, missing = render_local_windows(
            orientation,
            roi_list,
            vertices,
            sample,
            detector_list,
            range_map,
            exp_setup,
            simulator,
            window_size=window_size,
        )
        assert windows.shape == (len(roi_list), window_size, window_size)
        assert missing.shape == (len(roi_list),)
        assert windows.dtype == torch.float32

    def test_nominal_orientation_has_no_missing_peaks(self, roi_fixture):
        """Rendering at the exact orientation used to define the ROI set should
        reproduce every peak (that's the whole point of fixing peak identity)."""
        voxel, vertices, orientation, roi_list, ctx = roi_fixture
        sample, detector_list, range_map, exp_setup, simulator = ctx

        windows, missing = render_local_windows(
            orientation,
            roi_list,
            vertices,
            sample,
            detector_list,
            range_map,
            exp_setup,
            simulator,
            window_size=32,
        )
        assert not missing.any()
        # Every rendered window should have some nonzero signal near its center.
        assert (windows.sum(dim=(1, 2)) > 0).all()

    def test_sample_rotation_state_restored(self, roi_fixture):
        voxel, vertices, orientation, roi_list, ctx = roi_fixture
        sample, detector_list, range_map, exp_setup, simulator = ctx

        rotation_before = sample.sample_to_lab_matrix.clone()
        render_local_windows(
            orientation,
            roi_list,
            vertices,
            sample,
            detector_list,
            range_map,
            exp_setup,
            simulator,
            window_size=32,
        )
        assert torch.allclose(sample.sample_to_lab_matrix, rotation_before)

    def test_roi_stability_across_perturbation_range(self, roi_fixture):
        """Most ROI peaks should stay observable under small perturbations --
        if this fails often in practice, the plan calls for falling back to
        per-sample ROI re-derivation with padding/masking (see MIGRATION_HISTORY.md)."""
        voxel, vertices, orientation, roi_list, ctx = roi_fixture
        sample, detector_list, range_map, exp_setup, simulator = ctx

        rng = np.random.default_rng(0)
        matrices, _quats = sample_local_perturbations(orientation, 20, max_angle_deg=1.0, rng=rng)

        total_missing = 0
        total_peaks = 0
        for mat in matrices:
            perturbed = torch.from_numpy(mat).float()
            _windows, missing = render_local_windows(
                perturbed,
                roi_list,
                vertices,
                sample,
                detector_list,
                range_map,
                exp_setup,
                simulator,
                window_size=32,
            )
            total_missing += int(missing.sum().item())
            total_peaks += len(roi_list)

        dropout_rate = total_missing / total_peaks
        assert (
            dropout_rate < 0.5
        ), f"ROI set too unstable under 1deg perturbations: {dropout_rate:.1%} dropout"


# ============================================================================
# Perturbation sampling
# ============================================================================


class TestSampleLocalPerturbations:
    def test_returns_valid_rotation_matrices(self, roi_fixture):
        voxel, _vertices, orientation, _roi_list, _ctx = roi_fixture
        rng = np.random.default_rng(1)
        matrices, quats = sample_local_perturbations(orientation, 10, max_angle_deg=2.0, rng=rng)

        assert len(matrices) == 10
        assert len(quats) == 10
        for mat in matrices:
            assert mat.shape == (3, 3)
            # proper rotation: R^T R = I, det = 1
            assert np.allclose(mat.T @ mat, np.eye(3), atol=1e-5)
            assert np.isclose(np.linalg.det(mat), 1.0, atol=1e-4)

    def test_perturbations_stay_close_to_nominal(self, roi_fixture):
        from icenine.sampling import matrix_to_quaternion

        voxel, _vertices, orientation, _roi_list, _ctx = roi_fixture
        nominal_q = matrix_to_quaternion(orientation.numpy())

        rng = np.random.default_rng(2)
        _matrices, quats = sample_local_perturbations(orientation, 50, max_angle_deg=2.0, rng=rng)

        angles_deg = quat_misorientation_deg_batch(
            torch.from_numpy(nominal_q).float().unsqueeze(0).expand(50, -1),
            torch.from_numpy(np.stack(quats)).float(),
        )
        assert (angles_deg < 10.0).all()  # generous bound, radius formula is approximate


# ============================================================================
# Loss / metric
# ============================================================================


class TestLossAndMetric:
    def test_loss_zero_at_identical_quaternions(self):
        q = torch.tensor([[1.0, 0.0, 0.0, 0.0]])
        loss = quaternion_regression_loss(q, q)
        assert torch.isclose(loss, torch.tensor(0.0), atol=1e-6)

    def test_loss_handles_sign_ambiguity(self):
        """q and -q represent the same rotation; loss must be sign-invariant."""
        q1 = torch.tensor([[1.0, 0.0, 0.0, 0.0]])
        q2 = torch.tensor([[-1.0, 0.0, 0.0, 0.0]])
        loss = quaternion_regression_loss(q1, q2)
        assert torch.isclose(loss, torch.tensor(0.0), atol=1e-6)

    def test_loss_gradient_no_nan_near_optimum(self):
        """Regression check for the arccos-singularity this loss is designed to avoid."""
        q_true = torch.tensor([[1.0, 0.0, 0.0, 0.0]])
        q_pred = torch.tensor([[0.9999, 0.01, 0.0, 0.0]], requires_grad=True)
        q_pred_n = q_pred / q_pred.norm(dim=-1, keepdim=True)

        loss = quaternion_regression_loss(q_pred_n, q_true)
        loss.backward()

        assert torch.isfinite(loss)
        assert torch.isfinite(q_pred.grad).all()

    def test_misorientation_matches_quat_misorientation_deg(self):
        from icenine.orientation_search import _quat_misorientation_deg

        rng = np.random.default_rng(3)
        q1 = rng.normal(size=4)
        q1 /= np.linalg.norm(q1)
        q2 = rng.normal(size=4)
        q2 /= np.linalg.norm(q2)

        expected = _quat_misorientation_deg(q1, q2)
        actual = quat_misorientation_deg_batch(
            torch.from_numpy(q1).float().unsqueeze(0), torch.from_numpy(q2).float().unsqueeze(0)
        ).item()
        assert actual == pytest.approx(expected, abs=1e-3)


# ============================================================================
# Model sanity (no full training)
# ============================================================================


class TestToyOrientationNet:
    def test_output_is_unit_quaternion(self, roi_fixture):
        from icenine.toy_orientation_model import ToyOrientationNet

        voxel, vertices, orientation, roi_list, ctx = roi_fixture
        sample, detector_list, range_map, exp_setup, simulator = ctx
        window_size = 32

        windows, _missing = render_local_windows(
            orientation,
            roi_list,
            vertices,
            sample,
            detector_list,
            range_map,
            exp_setup,
            simulator,
            window_size=window_size,
        )
        model = ToyOrientationNet(n_peaks=len(roi_list), window_size=window_size)
        q_pred = model(windows.unsqueeze(0))

        assert q_pred.shape == (1, 4)
        assert torch.isclose(q_pred.norm(), torch.tensor(1.0), atol=1e-5)


def test_split_by_voxel_disjoint_deterministic_and_spans_r_perp():
    from icenine.orientation_nn import split_by_voxel

    r_perp = np.linspace(0.0, 500.0, 24)
    vid = np.repeat(np.arange(24), 5)
    tr1, va1, vv1 = split_by_voxel(vid, r_perp, 4, seed=3)
    tr2, va2, vv2 = split_by_voxel(vid, r_perp, 4, seed=3)
    assert np.array_equal(tr1, tr2) and np.array_equal(va1, va2) and vv1 == vv2
    assert len(vv1) == 4 and len(va1) == 20 and len(tr1) + len(va1) == len(vid)
    assert not set(vid[tr1]) & set(vid[va1])  # no voxel on both sides
    assert set(vid[va1]) == set(vv1)
    # one voxel per r_perp quartile
    assert [v // 6 for v in vv1] == [0, 1, 2, 3]
    with pytest.raises(ValueError):
        split_by_voxel(vid, r_perp, 24)
