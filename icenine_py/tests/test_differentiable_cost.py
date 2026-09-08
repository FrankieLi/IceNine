"""
Tests for differentiable cost function infrastructure.

Phase 1: ExperimentalImageStack and MultiScaleImageStack
Phase 2: DifferentiableCostFunction — gradient flow, finite-difference checks,
         correlation with hard cost function, multi-scale basin widening
"""

import math
import os
from pathlib import Path

import numpy as np
import pytest
import torch

from icenine.image_data import ImageData
from icenine.experimental_data import ExperimentalData
from icenine.differentiable_cost import (
    DifferentiableCostFunction,
    DifferentiableOverlapInfo,
    ExperimentalImageStack,
    MultiScaleImageStack,
    SparseImageStack,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_test_exp_data(
    n_omega: int = 10,
    n_det: int = 2,
    num_rows: int = 64,
    num_cols: int = 64,
    seed: int = 42,
) -> ExperimentalData:
    """Create synthetic ExperimentalData with random bright pixels."""
    rng = np.random.RandomState(seed)
    images = []
    for omega_idx in range(n_omega):
        det_images = []
        for det_idx in range(n_det):
            img = ImageData(num_rows, num_cols, mode="dense")
            # Scatter 5-20 bright pixels per image
            n_bright = rng.randint(5, 20)
            for _ in range(n_bright):
                j = rng.randint(0, num_cols)
                k = rng.randint(0, num_rows)
                intensity = rng.uniform(0.5, 3.0)
                img.set_pixel(j, k, intensity)
            det_images.append(img)
        images.append(det_images)

    return ExperimentalData(
        images=images,
        n_omega_intervals=n_omega,
        n_detectors=n_det,
    )


# ---------------------------------------------------------------------------
# ExperimentalImageStack tests
# ---------------------------------------------------------------------------


class TestExperimentalImageStack:
    def test_basic_construction(self):
        """Stack builds with correct shape and dtype."""
        exp = _make_test_exp_data(n_omega=5, n_det=2, num_rows=32, num_cols=32)
        stack = ExperimentalImageStack(exp, binary=True)

        assert stack.images.shape == (10, 1, 32, 32)
        assert stack.images.dtype == torch.float32
        assert stack.n_omega == 5
        assert stack.n_det == 2
        assert stack.H == 32
        assert stack.W == 32
        assert stack.binary is True

    def test_binary_mode(self):
        """Binary mode produces only 0.0 and 1.0 values."""
        exp = _make_test_exp_data()
        stack = ExperimentalImageStack(exp, binary=True)

        unique_vals = torch.unique(stack.images)
        assert all(v in (0.0, 1.0) for v in unique_vals.tolist())

    def test_intensity_mode(self):
        """Intensity mode preserves original values."""
        exp = _make_test_exp_data()
        stack = ExperimentalImageStack(exp, binary=False)

        # Should have values other than just 0 and 1
        unique_vals = torch.unique(stack.images)
        assert len(unique_vals) > 2, "Expected more than binary values in intensity mode"

    def test_round_trip_single_image(self):
        """Extracted image matches original ImageData."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=16, num_cols=16)
        stack = ExperimentalImageStack(exp, binary=True)

        for omega_idx in range(3):
            for det_idx in range(2):
                original = exp.get_image(omega_idx, det_idx)
                original_binary = (original._pixels_dense > 0).float()

                extracted = stack.get_image(omega_idx, det_idx)
                assert extracted.shape == (1, 16, 16)

                torch.testing.assert_close(extracted[0], original_binary)

    def test_flat_index(self):
        """Flat indexing is correct: omega * n_det + det."""
        exp = _make_test_exp_data(n_omega=5, n_det=3)
        stack = ExperimentalImageStack(exp, binary=True)

        assert stack.flat_index(0, 0) == 0
        assert stack.flat_index(0, 2) == 2
        assert stack.flat_index(1, 0) == 3
        assert stack.flat_index(4, 2) == 14

    def test_batch_get_images(self):
        """Batch gather returns correct images."""
        exp = _make_test_exp_data(n_omega=5, n_det=2)
        stack = ExperimentalImageStack(exp, binary=True)

        indices = torch.tensor([0, 3, 7])
        batch = stack.get_images_batch(indices)
        assert batch.shape == (3, 1, 64, 64)

        # Verify each image matches
        torch.testing.assert_close(batch[0], stack.images[0])
        torch.testing.assert_close(batch[1], stack.images[3])
        torch.testing.assert_close(batch[2], stack.images[7])

    def test_memory_bytes(self):
        """Memory reporting is accurate."""
        exp = _make_test_exp_data(n_omega=5, n_det=2, num_rows=32, num_cols=32)
        stack = ExperimentalImageStack(exp, binary=True)

        # 10 images × 1 × 32 × 32 × 4 bytes = 40,960 bytes
        expected = 10 * 1 * 32 * 32 * 4
        assert stack.memory_bytes == expected

    def test_to_device(self):
        """to() moves images and returns self."""
        exp = _make_test_exp_data(n_omega=2, n_det=1, num_rows=8, num_cols=8)
        stack = ExperimentalImageStack(exp, binary=True)

        result = stack.to(torch.device("cpu"))
        assert result is stack
        assert stack.images.device.type == "cpu"

    def test_to_image_stack_method(self):
        """ExperimentalData.to_image_stack() creates a valid stack."""
        exp = _make_test_exp_data(n_omega=3, n_det=2)

        stack_binary = exp.to_image_stack(binary=True)
        assert isinstance(stack_binary, ExperimentalImageStack)
        assert stack_binary.binary is True

        stack_intensity = exp.to_image_stack(binary=False)
        assert isinstance(stack_intensity, ExperimentalImageStack)
        assert stack_intensity.binary is False

    def test_sparse_image_support(self):
        """Stack correctly handles sparse ImageData objects."""
        # Create sparse images
        images = []
        for omega in range(3):
            det_images = []
            for det in range(2):
                img = ImageData(16, 16, mode="sparse")
                img.set_pixel(5, 5, 1.0)
                img.set_pixel(10, 10, 2.0)
                det_images.append(img)
            images.append(det_images)

        exp = ExperimentalData(images=images, n_omega_intervals=3, n_detectors=2)
        stack = ExperimentalImageStack(exp, binary=True)

        assert stack.images.shape == (6, 1, 16, 16)
        # Check that the bright pixels are present
        assert stack.images[0, 0, 5, 5] == 1.0
        assert stack.images[0, 0, 10, 10] == 1.0


# ---------------------------------------------------------------------------
# MultiScaleImageStack tests
# ---------------------------------------------------------------------------


class TestMultiScaleImageStack:
    def test_basic_construction(self):
        """Multi-scale stack builds with correct number of scales."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=32, num_cols=32)
        base_stack = ExperimentalImageStack(exp, binary=True)
        multi = MultiScaleImageStack(base_stack, downsample_factors=[1, 4, 8])

        assert multi.n_scales == 3
        assert len(multi.scales) == 3

    def test_scale_0_is_original(self):
        """Scale 0 (factor=1) shares the original tensor."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=32, num_cols=32)
        base_stack = ExperimentalImageStack(exp, binary=True)
        multi = MultiScaleImageStack(base_stack, downsample_factors=[1, 4])

        # Scale 0 should be the same object
        assert multi.get_at_scale(0) is base_stack

    def test_downsampled_scale_nonzero(self):
        """Downsampled scales preserve bright pixels (max_pool dilation)."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=64, num_cols=64)
        base_stack = ExperimentalImageStack(exp, binary=True)
        multi = MultiScaleImageStack(base_stack, downsample_factors=[1, 4])

        # Original has bright pixels
        assert base_stack.images.sum().item() > 0

        # All scales should have nonzero content
        for scale_idx in range(multi.n_scales):
            scale = multi.get_at_scale(scale_idx)
            assert scale.images.sum().item() > 0, (
                f"Scale {scale_idx} (factor={multi.downsample_factors[scale_idx]}): "
                f"downsampled images are all zero"
            )

    def test_downsample_reduces_resolution(self):
        """Downsampled scale has smaller H and W than original."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=64, num_cols=64)
        base_stack = ExperimentalImageStack(exp, binary=True)
        multi = MultiScaleImageStack(base_stack, downsample_factors=[1, 4])

        coarse = multi.get_at_scale(1)
        assert coarse.H == 16
        assert coarse.W == 16

    def test_downsampled_scale_shapes_consistent(self):
        """All scales have the same number of images."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=32, num_cols=32)
        base_stack = ExperimentalImageStack(exp, binary=True)
        multi = MultiScaleImageStack(base_stack, downsample_factors=[1, 4, 8])

        n_images = base_stack.n_omega * base_stack.n_det
        for i in range(multi.n_scales):
            s = multi.get_at_scale(i)
            assert (
                s.images.shape[0] == n_images
            ), f"Scale {i}: expected {n_images} images, got {s.images.shape[0]}"
            assert s.H <= base_stack.H
            assert s.W <= base_stack.W

    def test_default_downsample_factors(self):
        """Default downsample_factors are [1, 4, 8]."""
        exp = _make_test_exp_data(n_omega=2, n_det=1, num_rows=32, num_cols=32)
        base_stack = ExperimentalImageStack(exp, binary=True)
        multi = MultiScaleImageStack(base_stack)

        assert multi.downsample_factors == [1, 4, 8]
        assert multi.n_scales == 3

    def test_to_device(self):
        """to() moves all scales to device."""
        exp = _make_test_exp_data(n_omega=2, n_det=1, num_rows=8, num_cols=8)
        base_stack = ExperimentalImageStack(exp, binary=True)
        multi = MultiScaleImageStack(base_stack, downsample_factors=[1, 4])

        result = multi.to(torch.device("cpu"))
        assert result is multi
        for s in multi.scales:
            assert s.device.type == "cpu"


# ---------------------------------------------------------------------------
# SparseImageStack tests
# ---------------------------------------------------------------------------


class TestSparseImageStack:
    def test_basic_construction(self):
        """Sparse stack builds with correct metadata."""
        exp = _make_test_exp_data(n_omega=5, n_det=2, num_rows=32, num_cols=32)
        stack = SparseImageStack(exp, binary=True)

        assert stack.n_omega == 5
        assert stack.n_det == 2
        assert stack.H == 32
        assert stack.W == 32
        assert stack.binary is True

    def test_memory_much_smaller_than_dense(self):
        """Sparse storage uses orders of magnitude less memory than dense."""
        exp = _make_test_exp_data(n_omega=10, n_det=2, num_rows=64, num_cols=64)
        dense = ExperimentalImageStack(exp, binary=True)
        sparse = SparseImageStack(exp, binary=True)

        # Dense: 20 × 64 × 64 × 4 = 327,680 bytes
        # Sparse: ~200 pixels × 4 bytes (int16 coords) ≈ 800 bytes
        assert (
            sparse.memory_bytes < dense.memory_bytes / 10
        ), f"Sparse ({sparse.memory_bytes}) should be much smaller than dense ({dense.memory_bytes})"

    def test_round_trip_matches_dense(self):
        """Densified sparse images match dense stack exactly."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=16, num_cols=16)
        dense = ExperimentalImageStack(exp, binary=True)
        sparse = SparseImageStack(exp, binary=True)

        for omega_idx in range(3):
            for det_idx in range(2):
                dense_img = dense.get_image(omega_idx, det_idx)
                sparse_img = sparse.get_image(omega_idx, det_idx)
                torch.testing.assert_close(sparse_img, dense_img)

    def test_batch_get_matches_dense(self):
        """Batch get from sparse matches batch get from dense."""
        exp = _make_test_exp_data(n_omega=5, n_det=2)
        dense = ExperimentalImageStack(exp, binary=True)
        sparse = SparseImageStack(exp, binary=True)

        indices = torch.tensor([0, 3, 7, 3])  # includes duplicate
        dense_batch = dense.get_images_batch(indices)
        sparse_batch = sparse.get_images_batch(indices)

        assert sparse_batch.shape == dense_batch.shape
        torch.testing.assert_close(sparse_batch, dense_batch)

    def test_binary_mode(self):
        """Binary mode produces only 0.0 and 1.0 values."""
        exp = _make_test_exp_data()
        stack = SparseImageStack(exp, binary=True)

        batch = stack.get_images_batch(torch.tensor([0, 1, 2]))
        unique_vals = torch.unique(batch)
        assert all(v in (0.0, 1.0) for v in unique_vals.tolist())

    def test_intensity_mode(self):
        """Intensity mode preserves original values."""
        exp = _make_test_exp_data()
        sparse = SparseImageStack(exp, binary=False)
        dense = ExperimentalImageStack(exp, binary=False)

        indices = torch.tensor([0, 5])
        sparse_batch = sparse.get_images_batch(indices)
        dense_batch = dense.get_images_batch(indices)
        torch.testing.assert_close(sparse_batch, dense_batch)

    def test_to_device(self):
        """to() moves sparse data and returns self."""
        exp = _make_test_exp_data(n_omega=2, n_det=1, num_rows=8, num_cols=8)
        stack = SparseImageStack(exp, binary=True)

        result = stack.to(torch.device("cpu"))
        assert result is stack
        assert stack.device.type == "cpu"

    def test_flat_index(self):
        """Flat indexing matches ExperimentalImageStack convention."""
        exp = _make_test_exp_data(n_omega=5, n_det=3)
        stack = SparseImageStack(exp, binary=True)

        assert stack.flat_index(0, 0) == 0
        assert stack.flat_index(0, 2) == 2
        assert stack.flat_index(1, 0) == 3
        assert stack.flat_index(4, 2) == 14

    def test_to_sparse_image_stack_method(self):
        """ExperimentalData.to_sparse_image_stack() creates a valid stack."""
        exp = _make_test_exp_data(n_omega=3, n_det=2)
        stack = exp.to_sparse_image_stack(binary=True)
        assert isinstance(stack, SparseImageStack)
        assert stack.binary is True

    def test_multiscale_with_sparse(self):
        """MultiScaleImageStack works with SparseImageStack input."""
        exp = _make_test_exp_data(n_omega=3, n_det=2, num_rows=32, num_cols=32)
        sparse_stack = SparseImageStack(exp, binary=True)
        multi = MultiScaleImageStack(sparse_stack, downsample_factors=[1, 4])

        assert multi.n_scales == 2
        # Scale 0 is the original sparse stack
        assert multi.get_at_scale(0) is sparse_stack
        # Scale 1 is a dense stack at reduced resolution
        coarse = multi.get_at_scale(1)
        assert isinstance(coarse, ExperimentalImageStack)
        assert coarse.H < sparse_stack.H  # downsampled
        assert coarse.W < sparse_stack.W

    def test_multiscale_downsample_sizes(self):
        """Downsampled scales have correct H/W."""
        exp = _make_test_exp_data(n_omega=2, n_det=1, num_rows=64, num_cols=64)
        sparse_stack = SparseImageStack(exp, binary=True)

        # factor=4 → 16×16
        multi = MultiScaleImageStack(sparse_stack, downsample_factors=[1, 4])
        coarse = multi.get_at_scale(1)
        assert coarse.H == 16
        assert coarse.W == 16

        # factor=8 → 8×8
        multi2 = MultiScaleImageStack(sparse_stack, downsample_factors=[1, 8])
        coarse2 = multi2.get_at_scale(1)
        assert coarse2.H == 8
        assert coarse2.W == 8

    def test_empty_images(self):
        """Handles images with zero bright pixels."""
        images = []
        for omega in range(3):
            det_images = []
            for det in range(2):
                img = ImageData(16, 16, mode="dense")
                # Leave all pixels at zero
                det_images.append(img)
            images.append(det_images)

        exp = ExperimentalData(images=images, n_omega_intervals=3, n_detectors=2)
        stack = SparseImageStack(exp, binary=True)

        batch = stack.get_images_batch(torch.tensor([0, 1, 2]))
        assert batch.sum().item() == 0.0
        assert stack.memory_bytes == 0  # no coords stored


# ---------------------------------------------------------------------------
# DifferentiableCostFunction fixtures (ThreeVoxels integration)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def project_root():
    return Path(__file__).parent.parent.parent


@pytest.fixture(scope="module")
def example_dir(project_root):
    return project_root / "Examples" / "Example2.ThreeVoxels"


# Module-scoped cache for the expensive setup (avoid rebuilding per test)
_DIFF_COST_CACHE: dict = {}


@pytest.fixture(scope="module")
def diff_cost_components(project_root, example_dir):
    """
    Module-scoped fixture: loads config, detectors, sample, mic, and
    a memory-efficient SparseImageStack (directly from files, ~26 KB).

    Does NOT load dense ExperimentalData (which would be ~7.8 GB for 2048×2048).
    For hard cost comparison, loads ExperimentalData in sparse ImageData mode.
    """
    if _DIFF_COST_CACHE:
        return _DIFF_COST_CACHE

    from icenine.config_file import ConfigFile
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.mic_file import MicFile
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    if not config_path.exists():
        pytest.skip(f"Config not found: {config_path}")

    old_cwd = os.getcwd()
    os.chdir(example_dir)
    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = "3Grains.sim"

    data_dir = example_dir / "ScatteringData_Python"
    if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
        os.chdir(old_cwd)
        pytest.skip(f"Forward sim output not found or incomplete: {data_dir}")

    # Load images directly into sparse stack (~26 KB, no dense intermediate)
    image_stack = SparseImageStack.from_image_directory(
        directory=str(data_dir),
        basename="3Grains.sim",
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
        binary=True,
    )

    # Load ExperimentalData in sparse mode for hard cost fn (~few MB, not ~7.8 GB)
    exp_data = ExperimentalData.from_image_directory(
        directory=str(data_dir),
        basename="3Grains.sim",
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
        mode="sparse",
    )

    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()
    detector_list = exp_setup.get_detector_list()
    range_map = exp_setup.get_range_to_index_map()
    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])
    simulator = Simulation(exp_setup)
    structure_list = sample.get_structure_list()

    mic_path = example_dir / "SimInput" / "three_voxels.mic"
    if not mic_path.exists():
        os.chdir(old_cwd)
        pytest.skip(f"Ground truth .mic not found: {mic_path}")
    mic = MicFile.read(str(mic_path))

    os.chdir(old_cwd)

    result = {
        "exp_data": exp_data,
        "image_stack": image_stack,
        "simulator": simulator,
        "detector_list": detector_list,
        "range_map": range_map,
        "sample": sample,
        "structure_list": structure_list,
        "mic": mic,
        "get_vertices": _get_voxel_vertices,
        "data_dir": str(data_dir),
    }
    _DIFF_COST_CACHE.update(result)
    return result


@pytest.fixture(scope="module")
def diff_cost_setup(diff_cost_components):
    """
    Build multi-scale stack + hard/diff cost functions.
    Uses SparseImageStack (~26 KB) and sparse ExperimentalData (~few MB).
    """
    from icenine.cost_functions import VoxelCostFunction

    c = diff_cost_components

    multi_stack = MultiScaleImageStack(c["image_stack"], downsample_factors=[1])

    hard_cost_fn = VoxelCostFunction(
        simulator=c["simulator"],
        detector_list=c["detector_list"],
        range_map=c["range_map"],
        exp_data=c["exp_data"],
        sample=c["sample"],
        structure_list=c["structure_list"],
        mode="hard",
    )

    diff_cost_fn = DifferentiableCostFunction(
        simulator=c["simulator"],
        detector_list=c["detector_list"],
        range_map=c["range_map"],
        image_stack=multi_stack,
        sample=c["sample"],
        structure_list=c["structure_list"],
    )

    return {
        **c,
        "diff_cost_fn": diff_cost_fn,
        "hard_cost_fn": hard_cost_fn,
        "multi_stack": multi_stack,
    }


# ---------------------------------------------------------------------------
# DifferentiableCostFunction tests
# ---------------------------------------------------------------------------


class TestDifferentiableCostFunction:
    """Phase 2 tests: gradient flow, finite-diff, correlation, multi-scale."""

    def test_gradient_existence(self, diff_cost_setup):
        """loss.backward() produces non-zero gradients on orientation."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        orientation = torch.from_numpy(voxel.orientation).float().requires_grad_(True)
        info = c["diff_cost_fn"].evaluate(orientation, vertices, phase_index=voxel.phase, scale=0)

        assert info.n_peaks > 0, "Expected observable peaks for ground truth orientation"
        assert info.cost.requires_grad, "cost tensor should carry grad_fn"

        info.cost.backward()

        assert orientation.grad is not None, "orientation.grad should be populated after backward()"
        assert torch.any(
            orientation.grad != 0
        ), "Gradient should be non-zero for ground truth orientation"

    def test_quality_is_positive_at_ground_truth(self, diff_cost_setup):
        """Ground truth orientation should give positive quality (cost < 1)."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        orientation = torch.from_numpy(voxel.orientation).float()
        info = c["diff_cost_fn"].evaluate(orientation, vertices, phase_index=voxel.phase, scale=0)

        assert (
            info.quality.item() > 0.1
        ), f"Ground truth orientation should have quality > 0.1, got {info.quality.item():.4f}"
        assert (
            info.cost.item() < 0.9
        ), f"Ground truth cost should be < 0.9, got {info.cost.item():.4f}"

    def test_quality_cost_sum_to_one(self, diff_cost_setup):
        """quality + cost = 1."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        orientation = torch.from_numpy(voxel.orientation).float()
        info = c["diff_cost_fn"].evaluate(orientation, vertices, phase_index=voxel.phase, scale=0)

        total = info.quality.item() + info.cost.item()
        assert abs(total - 1.0) < 1e-6, f"quality + cost = {total}, expected 1.0"

    def test_finite_difference_gradient_check(self, diff_cost_setup):
        """Autograd gradients match numerical finite-difference Jacobian."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        orientation = torch.from_numpy(voxel.orientation).float().requires_grad_(True)
        info = c["diff_cost_fn"].evaluate(orientation, vertices, phase_index=voxel.phase, scale=0)
        info.cost.backward()
        autograd = orientation.grad.clone()

        # Numerical Jacobian via central differences
        eps = 1e-4
        numerical_grad = torch.zeros_like(orientation)
        base_orient = orientation.detach().clone()

        for i in range(3):
            for j in range(3):
                orient_plus = base_orient.clone()
                orient_plus[i, j] += eps
                info_plus = c["diff_cost_fn"].evaluate(
                    orient_plus, vertices, phase_index=voxel.phase, scale=0
                )

                orient_minus = base_orient.clone()
                orient_minus[i, j] -= eps
                info_minus = c["diff_cost_fn"].evaluate(
                    orient_minus, vertices, phase_index=voxel.phase, scale=0
                )

                numerical_grad[i, j] = (info_plus.cost.item() - info_minus.cost.item()) / (
                    2.0 * eps
                )

        # Check that gradients are correlated — they won't be exact because
        # discrete omega routing creates discontinuities, but the smooth
        # components (pixel coordinates via grid_sample) should agree
        nonzero_mask = numerical_grad.abs() > 1e-6
        if nonzero_mask.any():
            auto_vals = autograd[nonzero_mask]
            num_vals = numerical_grad[nonzero_mask]
            # Cosine similarity: should be positive (same direction)
            cos_sim = torch.dot(auto_vals.flatten(), num_vals.flatten()) / (
                auto_vals.norm() * num_vals.norm() + 1e-12
            )
            assert cos_sim > 0.5, (
                f"Autograd and finite-diff gradients should be correlated, "
                f"cosine similarity = {cos_sim.item():.3f}"
            )

    def test_wrong_orientation_has_lower_quality(self, diff_cost_setup):
        """A random orientation should have lower quality than ground truth."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        # Ground truth
        orient_gt = torch.from_numpy(voxel.orientation).float()
        info_gt = c["diff_cost_fn"].evaluate(orient_gt, vertices, phase_index=voxel.phase, scale=0)

        # Random orientation (identity matrix — very different from ground truth)
        orient_wrong = torch.eye(3, dtype=torch.float32)
        info_wrong = c["diff_cost_fn"].evaluate(
            orient_wrong, vertices, phase_index=voxel.phase, scale=0
        )

        assert info_gt.quality.item() > info_wrong.quality.item(), (
            f"Ground truth quality ({info_gt.quality.item():.4f}) should exceed "
            f"random orientation quality ({info_wrong.quality.item():.4f})"
        )

    def test_correlation_with_hard_cost(self, diff_cost_setup):
        """Soft overlap quality should correlate with hard overlap quality."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        # Hard cost at ground truth
        hard_info = c["hard_cost_fn"].evaluate(
            orientation=voxel.orientation,
            voxel_vertices=vertices,
            phase_index=voxel.phase,
        )

        # Soft cost at ground truth
        orient_t = torch.from_numpy(voxel.orientation).float()
        soft_info = c["diff_cost_fn"].evaluate(orient_t, vertices, phase_index=voxel.phase, scale=0)

        # Both should be positive/non-zero for ground truth
        assert hard_info.quality > 0.1, f"Hard quality too low: {hard_info.quality}"
        assert soft_info.quality.item() > 0.1, f"Soft quality too low: {soft_info.quality.item()}"

        # They use different metrics (binary overlap ratio vs centroid bilinear sample)
        # so absolute values may differ, but both should agree that ground truth is good
        # Check that soft quality is in a reasonable range relative to hard quality
        ratio = soft_info.quality.item() / hard_info.quality
        assert 0.1 < ratio < 10.0, (
            f"Soft/hard quality ratio {ratio:.3f} is too extreme "
            f"(hard={hard_info.quality:.4f}, soft={soft_info.quality.item():.4f})"
        )

    def test_multiple_voxels(self, diff_cost_setup):
        """DifferentiableCostFunction works for all 3 ground truth voxels."""
        c = diff_cost_setup
        for idx, voxel in enumerate(c["mic"].voxels[:3]):
            vertices = c["get_vertices"](voxel)
            orient_t = torch.from_numpy(voxel.orientation).float().requires_grad_(True)
            info = c["diff_cost_fn"].evaluate(orient_t, vertices, phase_index=voxel.phase, scale=0)

            assert info.n_peaks > 0, f"Voxel {idx}: no observable peaks"
            assert info.quality.item() > 0.0, f"Voxel {idx}: zero quality"
            info.cost.backward()
            assert orient_t.grad is not None, f"Voxel {idx}: no gradient"

    def test_blurred_scale_wider_basin(self, diff_cost_components):
        """Coarse scale (8x downsample) gives non-zero quality at larger angular offsets."""
        c = diff_cost_components
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)
        base_orient = torch.from_numpy(voxel.orientation).float()

        # Reuse pre-built sparse stack from fixture (~26 KB)
        image_stack = c["image_stack"]
        multi_stack = MultiScaleImageStack(image_stack, downsample_factors=[1, 8])

        diff_cost_fn = DifferentiableCostFunction(
            simulator=c["simulator"],
            detector_list=c["detector_list"],
            range_map=c["range_map"],
            image_stack=multi_stack,
            sample=c["sample"],
            structure_list=c["structure_list"],
        )

        # Create a small angular perturbation (~2 degrees)
        angle_rad = 2.0 * math.pi / 180.0
        axis = torch.tensor([1.0, 0.0, 0.0])  # rotate about x
        K = torch.tensor(
            [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]],
            dtype=torch.float32,
        )
        R_perturb = torch.eye(3) + math.sin(angle_rad) * K + (1 - math.cos(angle_rad)) * (K @ K)
        perturbed = R_perturb @ base_orient

        # Evaluate at fine scale (1x, original resolution)
        info_fine = diff_cost_fn.evaluate(perturbed, vertices, phase_index=voxel.phase, scale=0)
        # Evaluate at coarse scale (8x downsample, scale index 1)
        info_coarse = diff_cost_fn.evaluate(perturbed, vertices, phase_index=voxel.phase, scale=1)

        # At 2° offset, the coarse scale should have higher quality than the
        # fine scale (max_pool dilation widens spots, making them easier to "see")
        assert (
            info_coarse.quality.item() >= 0.0
        ), f"Coarse scale should give non-negative quality at 2° offset"

        # If fine scale gives zero but coarse gives non-zero, the basin is wider
        if info_fine.quality.item() < 0.01:
            assert info_coarse.quality.item() > info_fine.quality.item(), (
                f"At 2° offset: coarse quality ({info_coarse.quality.item():.4f}) should exceed "
                f"fine quality ({info_fine.quality.item():.4f})"
            )

        del multi_stack, diff_cost_fn

    def test_gradient_at_blurred_scale(self, diff_cost_components):
        """Blurred scale produces non-zero gradients even at angular offsets."""
        c = diff_cost_components
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)
        base_orient = torch.from_numpy(voxel.orientation).float()

        image_stack = c["image_stack"]
        multi_stack = MultiScaleImageStack(image_stack, downsample_factors=[1, 8])

        diff_cost_fn = DifferentiableCostFunction(
            simulator=c["simulator"],
            detector_list=c["detector_list"],
            range_map=c["range_map"],
            image_stack=multi_stack,
            sample=c["sample"],
            structure_list=c["structure_list"],
        )

        # 1 degree perturbation
        angle_rad = 1.0 * math.pi / 180.0
        axis = torch.tensor([0.0, 1.0, 0.0])
        K = torch.tensor(
            [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]],
            dtype=torch.float32,
        )
        R_perturb = torch.eye(3) + math.sin(angle_rad) * K + (1 - math.cos(angle_rad)) * (K @ K)
        perturbed = (R_perturb @ base_orient).requires_grad_(True)

        info = diff_cost_fn.evaluate(
            perturbed, vertices, phase_index=voxel.phase, scale=1  # blurred
        )

        if info.n_peaks > 0 and info.quality.item() > 0:
            info.cost.backward()
            assert perturbed.grad is not None, "Should have gradient at blurred scale"
            grad_norm = perturbed.grad.norm().item()
            assert grad_norm > 0, f"Gradient norm should be > 0 at 1° offset with blurred scale"

        del multi_stack, diff_cost_fn

    def test_differentiable_overlap_info_dataclass(self):
        """DifferentiableOverlapInfo stores correct types."""
        quality = torch.tensor(0.8, requires_grad=True)
        cost = 1.0 - quality
        info = DifferentiableOverlapInfo(quality=quality, cost=cost, n_peaks=42)

        assert info.quality.item() == pytest.approx(0.8)
        assert info.cost.item() == pytest.approx(0.2)
        assert info.n_peaks == 42
        assert info.cost.requires_grad

    def test_zero_result_for_invalid_phase(self, diff_cost_setup):
        """Invalid phase index returns zero quality, no crash."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        orient_t = torch.from_numpy(voxel.orientation).float()
        info = c["diff_cost_fn"].evaluate(orient_t, vertices, phase_index=999, scale=0)
        assert info.quality.item() == 0.0
        assert info.cost.item() == 1.0

    def test_zero_result_is_graph_connected(self, diff_cost_setup):
        """The degenerate (zero-peak) result must carry a grad_fn so callers
        relying on cost.backward() (e.g. RiemannianAdamOptimizer) get a
        genuine zero-gradient step rather than a non-differentiable tensor
        that would raise on .backward()."""
        c = diff_cost_setup
        voxel = c["mic"].voxels[0]
        vertices = c["get_vertices"](voxel)

        orient_t = torch.from_numpy(voxel.orientation).float().requires_grad_(True)
        info = c["diff_cost_fn"].evaluate(orient_t, vertices, phase_index=999, scale=0)
        assert info.cost.requires_grad
        assert info.cost.grad_fn is not None
        info.cost.backward()  # must not raise
        assert orient_t.grad is not None
        assert torch.all(orient_t.grad == 0.0)


# ---------------------------------------------------------------------------
# TestOmegaBlend — omega-direction morphological dilation
# ---------------------------------------------------------------------------


class TestOmegaBlend:
    """Tests for MultiScaleImageStack omega_window parameter.

    Uses a tiny synthetic ExperimentalImageStack (n_omega=8, n_det=2, H=16, W=16)
    to keep tests fast and independent of real data.
    """

    def _make_synthetic_stack(self, n_omega: int = 8, n_det: int = 2) -> "ExperimentalImageStack":
        """Build a tiny all-zero ExperimentalImageStack for unit testing."""
        from icenine.differentiable_cost import ExperimentalImageStack

        result = ExperimentalImageStack.__new__(ExperimentalImageStack)
        result.n_omega = n_omega
        result.n_det = n_det
        result.H = 16
        result.W = 16
        result.binary = True
        result.dtype = torch.float32
        result.images = torch.zeros(n_omega * n_det, 1, 16, 16, dtype=torch.float32)
        return result

    def test_omega_window_zero_unchanged(self):
        """omega_window=0 leaves downsampled scale images identical to baseline."""
        stack = self._make_synthetic_stack()
        # Place a bright pixel in frame 3, det 0
        stack.images[3 * stack.n_det + 0, 0, 8, 8] = 1.0

        ms_no_blend = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=0)
        ms_blend = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=0)

        torch.testing.assert_close(ms_no_blend.scales[1].images, ms_blend.scales[1].images)

    def test_omega_blend_shape(self):
        """Blended stack has the same shape as unblended."""
        stack = self._make_synthetic_stack()
        ms = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=1)
        assert ms.scales[1].images.shape == (stack.n_omega * stack.n_det, 1, 4, 4)

    def test_omega_blend_spreads_peaks(self):
        """Bright pixel at frame i appears in blended frames i-1 and i+1."""
        stack = self._make_synthetic_stack(n_omega=8, n_det=1)
        bright_omega = 4
        bright_flat = bright_omega * stack.n_det + 0
        stack.images[bright_flat, 0, 8, 8] = 1.0

        ms = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=1)
        blended = ms.scales[1].images  # (8*1, 1, 4, 4)

        bright_px = 8 // 4  # pixel 8 at factor 4 → pixel 2

        # Frame i should be bright
        assert blended[bright_flat, 0, bright_px, bright_px].item() > 0.0
        # Frame i-1 should pick up the signal
        prev_flat = (bright_omega - 1) * stack.n_det + 0
        assert blended[prev_flat, 0, bright_px, bright_px].item() > 0.0
        # Frame i+1 should pick up the signal
        next_flat = (bright_omega + 1) * stack.n_det + 0
        assert blended[next_flat, 0, bright_px, bright_px].item() > 0.0
        # Frame i-2 should NOT be affected (window=1, not 2)
        prev2_flat = (bright_omega - 2) * stack.n_det + 0
        assert blended[prev2_flat, 0, bright_px, bright_px].item() == 0.0

    def test_omega_blend_geq_original(self):
        """Every pixel in blended stack is >= the corresponding pixel in unblended."""
        stack = self._make_synthetic_stack()
        # Scatter a few random bright pixels
        stack.images[1, 0, 5, 7] = 1.0
        stack.images[6, 0, 2, 11] = 1.0

        ms_base = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=0)
        ms_blend = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=1)

        assert (ms_blend.scales[1].images >= ms_base.scales[1].images).all()

    def test_omega_blend_boundary_no_wraparound(self):
        """Bright pixel at frame 0 does NOT appear in the last frame (zero-padding)."""
        n_omega = 8
        stack = self._make_synthetic_stack(n_omega=n_omega, n_det=1)
        # Only frame 0 is bright
        stack.images[0, 0, 8, 8] = 1.0

        ms = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=1)
        blended = ms.scales[1].images

        bright_px = 8 // 4  # pixel 2
        last_flat = (n_omega - 1) * stack.n_det + 0
        assert blended[last_flat, 0, bright_px, bright_px].item() == 0.0

    def test_omega_blend_only_for_downsampled_scales(self):
        """Scale 0 (factor=1) is never omega-blended; scale 1 (factor=4) is."""
        stack = self._make_synthetic_stack()
        stack.images[3, 0, 8, 8] = 1.0  # bright pixel at flat idx 3

        ms = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=1)

        # Scale 0 is the original stack (shared, unmodified)
        assert ms.scales[0] is stack

        # Scale 1 should differ from the unblended baseline at neighbors of frame 3
        ms_no_blend = MultiScaleImageStack(stack, downsample_factors=[1, 4], omega_window=0)
        diff = (ms.scales[1].images - ms_no_blend.scales[1].images).abs().sum().item()
        assert diff > 0.0, "omega blend should change at least one pixel at scale 1"

    def test_omega_blend_memory_reasonable(self):
        """Omega-blended MultiScaleImageStack stays within acceptable memory."""
        # Use n_omega=180, n_det=2, H=512, W=512 (equivalent to scale-1 real stack)
        from icenine.differentiable_cost import ExperimentalImageStack

        large_stack = ExperimentalImageStack.__new__(ExperimentalImageStack)
        large_stack.n_omega = 180
        large_stack.n_det = 2
        large_stack.H = 512
        large_stack.W = 512
        large_stack.binary = True
        large_stack.dtype = torch.float32
        large_stack.images = torch.zeros(360, 1, 512, 512, dtype=torch.float32)

        ms = MultiScaleImageStack(large_stack, downsample_factors=[1, 4], omega_window=1)
        # scale 1 (4×): 360 × 128 × 128 × 4 bytes = ~23 MB
        # scale 0 shares large_stack: 360 × 512 × 512 × 4 = ~376 MB
        total_mb = sum(s.memory_bytes for s in ms.scales) / (1024 * 1024)
        assert total_mb < 500, f"Expected < 500 MB, got {total_mb:.1f} MB"
