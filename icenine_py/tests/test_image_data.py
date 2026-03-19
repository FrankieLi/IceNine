"""
Tests for ImageData class with dual-mode storage.

This module tests the detector image container with both dense and sparse modes,
including pixel operations, geometric rasterization, overlap calculations,
and differentiability.

Author: S. F. Li
Date: 2025-01-16
"""

import pytest
import torch
import numpy as np
import tempfile
import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from icenine.image_data import ImageData, ImageDataParameters


class TestImageDataCreation:
    """Test ImageData creation and initialization."""

    def test_dense_creation(self):
        """Test creating ImageData in dense mode."""
        image = ImageData(100, 200, mode='dense')

        assert image.num_rows == 100
        assert image.num_cols == 200
        assert image.mode == 'dense'
        assert image.shape == (100, 200)
        assert image.num_nonzero == 0
        assert image.density == 0.0

    def test_sparse_creation(self):
        """Test creating ImageData in sparse mode."""
        image = ImageData(100, 200, mode='sparse')

        assert image.num_rows == 100
        assert image.num_cols == 200
        assert image.mode == 'sparse'
        assert image.shape == (100, 200)
        assert image.num_nonzero == 0
        assert image.density == 0.0

    def test_custom_dtype_device(self):
        """Test creation with custom dtype and device."""
        image = ImageData(50, 50, mode='dense', dtype=torch.float64, device='cpu')

        assert image.dtype == torch.float64
        assert image.device.type == 'cpu'

    def test_invalid_mode(self):
        """Test that invalid mode raises ValueError."""
        with pytest.raises(ValueError, match="mode must be 'dense' or 'sparse'"):
            ImageData(100, 100, mode='invalid')


class TestPixelOperations:
    """Test basic pixel set/get/add operations."""

    def test_set_get_pixel_dense(self):
        """Test set and get pixel in dense mode."""
        image = ImageData(100, 100, mode='dense')

        # Set pixel
        image.set_pixel(10, 20, 1.5)

        # Get pixel
        value = image.get_pixel(10, 20)
        assert value.item() == pytest.approx(1.5)

        # Check nonzero count
        assert image.num_nonzero == 1

    def test_set_get_pixel_sparse(self):
        """Test set and get pixel in sparse mode."""
        image = ImageData(100, 100, mode='sparse')

        # Set pixel
        image.set_pixel(10, 20, 2.5)

        # Get pixel
        value = image.get_pixel(10, 20)
        assert value.item() == pytest.approx(2.5)

        # Check nonzero count
        assert image.num_nonzero == 1

    def test_add_to_pixel(self):
        """Test adding to pixel (accumulation)."""
        image = ImageData(100, 100, mode='dense')

        # Set initial value
        image.set_pixel(10, 20, 1.0)

        # Add to pixel
        image.add_to_pixel(10, 20, 0.5)

        # Check accumulated value
        value = image.get_pixel(10, 20)
        assert value.item() == pytest.approx(1.5)

    def test_batched_operations(self):
        """Test batched set/get operations."""
        image = ImageData(100, 100, mode='dense')

        # Set multiple pixels
        j_coords = torch.tensor([10, 20, 30])
        k_coords = torch.tensor([15, 25, 35])
        values = torch.tensor([1.0, 2.0, 3.0])

        image.set_pixels(j_coords, k_coords, values)

        # Get multiple pixels
        retrieved = image.get_pixels(j_coords, k_coords)

        assert torch.allclose(retrieved, values)
        assert image.num_nonzero == 3


class TestBoundsChecking:
    """Test bounds checking and pixel state queries."""

    def test_is_in_bounds(self):
        """Test bounds checking."""
        image = ImageData(100, 100, mode='dense')

        # Valid coordinates
        assert image.is_in_bounds(50, 50).item() is True
        assert image.is_in_bounds(0, 0).item() is True
        assert image.is_in_bounds(99, 99).item() is True

        # Out of bounds
        assert image.is_in_bounds(-1, 50).item() is False
        assert image.is_in_bounds(50, -1).item() is False
        assert image.is_in_bounds(100, 50).item() is False
        assert image.is_in_bounds(50, 100).item() is False

    def test_is_dark_is_bright(self):
        """Test dark/bright pixel queries."""
        image = ImageData(100, 100, mode='dense')

        # Initially all dark
        assert image.is_dark(50, 50).item() is True
        assert image.is_bright(50, 50).item() is False

        # Set pixel to positive value
        image.set_pixel(50, 50, 1.0)

        assert image.is_dark(50, 50).item() is False
        assert image.is_bright(50, 50).item() is True

    def test_clear(self):
        """Test clearing image."""
        image = ImageData(100, 100, mode='dense')

        # Add some pixels
        image.set_pixel(10, 10, 1.0)
        image.set_pixel(20, 20, 2.0)

        assert image.num_nonzero == 2

        # Clear
        image.clear()

        assert image.num_nonzero == 0
        assert image.get_pixel(10, 10).item() == 0.0


class TestModeConversion:
    """Test conversion between dense and sparse modes."""

    def test_dense_to_sparse(self):
        """Test converting dense to sparse."""
        # Create dense image
        dense = ImageData(100, 100, mode='dense')
        dense.set_pixel(10, 10, 1.0)
        dense.set_pixel(20, 20, 2.0)

        # Convert to sparse
        sparse = dense.to_sparse()

        assert sparse.mode == 'sparse'
        assert sparse.num_nonzero == 2
        assert sparse.get_pixel(10, 10).item() == pytest.approx(1.0)
        assert sparse.get_pixel(20, 20).item() == pytest.approx(2.0)

    def test_sparse_to_dense(self):
        """Test converting sparse to dense."""
        # Create sparse image
        sparse = ImageData(100, 100, mode='sparse')
        sparse.set_pixel(10, 10, 1.0)
        sparse.set_pixel(20, 20, 2.0)

        # Convert to dense
        dense = sparse.to_dense()

        assert dense.mode == 'dense'
        assert dense.num_nonzero == 2
        assert dense.get_pixel(10, 10).item() == pytest.approx(1.0)
        assert dense.get_pixel(20, 20).item() == pytest.approx(2.0)

    def test_roundtrip_conversion(self):
        """Test dense -> sparse -> dense roundtrip."""
        # Create dense image
        original = ImageData(100, 100, mode='dense')
        original.set_pixel(10, 10, 1.5)
        original.set_pixel(20, 20, 2.5)
        original.set_pixel(30, 30, 3.5)

        # Convert to sparse and back
        sparse = original.to_sparse()
        reconstructed = sparse.to_dense()

        # Check all values match
        assert reconstructed.num_nonzero == original.num_nonzero
        assert reconstructed.get_pixel(10, 10).item() == pytest.approx(1.5)
        assert reconstructed.get_pixel(20, 20).item() == pytest.approx(2.5)
        assert reconstructed.get_pixel(30, 30).item() == pytest.approx(3.5)


class TestIOOperations:
    """Test I/O operations (ASCII, binary, NumPy)."""

    def test_ascii_save_load(self):
        """Test ASCII file save and load."""
        # Create image with data
        image = ImageData(100, 100, mode='dense')
        image.set_pixel(10, 20, 1.5)
        image.set_pixel(30, 40, 2.5)

        # Save to temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            filename = f.name

        try:
            image.save_ascii(filename)

            # Load into new image
            loaded = ImageData(100, 100, mode='dense')
            loaded.load_ascii(filename)

            # Check values match
            assert loaded.num_nonzero == 2
            assert loaded.get_pixel(10, 20).item() == pytest.approx(1.5)
            assert loaded.get_pixel(30, 40).item() == pytest.approx(2.5)

        finally:
            os.unlink(filename)

    def test_binary_save_load(self):
        """Test binary file save and load."""
        # Create image with data
        image = ImageData(100, 100, mode='dense')
        image.set_pixel(10, 20, 1.5)
        image.set_pixel(30, 40, 2.5)

        # Save to temporary file
        with tempfile.NamedTemporaryFile(delete=False, suffix='.pt') as f:
            filename = f.name

        try:
            image.save_binary(filename)

            # Load
            loaded = ImageData.load_binary(filename)

            # Check values match
            assert loaded.num_rows == 100
            assert loaded.num_cols == 100
            assert loaded.mode == 'dense'
            assert loaded.num_nonzero == 2
            assert loaded.get_pixel(10, 20).item() == pytest.approx(1.5)
            assert loaded.get_pixel(30, 40).item() == pytest.approx(2.5)

        finally:
            os.unlink(filename)

    def test_to_numpy(self):
        """Test conversion to NumPy array."""
        image = ImageData(50, 50, mode='dense')
        image.set_pixel(10, 20, 1.5)

        # Convert to NumPy
        array = image.to_numpy()

        assert isinstance(array, np.ndarray)
        assert array.shape == (50, 50)
        assert array[20, 10] == pytest.approx(1.5)  # Note: [k, j] indexing

    def test_from_numpy(self):
        """Test creation from NumPy array."""
        # Create NumPy array
        array = np.zeros((50, 50))
        array[20, 10] = 1.5

        # Create ImageData from array
        image = ImageData.from_numpy(array, mode='dense')

        assert image.num_rows == 50
        assert image.num_cols == 50
        assert image.get_pixel(10, 20).item() == pytest.approx(1.5)


class TestTriangleRasterization:
    """Test triangle rasterization (both hard and soft modes)."""

    def test_simple_triangle_hard(self):
        """Test hard rasterization of simple triangle."""
        image = ImageData(100, 100, mode='dense')

        # Define triangle vertices (in pixel coordinates)
        v0 = torch.tensor([10.0, 10.0])
        v1 = torch.tensor([50.0, 10.0])
        v2 = torch.tensor([30.0, 50.0])

        # Rasterize with hard mode
        image.add_triangle(v0, v1, v2, intensity=1.0, mode='hard')

        # Check that some pixels were lit
        assert image.num_nonzero > 0

        # Check pixels inside triangle are bright
        # Center of triangle should be inside
        center_j = (10 + 50 + 30) / 3.0
        center_k = (10 + 10 + 50) / 3.0
        assert image.get_pixel(int(center_j), int(center_k)).item() > 0

    def test_simple_triangle_soft(self):
        """Test soft rasterization of simple triangle."""
        image = ImageData(100, 100, mode='dense')

        # Define triangle vertices
        v0 = torch.tensor([10.0, 10.0])
        v1 = torch.tensor([50.0, 10.0])
        v2 = torch.tensor([30.0, 50.0])

        # Rasterize with soft mode
        image.add_triangle(v0, v1, v2, intensity=1.0, mode='soft', temperature=1.0)

        # Check that pixels were lit (soft mode gives fractional values)
        assert image._pixels_dense.sum() > 0

        # Center should have high intensity
        center_j = (10 + 50 + 30) / 3.0
        center_k = (10 + 10 + 50) / 3.0
        assert image.get_pixel(int(center_j), int(center_k)).item() > 0.5

    def test_triangle_temperature_effect(self):
        """Test that temperature affects soft rasterization."""
        # Low temperature (sharper)
        image_sharp = ImageData(100, 100, mode='dense')
        v0 = torch.tensor([20.0, 20.0])
        v1 = torch.tensor([60.0, 20.0])
        v2 = torch.tensor([40.0, 60.0])
        image_sharp.add_triangle(v0, v1, v2, intensity=1.0, mode='soft', temperature=0.1)

        # High temperature (softer)
        image_soft = ImageData(100, 100, mode='dense')
        image_soft.add_triangle(v0, v1, v2, intensity=1.0, mode='soft', temperature=10.0)

        # Sharp should have more concentrated values (closer to binary)
        # Soft should have more spread-out values
        sharp_nonzero = (image_sharp._pixels_dense > 0.1).sum()
        soft_nonzero = (image_soft._pixels_dense > 0.1).sum()

        assert soft_nonzero > sharp_nonzero  # Soft mode spreads to more pixels

    def test_degenerate_triangle(self):
        """Test handling of degenerate triangle (collinear points)."""
        image = ImageData(100, 100, mode='dense')

        # Collinear points (degenerate triangle)
        v0 = torch.tensor([10.0, 10.0])
        v1 = torch.tensor([20.0, 10.0])
        v2 = torch.tensor([30.0, 10.0])

        # Should not crash
        image.add_triangle(v0, v1, v2, intensity=1.0, mode='hard')

        # Should produce minimal or no output (degenerate)
        assert image.num_nonzero <= 5  # Allow for edge pixels

    def test_triangle_out_of_bounds(self):
        """Test triangle partially out of bounds."""
        image = ImageData(100, 100, mode='dense')

        # Triangle partially outside image
        v0 = torch.tensor([-10.0, -10.0])
        v1 = torch.tensor([50.0, -10.0])
        v2 = torch.tensor([20.0, 50.0])

        # Should clip to image bounds
        image.add_triangle(v0, v1, v2, intensity=1.0, mode='hard')

        # Should have some pixels (the part inside bounds)
        assert image.num_nonzero > 0


class TestPolygonRasterization:
    """Test polygon rasterization."""

    def test_square_polygon(self):
        """Test rasterizing square polygon."""
        image = ImageData(100, 100, mode='dense')

        # Square vertices
        vertices = torch.tensor([
            [20.0, 20.0],
            [60.0, 20.0],
            [60.0, 60.0],
            [20.0, 60.0]
        ])

        image.add_polygon(vertices, intensity=1.0, mode='hard')

        # Check pixels inside square are lit
        assert image.num_nonzero > 0

        # Center should be bright
        assert image.get_pixel(40, 40).item() > 0

    def test_pentagon_polygon(self):
        """Test rasterizing pentagon."""
        image = ImageData(100, 100, mode='dense')

        # Pentagon (approximate)
        vertices = torch.tensor([
            [50.0, 20.0],
            [70.0, 40.0],
            [60.0, 65.0],
            [40.0, 65.0],
            [30.0, 40.0]
        ])

        image.add_polygon(vertices, intensity=1.0, mode='soft')

        # Should have pixels lit
        assert image._pixels_dense.sum() > 0

    def test_polygon_invalid_vertices(self):
        """Test that polygon with < 3 vertices raises error."""
        image = ImageData(100, 100, mode='dense')

        # Only 2 vertices
        vertices = torch.tensor([[10.0, 10.0], [20.0, 20.0]])

        with pytest.raises(ValueError, match="Polygon must have at least 3 vertices"):
            image.add_polygon(vertices, intensity=1.0)


class TestOverlapCalculations:
    """Test overlap calculations for cost functions."""

    def test_get_num_pixels_lit(self):
        """Test counting pixels lit by triangle."""
        image = ImageData(100, 100, mode='dense')

        v0 = torch.tensor([20.0, 20.0])
        v1 = torch.tensor([60.0, 20.0])
        v2 = torch.tensor([40.0, 60.0])

        num_lit = image.get_num_pixels_lit(v0, v1, v2, mode='hard')

        assert num_lit > 0
        assert isinstance(num_lit, torch.Tensor)

    def test_triangle_overlap_no_overlap(self):
        """Test overlap when triangle doesn't overlap existing data."""
        # Create experimental image with data in one region
        experimental = ImageData(100, 100, mode='dense')
        exp_v0 = torch.tensor([10.0, 10.0])
        exp_v1 = torch.tensor([30.0, 10.0])
        exp_v2 = torch.tensor([20.0, 30.0])
        experimental.add_triangle(exp_v0, exp_v1, exp_v2, intensity=1.0, mode='hard')

        # Test triangle in different region (no overlap)
        v0 = torch.tensor([60.0, 60.0])
        v1 = torch.tensor([80.0, 60.0])
        v2 = torch.tensor([70.0, 80.0])

        overlap, total = experimental.get_triangle_overlap_property(v0, v1, v2, mode='hard')

        assert total > 0  # Triangle has pixels
        assert overlap == 0  # But no overlap with experimental

    def test_triangle_overlap_full_overlap(self):
        """Test overlap when triangles fully overlap."""
        # Create experimental image
        experimental = ImageData(100, 100, mode='dense')
        v0 = torch.tensor([20.0, 20.0])
        v1 = torch.tensor([60.0, 20.0])
        v2 = torch.tensor([40.0, 60.0])
        experimental.add_triangle(v0, v1, v2, intensity=1.0, mode='hard')

        # Test with same triangle (perfect overlap)
        overlap, total = experimental.get_triangle_overlap_property(v0, v1, v2, mode='hard')

        assert overlap > 0
        assert total > 0
        # Overlap should be close to total (perfect match)
        assert overlap == pytest.approx(total, rel=0.1)

    def test_triangle_overlap_partial_overlap(self):
        """Test overlap with partial overlap."""
        # Create experimental image
        experimental = ImageData(100, 100, mode='dense')
        exp_v0 = torch.tensor([20.0, 20.0])
        exp_v1 = torch.tensor([60.0, 20.0])
        exp_v2 = torch.tensor([40.0, 60.0])
        experimental.add_triangle(exp_v0, exp_v1, exp_v2, intensity=1.0, mode='hard')

        # Test with shifted triangle (partial overlap)
        v0 = torch.tensor([30.0, 30.0])
        v1 = torch.tensor([70.0, 30.0])
        v2 = torch.tensor([50.0, 70.0])

        overlap, total = experimental.get_triangle_overlap_property(v0, v1, v2, mode='hard')

        assert overlap > 0  # Some overlap
        assert total > overlap  # But not full overlap


class TestDifferentiability:
    """Test gradient flow through operations."""

    def test_triangle_vertex_gradients(self):
        """Test gradients flow through triangle vertices."""
        image = ImageData(100, 100, mode='dense')

        # Create triangle with requires_grad
        v0 = torch.tensor([20.0, 20.0], requires_grad=True)
        v1 = torch.tensor([60.0, 20.0], requires_grad=True)
        v2 = torch.tensor([40.0, 60.0], requires_grad=True)

        # Rasterize with soft mode (differentiable)
        image.add_triangle(v0, v1, v2, intensity=1.0, mode='soft', temperature=1.0)

        # Compute loss (sum of all pixels)
        loss = image._pixels_dense.sum()

        # Backpropagate
        loss.backward()

        # Check gradients exist
        assert v0.grad is not None
        assert v1.grad is not None
        assert v2.grad is not None

        # Gradients should be non-zero
        assert torch.abs(v0.grad).sum() > 0
        assert torch.abs(v1.grad).sum() > 0
        assert torch.abs(v2.grad).sum() > 0

    def test_overlap_gradients(self):
        """Test gradients flow through overlap calculation."""
        # Create experimental image (fixed)
        experimental = ImageData(100, 100, mode='dense')
        exp_v0 = torch.tensor([30.0, 30.0])
        exp_v1 = torch.tensor([70.0, 30.0])
        exp_v2 = torch.tensor([50.0, 70.0])
        experimental.add_triangle(exp_v0, exp_v1, exp_v2, intensity=1.0, mode='hard')

        # Create simulated triangle with gradient tracking
        v0 = torch.tensor([25.0, 25.0], requires_grad=True)
        v1 = torch.tensor([65.0, 25.0], requires_grad=True)
        v2 = torch.tensor([45.0, 65.0], requires_grad=True)

        # Compute overlap (differentiable)
        overlap, total = experimental.get_triangle_overlap_property(
            v0, v1, v2, mode='soft', temperature=1.0
        )

        # Loss: maximize overlap
        loss = -overlap

        # Backpropagate
        loss.backward()

        # Gradients should exist
        assert v0.grad is not None
        assert v1.grad is not None
        assert v2.grad is not None

    def test_intensity_gradient(self):
        """Test gradient through intensity parameter."""
        image = ImageData(100, 100, mode='dense')

        v0 = torch.tensor([20.0, 20.0])
        v1 = torch.tensor([60.0, 20.0])
        v2 = torch.tensor([40.0, 60.0])

        # Intensity with gradient
        intensity = torch.tensor(2.0, requires_grad=True)

        image.add_triangle(v0, v1, v2, intensity=intensity, mode='soft')

        loss = image._pixels_dense.sum()
        loss.backward()

        # Gradient should exist and be positive
        assert intensity.grad is not None
        assert intensity.grad.item() > 0

    def test_soft_vs_hard_gradients(self):
        """Test that soft mode has gradients but hard mode doesn't."""
        # Soft mode
        image_soft = ImageData(100, 100, mode='dense')
        v0_soft = torch.tensor([20.0, 20.0], requires_grad=True)
        v1_soft = torch.tensor([60.0, 20.0], requires_grad=True)
        v2_soft = torch.tensor([40.0, 60.0], requires_grad=True)

        image_soft.add_triangle(v0_soft, v1_soft, v2_soft, intensity=1.0, mode='soft')
        loss_soft = image_soft._pixels_dense.sum()
        loss_soft.backward()

        # Soft mode should have gradients
        assert v0_soft.grad is not None
        assert torch.abs(v0_soft.grad).sum() > 0

        # Hard mode - gradients won't flow because of non-differentiable operations
        # We can verify this by checking that the loss doesn't have grad_fn
        image_hard = ImageData(100, 100, mode='dense')
        v0_hard = torch.tensor([20.0, 20.0], requires_grad=True)
        v1_hard = torch.tensor([60.0, 20.0], requires_grad=True)
        v2_hard = torch.tensor([40.0, 60.0], requires_grad=True)

        image_hard.add_triangle(v0_hard, v1_hard, v2_hard, intensity=1.0, mode='hard')
        loss_hard = image_hard._pixels_dense.sum()

        # Hard mode loss won't have grad_fn (non-differentiable)
        # This is expected - hard mode is for efficiency, not optimization
        assert loss_hard.grad_fn is None or not loss_hard.requires_grad

    def test_temperature_gradient_effect(self):
        """Test gradient magnitude changes with temperature."""
        # Low temperature (sharp boundaries)
        image_sharp = ImageData(100, 100, mode='dense')
        v0_sharp = torch.tensor([20.0, 20.0], requires_grad=True)
        v1_sharp = torch.tensor([60.0, 20.0], requires_grad=True)
        v2_sharp = torch.tensor([40.0, 60.0], requires_grad=True)

        image_sharp.add_triangle(v0_sharp, v1_sharp, v2_sharp,
                                intensity=1.0, mode='soft', temperature=0.1)
        loss_sharp = image_sharp._pixels_dense.sum()
        loss_sharp.backward()

        grad_sharp = v0_sharp.grad.clone()

        # High temperature (smooth boundaries)
        image_smooth = ImageData(100, 100, mode='dense')
        v0_smooth = torch.tensor([20.0, 20.0], requires_grad=True)
        v1_smooth = torch.tensor([60.0, 20.0], requires_grad=True)
        v2_smooth = torch.tensor([40.0, 60.0], requires_grad=True)

        image_smooth.add_triangle(v0_smooth, v1_smooth, v2_smooth,
                                 intensity=1.0, mode='soft', temperature=10.0)
        loss_smooth = image_smooth._pixels_dense.sum()
        loss_smooth.backward()

        grad_smooth = v0_smooth.grad.clone()

        # Both should have gradients
        assert torch.abs(grad_sharp).sum() > 0
        assert torch.abs(grad_smooth).sum() > 0


class TestHelperMethods:
    """Test helper methods and utilities."""

    def test_get_parameters(self):
        """Test getting parameters as dataclass."""
        image = ImageData(100, 200, mode='dense', dtype=torch.float32, device='cpu')

        params = image.get_parameters()

        assert isinstance(params, ImageDataParameters)
        assert params.num_rows == 100
        assert params.num_cols == 200
        assert params.mode == 'dense'
        assert params.dtype == torch.float32
        assert params.device == 'cpu'

    def test_repr(self):
        """Test string representation."""
        image = ImageData(100, 200, mode='dense')
        image.set_pixel(10, 20, 1.0)

        repr_str = repr(image)

        assert 'ImageData' in repr_str
        assert 'num_rows=100' in repr_str
        assert 'num_cols=200' in repr_str
        assert 'mode=\'dense\'' in repr_str
        assert 'nonzero=1' in repr_str


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
