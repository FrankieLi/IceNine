"""
Tests for detector geometry and coordinate transformations.

These tests verify the Detector class implementation against expected behavior
and ensure consistency with the C++ implementation.
"""

import pytest
import torch
import numpy as np

from icenine.detector import Detector, DetectorParameters
from icenine.geometry import Ray, Plane


class TestDetectorCreation:
    """Test detector creation and initialization."""

    def test_basic_creation(self):
        """Test creating a detector with default parameters."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,  # pixels
            beam_center_k=1024.0,  # pixels
            dtype=torch.float32,
            device='cpu'
        )

        assert detector.num_rows == 2048
        assert detector.num_cols == 2048
        assert detector.pixel_height == 0.2
        assert detector.pixel_width == 0.2
        assert detector.beam_center_j == 1024.0
        assert detector.beam_center_k == 1024.0

        # Check default position (origin)
        assert torch.allclose(detector.position, torch.zeros(3))

        # Check default orientation (identity)
        assert torch.allclose(detector.orientation, torch.eye(3))

    def test_custom_position_orientation(self):
        """Test creating detector with custom position and orientation."""
        position = torch.tensor([100.0, 0.0, 0.0])
        orientation = torch.eye(3)  # Identity for now

        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            pixel_height=0.1,
            pixel_width=0.1,
            beam_center_j=512.0,  # pixels
            beam_center_k=512.0,  # pixels
            position=position,
            orientation=orientation
        )

        assert torch.allclose(detector.position, position)
        assert torch.allclose(detector.orientation, orientation)

    def test_from_beam_center_factory(self):
        """Test factory method for creating detector."""
        position = torch.tensor([100.0, 0.0, 0.0])

        detector = Detector.from_beam_center(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=position
        )

        assert detector.num_rows == 2048
        assert torch.allclose(detector.position, position)

    def test_detector_properties(self):
        """Test detector property accessors."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0
        )

        # Check computed properties
        assert detector.detector_width == 2048 * 0.2
        assert detector.detector_height == 2048 * 0.2

        # Check basis vectors (should be Y and Z for unrotated detector)
        j_basis, k_basis = detector.basis_vectors
        assert torch.allclose(j_basis, torch.tensor([0.0, 1.0, 0.0]))
        assert torch.allclose(k_basis, torch.tensor([0.0, 0.0, 1.0]))


class TestCoordinateTransformations:
    """Test coordinate transformations between different frames."""

    def test_lab_to_detector_coordinate_origin(self):
        """Test lab to detector coordinate transformation at detector center."""
        position = torch.tensor([100.0, 0.0, 0.0])
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=position
        )

        # Point at detector center (rotation center) should map to beam center in mm
        j, k = detector.lab_to_detector_coordinate(position)

        # At detector center, we should be at beam center coordinates in mm
        # beam_center is 1024 pixels * 0.2 mm/pixel = 204.8 mm
        assert j.item() == pytest.approx(204.8, abs=1e-4)
        assert k.item() == pytest.approx(204.8, abs=1e-4)

    def test_detector_to_lab_coordinate_roundtrip(self):
        """Test roundtrip: lab -> detector -> lab."""
        position = torch.tensor([100.0, 0.0, 0.0])
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=position
        )

        # Test point on detector plane
        lab_point = torch.tensor([100.0, 10.0, 5.0])

        # Lab -> detector -> lab
        j, k = detector.lab_to_detector_coordinate(lab_point)
        lab_point_reconstructed = detector.detector_to_lab_coordinate(j, k)

        assert torch.allclose(lab_point, lab_point_reconstructed, atol=1e-5)

    def test_pixel_to_lab_coordinate_roundtrip(self):
        """Test roundtrip: pixel -> lab -> pixel."""
        position = torch.tensor([100.0, 0.0, 0.0])
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=position
        )

        # Test pixel coordinates
        col_original = torch.tensor(1024.0)
        row_original = torch.tensor(1024.0)

        # Pixel -> lab -> pixel
        lab_point = detector.pixel_to_lab_coordinate(col_original, row_original)
        row_reconstructed, col_reconstructed = detector.lab_to_pixel(lab_point)

        # Note: There's a 0.5 pixel offset due to the half-pixel adjustment in lab_to_pixel
        # pixel_to_lab returns the pixel edge, while lab_to_pixel adds 0.5 to find the pixel center
        assert torch.isclose(col_original, col_reconstructed, atol=0.6)
        assert torch.isclose(row_original, row_reconstructed, atol=0.6)

    def test_beam_center_pixel_coordinates(self):
        """Test that beam center maps to expected pixel coordinates."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,  # 1024 pixels * 0.2 mm/pixel
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # The beam center (at detector position) should map to beam center pixels
        beam_center_lab = detector.position
        row, col = detector.lab_to_pixel(beam_center_lab)

        # Beam center should be at (1024, 1024) pixels
        assert row.item() == pytest.approx(1024.0, abs=0.5)
        assert col.item() == pytest.approx(1024.0, abs=0.5)

    def test_batched_coordinate_transformations(self):
        """Test batched coordinate transformations."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Create batch of lab points
        batch_size = 10
        lab_points = torch.randn(batch_size, 3) * 10.0
        lab_points[:, 0] = 100.0  # All points on detector plane

        # Transform to detector coordinates
        j, k = detector.lab_to_detector_coordinate(lab_points)

        assert j.shape == (batch_size,)
        assert k.shape == (batch_size,)

        # Roundtrip test
        lab_points_reconstructed = detector.detector_to_lab_coordinate(j, k)
        assert torch.allclose(lab_points, lab_points_reconstructed, atol=1e-4)


class TestDetectorTransformations:
    """Test detector position and orientation transformations."""

    def test_set_position(self):
        """Test setting detector position."""
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=512.0,
            beam_center_k=512.0
        )

        # Set new position
        new_position = torch.tensor([150.0, 10.0, -5.0])
        detector.set_position(new_position)

        assert torch.allclose(detector.position, new_position)

    def test_translate(self):
        """Test translating detector."""
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=512.0,
            beam_center_k=512.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Translate detector
        translation = torch.tensor([10.0, 5.0, -3.0])
        detector.translate(translation)

        expected_position = torch.tensor([110.0, 5.0, -3.0])
        assert torch.allclose(detector.position, expected_position)

    def test_set_orientation_euler(self):
        """Test setting orientation using Euler angles."""
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=512.0,
            beam_center_k=512.0
        )

        # Set orientation (small rotation)
        detector.set_orientation_euler(phi=10.0, theta=5.0, psi=0.0)

        # Orientation should not be identity
        assert not torch.allclose(detector.orientation, torch.eye(3))

        # Orientation should be a valid rotation matrix (orthogonal, det=1)
        orientation = detector.orientation
        assert torch.allclose(
            orientation @ orientation.T,
            torch.eye(3),
            atol=1e-5
        )
        assert torch.isclose(torch.det(orientation), torch.tensor(1.0), atol=1e-5)

    def test_rotate(self):
        """Test rotating detector by additional Euler angles."""
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=512.0,
            beam_center_k=512.0
        )

        # Initial orientation
        initial_orientation = detector.orientation.clone()

        # Rotate detector
        detector.rotate(phi=5.0, theta=2.0, psi=0.0)

        # Orientation should have changed
        assert not torch.allclose(detector.orientation, initial_orientation)

        # Still a valid rotation matrix
        orientation = detector.orientation
        assert torch.allclose(
            orientation @ orientation.T,
            torch.eye(3),
            atol=1e-5
        )

    def test_orientation_affects_basis_vectors(self):
        """Test that orientation changes affect basis vectors."""
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=512.0,
            beam_center_k=512.0
        )

        # Initial basis vectors (Y and Z for unrotated)
        j_basis_initial, k_basis_initial = detector.basis_vectors

        # Rotate detector
        detector.set_orientation_euler(phi=45.0, theta=0.0, psi=0.0)

        # Basis vectors should have changed
        j_basis_rotated, k_basis_rotated = detector.basis_vectors

        assert not torch.allclose(j_basis_initial, j_basis_rotated)
        # Note: k_basis might not change for rotation around Z


class TestRayIntersection:
    """Test ray-detector intersection calculations."""

    def test_perpendicular_ray_intersection(self):
        """Test intersection with ray perpendicular to detector."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Ray traveling along +X towards detector
        ray = Ray(
            origin=torch.tensor([0.0, 0.0, 0.0]),
            direction=torch.tensor([1.0, 0.0, 0.0])
        )

        intersects, t = detector.intersect_ray(ray)

        assert intersects.item() is True
        assert t.item() == pytest.approx(100.0, abs=1e-5)

        # Verify intersection point
        intersection_point = ray.at(t)
        assert torch.allclose(intersection_point, detector.position, atol=1e-5)

    def test_angled_ray_intersection(self):
        """Test intersection with angled ray."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Ray at angle
        direction = torch.tensor([1.0, 0.5, 0.3])
        direction = direction / torch.norm(direction)

        ray = Ray(
            origin=torch.tensor([0.0, 0.0, 0.0]),
            direction=direction
        )

        intersects, t = detector.intersect_ray(ray)

        assert intersects.item() is True
        assert t.item() > 0

        # Intersection point should be on detector plane
        intersection_point = ray.at(t)
        assert intersection_point[0].item() == pytest.approx(100.0, abs=1e-4)

    def test_parallel_ray_no_intersection(self):
        """Test that parallel ray doesn't intersect."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Ray parallel to detector plane (along Y)
        ray = Ray(
            origin=torch.tensor([0.0, 0.0, 0.0]),
            direction=torch.tensor([0.0, 1.0, 0.0])
        )

        intersects, t = detector.intersect_ray(ray)

        assert intersects.item() is False


class TestInRange:
    """Test pixel bounds checking."""

    def test_in_range_center(self):
        """Test that center pixel is in range."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0
        )

        col = torch.tensor(1024.0)
        row = torch.tensor(1024.0)

        assert detector.in_range(col, row).item() is True

    def test_in_range_corners(self):
        """Test corner pixels."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0
        )

        # Top-left corner
        assert detector.in_range(torch.tensor(0.0), torch.tensor(0.0)).item() is True

        # Bottom-right corner (last valid pixel is 2047, 2047)
        assert detector.in_range(torch.tensor(2047.0), torch.tensor(2047.0)).item() is True

    def test_out_of_range(self):
        """Test pixels outside detector bounds."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0
        )

        # Negative pixels
        assert detector.in_range(torch.tensor(-1.0), torch.tensor(100.0)).item() is False
        assert detector.in_range(torch.tensor(100.0), torch.tensor(-1.0)).item() is False

        # Beyond bounds
        assert detector.in_range(torch.tensor(2048.0), torch.tensor(100.0)).item() is False
        assert detector.in_range(torch.tensor(100.0), torch.tensor(2048.0)).item() is False


class TestDifferentiability:
    """Test that operations are differentiable."""

    def test_lab_to_pixel_gradient(self):
        """Test gradient flow through lab_to_pixel."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Lab point with gradient tracking
        lab_point = torch.tensor([100.0, 10.0, 5.0], requires_grad=True)

        # Transform to pixel coordinates
        row, col = detector.lab_to_pixel(lab_point)

        # Compute loss (e.g., distance from center)
        loss = (row - 1024.0)**2 + (col - 1024.0)**2

        # Backpropagate
        loss.backward()

        # Gradient should exist
        assert lab_point.grad is not None
        assert not torch.all(lab_point.grad == 0)

    def test_orientation_euler_gradient(self):
        """Test gradient flow through Euler angle setting."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Note: set_orientation_euler creates new tensors internally,
        # so we test the underlying euler_to_matrix_torch directly
        from icenine.geometry import euler_to_matrix_torch

        phi = torch.tensor(10.0, requires_grad=True)
        theta = torch.tensor(5.0, requires_grad=True)
        psi = torch.tensor(0.0, requires_grad=True)

        # Build orientation matrix
        orientation = euler_to_matrix_torch(phi, theta, psi)

        # Compute loss (e.g., difference from identity)
        loss = torch.sum((orientation - torch.eye(3))**2)

        # Backpropagate
        loss.backward()

        # Gradients should exist
        assert phi.grad is not None
        assert theta.grad is not None
        assert psi.grad is not None

    def test_ray_intersection_gradient(self):
        """Test gradient flow through ray intersection."""
        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=torch.tensor([100.0, 0.0, 0.0])
        )

        # Ray with gradient tracking
        origin = torch.tensor([0.0, 0.0, 0.0], requires_grad=True)
        direction = torch.tensor([1.0, 0.1, 0.05])
        direction = direction / torch.norm(direction)

        ray = Ray(origin=origin, direction=direction)

        # Intersect ray with detector
        intersects, t = detector.intersect_ray(ray)

        # Backpropagate through t
        t.backward()

        # Gradient should exist for origin
        assert origin.grad is not None


class TestDetectorParameters:
    """Test detector parameter serialization."""

    def test_get_parameters(self):
        """Test getting detector parameters."""
        position = torch.tensor([100.0, 0.0, 0.0])
        orientation = torch.eye(3)

        detector = Detector(
            num_rows=2048,
            num_cols=2048,
            pixel_height=0.2,
            pixel_width=0.2,
            beam_center_j=1024.0,
            beam_center_k=1024.0,
            position=position,
            orientation=orientation
        )

        params = detector.get_parameters()

        assert isinstance(params, DetectorParameters)
        assert params.num_rows == 2048
        assert params.num_cols == 2048
        assert params.pixel_height == 0.2
        assert params.pixel_width == 0.2
        assert params.beam_center_j == 1024.0  # pixels
        assert params.beam_center_k == 1024.0  # pixels
        assert torch.allclose(params.position, position)
        assert torch.allclose(params.orientation, orientation)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
