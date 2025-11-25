"""
Unit tests for ray projection and reflection functions in diffraction_core.

Tests the geometric optics functions used for projecting diffraction peaks
onto detectors.
"""

import pytest
import torch
import numpy as np

from icenine.diffraction_core import (
    get_reflection_vector,
    get_reflected_ray_dir,
    build_reflected_ray,
    get_reflected_ray,
    get_illuminated_pixel
)
from icenine.sample import Sample
from icenine.detector import Detector
from icenine.geometry import Ray


class TestReflectionVector:
    """Test get_reflection_vector function."""

    def test_reflection_downward_ray(self):
        """Test reflection of downward ray off horizontal surface."""
        r_in = torch.tensor([0., 0., -1.])
        normal = torch.tensor([0., 0., 1.])

        r_out = get_reflection_vector(r_in, normal)

        expected = torch.tensor([0., 0., 1.])
        assert torch.allclose(r_out, expected, atol=1e-6)

    def test_reflection_diagonal_ray(self):
        """Test reflection at 45 degrees."""
        # Ray coming down and right at 45°
        r_in = torch.tensor([1., 0., -1.]) / np.sqrt(2)
        normal = torch.tensor([0., 0., 1.])

        r_out = get_reflection_vector(r_in, normal)

        # Should reflect up and right at 45°
        expected = torch.tensor([1., 0., 1.]) / np.sqrt(2)
        assert torch.allclose(r_out, expected, atol=1e-6)

    def test_reflection_perpendicular(self):
        """Test perpendicular incidence (straight bounce back)."""
        r_in = torch.tensor([0., 0., -1.])
        normal = torch.tensor([0., 0., 1.])

        r_out = get_reflection_vector(r_in, normal)

        expected = torch.tensor([0., 0., 1.])
        assert torch.allclose(r_out, expected, atol=1e-6)

    def test_reflection_preserves_magnitude(self):
        """Test that reflection preserves vector magnitude."""
        r_in = torch.tensor([1., 2., -3.])
        normal = torch.tensor([0., 0., 1.])
        normal = normal / torch.norm(normal)

        r_out = get_reflection_vector(r_in, normal)

        assert torch.isclose(torch.norm(r_out), torch.norm(r_in), atol=1e-6)

    def test_reflection_batched(self):
        """Test batched reflection of multiple rays."""
        r_in = torch.tensor([[0., 0., -1.],
                            [1., 0., -1.],
                            [0., 1., -1.]])
        normal = torch.tensor([0., 0., 1.])

        r_out = get_reflection_vector(r_in, normal)

        expected = torch.tensor([[0., 0., 1.],
                                [1., 0., 1.],
                                [0., 1., 1.]])
        assert torch.allclose(r_out, expected, atol=1e-6)


class TestReflectedRayDir:
    """Test get_reflected_ray_dir function."""

    def test_with_identity_orientation(self):
        """Test with sample at identity orientation."""
        sample = Sample()
        sample.set_location(np.array([0., 0., 0.]))
        sample.set_orientation_matrix(np.eye(3))

        normal = torch.tensor([0., 0., 1.])  # Sample frame
        beam_dir = torch.tensor([0., 0., 1.])  # Lab frame

        reflected_dir = get_reflected_ray_dir(sample, normal, beam_dir)

        # Normal reflection: beam along +Z reflects along -Z
        expected = torch.tensor([0., 0., -1.])
        assert torch.allclose(reflected_dir, expected, atol=1e-5)

    def test_with_rotated_sample(self):
        """Test with rotated sample orientation."""
        sample = Sample()
        sample.set_location(np.array([0., 0., 0.]))

        # Rotate 90° around Y-axis: Z -> -X
        orientation = np.array([[0., 0., -1.],
                               [0., 1., 0.],
                               [1., 0., 0.]], dtype=np.float32)
        sample.set_orientation_matrix(orientation)

        normal = torch.tensor([0., 0., 1.])  # Sample frame (will become -X in lab)
        beam_dir = torch.tensor([0., 0., 1.])  # Lab frame

        reflected_dir = get_reflected_ray_dir(sample, normal, beam_dir)

        # Check magnitude is preserved
        assert torch.isclose(torch.norm(reflected_dir), torch.tensor(1.0), atol=1e-5)


class TestBuildReflectedRay:
    """Test build_reflected_ray function."""

    def test_ray_construction(self):
        """Test that ray is correctly constructed."""
        sample = Sample()
        sample.set_location(np.array([1., 0., 0.]))  # Sample offset
        sample.set_orientation_matrix(np.eye(3))

        vertex = torch.tensor([0., 0., 0.])  # Sample frame origin
        ref_dir = torch.tensor([0., 0., 1.])  # Lab frame direction

        ray = build_reflected_ray(sample, vertex, ref_dir)

        # Vertex should be transformed to lab frame
        expected_origin = torch.tensor([1., 0., 0.])
        assert torch.allclose(ray.origin, expected_origin, atol=1e-5)
        assert torch.allclose(ray.direction, ref_dir, atol=1e-5)

    def test_with_rotated_sample(self):
        """Test with rotated sample."""
        sample = Sample()
        sample.set_location(np.array([0., 0., 0.]))

        # Rotate 90° around Z-axis: X -> Y
        orientation = np.array([[0., -1., 0.],
                               [1., 0., 0.],
                               [0., 0., 1.]], dtype=np.float32)
        sample.set_orientation_matrix(orientation)

        vertex = torch.tensor([1., 0., 0.])  # Sample frame
        ref_dir = torch.tensor([0., 0., 1.])

        ray = build_reflected_ray(sample, vertex, ref_dir)

        # Vertex [1,0,0] in sample frame -> [0,1,0] in lab frame
        expected_origin = torch.tensor([0., 1., 0.])
        assert torch.allclose(ray.origin, expected_origin, atol=1e-5)


class TestGetReflectedRay:
    """Test get_reflected_ray convenience wrapper."""

    def test_complete_workflow(self):
        """Test complete ray generation workflow."""
        sample = Sample()
        sample.set_location(np.array([0., 0., 0.]))
        sample.set_orientation_matrix(np.eye(3))

        vertex = torch.tensor([0., 0., 0.])
        normal = torch.tensor([0., 0., 1.])
        beam_dir = torch.tensor([0., 0., 1.])

        ray = get_reflected_ray(sample, vertex, normal, beam_dir)

        # Should be Ray object
        assert isinstance(ray, Ray)
        assert ray.origin is not None
        assert ray.direction is not None


class TestGetIlluminatedPixel:
    """Test get_illuminated_pixel function."""

    def test_ray_hitting_detector_center(self):
        """Test ray hitting detector at beam center."""
        # Create detector
        # Default detector orientation: surface perpendicular to X-axis (in YZ plane)
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            beam_center_j=512.0,
            beam_center_k=512.0,
            pixel_width=0.004,
            pixel_height=0.004,
            position=torch.tensor([1.0, 0., 0.]),  # 1m away along X
            orientation=torch.eye(3)
        )

        # Ray from origin toward detector along +X
        ray = Ray(
            origin=torch.tensor([0., 0., 0.]),
            direction=torch.tensor([1., 0., 0.])
        )

        hit, pixel_col, pixel_row = get_illuminated_pixel(detector, ray)

        assert hit.item() == True
        # Should hit near beam center
        assert abs(pixel_col.item() - 512.0) < 10
        assert abs(pixel_row.item() - 512.0) < 10

    def test_ray_missing_detector(self):
        """Test ray that doesn't intersect detector."""
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            beam_center_j=512.0,
            beam_center_k=512.0,
            pixel_width=0.004,
            pixel_height=0.004,
            position=torch.tensor([1.0, 0., 0.]),
            orientation=torch.eye(3)
        )

        # Ray pointing away from detector
        ray = Ray(
            origin=torch.tensor([0., 0., 0.]),
            direction=torch.tensor([-1., 0., 0.])  # Wrong direction
        )

        hit, pixel_col, pixel_row = get_illuminated_pixel(detector, ray)

        assert hit.item() == False
        assert pixel_col.item() == 0.0
        assert pixel_row.item() == 0.0

    def test_ray_at_angle(self):
        """Test ray hitting detector at an angle."""
        detector = Detector(
            num_rows=1024,
            num_cols=1024,
            beam_center_j=512.0,
            beam_center_k=512.0,
            pixel_width=0.004,
            pixel_height=0.004,
            position=torch.tensor([1.0, 0., 0.]),
            orientation=torch.eye(3)
        )

        # Ray at slight angle in Y direction
        ray = Ray(
            origin=torch.tensor([0., 0., 0.]),
            direction=torch.tensor([1., 0.1, 0.]).float()  # Slight angle
        )

        hit, pixel_col, pixel_row = get_illuminated_pixel(detector, ray)

        assert hit.item() == True
        # Should be offset from center
        assert pixel_col.item() > 512.0  # Offset in +J direction


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
