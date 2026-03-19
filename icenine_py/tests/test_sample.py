"""
Tests for Sample class.

Validates coordinate transformations, file I/O, and crystal structure management
against C++ reference implementation.
"""

import pytest
import numpy as np
import torch
from pathlib import Path

from icenine.sample import Sample
from icenine.crystal_structure import CrystalStructure
from icenine.symmetry import CrystalSymmetry


class TestSample:
    """Test suite for Sample class."""

    def test_construction(self):
        """Test default construction creates identity transformation."""
        sample = Sample()

        # Should have identity transformation matrix
        expected_matrix = torch.eye(4, dtype=torch.float32)
        assert torch.allclose(sample.sample_to_lab_matrix, expected_matrix)

        # Should have zero location and orientation
        assert torch.allclose(sample.location, torch.zeros(3))
        assert torch.allclose(sample.orientation_euler, torch.zeros(3))

        # Should have no mic file or crystal structures
        assert sample.mic_file is None
        assert len(sample.crystal_structures) == 0

    def test_set_location(self):
        """Test setting sample location."""
        sample = Sample()
        location = np.array([1.0, 2.0, 3.0])

        sample.set_location(location)

        # Check location stored correctly
        assert np.allclose(sample.get_location(), location)

        # Check translation part of matrix updated
        assert torch.allclose(
            sample.sample_to_lab_matrix[:3, 3],
            torch.tensor(location, dtype=torch.float32)
        )

    def test_translate(self):
        """Test applying translation."""
        sample = Sample()

        # Set initial location
        sample.set_location(np.array([1.0, 0.0, 0.0]))

        # Apply translation
        sample.translate(np.array([0.0, 2.0, 0.0]))

        # Should be at [1, 2, 0]
        expected = np.array([1.0, 2.0, 0.0])
        assert np.allclose(sample.get_location(), expected)

    def test_set_orientation_90deg_z_rotation(self):
        """Test 90° rotation around Z-axis."""
        sample = Sample()

        # Rotate 90° around Z (phi=90, theta=0, psi=0)
        sample.set_orientation(90, 0, 0)

        # Transform vector [1, 0, 0] → should become [0, -1, 0]
        # This is per C++ SetPassiveEulerMatrix convention
        v_sample = torch.tensor([1.0, 0.0, 0.0])
        v_lab = sample.to_lab_frame(v_sample)

        expected = torch.tensor([0.0, -1.0, 0.0])
        assert torch.allclose(v_lab, expected, atol=1e-6), \
            f"Expected {expected}, got {v_lab}"

    def test_set_orientation_180deg_z_rotation(self):
        """Test 180° rotation around Z-axis."""
        sample = Sample()

        # Rotate 180° around Z
        sample.set_orientation(180, 0, 0)

        # Transform vector [1, 0, 0] → should become [-1, 0, 0]
        v_sample = torch.tensor([1.0, 0.0, 0.0])
        v_lab = sample.to_lab_frame(v_sample)

        expected = torch.tensor([-1.0, 0.0, 0.0])
        assert torch.allclose(v_lab, expected, atol=1e-6)

    def test_rotate_composition(self):
        """Test rotation composition.

        C++ SetOrientation uses PassiveEuler, Rotate uses ActiveEuler.
        For pure Z rotation: Active(phi) @ Passive(phi) = Identity,
        so we test with a non-trivial composition instead.

        Verify that rotate_z (active, radians) composed with
        set_orientation (passive) gives correct combined result.
        """
        sample = Sample()

        # Set orientation with passive 45° Z rotation
        sample.set_orientation(45, 0, 0)

        # Then compose with active 90° Z rotation via rotate()
        sample.rotate(90, 0, 0)

        # Active(90) @ Passive(45):
        # Active(90) = [[0, -1], [1, 0]]
        # Passive(45) = [[cos45, sin45], [-sin45, cos45]]
        # Product = [[sin45, -cos45], [cos45, sin45]]
        # Applied to [1,0,0] → [sin45, cos45, 0]
        v_sample = torch.tensor([1.0, 0.0, 0.0])
        v_lab = sample.to_lab_frame(v_sample)

        sin45 = np.sin(np.deg2rad(45))
        cos45 = np.cos(np.deg2rad(45))
        expected = torch.tensor([sin45, cos45, 0.0], dtype=torch.float32)
        assert torch.allclose(v_lab, expected, atol=1e-6)

    def test_to_lab_frame_batched(self):
        """Test batched vector transformation."""
        sample = Sample()
        sample.set_orientation(90, 0, 0)

        # Transform multiple vectors at once
        vectors = torch.tensor([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ])

        transformed = sample.to_lab_frame(vectors)

        expected = torch.tensor([
            [0.0, -1.0, 0.0],   # [1,0,0] → [0,-1,0]
            [1.0, 0.0, 0.0],    # [0,1,0] → [1,0,0]
            [0.0, 0.0, 1.0]     # [0,0,1] → [0,0,1] (Z unchanged)
        ])

        assert torch.allclose(transformed, expected, atol=1e-6)

    def test_to_lab_frame_identity(self):
        """Test that identity transformation doesn't change vectors."""
        sample = Sample()  # Identity by default

        v = torch.tensor([1.0, 2.0, 3.0])
        v_transformed = sample.to_lab_frame(v)

        assert torch.allclose(v, v_transformed)

    def test_rotate_axis_angle(self):
        """Test axis-angle rotation.

        C++ uses ACTIVE rotation (Sample.cpp:138: BuildRotationAboutAxis).
        Active 90° Z rotation: [1,0,0] → [0,1,0] (counterclockwise).
        """
        sample = Sample()

        # Rotate 90° around Z-axis using axis-angle
        axis = np.array([0.0, 0.0, 1.0])
        sample.rotate_axis_angle(axis, 90)

        # Active rotation: [1, 0, 0] → [0, 1, 0]
        v = torch.tensor([1.0, 0.0, 0.0])
        v_lab = sample.to_lab_frame(v)

        expected = torch.tensor([0.0, 1.0, 0.0])
        assert torch.allclose(v_lab, expected, atol=1e-6)

    def test_rotate_z_optimized(self):
        """Test optimized Z rotation.

        rotate_z takes RADIANS (matching C++ RotateZ which calls cos/sin directly).
        C++ uses ACTIVE rotation (Sample.cpp:162 comment).
        """
        sample = Sample()
        sample.rotate_z(np.pi / 2)  # 90° in radians

        # Active rotation: [1, 0, 0] → [0, 1, 0]
        v = torch.tensor([1.0, 0.0, 0.0])
        v_lab = sample.to_lab_frame(v)

        expected = torch.tensor([0.0, 1.0, 0.0])
        assert torch.allclose(v_lab, expected, atol=1e-6)

    def test_get_orientation_matrix(self):
        """Test retrieving orientation as rotation matrix."""
        sample = Sample()
        sample.set_orientation(90, 0, 0)

        R = sample.get_orientation_matrix()

        # Should be 3x3 rotation matrix
        assert R.shape == (3, 3)

        # Should rotate [1,0,0] → [0,-1,0]
        v = np.array([1.0, 0.0, 0.0])
        v_rotated = R @ v

        expected = np.array([0.0, -1.0, 0.0])
        assert np.allclose(v_rotated, expected, atol=1e-6)

    def test_set_orientation_matrix(self):
        """Test setting orientation from rotation matrix."""
        sample = Sample()

        # 90° passive rotation matrix
        R = np.array([
            [0, 1, 0],
            [-1, 0, 0],
            [0, 0, 1]
        ], dtype=np.float32)

        sample.set_orientation_matrix(R)

        # Transform [1, 0, 0]
        v = torch.tensor([1.0, 0.0, 0.0])
        v_lab = sample.to_lab_frame(v)

        expected = torch.tensor([0.0, -1.0, 0.0])
        assert torch.allclose(v_lab, expected, atol=1e-6)

    def test_add_crystal_structure(self):
        """Test adding crystal structures."""
        sample = Sample()

        # Create FCC gold structure
        gold = CrystalStructure.create_fcc("Au", 4.0782)
        sample.add_crystal_structure(gold)

        assert len(sample.get_structure_list()) == 1
        # Check lattice parameter from pymatgen structure
        assert abs(sample.get_structure_list()[0].structure.lattice.a - 4.0782) < 0.001

    def test_multiple_crystal_structures(self):
        """Test multi-phase material support."""
        sample = Sample()

        # Add multiple phases
        gold = CrystalStructure.create_fcc("Au", 4.0782)
        copper = CrystalStructure.create_fcc("Cu", 3.615)

        sample.add_crystal_structure(gold)
        sample.add_crystal_structure(copper)

        structures = sample.get_structure_list()
        assert len(structures) == 2
        assert abs(structures[0].structure.lattice.a - 4.0782) < 0.001
        assert abs(structures[1].structure.lattice.a - 3.615) < 0.001

    def test_set_sample_symmetry(self):
        """Test setting crystal symmetry."""
        sample = Sample()

        # Create symmetry from a crystal structure
        gold = CrystalStructure.create_fcc("Au", 4.0782)
        symmetry = CrystalSymmetry(gold.structure)
        sample.set_sample_symmetry(symmetry)

        assert sample.get_sample_symmetry() is not None
        assert sample.get_sample_symmetry().point_group == "m-3m"  # Cubic

    def test_num_voxels_no_mic(self):
        """Test voxel count with no mic file loaded."""
        sample = Sample()
        assert sample.num_voxels() == 0

    def test_repr(self):
        """Test string representation."""
        sample = Sample()
        sample.set_location(np.array([1.0, 2.0, 3.0]))
        sample.set_orientation(90, 0, 0)

        repr_str = repr(sample)
        assert "Sample" in repr_str
        assert "location" in repr_str
        assert "orientation" in repr_str


class TestSampleFileIO:
    """Test file I/O operations."""

    def test_load_sample_nonexistent_file(self):
        """Test loading nonexistent file fails gracefully."""
        sample = Sample()
        success = sample.load_sample("nonexistent_file.mic")

        assert success is False
        assert sample.mic_file is None

    def test_load_sample_success(self):
        """Test loading actual .mic file."""
        sample = Sample()

        # Try to load Au1007_small.mic if it exists
        mic_path = Path(__file__).parent.parent.parent.parent / "DataFiles" / "Au1007_small.mic"

        if mic_path.exists():
            success = sample.load_sample(str(mic_path))
            assert success is True
            assert sample.mic_file is not None
            assert sample.num_voxels() > 0
            print(f"Loaded {sample.num_voxels()} voxels")
        else:
            pytest.skip(f"Test file not found: {mic_path}")

    def test_get_mic(self):
        """Test accessing mic file."""
        sample = Sample()

        # No mic file initially
        assert sample.get_mic() is None

        # Load mic file
        mic_path = Path(__file__).parent.parent.parent.parent / "DataFiles" / "Au1007_small.mic"

        if mic_path.exists():
            sample.load_sample(str(mic_path))
            mic = sample.get_mic()

            assert mic is not None
            assert hasattr(mic, 'voxels')
            assert len(mic.voxels) > 0
        else:
            pytest.skip(f"Test file not found: {mic_path}")


class TestSampleTransformationsAdvanced:
    """Advanced coordinate transformation tests."""

    def test_combined_rotation_translation(self):
        """Test combined rotation and translation."""
        sample = Sample()

        # Set orientation and location
        sample.set_orientation(90, 0, 0)  # 90° around Z
        sample.set_location(np.array([10.0, 0.0, 0.0]))

        # Get full transformation matrix
        T = sample.get_transformation_matrix()

        # Check rotation part
        R = T[:3, :3]
        v = np.array([1.0, 0.0, 0.0])
        v_rotated = R @ v
        assert np.allclose(v_rotated, [0.0, -1.0, 0.0], atol=1e-6)

        # Check translation part
        t = T[:3, 3]
        assert np.allclose(t, [10.0, 0.0, 0.0])

    def test_euler_angle_roundtrip(self):
        """Test Euler angle extraction after setting orientation."""
        sample = Sample()

        # Set specific Euler angles
        phi, theta, psi = 45.0, 30.0, 60.0
        sample.set_orientation(phi, theta, psi)

        # Get orientation back
        retrieved = sample.get_orientation()

        # Convert back to degrees
        retrieved_deg = np.rad2deg(retrieved)

        # Should match (within numerical precision)
        # Note: Euler angles may have multiple representations
        # So we check the rotation matrix instead
        R_original = sample.get_orientation_matrix()

        sample2 = Sample()
        sample2.set_orientation(
            retrieved_deg[0],
            retrieved_deg[1],
            retrieved_deg[2]
        )
        R_retrieved = sample2.get_orientation_matrix()

        assert np.allclose(R_original, R_retrieved, atol=1e-6)

    def test_orthogonality_of_rotation_matrix(self):
        """Test that rotation matrices are orthogonal."""
        sample = Sample()

        # Set arbitrary orientation
        sample.set_orientation(45, 30, 60)

        R = sample.get_orientation_matrix()

        # Check R^T @ R = I (orthogonality)
        I = R.T @ R
        expected_I = np.eye(3)
        assert np.allclose(I, expected_I, atol=1e-6)

        # Check det(R) = 1 (proper rotation, not reflection)
        det = np.linalg.det(R)
        assert np.allclose(det, 1.0, atol=1e-6)

    def test_inverse_transformation(self):
        """Test that lab→sample is inverse of sample→lab."""
        sample = Sample()
        sample.set_orientation(45, 30, 60)

        # Transform sample → lab → sample should give identity
        v_sample = torch.tensor([1.0, 2.0, 3.0])
        v_lab = sample.to_lab_frame(v_sample)

        # Inverse transformation (lab → sample)
        R = sample.get_orientation_matrix()
        R_inv = np.linalg.inv(R)
        v_back = R_inv @ v_lab.numpy()

        assert np.allclose(v_back, v_sample.numpy(), atol=1e-6)


class TestSamplePerformance:
    """Performance-related tests."""

    def test_batched_transformation_shape(self):
        """Test that batched transformation preserves shape."""
        sample = Sample()
        sample.set_orientation(90, 0, 0)

        # Test various batch sizes
        for n in [1, 10, 100, 1000]:
            vectors = torch.randn(n, 3)
            transformed = sample.to_lab_frame(vectors)

            assert transformed.shape == (n, 3)

    def test_single_vs_batched_consistency(self):
        """Test that single and batched transformations give same result."""
        sample = Sample()
        sample.set_orientation(45, 30, 60)

        # Create test vectors
        vectors = torch.randn(10, 3)

        # Transform one at a time
        transformed_single = torch.stack([
            sample.to_lab_frame(v) for v in vectors
        ])

        # Transform as batch
        transformed_batch = sample.to_lab_frame(vectors)

        assert torch.allclose(transformed_single, transformed_batch, atol=1e-6)

    def test_gradient_preservation(self):
        """Test that PyTorch gradients are preserved."""
        sample = Sample()
        sample.set_orientation(90, 0, 0)

        # Create vector with gradient tracking
        v = torch.tensor([1.0, 0.0, 0.0], requires_grad=True)

        # Transform
        v_lab = sample.to_lab_frame(v)

        # Compute some scalar loss
        loss = v_lab.sum()

        # Backpropagate
        loss.backward()

        # Check gradient exists
        assert v.grad is not None
        assert v.grad.shape == (3,)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
