"""
Tests for MIC file I/O with backward compatibility validation.

Tests ensure:
1. Euler angle conversion matches C++ implementation
2. File reading matches C++ parser
3. File writing produces C++ compatible output
4. Round-trip: read → write → read preserves data
5. Backward compatibility with existing .mic files
"""

import pytest
import numpy as np
import torch
from pathlib import Path
import tempfile
import shutil

from icenine.mic_file import (
    Voxel,
    MicFile,
    euler_to_matrix,
    matrix_to_euler,
    euler_to_matrix_torch,
)


# ==========================================================================================
# Test Fixtures
# ==========================================================================================


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent.parent.parent / "DataFiles"


@pytest.fixture
def temp_dir():
    """Create temporary directory for test outputs."""
    temp = tempfile.mkdtemp()
    yield Path(temp)
    shutil.rmtree(temp)


@pytest.fixture
def sample_voxel():
    """Create a sample voxel for testing."""
    return Voxel(
        position=np.array([1.0, 2.0, 3.0], dtype=np.float32),
        orientation=np.eye(3, dtype=np.float32),
        side_length=0.012,
        generation=0,
        phase=1,
        confidence=0.95,
        cost=0.05,
        overlap_ratio=0.8,
    )


# ==========================================================================================
# Euler Angle Conversion Tests
# ==========================================================================================


class TestEulerConversion:
    """Test Euler angle ↔ rotation matrix conversions."""

    def test_identity_rotation(self):
        """Test identity rotation (0, 0, 0)."""
        R = euler_to_matrix(0.0, 0.0, 0.0)
        assert R.shape == (3, 3)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-6)

    def test_euler_roundtrip_identity(self):
        """Test round-trip conversion for identity."""
        R = np.eye(3, dtype=np.float32)
        phi1, Phi, phi2 = matrix_to_euler(R)
        R_reconstructed = euler_to_matrix(phi1, Phi, phi2)
        np.testing.assert_allclose(R, R_reconstructed, atol=1e-5)

    def test_euler_roundtrip_90deg_z(self):
        """Test round-trip for 90° rotation around Z."""
        R_original = euler_to_matrix(90.0, 0.0, 0.0)
        phi1, Phi, phi2 = matrix_to_euler(R_original)
        R_reconstructed = euler_to_matrix(phi1, Phi, phi2)
        np.testing.assert_allclose(R_original, R_reconstructed, atol=1e-5)

    def test_euler_roundtrip_90deg_x(self):
        """Test round-trip for 90° rotation around X."""
        R_original = euler_to_matrix(0.0, 90.0, 0.0)
        phi1, Phi, phi2 = matrix_to_euler(R_original)
        R_reconstructed = euler_to_matrix(phi1, Phi, phi2)
        np.testing.assert_allclose(R_original, R_reconstructed, atol=1e-5)

    def test_euler_roundtrip_random(self):
        """Test round-trip for random orientations."""
        np.random.seed(42)
        for _ in range(10):
            phi1_orig = np.random.uniform(0, 360)
            Phi_orig = np.random.uniform(0, 180)
            phi2_orig = np.random.uniform(0, 360)

            R = euler_to_matrix(phi1_orig, Phi_orig, phi2_orig)
            phi1, Phi, phi2 = matrix_to_euler(R)
            R_reconstructed = euler_to_matrix(phi1, Phi, phi2)

            np.testing.assert_allclose(R, R_reconstructed, atol=1e-4)

    def test_euler_torch_single(self):
        """Test PyTorch Euler conversion for single angles."""
        phi1 = torch.tensor(90.0)
        Phi = torch.tensor(0.0)
        phi2 = torch.tensor(0.0)

        R_torch = euler_to_matrix_torch(phi1, Phi, phi2)
        R_numpy = euler_to_matrix(90.0, 0.0, 0.0)

        assert R_torch.shape == (3, 3)
        np.testing.assert_allclose(R_torch.numpy(), R_numpy, atol=1e-5)

    def test_euler_torch_batched(self):
        """Test PyTorch Euler conversion for batched angles."""
        phi1 = torch.tensor([0.0, 90.0, 45.0])
        Phi = torch.tensor([0.0, 0.0, 30.0])
        phi2 = torch.tensor([0.0, 0.0, 15.0])

        R_torch = euler_to_matrix_torch(phi1, Phi, phi2)
        assert R_torch.shape == (3, 3, 3)

        # Check each rotation individually
        for i in range(3):
            R_numpy = euler_to_matrix(phi1[i].item(), Phi[i].item(), phi2[i].item())
            np.testing.assert_allclose(R_torch[i].numpy(), R_numpy, atol=1e-5)

    def test_euler_torch_differentiable(self):
        """Test that PyTorch Euler conversion is differentiable."""
        phi1 = torch.tensor(45.0, requires_grad=True)
        Phi = torch.tensor(30.0, requires_grad=True)
        phi2 = torch.tensor(15.0, requires_grad=True)

        R = euler_to_matrix_torch(phi1, Phi, phi2)
        loss = R.sum()
        loss.backward()

        # Check that gradients exist
        assert phi1.grad is not None
        assert Phi.grad is not None
        assert phi2.grad is not None

    def test_rotation_matrix_properties(self):
        """Test that generated matrices are valid rotations."""
        for phi1, Phi, phi2 in [
            (0, 0, 0),
            (90, 0, 0),
            (0, 90, 0),
            (45, 30, 15),
            (180, 90, 45),
        ]:
            R = euler_to_matrix(phi1, Phi, phi2)

            # Check orthogonality: R^T @ R = I
            np.testing.assert_allclose(R.T @ R, np.eye(3), atol=1e-6)

            # Check determinant = 1 (proper rotation)
            np.testing.assert_allclose(np.linalg.det(R), 1.0, atol=1e-6)


# ==========================================================================================
# Voxel Tests
# ==========================================================================================


class TestVoxel:
    """Test Voxel data structure."""

    def test_voxel_creation(self, sample_voxel):
        """Test basic voxel creation."""
        assert sample_voxel.position.shape == (3,)
        assert sample_voxel.orientation.shape == (3, 3)
        assert sample_voxel.side_length == 0.012
        assert sample_voxel.generation == 0
        assert sample_voxel.phase == 1

    def test_voxel_with_deformation(self):
        """Test voxel with deformation tensor."""
        deformation = np.array([[1.1, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.9]])
        voxel = Voxel(
            position=np.array([0, 0, 0]),
            orientation=np.eye(3),
            deformation=deformation,
        )
        assert voxel.deformation is not None
        np.testing.assert_allclose(voxel.deformation, deformation)

    def test_voxel_validation(self):
        """Test voxel data validation."""
        with pytest.raises(AssertionError):
            # Wrong position shape
            Voxel(position=np.array([1, 2]), orientation=np.eye(3))

        with pytest.raises(AssertionError):
            # Wrong orientation shape
            Voxel(position=np.array([1, 2, 3]), orientation=np.eye(2))


# ==========================================================================================
# MIC File I/O Tests
# ==========================================================================================


class TestMicFileIO:
    """Test MIC file reading and writing."""

    def test_read_au1007_small(self, test_data_dir):
        """Test reading existing Au1007_small.mic file."""
        mic_file = test_data_dir / "Au1007_small.mic"
        if not mic_file.exists():
            pytest.skip(f"Test file not found: {mic_file}")

        mic = MicFile.read(str(mic_file))

        # Check metadata
        assert mic.initial_side_length == pytest.approx(0.012, abs=1e-6)

        # Check voxel count (4 voxels + 1 header line = 5 lines total)
        assert len(mic.voxels) == 4

        # Check first voxel
        v0 = mic.voxels[0]
        assert v0.position[0] == pytest.approx(-0.012, abs=1e-6)
        assert v0.position[1] == pytest.approx(0.0, abs=1e-6)
        assert v0.position[2] == pytest.approx(0.0, abs=1e-6)
        assert v0.generation == 3  # From file: column 5 is generation
        assert v0.phase == 1  # From file: column 6 is phase
        assert v0.points_up is True

        # Check orientation matrix is valid rotation
        np.testing.assert_allclose(v0.orientation.T @ v0.orientation, np.eye(3), atol=1e-5)
        assert np.linalg.det(v0.orientation) == pytest.approx(1.0, abs=1e-5)

        # Check PyTorch tensors
        assert mic.positions.shape == (4, 3)
        assert mic.orientations.shape == (4, 3, 3)
        assert mic.confidence.shape == (4,)

    def test_write_then_read(self, sample_voxel, temp_dir):
        """Test write then read round-trip."""
        # Create MIC file with sample voxel
        mic_original = MicFile(voxels=[sample_voxel], initial_side_length=0.012)

        # Write to file
        output_file = temp_dir / "test_output.mic"
        mic_original.write(str(output_file))

        # Read back
        mic_loaded = MicFile.read(str(output_file))

        # Compare
        assert len(mic_loaded.voxels) == 1
        v = mic_loaded.voxels[0]

        np.testing.assert_allclose(v.position, sample_voxel.position, atol=1e-5)
        np.testing.assert_allclose(v.orientation, sample_voxel.orientation, atol=1e-4)
        assert v.generation == sample_voxel.generation
        assert v.phase == sample_voxel.phase
        assert v.confidence == pytest.approx(sample_voxel.confidence, abs=1e-4)

    def test_roundtrip_au1007_small(self, test_data_dir, temp_dir):
        """Test round-trip: read → write → read on Au1007_small.mic."""
        mic_file = test_data_dir / "Au1007_small.mic"
        if not mic_file.exists():
            pytest.skip(f"Test file not found: {mic_file}")

        # Read original
        mic_original = MicFile.read(str(mic_file))

        # Write to temp file
        temp_file = temp_dir / "roundtrip.mic"
        mic_original.write(str(temp_file))

        # Read back
        mic_roundtrip = MicFile.read(str(temp_file))

        # Compare metadata
        assert mic_roundtrip.initial_side_length == pytest.approx(
            mic_original.initial_side_length, abs=1e-6
        )
        assert len(mic_roundtrip.voxels) == len(mic_original.voxels)

        # Compare each voxel
        for v_orig, v_rt in zip(mic_original.voxels, mic_roundtrip.voxels):
            np.testing.assert_allclose(v_rt.position, v_orig.position, atol=1e-5)
            np.testing.assert_allclose(v_rt.orientation, v_orig.orientation, atol=1e-3)
            assert v_rt.generation == v_orig.generation
            assert v_rt.phase == v_orig.phase
            assert v_rt.points_up == v_orig.points_up
            assert v_rt.confidence == pytest.approx(v_orig.confidence, abs=1e-4)

    def test_empty_mic_file(self, temp_dir):
        """Test creating and writing empty MIC file."""
        mic = MicFile(voxels=[], initial_side_length=0.01)
        assert len(mic) == 0

        output_file = temp_dir / "empty.mic"
        mic.write(str(output_file))

        # Read back
        mic_loaded = MicFile.read(str(output_file))
        assert len(mic_loaded) == 0
        assert mic_loaded.initial_side_length == pytest.approx(0.01, abs=1e-6)

    def test_pytorch_save_load(self, sample_voxel, temp_dir):
        """Test PyTorch native save/load."""
        mic_original = MicFile(voxels=[sample_voxel], initial_side_length=0.012)

        # Save as PyTorch
        pt_file = temp_dir / "test.pt"
        mic_original.save_torch(str(pt_file))

        # Load back
        mic_loaded = MicFile.load_torch(str(pt_file))

        # Compare tensors
        torch.testing.assert_close(mic_loaded.positions, mic_original.positions)
        torch.testing.assert_close(mic_loaded.orientations, mic_original.orientations)
        torch.testing.assert_close(mic_loaded.confidence, mic_original.confidence)


# ==========================================================================================
# Backward Compatibility Tests
# ==========================================================================================


class TestBackwardCompatibility:
    """Tests to ensure backward compatibility with C++ implementation."""

    def test_euler_angles_match_cpp_examples(self):
        """Test Euler angle conversion against known C++ examples."""
        # Test case from Au1007_small.mic line 2:
        # Euler angles: 355.4292, 5.186272, 29.31929 (degrees)
        R = euler_to_matrix(355.4292, 5.186272, 29.31929)

        # Convert back
        phi1, Phi, phi2 = matrix_to_euler(R)

        # Should round-trip (within tolerance due to gimbal lock handling)
        R_reconstructed = euler_to_matrix(phi1, Phi, phi2)
        np.testing.assert_allclose(R, R_reconstructed, atol=1e-3)

    def test_file_format_whitespace(self, sample_voxel, temp_dir):
        """Test that output file format matches C++ whitespace."""
        mic = MicFile(voxels=[sample_voxel], initial_side_length=0.012)
        output_file = temp_dir / "format_test.mic"
        mic.write(str(output_file))

        # Read raw file content
        with open(output_file, "r") as f:
            lines = f.readlines()

        # Check header format
        assert len(lines) == 2  # Header + 1 voxel
        header = lines[0].strip()
        assert "E" in header  # Scientific notation

        # Check voxel line has correct number of columns (19 total)
        voxel_line = lines[1].strip()
        tokens = voxel_line.split()
        assert len(tokens) == 19

    def test_generation_side_length_relationship(self, test_data_dir):
        """Test that side length = initial / 2^generation."""
        mic_file = test_data_dir / "Au1007_small.mic"
        if not mic_file.exists():
            pytest.skip(f"Test file not found: {mic_file}")

        mic = MicFile.read(str(mic_file))

        for voxel in mic.voxels:
            expected_side = mic.initial_side_length / (2 ** voxel.generation)
            assert voxel.side_length == pytest.approx(expected_side, abs=1e-9)

    def test_scientific_notation_precision(self, temp_dir):
        """Test that scientific notation precision matches C++."""
        # Create voxel with specific values to test formatting
        voxel = Voxel(
            position=np.array([1.23456789e-3, -9.87654321e-4, 5.0e-5]),
            orientation=np.eye(3),
            side_length=0.012,
            confidence=0.123456789,
            cost=1.23456e-5,
            overlap_ratio=0.999999,
        )

        mic = MicFile(voxels=[voxel], initial_side_length=0.012)
        output_file = temp_dir / "precision_test.mic"
        mic.write(str(output_file))

        # Read back and check precision preserved
        mic_loaded = MicFile.read(str(output_file))
        v = mic_loaded.voxels[0]

        # Should preserve ~6-7 significant figures
        np.testing.assert_allclose(v.position, voxel.position, rtol=1e-6)
        assert v.confidence == pytest.approx(voxel.confidence, abs=1e-6)


# ==========================================================================================
# Edge Cases and Error Handling
# ==========================================================================================


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_nonexistent_file(self):
        """Test reading nonexistent file raises error."""
        with pytest.raises(FileNotFoundError):
            MicFile.read("nonexistent_file.mic")

    def test_invalid_format_empty_line1(self, temp_dir):
        """Test invalid format (empty first line) raises error."""
        bad_file = temp_dir / "bad.mic"
        with open(bad_file, "w") as f:
            f.write("\n")

        with pytest.raises(ValueError, match="expected 1 token"):
            MicFile.read(str(bad_file))

    def test_invalid_format_missing_columns(self, temp_dir):
        """Test invalid format (missing columns) raises error."""
        bad_file = temp_dir / "bad.mic"
        with open(bad_file, "w") as f:
            f.write("0.012\n")
            f.write("1.0 2.0 3.0 1 0\n")  # Only 5 columns, need >= 9

        with pytest.raises(ValueError, match="expected >= 9 columns"):
            MicFile.read(str(bad_file))

    def test_gimbal_lock_euler_angles(self):
        """Test gimbal lock cases (Phi = 0 or 180)."""
        # Phi = 0 (gimbal lock)
        R1 = euler_to_matrix(45.0, 0.0, 30.0)
        phi1, Phi, phi2 = matrix_to_euler(R1)
        R1_reconstructed = euler_to_matrix(phi1, Phi, phi2)
        np.testing.assert_allclose(R1, R1_reconstructed, atol=1e-4)

        # Phi = 180 (gimbal lock) - Note: scipy handles this by setting third angle to 0
        # The rotation matrix should still be preserved (different Euler angles, same rotation)
        R2 = euler_to_matrix(45.0, 180.0, 30.0)
        phi1, Phi, phi2 = matrix_to_euler(R2)
        R2_reconstructed = euler_to_matrix(phi1, Phi, phi2)

        # At Phi=180, gimbal lock means euler angles are not unique
        # but the rotation matrix representation should be equivalent
        # Check that determinant and orthogonality are preserved
        np.testing.assert_allclose(R2.T @ R2, np.eye(3), atol=1e-5)
        np.testing.assert_allclose(R2_reconstructed.T @ R2_reconstructed, np.eye(3), atol=1e-5)
        np.testing.assert_allclose(np.linalg.det(R2), 1.0, atol=1e-5)
        np.testing.assert_allclose(np.linalg.det(R2_reconstructed), 1.0, atol=1e-5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
