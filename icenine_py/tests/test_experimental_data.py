"""
Tests for experimental data loading module.

Tests verify:
1. Loading from forward simulation output (direct in-memory)
2. Loading from ASCII files on disk
3. Image dimensions and pixel counts are correct
4. is_bright queries work on loaded data
"""

from pathlib import Path

import pytest
import torch

from icenine.experimental_data import ExperimentalData
from icenine.image_data import ImageData


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def project_root():
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


@pytest.fixture(autouse=True)
def chdir_to_project_root(project_root, monkeypatch):
    """Change to project root so relative paths in config files resolve."""
    monkeypatch.chdir(project_root)


@pytest.fixture
def three_voxels_data_dir(project_root):
    """Path to ThreeVoxels forward simulation output."""
    return project_root / "Examples" / "Example2.ThreeVoxels" / "ScatteringData_Python"


@pytest.fixture
def three_voxels_config_path(project_root):
    """Path to ThreeVoxels simulation config."""
    return project_root / "Examples" / "Example2.ThreeVoxels" / "ConfigFiles" / "Example2.Simulation.config"


# ============================================================================
# Test: from synthetic images (in-memory)
# ============================================================================

class TestFromForwardSimulation:
    """Tests for loading from ForwardSimulation output."""

    def test_from_mock_forward_sim(self):
        """ExperimentalData wraps forward sim images correctly."""
        # Create mock forward simulation with images
        n_omega, n_det = 3, 2
        num_rows, num_cols = 64, 64

        mock_images = []
        for i in range(n_omega):
            det_images = []
            for j in range(n_det):
                img = ImageData(num_rows, num_cols)
                # Set a pixel so it's not empty
                img.set_pixel(10 + i, 20 + j, 1.0)
                det_images.append(img)
            mock_images.append(det_images)

        class MockForwardSim:
            images = mock_images

        exp_data = ExperimentalData.from_forward_simulation(MockForwardSim())

        assert exp_data.n_omega_intervals == n_omega
        assert exp_data.n_detectors == n_det

        # Verify pixel is preserved: omega=1, det=0 → set_pixel(j=11, k=20)
        img = exp_data.get_image(1, 0)
        assert img.is_bright(11, 20).item()
        assert not img.is_bright(0, 0).item()

    def test_from_empty_forward_sim_raises(self):
        """Should raise if forward sim has no images."""
        class MockForwardSim:
            images = []

        with pytest.raises(ValueError, match="no images"):
            ExperimentalData.from_forward_simulation(MockForwardSim())


# ============================================================================
# Test: loading from disk
# ============================================================================

class TestFromImageDirectory:
    """Tests for loading from ASCII files on disk."""

    @pytest.fixture
    def data_dir(self, three_voxels_data_dir):
        if not three_voxels_data_dir.exists():
            pytest.skip(f"Forward sim output not found: {three_voxels_data_dir}")
        return three_voxels_data_dir

    def test_load_single_omega(self, data_dir):
        """Load first omega step (2 detectors) from disk."""
        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=1,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
            file_start=0,
            det_offset=0,
        )

        assert exp_data.n_omega_intervals == 1
        assert exp_data.n_detectors == 2

        img0 = exp_data.get_image(0, 0)
        img1 = exp_data.get_image(0, 1)
        assert img0.num_rows == 2048
        assert img0.num_cols == 2048
        assert img1.num_rows == 2048

    def test_load_multiple_omegas(self, data_dir):
        """Load first 5 omega steps."""
        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=5,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
        )

        assert exp_data.n_omega_intervals == 5
        assert exp_data.n_detectors == 2

    def test_bright_pixels_present(self, data_dir):
        """Loaded data should have some bright pixels (non-empty images)."""
        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
        )

        total_bright = exp_data.count_bright_pixels()
        # ThreeVoxels forward sim produces ~3200 bright pixels total
        assert total_bright > 1000, f"Expected >1000 bright pixels, got {total_bright}"

    def test_missing_file_raises(self, data_dir):
        """Should raise FileNotFoundError for missing files."""
        with pytest.raises(FileNotFoundError):
            ExperimentalData.from_image_directory(
                directory=str(data_dir),
                basename="3Grains.sim",
                ext="d",
                serial_length=5,
                n_omega=200,  # Only 180 exist
                n_detectors=2,
                num_rows=2048,
                num_cols=2048,
            )


# ============================================================================
# Test: repr and metadata
# ============================================================================

class TestExperimentalDataMeta:
    """Tests for metadata and repr."""

    def test_repr(self):
        """repr should show dimensions."""
        images = [[ImageData(64, 64) for _ in range(2)] for _ in range(3)]
        exp_data = ExperimentalData(images, n_omega_intervals=3, n_detectors=2)
        assert "omega_intervals=3" in repr(exp_data)
        assert "detectors=2" in repr(exp_data)

    def test_count_bright_empty(self):
        """Empty images should have 0 bright pixels."""
        images = [[ImageData(64, 64) for _ in range(2)] for _ in range(3)]
        exp_data = ExperimentalData(images, n_omega_intervals=3, n_detectors=2)
        assert exp_data.count_bright_pixels() == 0
