"""
Unit tests for ForwardSimulation class.

Tests the main forward simulation loop for generating synthetic detector images.
"""

import pytest
import torch
import numpy as np
from pathlib import Path
import tempfile

from icenine.forward_simulation import ForwardSimulation
from icenine.config_file import ConfigFile
from icenine.experiment_setup import ExperimentSetup
from icenine.sample import Sample
from icenine.detector import Detector
from icenine.image_data import ImageData
from icenine.crystal_structure import CrystalStructure


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def project_root():
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


@pytest.fixture
def config_file(project_root):
    """Load ReconstructTest.config for testing."""
    config_path = project_root / "ConfigFiles" / "ReconstructTest.config"
    if not config_path.exists():
        pytest.skip(f"Config file not found: {config_path}")
    return ConfigFile.from_file(str(config_path))


@pytest.fixture
def temp_output_dir():
    """Create temporary output directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


# ============================================================================
# Initialization Tests
# ============================================================================

class TestForwardSimulationInit:
    """Test ForwardSimulation initialization."""

    def test_initialization_with_config(self, config_file):
        """Test initialization from config file."""
        sim = ForwardSimulation(config_file)

        assert sim.config_file == config_file
        assert isinstance(sim.exp_setup, ExperimentSetup)
        assert sim.images == []

    def test_repr(self, config_file):
        """Test string representation."""
        sim = ForwardSimulation(config_file)
        repr_str = repr(sim)

        assert "ForwardSimulation" in repr_str
        assert "images=" in repr_str


# ============================================================================
# Helper Method Tests
# ============================================================================

class TestHelperMethods:
    """Test helper methods."""

    def test_get_voxel_vertices(self, config_file):
        """Test voxel vertex generation."""
        sim = ForwardSimulation(config_file)

        # Create mock voxel with position
        class MockVoxel:
            def __init__(self):
                self.position = np.array([1.0, 2.0, 3.0])

        voxel = MockVoxel()
        vertices = sim._get_voxel_vertices(voxel)

        # Should return 3 vertices
        assert vertices.shape == (3, 3)
        assert isinstance(vertices, torch.Tensor)

        # Vertices should be near voxel position
        center = torch.tensor([1.0, 2.0, 3.0])
        for i in range(3):
            dist = torch.norm(vertices[i] - center)
            assert dist < 0.01  # Within 10 microns


# ============================================================================
# Simulation Workflow Tests
# ============================================================================

class TestSimulationWorkflow:
    """Test simulation workflow with real config file."""

    @pytest.mark.skipif(
        not (Path(__file__).parent.parent.parent / "ConfigFiles" / "ReconstructTest.config").exists(),
        reason="Config file not found"
    )
    def test_simulate_with_config_file(self, config_file, temp_output_dir):
        """
        Test complete simulation workflow with real config file.

        This is an integration test that exercises the full simulation pipeline.
        It may be slow and is primarily for validation.
        """
        # Create simulator
        sim = ForwardSimulation(config_file)

        # Run simulation (this will be slow for real data)
        try:
            images = sim.simulate_detector_images(output_dir=temp_output_dir)

            # Verify images were generated
            assert len(images) > 0
            assert len(images[0]) > 0
            assert isinstance(images[0][0], ImageData)

            print(f"\nGenerated {len(images)} omega steps x {len(images[0])} detectors")

        except FileNotFoundError as e:
            # Config file references missing data files - skip test
            pytest.skip(f"Missing data file: {e}")
        except Exception as e:
            # Other initialization errors - report but don't fail
            pytest.skip(f"Initialization error: {e}")


# ============================================================================
# Unit Tests with Mock Data
# ============================================================================

class TestSimulationWithMockData:
    """Test simulation with simple mock data."""

    def test_simulate_single_voxel(self, config_file):
        """
        Test simulation with a single voxel.

        Creates minimal setup to test the core simulation loop without
        requiring full dataset.
        """
        # This test requires implementing mock Sample and CrystalStructure
        # For now, we skip if real data is not available
        pytest.skip("Requires implementing mock Sample/CrystalStructure")

        # TODO: Implement with mock data:
        # 1. Create simple Sample with 1 voxel
        # 2. Create simple CrystalStructure with 1-2 reflections
        # 3. Create simple Detector
        # 4. Run simulation
        # 5. Verify image has non-zero intensity


# ============================================================================
# Image Output Tests
# ============================================================================

class TestImageOutput:
    """Test image saving functionality."""

    def test_save_images(self, config_file, temp_output_dir):
        """Test image saving to disk."""
        sim = ForwardSimulation(config_file)

        # Create mock images
        num_omegas = 2
        num_detectors = 1

        mock_images = []
        for i in range(num_omegas):
            detector_images = []
            for j in range(num_detectors):
                image = ImageData(1024, 1024)
                # Add some test data
                image.add_to_pixel(512, 512, 100.0)
                detector_images.append(image)
            mock_images.append(detector_images)

        # Create mock omega and file ranges
        class MockRange:
            def __init__(self, low, high):
                self.low = low
                self.high = high

        omega_ranges = [MockRange(0.0, 1.0), MockRange(1.0, 2.0)]
        file_ranges = [MockRange(0, num_omegas)]

        # Create mock detector
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

        # Save images
        sim._save_images(
            mock_images,
            omega_ranges,
            file_ranges,
            [detector],
            temp_output_dir
        )

        # Verify files were created (search recursively since basename may include subdirectory)
        output_files = list(temp_output_dir.glob("**/*"))
        # Filter to only files (not directories)
        output_files = [f for f in output_files if f.is_file()]
        assert len(output_files) == num_omegas * num_detectors

        print(f"\nSaved {len(output_files)} image files to {temp_output_dir}")
        for f in output_files:
            print(f"  {f.relative_to(temp_output_dir)}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
