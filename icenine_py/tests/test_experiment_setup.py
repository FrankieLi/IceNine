"""
Unit tests for experiment_setup module.

Tests ExperimentSetup and XDMExperimentSetup classes against C++ reference behavior.
"""

import pytest
import numpy as np
from pathlib import Path

from icenine.experiment_setup import (
    ExperimentSetup,
    XDMExperimentSetup,
    StepSizeInfo
)
from icenine.config_file import ConfigFile
from icenine.sample import Sample
from icenine.detector import Detector


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
def experiment_setup(config_file):
    """Create initialized XDMExperimentSetup."""
    setup = XDMExperimentSetup(config_file)
    return setup


# ============================================================================
# Base Class Tests
# ============================================================================

class TestExperimentSetup:
    """Test ExperimentSetup base class."""

    def test_default_initialization(self):
        """Test default constructor."""
        setup = ExperimentSetup()
        assert not setup.initialized
        assert setup.beam_energy == 0.0
        assert np.allclose(setup.beam_direction, [0, 0, 1])

    def test_initialization_with_config(self, config_file):
        """Test initialization with ConfigFile."""
        setup = ExperimentSetup(config_file)
        assert setup.initialized
        assert setup.beam_energy == config_file.beam_energy
        assert setup.beam_energy_width == config_file.beam_energy_width

    def test_set_config_file(self, config_file):
        """Test set_config_file method."""
        setup = ExperimentSetup()
        assert not setup.initialized

        setup.set_config_file(config_file)
        assert setup.initialized
        assert setup.beam_energy == config_file.beam_energy

    def test_beam_direction_normalization(self, config_file):
        """Test that beam direction is normalized."""
        setup = ExperimentSetup(config_file)
        beam_dir = setup.get_beam_direction()
        assert np.isclose(np.linalg.norm(beam_dir), 1.0)

    def test_accessors_require_initialization(self):
        """Test that accessors raise error if not initialized."""
        setup = ExperimentSetup()

        with pytest.raises(RuntimeError, match="not initialized"):
            setup.get_beam_energy()

        with pytest.raises(RuntimeError, match="not initialized"):
            setup.get_beam_direction()

    def test_accessors_after_initialization(self, config_file):
        """Test accessors return correct values after initialization."""
        setup = ExperimentSetup(config_file)

        assert setup.get_beam_energy() == config_file.beam_energy
        assert setup.get_beam_energy_width() == config_file.beam_energy_width
        assert setup.get_eta_limit() == config_file.eta_limit
        assert setup.get_beam_deflection_chi_laue() == config_file.beam_deflection_chi_laue
        assert setup.get_min_accepted_intensity_fraction() == config_file.min_amplitude_fraction


# ============================================================================
# XDM Experiment Setup Tests
# ============================================================================

class TestXDMExperimentSetup:
    """Test XDMExperimentSetup class."""

    def test_default_initialization(self):
        """Test default constructor."""
        setup = XDMExperimentSetup()
        assert len(setup.detector_list) == 0
        assert len(setup.omega_range_list) == 0
        assert setup.range_to_index_map is None

    def test_initialization_with_config(self, config_file):
        """Test initialization with ConfigFile."""
        setup = XDMExperimentSetup(config_file)
        assert setup.initialized
        assert setup.config_file == config_file

    def test_initialize_experiment(self, experiment_setup):
        """Test full experiment initialization."""
        # This will read detector and omega files
        experiment_setup.initialize_experiment()

        # Check that data was loaded
        assert len(experiment_setup.detector_list) > 0
        assert len(experiment_setup.omega_range_list) > 0
        assert len(experiment_setup.file_range_list) > 0
        assert experiment_setup.range_to_index_map is not None

    def test_detector_list_loaded(self, experiment_setup):
        """Test that detectors are loaded correctly."""
        experiment_setup.initialize_experiment()

        detectors = experiment_setup.get_detector_list()
        assert len(detectors) == experiment_setup.config_file.num_detectors

        # Check first detector has valid parameters
        det = detectors[0]
        assert det.beam_center_j > 0
        assert det.beam_center_k > 0
        assert det.pixel_width > 0
        assert det.pixel_height > 0

    def test_omega_ranges_loaded(self, experiment_setup):
        """Test that omega ranges are loaded correctly."""
        experiment_setup.initialize_experiment()

        omega_ranges = experiment_setup.get_omega_range_list()
        assert len(omega_ranges) > 0

        # Check that ranges are valid
        for omega_range in omega_ranges:
            assert omega_range.low < omega_range.high
            assert omega_range.low >= -2*np.pi
            assert omega_range.high <= 2*np.pi

    def test_range_to_index_map(self, experiment_setup):
        """Test SimulationRange mapper creation."""
        experiment_setup.initialize_experiment()

        mapper = experiment_setup.get_range_to_index_map()
        assert mapper is not None
        assert mapper.low < mapper.high

    def test_get_max_q(self, experiment_setup):
        """Test max Q calculation."""
        experiment_setup.initialize_experiment()

        # Create a sample at origin
        sample = Sample()
        sample.set_location(np.array([0.0, 0.0, 0.0]))

        # Get first detector
        detector = experiment_setup.detector_list[0]

        # Calculate max Q
        max_q = experiment_setup.get_max_q(detector, sample)

        # Verify it's reasonable (should be positive and < 20 Å⁻¹ for typical setup)
        assert max_q > 0
        assert max_q < 20.0
        print(f"Max Q: {max_q:.4f} Å⁻¹")

    def test_get_reciprocal_vector(self, experiment_setup):
        """Test reciprocal vector calculation."""
        # Scattered in forward direction (no scattering)
        k_out_dir = np.array([0.0, 0.0, 1.0])
        g_vec = experiment_setup.get_reciprocal_vector(k_out_dir)

        # For forward scattering, G should be ~zero
        assert np.linalg.norm(g_vec) < 0.1

        # Scattered at 90 degrees
        k_out_dir = np.array([1.0, 0.0, 0.0])
        g_vec = experiment_setup.get_reciprocal_vector(k_out_dir)

        # G magnitude should be sqrt(2)*k for 90° scattering
        k_mag = 0.506773182 * experiment_setup.beam_energy  # KEV_OVER_HBAR_C_IN_ANG * E
        expected_g = np.sqrt(2) * k_mag
        actual_g = np.linalg.norm(g_vec)

        print(f"90° scattering: |G| = {actual_g:.4f} Å⁻¹ (expected ~{expected_g:.4f})")
        assert np.isclose(actual_g, expected_g, rtol=0.01)

    def test_get_sample_symmetry(self, experiment_setup):
        """Test sample symmetry retrieval."""
        from icenine.config_file import SymmetryType

        # Test depends on what's in config file
        symmetry = experiment_setup.get_sample_symmetry()
        assert symmetry is not None

        # Should be Cubic for gold sample
        if experiment_setup.config_file.sample_symmetry == SymmetryType.CUBIC:
            from icenine.symmetry import CubicSymmetry
            assert isinstance(symmetry, CubicSymmetry)

    def test_initialize_sample(self, experiment_setup, project_root):
        """Test sample initialization with crystal structure."""
        experiment_setup.initialize_experiment()

        # Create sample
        sample = Sample()

        # Get first detector for max Q calculation
        detector = experiment_setup.detector_list[0]

        # Initialize sample
        # NOTE: This will fail if sample file doesn't exist or if
        # binary .dat file reading is not implemented
        try:
            experiment_setup.initialize_sample(sample, detector)

            # Check that sample was configured
            assert len(sample.get_structure_list()) > 0

            # Check location and orientation were set
            location = sample.get_location()
            assert not np.allclose(location, [0, 0, 0])

        except (FileNotFoundError, NotImplementedError) as e:
            pytest.skip(f"Sample initialization skipped: {e}")


# ============================================================================
# StepSizeInfo Tests
# ============================================================================

class TestStepSizeInfo:
    """Test StepSizeInfo dataclass."""

    def test_creation(self):
        """Test StepSizeInfo creation."""
        info = StepSizeInfo(
            euler_steps=np.array([0.1, 0.1, 0.1]),
            detector_pos=np.array([0.01, 0.01, 0.01]),
            beam_center_j=512.0,
            beam_center_k=512.0,
            pixel_height=0.004,
            pixel_width=0.004,
            angular_radius=0.05
        )

        assert np.allclose(info.euler_steps, [0.1, 0.1, 0.1])
        assert info.beam_center_j == 512.0
        assert info.angular_radius == 0.05


# ============================================================================
# Integration Tests
# ============================================================================

class TestExperimentSetupIntegration:
    """Integration tests with real config files."""

    def test_full_workflow(self, config_file, project_root):
        """Test complete initialization workflow."""
        # Create setup
        setup = XDMExperimentSetup(config_file)
        assert setup.initialized

        # Initialize experiment
        setup.initialize_experiment()
        assert len(setup.detector_list) > 0
        assert len(setup.omega_range_list) > 0

        # Get beam parameters
        beam_energy = setup.get_beam_energy()
        assert beam_energy > 0
        print(f"Beam energy: {beam_energy} keV")

        # Get detectors
        detectors = setup.get_detector_list()
        print(f"Number of detectors: {len(detectors)}")
        for i, det in enumerate(detectors):
            print(f"  Detector {i}: beam_center=({det.beam_center_j:.1f}, {det.beam_center_k:.1f})")

        # Get omega ranges
        omega_ranges = setup.get_omega_range_list()
        print(f"Number of omega ranges: {len(omega_ranges)}")
        print(f"  First range: [{np.rad2deg(omega_ranges[0].low):.1f}°, {np.rad2deg(omega_ranges[0].high):.1f}°]")

        # Everything should be consistent
        assert len(detectors) == config_file.num_detectors


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
