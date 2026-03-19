"""
Tests for ConfigFile parser.

Validates configuration file parsing against C++ reference implementation.
"""

import pytest
import numpy as np
from pathlib import Path
import tempfile

from icenine.config_file import (
    ConfigFile,
    FileType,
    SO3SearchMethod,
    SymmetryType,
    GridType
)


class TestConfigFileBasics:
    """Basic ConfigFile functionality tests."""

    def test_construction(self):
        """Test default construction."""
        config = ConfigFile()

        # Check defaults
        assert config.beam_energy == 0.0
        assert config.num_detectors == -1
        assert config.in_file_type == FileType.ASCII
        assert config.sample_symmetry == SymmetryType.CUBIC
        assert config.mic_grid_type == GridType.TRIANGULAR

        # Check numpy arrays initialized
        assert isinstance(config.beam_direction, np.ndarray)
        assert config.beam_direction.shape == (3,)

    def test_enums(self):
        """Test enum types."""
        # FileType
        assert FileType.BIN.value == 0
        assert FileType.TIF.value == 1
        assert FileType.ASCII.value == 2

        # SO3SearchMethod
        assert SO3SearchMethod.CONSTRAINED_EULER.value == 0
        assert SO3SearchMethod.UNIFORM_QUATERNION.value == 1

        # SymmetryType
        assert SymmetryType.CUBIC.value == 0
        assert SymmetryType.HEXAGONAL.value == 1

        # GridType
        assert GridType.TRIANGULAR.value == 0
        assert GridType.SQUARE.value == 1


class TestConfigFileParsing:
    """Test parsing of different data types."""

    def test_parse_string(self):
        """Test string parameter parsing."""
        config_text = """
InfileBasename Test_Input_
OutfileBasename Test_Output_
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.config', delete=False) as f:
            f.write(config_text)
            config_file = f.name

        try:
            config = ConfigFile()
            lines = config._tokenize(config_text)
            config._parse_lines(lines, config_file)

            assert config.in_file_basename == "Test_Input_"
            assert config.out_file_basename == "Test_Output_"
        finally:
            Path(config_file).unlink()

    def test_parse_float(self):
        """Test float parameter parsing."""
        config_text = """
BeamEnergy 50.02099
BeamEnergyWidth 0.05
MaxQ 8.0
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert config.beam_energy == 50.02099
        assert config.beam_energy_width == 0.05
        assert config.max_q == 8.0

    def test_parse_int(self):
        """Test integer parameter parsing."""
        config_text = """
MaxDiscreteCandidates 100
NumDetectors 2
MinLocalResolution 0
MaxLocalResolution 3
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert config.max_discrete_candidates == 100
        assert config.num_detectors == 2
        assert config.min_local_resolution == 0
        assert config.max_local_resolution == 3

    def test_parse_bool(self):
        """Test boolean parameter parsing."""
        config_text = """
EnableStrain 1
ConstrainedOptimization 0
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert config.strain_enabled is True
        assert config.constrained_param_mc is False

    def test_parse_flag(self):
        """Test flag parameter parsing (presence = True)."""
        config_text = """
IntensityDecomposition
LazyBFS
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert config.intensity_decomposition is True
        assert config.lazy_bfs is True
        assert config.lazy_strain is False  # Not present

    def test_parse_angle_conversion(self):
        """Test angle conversion from degrees to radians."""
        config_text = """
EtaLimit 81
LocalOrientationGridRadius 5
ConsistencyError 0.5
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert np.isclose(config.eta_limit, np.deg2rad(81))
        assert np.isclose(config.local_orientation_grid_radius, np.deg2rad(5))
        assert np.isclose(config.consistency_error, np.deg2rad(0.5))

    def test_parse_vector3(self):
        """Test 3D vector parsing."""
        config_text = """
BeamDirection 1 0 0
SampleLocation 0 0 0.1
SampleCenter 0.5 0.5 0.5
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert np.allclose(config.beam_direction, [1, 0, 0])
        assert np.allclose(config.sample_location, [0, 0, 0.1])
        assert np.allclose(config.sample_center, [0.5, 0.5, 0.5])

    def test_parse_vector3_with_angle_conversion(self):
        """Test 3D vector parsing with angle conversion."""
        config_text = """
SampleOrientation 45 30 60
DetectorOrientationDeviationInEuler 0.5 0.5 0.5
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        expected_orientation = np.deg2rad([45, 30, 60])
        expected_deviation = np.deg2rad([0.5, 0.5, 0.5])

        assert np.allclose(config.sample_orientation, expected_orientation)
        assert np.allclose(config.detector_orient_deviation_euler, expected_deviation)

    def test_parse_file_type(self):
        """Test file type enum parsing."""
        for file_type_str, file_type_enum in [
            ('bin', FileType.BIN),
            ('tif', FileType.TIF),
            ('ascii', FileType.ASCII),
            ('BIN', FileType.BIN),  # Case insensitive
            ('ASCII', FileType.ASCII),
        ]:
            config_text = f"InFileType {file_type_str}"
            config = ConfigFile()
            lines = config._tokenize(config_text)
            config._parse_lines(lines, "test.config")
            assert config.in_file_type == file_type_enum

    def test_parse_symmetry(self):
        """Test symmetry enum parsing."""
        for sym_str, sym_enum in [
            ('Cubic', SymmetryType.CUBIC),
            ('Hexagonal', SymmetryType.HEXAGONAL),
            ('Tetragonal', SymmetryType.TETRAGONAL),
            ('None', SymmetryType.NONE),
        ]:
            config_text = f"SampleSymmetry {sym_str}"
            config = ConfigFile()
            lines = config._tokenize(config_text)
            config._parse_lines(lines, "test.config")
            assert config.sample_symmetry == sym_enum

    def test_parse_search_method(self):
        """Test search method enum parsing."""
        for method_str, method_enum in [
            ('ConstrainedEuler', SO3SearchMethod.CONSTRAINED_EULER),
            ('UniformSO3', SO3SearchMethod.UNIFORM_QUATERNION),
            ('UniformQuaternion', SO3SearchMethod.UNIFORM_QUATERNION),
        ]:
            config_text = f"OrientationSearchMethod {method_str}"
            config = ConfigFile()
            lines = config._tokenize(config_text)
            config._parse_lines(lines, "test.config")
            assert config.orientation_search_method == method_enum

    def test_parse_grid_type(self):
        """Test grid type enum parsing."""
        for grid_str, grid_enum in [
            ('Triangular', GridType.TRIANGULAR),
            ('Square', GridType.SQUARE),
        ]:
            config_text = f"GridType {grid_str}"
            config = ConfigFile()
            lines = config._tokenize(config_text)
            config._parse_lines(lines, "test.config")
            assert config.mic_grid_type == grid_enum

    def test_parse_detector_spacing(self):
        """Test detector spacing multi-value parsing."""
        config_text = """
DetectorSpacing 0 2.0
DetectorSpacing 1 3.5
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert len(config.detector_spacing) == 2
        assert config.detector_spacing[0] == 2.0
        assert config.detector_spacing[1] == 3.5

    def test_parse_boundary_voxels(self):
        """Test boundary voxel selection parsing."""
        config_text = """
SelectBoundaryVoxels 0.90 5.0 0.5
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        assert config.select_boundary_voxels is not None
        max_cost, angle_rad, radius = config.select_boundary_voxels
        assert max_cost == 0.90
        assert np.isclose(angle_rad, np.deg2rad(5.0))
        assert radius == 0.5


class TestConfigFileComments:
    """Test comment and empty line handling."""

    def test_skip_comments(self):
        """Test that comment lines are skipped."""
        config_text = """
# This is a comment
BeamEnergy 50.0
# Another comment
MaxQ 8.0
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        # Should only have 2 non-comment lines
        assert len(lines) == 2
        assert lines[0][1][0] == 'BeamEnergy'
        assert lines[1][1][0] == 'MaxQ'

    def test_skip_empty_lines(self):
        """Test that empty lines are skipped."""
        config_text = """

BeamEnergy 50.0

MaxQ 8.0

        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        assert len(lines) == 2


class TestConfigFileValidation:
    """Test validation logic."""

    def test_missing_required_parameter(self):
        """Test error on missing required parameter."""
        config_text = """
BeamEnergy 50.0
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)
        config._parse_lines(lines, "test.config")

        with pytest.raises(ValueError, match="Missing required parameters"):
            config._validate()

    def test_invalid_beam_energy(self):
        """Test error on invalid beam energy."""
        config = ConfigFile()
        config.beam_energy = -1.0

        with pytest.raises(ValueError, match="BeamEnergy must be > 0"):
            config._validate_ranges()

    def test_invalid_max_q(self):
        """Test error on invalid max Q."""
        config = ConfigFile()
        config.beam_energy = 50.0  # Valid value for earlier check
        config.max_q = 0.0  # Invalid value being tested

        with pytest.raises(ValueError, match="MaxQ must be > 0"):
            config._validate_ranges()

    def test_invalid_amplitude_fraction(self):
        """Test error on invalid amplitude fraction."""
        config = ConfigFile()
        config.beam_energy = 50.0  # Valid value for earlier check
        config.max_q = 5.0  # Valid value for earlier check
        config.min_amplitude_fraction = 1.5  # Invalid value being tested

        with pytest.raises(ValueError, match="MinAmplitudeFraction must be in"):
            config._validate_ranges()

    def test_detector_spacing_consistency(self):
        """Test detector spacing count validation."""
        config = ConfigFile()
        config.num_detectors = 2
        config.detector_spacing = [1.0, 2.0, 3.0]  # Too many!
        config._initialized.add('DetectorSpacing')

        with pytest.raises(ValueError, match="DetectorSpacing requires 1 entries"):
            config._validate_consistency()

    def test_strain_requires_config_file(self):
        """Test that strain enabled requires config filename."""
        config = ConfigFile()
        config.strain_enabled = True
        config.strain_opt_config_filename = ""

        with pytest.raises(ValueError, match="StrainOptConfigFilename required"):
            config._validate_consistency()

    def test_resolution_range(self):
        """Test resolution range validation."""
        config = ConfigFile()
        config.beam_energy = 50.0  # Valid value for earlier check
        config.max_q = 5.0  # Valid value for earlier check
        config.min_amplitude_fraction = 0.5  # Valid value for earlier check
        config.num_detectors = 2  # Valid value for earlier check
        config.min_local_resolution = 5  # Valid value, but...
        config.max_local_resolution = 3  # ...max < min (invalid, being tested)

        with pytest.raises(ValueError, match="MaxLocalResolution must be >= MinLocalResolution"):
            config._validate_ranges()


class TestConfigFileErrors:
    """Test error handling."""

    def test_unknown_keyword(self):
        """Test error on unknown keyword."""
        config_text = """
InvalidKeyword 123
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        with pytest.raises(ValueError, match="Unknown keyword"):
            config._parse_lines(lines, "test.config")

    def test_missing_value(self):
        """Test error when value is missing."""
        config_text = """
BeamEnergy
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        with pytest.raises(ValueError, match="Missing value"):
            config._parse_lines(lines, "test.config")

    def test_invalid_file_type(self):
        """Test error on invalid file type."""
        config_text = """
InFileType invalid
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        with pytest.raises(ValueError, match="Unknown file type"):
            config._parse_lines(lines, "test.config")

    def test_invalid_symmetry(self):
        """Test error on invalid symmetry."""
        config_text = """
SampleSymmetry Invalid
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        with pytest.raises(ValueError, match="Unknown symmetry"):
            config._parse_lines(lines, "test.config")

    def test_vector3_missing_values(self):
        """Test error when vector3 has insufficient values."""
        config_text = """
BeamDirection 1 0
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        with pytest.raises(ValueError, match="Vector3 requires 3 values"):
            config._parse_lines(lines, "test.config")

    def test_bool_invalid_value(self):
        """Test error on invalid bool value."""
        config_text = """
EnableStrain 2
        """
        config = ConfigFile()
        lines = config._tokenize(config_text)

        with pytest.raises(ValueError, match="must be 0 or 1"):
            config._parse_lines(lines, "test.config")

    def test_file_not_found(self):
        """Test error when file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            ConfigFile.from_file("nonexistent_file.config")


class TestConfigFileIntegration:
    """Integration tests with real config files."""

    def test_parse_minimal_config(self):
        """Test parsing minimal valid configuration."""
        # Create minimal config with all required parameters
        config_text = """
# Minimal valid configuration
InfileBasename Test_
InfileExtension d
InFileType ascii
InfileSerialLength 3
OutfileBasename Out_
OutfileExtension d
OutfileSerialLength 3
OutStructureBasename Struct_
FileNumStart 0
FileNumEnd 0

BeamEnergy 50.0
BeamEnergyWidth 0.05
BeamHeight 0.001
BeamDirection 1 0 0
BeamDeflectionChiLaue 0

DetectorFilename detector.txt
EtaLimit 90

SampleLocation 0 0 0
SampleRadius 1.0
SampleCenter 0 0 0
SampleOrientation 0 0 0
SampleFilename sample.mic
StructureFilename structure.dat
FundamentalZoneFilename fz.dat
SampleSymmetry Cubic
MaxInitSideLength 0.01
MinSideLength 0.01

MinAmplitudeFraction 0.1
MaxQ 5.0

RotationRangeFilename omega.dat

LocalOrientationGridRadius 5
MinLocalResolution 0
MaxLocalResolution 3
MaxDiscreteCandidates 100
MaxAcceptedCost 0.9
MaxConvergenceCost 0.01
MaxDeepeningHitRatio 0.8
MaxMCSteps 1000
MCRadiusScaleFactor 0.5
SuccessiveRestarts 2
MinAccelerationThreshold 0.85
SecondsBetweenSave 3600

OptimizationFilename opt.txt
OptimizationConstrainFilename constrain.txt
DetectionLimitFilename limit.txt
NumParameterOptimizationSteps 10
NumElementToOptimizePerPE 20
ParameterMCInitTemperature 0.0
OrientationSearchMethod ConstrainedEuler
CoolingFraction 0.001
ThermalizeFraction 0.001
ParameterRefinements 100
NumDetectors 2
DetectorSpacing 0 2.0
DetectorSpacingDeviation 0.01
DetectorOrientationDeviationInEuler 0.5 0.5 0.5
DetectorOrientationDeviationInSO3 1.0
ParamMCMaxLocalRestarts 3
ParamMCMaxGlobalRestarts 3
ParamMCNumGlobalSearchElements 10
ConstrainedOptimization 1
SearchVolumeReductionFactor 4

ConsistencyError 0
BraggFilterTolerance 0
        """

        with tempfile.NamedTemporaryFile(mode='w', suffix='.config', delete=False) as f:
            f.write(config_text)
            config_file = f.name

        try:
            config = ConfigFile.from_file(config_file)

            # Verify key parameters
            assert config.beam_energy == 50.0
            assert np.allclose(config.beam_direction, [1, 0, 0])
            assert config.sample_symmetry == SymmetryType.CUBIC
            assert config.num_detectors == 2
            assert len(config.detector_spacing) == 1
            assert config.detector_spacing[0] == 2.0
            assert config.constrained_param_mc is True

        finally:
            Path(config_file).unlink()

    def test_parse_reconstruct_test_config_if_exists(self):
        """Test parsing actual ReconstructTest.config if it exists."""
        # Path from test file: icenine_py/tests/test_config_file.py -> IceNine/ConfigFiles/
        config_path = Path(__file__).parent.parent.parent / "ConfigFiles" / "ReconstructTest.config"

        if not config_path.exists():
            pytest.skip(f"Config file not found: {config_path}")

        config = ConfigFile.from_file(str(config_path))

        # Validate key parameters from ReconstructTest.config
        assert config.beam_energy == pytest.approx(50.02099)
        assert np.allclose(config.beam_direction, [1, 0, 0])
        assert config.sample_symmetry == SymmetryType.CUBIC
        assert config.num_detectors == 2
        assert config.in_file_type == FileType.ASCII
        assert config.max_discrete_candidates == 100
        assert config.orientation_search_method == SO3SearchMethod.CONSTRAINED_EULER


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
