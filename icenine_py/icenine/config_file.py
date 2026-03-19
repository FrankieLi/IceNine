"""
Configuration file parser for IceNine experiments.

Python port of CConfigFile from Src/ConfigFile.h and Src/ConfigFile.cpp.

The ConfigFile class parses text-based .config files with keyword-value pairs,
handling 80+ parameters across multiple categories (beam, sample, detector,
search, optimization, etc.).

File format:
    # Comments start with #
    Keyword value1 value2 ...
    Keyword value

Example:
    >>> config = ConfigFile.from_file('ReconstructTest.config')
    >>> print(config.beam_energy)
    50.02099
    >>> print(config.num_detectors)
    2

Author: S. F. Li
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Tuple, Set, Callable, Any, Dict, Type, ClassVar
from functools import partial
import numpy as np


# ============================================================================
# Enumerations
# ============================================================================

class FileType(Enum):
    """Input/output file type.

    C++ Reference: ConfigFile.h EFileType
    """
    BIN = 0
    TIF = 1
    ASCII = 2


class SO3SearchMethod(Enum):
    """Orientation search method.

    C++ Reference: ConfigFile.h SO3SearchMethod
    """
    CONSTRAINED_EULER = 0
    UNIFORM_QUATERNION = 1


class SymmetryType(Enum):
    """Crystal symmetry type.

    C++ Reference: Symmetry.h ESymmetryT
    """
    CUBIC = 0
    HEXAGONAL = 1
    TETRAGONAL = 2
    NONE = 3


class GridType(Enum):
    """Microstructure grid type.

    C++ Reference: MicIO.h GridType
    """
    TRIANGULAR = 0
    SQUARE = 1


# ============================================================================
# ConfigFile Class
# ============================================================================

@dataclass
class ConfigFile:
    """
    Configuration file parser for IceNine experiments.

    Python port of C++ CConfigFile class (Src/ConfigFile.h, Src/ConfigFile.cpp).

    Parses text-based .config files with 80+ parameters organized into categories:
    - File I/O (11 parameters)
    - Beam (5 parameters)
    - Detector (2 parameters)
    - Sample (9 parameters)
    - Simulation (2 parameters)
    - Search Algorithm (12 parameters)
    - Parameter Optimization (18 parameters)
    - Scattering Vector (2 parameters)
    - Pipeline Commands (8 flags)
    - Strain (2 parameters)
    - Partial Results (3 parameters)
    - Grid Type (1 parameter)

    Attributes are named using snake_case (Python convention) instead of C++ CamelCase.

    Example:
        >>> config = ConfigFile.from_file('ReconstructTest.config')
        >>> print(f"Beam energy: {config.beam_energy} keV")
        >>> print(f"Sample file: {config.sample_filename}")
    """

    # ========================================================================
    # Enum mapping constants (for parser efficiency)
    # ========================================================================
    _SYMMETRY_MAP: ClassVar[Dict[str, SymmetryType]] = {
        'Cubic': SymmetryType.CUBIC,
        'Hexagonal': SymmetryType.HEXAGONAL,
        'Tetragonal': SymmetryType.TETRAGONAL,
        'None': SymmetryType.NONE
    }

    _FILE_TYPE_MAP: ClassVar[Dict[str, FileType]] = {
        'bin': FileType.BIN,
        'tif': FileType.TIF,
        'ascii': FileType.ASCII
    }

    _SEARCH_METHOD_MAP: ClassVar[Dict[str, SO3SearchMethod]] = {
        'ConstrainedEuler': SO3SearchMethod.CONSTRAINED_EULER,
        'UniformSO3': SO3SearchMethod.UNIFORM_QUATERNION,
        'UniformQuaternion': SO3SearchMethod.UNIFORM_QUATERNION,
    }

    _GRID_TYPE_MAP: ClassVar[Dict[str, GridType]] = {
        'Triangular': GridType.TRIANGULAR,
        'Square': GridType.SQUARE
    }

    # ========================================================================
    # File I/O (11 parameters)
    # ========================================================================
    file_start_num: int = 0
    file_end_num: int = 0
    in_file_basename: str = ""
    in_file_ext: str = ""
    in_file_type: FileType = FileType.ASCII
    in_file_serial_length: int = 0
    out_file_basename: str = ""
    out_file_ext: str = ""
    out_file_serial_length: int = 0
    out_structure_basename: str = ""
    bc_peak_detector_offset: int = 0

    # ========================================================================
    # Beam Information (5 parameters)
    # ========================================================================
    beam_energy: float = 0.0              # keV
    beam_energy_width: float = 0.0        # fraction (e.g., 0.05)
    beam_height: float = 0.0              # mm
    beam_direction: np.ndarray = None     # [x, y, z] unit vector
    beam_deflection_chi_laue: float = 0.0 # radians

    # ========================================================================
    # Detector (2 parameters)
    # ========================================================================
    detector_filename: str = ""
    eta_limit: float = 0.0  # radians (input as degrees)

    # ========================================================================
    # Sample (9 parameters)
    # ========================================================================
    sample_location: np.ndarray = None       # [x, y, z] meters
    sample_radius: float = 0.0               # mm
    sample_center: np.ndarray = None         # [x, y, z]
    sample_orientation: np.ndarray = None    # Euler angles in radians (input as degrees)
    sample_filename: str = ""
    structure_filename: str = ""
    fundamental_zone_filename: str = ""
    sample_symmetry: SymmetryType = SymmetryType.CUBIC
    max_init_side_length: float = 0.0
    min_side_length: float = 0.0

    # ========================================================================
    # Simulation (2 parameters)
    # ========================================================================
    min_amplitude_fraction: float = 0.0  # [0, 1]
    max_q: float = 0.0                   # Å⁻¹

    # ========================================================================
    # Initialization Files (1 parameter)
    # ========================================================================
    rotation_range_filename: str = ""  # omega file

    # ========================================================================
    # Search Algorithm (12 parameters)
    # ========================================================================
    local_orientation_grid_radius: float = 0.0  # radians (input as degrees)
    min_local_resolution: int = 0
    max_local_resolution: int = 0
    max_discrete_candidates: int = 0
    max_accepted_cost: float = 0.0
    max_convergence_cost: float = 0.0
    max_deepening_hit_ratio: float = 0.0
    max_mc_steps: int = 0
    mc_radius_scale_factor: float = 0.0
    successive_restarts: int = 0
    min_acceleration_threshold: float = 0.0
    seconds_between_save: int = 0

    # ========================================================================
    # Parameter Optimization (18 parameters)
    # ========================================================================
    optimization_filename: str = ""
    optimization_constrain_filename: str = ""
    detection_limit_filename: str = ""
    num_param_opt_steps: int = 0
    num_element_per_pe: int = 0
    param_mc_temperature: float = 0.0
    orientation_search_method: SO3SearchMethod = SO3SearchMethod.CONSTRAINED_EULER
    cooling_fraction: float = 0.0
    thermalize_fraction: float = 0.0
    parameter_refinements: int = 0
    num_detectors: int = -1
    detector_spacing: List[float] = None
    detector_dist_deviation: float = 0.0
    detector_orient_deviation_euler: np.ndarray = None  # radians (input as degrees)
    detector_orient_deviation_so3: float = 0.0  # radians (input as degrees)
    max_param_mc_local_restarts: int = 0
    max_param_mc_global_restarts: int = 0
    param_mc_global_search_elements: int = 0
    constrained_param_mc: bool = False
    search_vol_reduction_factor: float = 0.0

    # ========================================================================
    # Scattering Vector (2 parameters)
    # ========================================================================
    consistency_error: float = 0.0       # radians (input as degrees)
    bragg_filter_tolerance: float = 0.0  # radians (input as degrees)

    # ========================================================================
    # Pipeline Commands (8 flags)
    # ========================================================================
    run_param_optimization: bool = False
    run_reconstruction: bool = False
    run_adp_reconstruction: bool = False
    select_boundary_voxels: Optional[Tuple[float, float, float]] = None  # (cost, angle_rad, radius)
    intensity_decomposition: bool = False
    lazy_bfs: bool = False
    lazy_strain: bool = False
    local_orientation_optimization: bool = False

    # ========================================================================
    # Strain Optimization (2 parameters)
    # ========================================================================
    strain_enabled: bool = False
    strain_opt_config_filename: str = ""

    # ========================================================================
    # Partial Results (3 parameters)
    # ========================================================================
    use_partial_result: bool = False
    partial_result_acceptance_conf: float = 0.0
    partial_result_filename: str = ""

    # ========================================================================
    # Grid Type (1 parameter)
    # ========================================================================
    mic_grid_type: GridType = GridType.TRIANGULAR

    # ========================================================================
    # Internal tracking
    # ========================================================================
    _initialized: Set[str] = field(default_factory=set, init=False, repr=False)

    def __post_init__(self):
        """Initialize numpy arrays that can't be in field defaults."""
        if self.beam_direction is None:
            self.beam_direction = np.zeros(3, dtype=np.float32)
        if self.sample_location is None:
            self.sample_location = np.zeros(3, dtype=np.float32)
        if self.sample_center is None:
            self.sample_center = np.zeros(3, dtype=np.float32)
        if self.sample_orientation is None:
            self.sample_orientation = np.zeros(3, dtype=np.float32)
        if self.detector_orient_deviation_euler is None:
            self.detector_orient_deviation_euler = np.zeros(3, dtype=np.float32)
        if self.detector_spacing is None:
            self.detector_spacing = []

    # ========================================================================
    # Main parsing interface
    # ========================================================================

    @classmethod
    def from_file(cls, filename: str) -> 'ConfigFile':
        """
        Load configuration from .config file.

        C++ Reference: ConfigFile.cpp InputConfigParameters()

        Args:
            filename: Path to .config file

        Returns:
            ConfigFile instance with parsed parameters

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If parsing fails or required parameters missing

        Example:
            >>> config = ConfigFile.from_file('ReconstructTest.config')
            >>> print(config.beam_energy)
            50.02099
        """
        config = cls()

        # Read file
        try:
            with open(filename, 'r') as f:
                content = f.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"Config file not found: {filename}")

        # Tokenize
        lines = config._tokenize(content)

        # Parse
        config._parse_lines(lines, filename)

        # Validate
        config._validate()

        return config

    def _tokenize(self, content: str) -> List[Tuple[int, List[str]]]:
        """
        Tokenize config file content.

        Rules:
        - Lines starting with '#' are comments (skip)
        - Empty lines are skipped
        - Split by whitespace (space, tab)

        Returns:
            List of (line_number, tokens) tuples
        """
        lines = []
        for line_num, line in enumerate(content.split('\n'), 1):
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Split into tokens
            tokens = line.split()
            if tokens:
                lines.append((line_num, tokens))

        return lines

    def _parse_lines(self, lines: List[Tuple[int, List[str]]], filename: str):
        """
        Parse tokenized lines using keyword dispatch.

        C++ Reference: ConfigFile.cpp Parse()
        """
        # Keyword parser dictionary
        # Maps config file keywords to parsing functions
        parsers = {
            # File I/O
            'InfileBasename': self._parse_string('in_file_basename'),
            'InfileExtension': self._parse_string('in_file_ext'),
            'InFileType': self._parse_file_type,
            'InfileSerialLength': self._parse_int('in_file_serial_length'),
            'OutfileBasename': self._parse_string('out_file_basename'),
            'OutfileExtension': self._parse_string('out_file_ext'),
            'OutfileSerialLength': self._parse_int('out_file_serial_length'),
            'OutStructureBasename': self._parse_string('out_structure_basename'),
            'FileNumStart': self._parse_int('file_start_num'),
            'FileNumEnd': self._parse_int('file_end_num'),
            'BCPeakDetectorOffset': self._parse_int('bc_peak_detector_offset'),

            # Beam
            'BeamEnergy': self._parse_float('beam_energy'),
            'BeamEnergyWidth': self._parse_float('beam_energy_width'),
            'BeamHeight': self._parse_float('beam_height'),
            'BeamDirection': self._parse_vector3('beam_direction'),
            'BeamDeflectionChiLaue': self._parse_float('beam_deflection_chi_laue'),

            # Detector
            'DetectorFilename': self._parse_string('detector_filename'),
            'EtaLimit': self._parse_angle('eta_limit'),

            # Sample
            'SampleLocation': self._parse_vector3('sample_location'),
            'SampleRadius': self._parse_float('sample_radius'),
            'SampleCenter': self._parse_vector3('sample_center'),
            'SampleOrientation': self._parse_vector3('sample_orientation', convert_angles=True),
            'SampleFilename': self._parse_string('sample_filename'),
            'StructureFilename': self._parse_string('structure_filename'),
            'FundamentalZoneFilename': self._parse_string('fundamental_zone_filename'),
            'SampleSymmetry': self._parse_symmetry,
            'MaxInitSideLength': self._parse_float('max_init_side_length'),
            'MinSideLength': self._parse_float('min_side_length'),

            # Simulation
            'MinAmplitudeFraction': self._parse_float('min_amplitude_fraction'),
            'MaxQ': self._parse_float('max_q'),

            # Initialization Files
            'RotationRangeFilename': self._parse_string('rotation_range_filename'),

            # Search Algorithm
            'LocalOrientationGridRadius': self._parse_angle('local_orientation_grid_radius'),
            'MinLocalResolution': self._parse_int('min_local_resolution'),
            'MaxLocalResolution': self._parse_int('max_local_resolution'),
            'MaxDiscreteCandidates': self._parse_int('max_discrete_candidates'),
            'MaxAcceptedCost': self._parse_float('max_accepted_cost'),
            'MaxConvergenceCost': self._parse_float('max_convergence_cost'),
            'MaxDeepeningHitRatio': self._parse_float('max_deepening_hit_ratio'),
            'MaxMCSteps': self._parse_int('max_mc_steps'),
            'MCRadiusScaleFactor': self._parse_float('mc_radius_scale_factor'),
            'SuccessiveRestarts': self._parse_int('successive_restarts'),
            'MinAccelerationThreshold': self._parse_float('min_acceleration_threshold'),
            'SecondsBetweenSave': self._parse_int('seconds_between_save'),

            # Parameter Optimization
            'OptimizationFilename': self._parse_string('optimization_filename'),
            'OptimizationConstrainFilename': self._parse_string('optimization_constrain_filename'),
            'DetectionLimitFilename': self._parse_string('detection_limit_filename'),
            'NumParameterOptimizationSteps': self._parse_int('num_param_opt_steps'),
            'NumElementToOptimizePerPE': self._parse_int('num_element_per_pe'),
            'ParameterMCInitTemperature': self._parse_float('param_mc_temperature'),
            'OrientationSearchMethod': self._parse_search_method,
            'CoolingFraction': self._parse_float('cooling_fraction'),
            'ThermalizeFraction': self._parse_float('thermalize_fraction'),
            'ParameterRefinements': self._parse_int('parameter_refinements'),
            'NumDetectors': self._parse_int('num_detectors'),
            'DetectorSpacing': self._parse_detector_spacing,
            'DetectorSpacingDeviation': self._parse_float('detector_dist_deviation'),
            'DetectorOrientationDeviationInEuler': self._parse_vector3('detector_orient_deviation_euler', convert_angles=True),
            'DetectorOrientationDeviationInSO3': self._parse_angle('detector_orient_deviation_so3'),
            'ParamMCMaxLocalRestarts': self._parse_int('max_param_mc_local_restarts'),
            'ParamMCMaxGlobalRestarts': self._parse_int('max_param_mc_global_restarts'),
            'ParamMCNumGlobalSearchElements': self._parse_int('param_mc_global_search_elements'),
            'ConstrainedOptimization': self._parse_bool('constrained_param_mc'),
            'SearchVolumeReductionFactor': self._parse_float('search_vol_reduction_factor'),

            # Scattering Vector
            'ConsistencyError': self._parse_angle('consistency_error'),
            'BraggFilterTolerance': self._parse_angle('bragg_filter_tolerance'),

            # Pipeline Commands
            'RunParameterOptimization': self._parse_flag('run_param_optimization'),
            'RunReconstruction': self._parse_flag('run_reconstruction'),
            'RunAdpReconstruction': self._parse_flag('run_adp_reconstruction'),
            'SelectBoundaryVoxels': self._parse_boundary_voxels,
            'IntensityDecomposition': self._parse_flag('intensity_decomposition'),
            'LazyBFS': self._parse_flag('lazy_bfs'),
            'LazyStrain': self._parse_flag('lazy_strain'),
            'LocalOrientationOptimization': self._parse_flag('local_orientation_optimization'),

            # Strain
            'EnableStrain': self._parse_bool('strain_enabled'),
            'StrainOptConfigFilename': self._parse_string('strain_opt_config_filename'),

            # Partial Results
            'PartialResultFilename': self._parse_string('partial_result_filename'),
            'PartialResultAcceptanceConfidence': self._parse_float('partial_result_acceptance_conf'),

            # Grid Type
            'GridType': self._parse_grid_type,
        }

        for line_num, tokens in lines:
            keyword = tokens[0]

            try:
                if keyword in parsers:
                    parsers[keyword](tokens)
                    self._initialized.add(keyword)
                else:
                    raise ValueError(f"Unknown keyword: {keyword}")
            except Exception as e:
                raise ValueError(f"Error in {filename} at line {line_num}: {e}")

        # Handle derived flags
        if self.partial_result_filename:
            self.use_partial_result = True

    # ========================================================================
    # Generic parser factory
    # ========================================================================

    def _make_parser(self, attr_name: str, converter: Callable = None,
                     min_tokens: int = 2, num_values: int = 1,
                     validator: Callable = None) -> Callable:
        """
        Generic parser factory that eliminates repetitive parser code.

        Args:
            attr_name: Attribute name to set
            converter: Function to convert string to target type (default: str)
            min_tokens: Minimum number of tokens required
            num_values: Number of values to parse (1 = scalar, 3 = vector3, etc.)
            validator: Optional validation function for converted value

        Returns:
            Parser function that takes token list
        """
        def parser(tokens):
            if len(tokens) < min_tokens:
                raise ValueError(f"Missing value for {tokens[0]}")

            if num_values == 1:
                # Scalar value
                value = tokens[1] if converter is None else converter(tokens[1])
                if validator:
                    validator(value, tokens[0])
                setattr(self, attr_name, value)
            else:
                # Vector/multiple values
                values = [converter(tokens[i+1]) for i in range(num_values)]
                value = np.array(values, dtype=np.float32) if converter == float else values
                setattr(self, attr_name, value)

        return parser

    # ========================================================================
    # Specialized parser factories using generic _make_parser
    # ========================================================================

    def _parse_string(self, attr_name: str):
        """Create parser for string parameter."""
        return self._make_parser(attr_name, converter=None)

    def _parse_float(self, attr_name: str):
        """Create parser for float parameter."""
        return self._make_parser(attr_name, converter=float)

    def _parse_int(self, attr_name: str):
        """Create parser for int parameter."""
        return self._make_parser(attr_name, converter=int)

    def _parse_bool(self, attr_name: str):
        """Create parser for bool parameter (0/1)."""
        def parser(tokens):
            if len(tokens) < 2:
                raise ValueError(f"Missing value for {tokens[0]}")
            value = int(tokens[1])
            if value not in (0, 1):
                raise ValueError(f"{tokens[0]} must be 0 or 1, got {value}")
            setattr(self, attr_name, bool(value))
        return parser

    def _parse_flag(self, attr_name: str):
        """Create parser for boolean flag (presence = True)."""
        def parser(tokens):
            setattr(self, attr_name, True)
        return parser

    def _parse_angle(self, attr_name: str):
        """Create parser for angle (degrees → radians)."""
        return self._make_parser(attr_name, converter=lambda x: np.deg2rad(float(x)))

    def _parse_vector3(self, attr_name: str, convert_angles: bool = False):
        """Create parser for 3D vector."""
        def parser(tokens):
            if len(tokens) < 4:
                raise ValueError(f"Vector3 requires 3 values for {tokens[0]}")
            converter = lambda x: np.deg2rad(float(x)) if convert_angles else float(x)
            vec = np.array([converter(tokens[1]), converter(tokens[2]), converter(tokens[3])],
                          dtype=np.float32)
            setattr(self, attr_name, vec)
        return parser

    def _make_enum_parser(self, attr_name: str, string_to_enum: Dict[str, Enum],
                         case_sensitive: bool = True, error_prefix: str = None):
        """
        Generic enum parser factory.

        Args:
            attr_name: Attribute name to set
            string_to_enum: Dictionary mapping string values to enum values
            case_sensitive: Whether string matching is case-sensitive
            error_prefix: Custom prefix for error messages (e.g., "symmetry", "file type")

        Returns:
            Parser function for the enum
        """
        def parser(tokens):
            if len(tokens) < 2:
                raise ValueError(f"Missing value for {tokens[0]}")

            value = tokens[1] if case_sensitive else tokens[1].lower()

            if value not in string_to_enum:
                valid_values = ', '.join(string_to_enum.keys())
                # Use custom error prefix if provided, otherwise use keyword name
                prefix = error_prefix if error_prefix else tokens[0]
                raise ValueError(f"Unknown {prefix}: {tokens[1]}. Must be one of: {valid_values}")

            setattr(self, attr_name, string_to_enum[value])

        return parser

    # ========================================================================
    # Enum-specific parsers using generic _make_enum_parser
    # ========================================================================

    def _parse_symmetry(self, tokens):
        """Parse symmetry enum."""
        parser = self._make_enum_parser('sample_symmetry', self._SYMMETRY_MAP, error_prefix='symmetry')
        parser(tokens)

    def _parse_file_type(self, tokens):
        """Parse file type enum."""
        parser = self._make_enum_parser('in_file_type', self._FILE_TYPE_MAP,
                                       case_sensitive=False, error_prefix='file type')
        parser(tokens)

    def _parse_search_method(self, tokens):
        """Parse search method enum."""
        parser = self._make_enum_parser('orientation_search_method', self._SEARCH_METHOD_MAP,
                                       error_prefix='search method')
        parser(tokens)

    def _parse_grid_type(self, tokens):
        """Parse grid type enum."""
        parser = self._make_enum_parser('mic_grid_type', self._GRID_TYPE_MAP,
                                       error_prefix='grid type')
        parser(tokens)

    def _parse_detector_spacing(self, tokens):
        """Parse detector spacing (index value pairs)."""
        if len(tokens) < 3:
            raise ValueError("DetectorSpacing requires index and value")

        index = int(tokens[1])
        value = float(tokens[2])

        # Ensure list is large enough
        while len(self.detector_spacing) <= index:
            self.detector_spacing.append(0.0)

        self.detector_spacing[index] = value

    def _parse_boundary_voxels(self, tokens):
        """Parse boundary voxel selection parameters."""
        if len(tokens) < 4:
            raise ValueError("SelectBoundaryVoxels requires MaxCost Angle Radius")

        max_cost = float(tokens[1])
        angle = np.deg2rad(float(tokens[2]))  # degrees → radians
        radius = float(tokens[3])

        self.select_boundary_voxels = (max_cost, angle, radius)

    # ========================================================================
    # Validation
    # ========================================================================

    def _validate(self):
        """
        Validate configuration.

        Checks:
        1. All required parameters are set
        2. Cross-parameter consistency
        3. Value ranges are valid

        C++ Reference: ConfigFile.cpp Parse() (lines 299-313 for required params)

        Raises:
            ValueError: If validation fails
        """
        # Required keywords (from C++ ConfigFile.cpp lines 299-313)
        required = {
            # File I/O
            'InfileBasename', 'InfileExtension', 'InFileType',
            'InfileSerialLength', 'OutfileBasename', 'OutfileExtension',
            'OutfileSerialLength', 'OutStructureBasename',
            'FileNumStart', 'FileNumEnd',  # Obsolete but still required

            # Beam
            'BeamEnergy', 'BeamEnergyWidth', 'BeamHeight',
            'BeamDirection', 'BeamDeflectionChiLaue',

            # Detector
            'DetectorFilename', 'EtaLimit',

            # Sample
            'SampleLocation', 'SampleRadius', 'SampleCenter',
            'SampleOrientation', 'SampleFilename', 'StructureFilename',
            'FundamentalZoneFilename', 'SampleSymmetry',
            'MaxInitSideLength', 'MinSideLength',

            # Simulation
            'MinAmplitudeFraction', 'MaxQ',

            # Initialization
            'RotationRangeFilename',

            # Search
            'LocalOrientationGridRadius', 'MinLocalResolution',
            'MaxLocalResolution', 'MaxDiscreteCandidates',
            'MaxAcceptedCost', 'MaxConvergenceCost',
            'MaxDeepeningHitRatio', 'MaxMCSteps',
            'MCRadiusScaleFactor', 'SuccessiveRestarts',
            'MinAccelerationThreshold', 'SecondsBetweenSave',

            # Parameter Optimization
            'OptimizationFilename', 'OptimizationConstrainFilename',
            'DetectionLimitFilename', 'NumParameterOptimizationSteps',
            'NumElementToOptimizePerPE', 'ParameterMCInitTemperature',
            'OrientationSearchMethod', 'CoolingFraction',
            'ThermalizeFraction', 'ParameterRefinements',
            'NumDetectors', 'DetectorSpacingDeviation',
            'DetectorOrientationDeviationInEuler', 'DetectorOrientationDeviationInSO3',
            'ParamMCMaxLocalRestarts', 'ParamMCMaxGlobalRestarts',
            'ParamMCNumGlobalSearchElements', 'ConstrainedOptimization',
            'SearchVolumeReductionFactor',

            # Scattering Vector
            'ConsistencyError', 'BraggFilterTolerance',
        }

        missing = required - self._initialized
        if missing:
            raise ValueError(f"Missing required parameters:\n  " + '\n  '.join(sorted(missing)))

        # Cross-parameter validation
        self._validate_ranges()
        self._validate_consistency()

    def _validate_ranges(self):
        """Validate parameter value ranges."""
        if self.beam_energy <= 0:
            raise ValueError("BeamEnergy must be > 0")

        if self.max_q <= 0:
            raise ValueError("MaxQ must be > 0")

        if not (0 <= self.min_amplitude_fraction <= 1):
            raise ValueError("MinAmplitudeFraction must be in [0, 1]")

        if self.num_detectors <= 0:
            raise ValueError("NumDetectors must be > 0")

        if self.min_local_resolution < 0:
            raise ValueError("MinLocalResolution must be >= 0")

        if self.max_local_resolution < self.min_local_resolution:
            raise ValueError("MaxLocalResolution must be >= MinLocalResolution")

    def _validate_consistency(self):
        """Validate cross-parameter consistency."""
        # Detector spacing must have num_detectors - 1 entries
        if 'DetectorSpacing' in self._initialized:
            expected = self.num_detectors - 1
            if len(self.detector_spacing) != expected:
                raise ValueError(
                    f"DetectorSpacing requires {expected} entries for {self.num_detectors} detectors, "
                    f"got {len(self.detector_spacing)}"
                )

        # If strain enabled, must have strain config file
        if self.strain_enabled and not self.strain_opt_config_filename:
            raise ValueError("StrainOptConfigFilename required when EnableStrain is 1")
