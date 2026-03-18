# ConfigFile Parser Port Plan

**Status**: Planning Phase
**Created**: 2025-11-16
**Estimated Complexity**: MEDIUM
**Estimated Timeline**: 3-4 days

## Executive Summary

The `CConfigFile` class is a comprehensive configuration parser that handles 80+ parameters for IceNine experiments. It parses text-based `.config` files using keyword-value pairs with automatic type conversion and validation.

**Key Features**:
- Text-based format with `keyword value1 value2 ...` syntax
- 80+ parameters across 10 categories
- Type conversion (string, int, float, bool, vector, enum)
- Automatic degree-to-radian conversion
- Required vs optional parameter tracking
- Cross-parameter validation

**No blockers** - can start implementation immediately.

## Architecture Overview

### File Locations

**C++ Source**:
- Header: `Src/ConfigFile.h` (289 lines)
- Implementation: `Src/ConfigFile.cpp` (945 lines)

**Python Target**:
- `icenine_py/icenine/config_file.py` (new file, ~500 lines estimated)
- `icenine_py/tests/test_config_file.py` (new file, ~300 lines estimated)

**Test Data**:
- `ConfigFiles/ReconstructTest.config` - Full reconstruction config
- `ConfigFiles/OneGrain.Config` - Simple single grain config
- `ConfigFiles/StrainOpt.config` - Strain optimization config

### Class Structure

```python
@dataclass
class ConfigFile:
    """Configuration file parser for IceNine experiments.

    Python port of CConfigFile from Src/ConfigFile.h.

    Attributes:
        80+ configuration parameters organized by category
    """

    # File I/O (11 parameters)
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

    # Beam (5 parameters)
    beam_energy: float = 0.0              # keV
    beam_energy_width: float = 0.0        # fraction
    beam_height: float = 0.0              # mm
    beam_direction: np.ndarray = None     # [x, y, z]
    beam_deflection_chi_laue: float = 0.0 # radians

    # Detector (2 parameters)
    detector_filename: str = ""
    eta_limit: float = 0.0  # degrees → radians

    # Sample (9 parameters)
    sample_location: np.ndarray = None       # [x, y, z] meters
    sample_radius: float = 0.0
    sample_center: np.ndarray = None         # [x, y, z]
    sample_orientation: np.ndarray = None    # Euler angles, degrees → radians
    sample_filename: str = ""
    structure_filename: str = ""
    fundamental_zone_filename: str = ""
    sample_symmetry: SymmetryType = SymmetryType.CUBIC
    max_init_side_length: float = 0.0
    min_side_length: float = 0.0

    # Simulation (2 parameters)
    min_amplitude_fraction: float = 0.0  # [0, 1]
    max_q: float = 0.0

    # Initialization Files (1 parameter)
    rotation_range_filename: str = ""  # omega file

    # Search Algorithm (12 parameters)
    local_orientation_grid_radius: float = 0.0  # degrees → radians
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

    # Parameter Optimization (18 parameters)
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
    detector_spacing: List[float] = None  # size = num_detectors - 1
    detector_dist_deviation: float = 0.0
    detector_orient_deviation_euler: np.ndarray = None  # degrees → radians
    detector_orient_deviation_so3: float = 0.0  # degrees → radians
    max_param_mc_local_restarts: int = 0
    max_param_mc_global_restarts: int = 0
    param_mc_global_search_elements: int = 0
    constrained_param_mc: bool = False
    search_vol_reduction_factor: float = 0.0

    # Scattering Vector (2 parameters)
    consistency_error: float = 0.0      # degrees → radians
    bragg_filter_tolerance: float = 0.0 # degrees → radians

    # Pipeline Commands (8 flags)
    run_param_optimization: bool = False
    run_reconstruction: bool = False
    run_adp_reconstruction: bool = False
    select_boundary_voxels: Optional[Tuple[float, float, float]] = None  # (cost, angle, radius)
    intensity_decomposition: bool = False
    lazy_bfs: bool = False
    lazy_strain: bool = False
    local_orientation_optimization: bool = False

    # Strain (2 parameters)
    strain_enabled: bool = False
    strain_opt_config_filename: str = ""

    # Partial Results (2 parameters)
    use_partial_result: bool = False
    partial_result_acceptance_conf: float = 0.0
    partial_result_filename: str = ""

    # Grid Type (1 parameter)
    mic_grid_type: GridType = GridType.TRIANGULAR
```

## Configuration File Format

### Syntax

```
# Comments start with '#'
Keyword value1 value2 ...
Keyword value
```

**Rules**:
1. One keyword per line
2. Case-sensitive keywords
3. Space/tab delimited values
4. Empty lines and comment lines (`#`) ignored
5. No inline comments
6. Vectors: space-separated values (e.g., `BeamDirection 1 0 0`)
7. Angles: input in degrees, stored as radians (automatic conversion)

### Example

```config
# Beam configuration
BeamEnergy                50.02099
BeamEnergyWidth           0.05
BeamDirection             1 0 0
BeamHeight                0.0012

# Sample
SampleFilename            DataFiles/Au1007_small.mic
SampleOrientation         0 0 0
SampleLocation            0 0 0

# Search parameters
MaxDiscreteCandidates     100
LocalOrientationGridRadius 5
MaxAcceptedCost           0.90

# Detector
DetectorFilename          ConfigFiles/DetectorFile.txt
NumDetectors              2
DetectorSpacing           0    2.0

# Pipeline
IntensityDecomposition
```

### Type Parsing

| Type | C++ Conversion | Python Equivalent | Example |
|------|----------------|-------------------|---------|
| String | Direct copy | `str` | `InfileBasename TestInput/data_` |
| Float | `atof()` | `float()` | `BeamEnergy 50.02099` |
| Int | `atoi()` | `int()` | `MaxDiscreteCandidates 100` |
| Bool | `atoi()` as 0/1 | `bool(int())` | `EnableStrain 0` |
| Vector3 | `ExtractVector()` | `np.array([float, float, float])` | `BeamDirection 1 0 0` |
| Enum | String comparison | `Enum['STRING']` | `SampleSymmetry Cubic` |
| Angle | `atof() * DEG_TO_RAD` | `np.deg2rad(float())` | `EtaLimit 81` → 1.414 rad |

### Degree → Radian Conversion

**Automatically converted at parse time** (15 parameters):
- `local_orientation_grid_radius`
- `eta_limit`
- `consistency_error`
- `bragg_filter_tolerance`
- `detector_orient_deviation_euler` (all 3 components)
- `detector_orient_deviation_so3`
- `sample_orientation` (all 3 Euler angles)
- `select_boundary_voxels` (angle component)

**Already in radians** (3 parameters):
- `beam_deflection_chi_laue`

## Parameter Categories

### Complete Parameter List (80+ total)

#### File I/O (11 parameters)
```python
file_start_num: int              # Starting file number (obsolete, default 0)
file_end_num: int                # Ending file number (obsolete, default 0)
in_file_basename: str            # REQUIRED - Input file base path
in_file_ext: str                 # REQUIRED - Input file extension
in_file_type: FileType           # REQUIRED - bin/tif/ascii
in_file_serial_length: int       # REQUIRED - Serial number padding
out_file_basename: str           # REQUIRED - Output file base path
out_file_ext: str                # REQUIRED - Output extension
out_file_serial_length: int      # REQUIRED - Output serial padding
out_structure_basename: str      # REQUIRED - Structure file base
bc_peak_detector_offset: int     # OPTIONAL - Default 0 (IceNine), 1 (xdmmpi)
```

**Config keywords**:
- `InfileBasename`, `InfileExtension`, `InFileType`
- `InfileSerialLength`, `OutfileBasename`, `OutfileExtension`
- `OutfileSerialLength`, `OutStructureBasename`
- `FileNumStart`, `FileNumEnd` (obsolete)
- `BCPeakDetectorOffset`

#### Beam Information (5 parameters)
```python
beam_energy: float               # REQUIRED - Energy in keV
beam_energy_width: float         # REQUIRED - Spread fraction (e.g., 0.05)
beam_height: float               # REQUIRED - Height in mm
beam_direction: np.ndarray       # REQUIRED - Unit vector [x, y, z]
beam_deflection_chi_laue: float  # REQUIRED - Deflection in radians
```

**Config keywords**:
- `BeamEnergy`, `BeamEnergyWidth`, `BeamHeight`
- `BeamDirection`, `BeamDeflectionChiLaue`

#### Detector (2 parameters)
```python
detector_filename: str           # REQUIRED - Detector geometry file
eta_limit: float                 # REQUIRED - Azimuth limit (deg → rad)
```

**Config keywords**:
- `DetectorFilename`, `EtaLimit`

#### Sample (9 parameters)
```python
sample_location: np.ndarray      # REQUIRED - Position [x, y, z] meters
sample_radius: float             # REQUIRED - Radius (mm)
sample_center: np.ndarray        # REQUIRED - Center [x, y, z]
sample_orientation: np.ndarray   # REQUIRED - Euler [phi, theta, psi] (deg → rad)
sample_filename: str             # REQUIRED - .mic file path
structure_filename: str          # REQUIRED - Crystal .dat file
fundamental_zone_filename: str   # REQUIRED - FZ orientation file
sample_symmetry: SymmetryType    # REQUIRED - Cubic/Hexagonal/Tetragonal
max_init_side_length: float      # REQUIRED - Max voxel size
min_side_length: float           # REQUIRED - Min voxel size
```

**Config keywords**:
- `SampleLocation`, `SampleRadius`, `SampleCenter`
- `SampleOrientation`, `SampleFilename`, `StructureFilename`
- `FundamentalZoneFilename`, `SampleSymmetry`
- `MaxInitSideLength`, `MinSideLength`

#### Simulation (2 parameters)
```python
min_amplitude_fraction: float    # REQUIRED - Min peak intensity [0, 1]
max_q: float                     # REQUIRED - Max scattering vector (Å⁻¹)
```

**Config keywords**:
- `MinAmplitudeFraction`, `MaxQ`

#### Initialization Files (1 parameter)
```python
rotation_range_filename: str     # REQUIRED - Omega range file
```

**Config keywords**:
- `RotationRangeFilename`

#### Search Algorithm (12 parameters)
```python
local_orientation_grid_radius: float  # REQUIRED - Search radius (deg → rad)
min_local_resolution: int             # REQUIRED - Min resolution level
max_local_resolution: int             # REQUIRED - Max resolution level
max_discrete_candidates: int          # REQUIRED - Max candidates to keep
max_accepted_cost: float              # REQUIRED - 1 - min confidence
max_convergence_cost: float           # REQUIRED - Early exit threshold
max_deepening_hit_ratio: float        # REQUIRED - Deepening threshold
max_mc_steps: int                     # REQUIRED - Max Monte Carlo steps
mc_radius_scale_factor: float         # REQUIRED - MC step scaling
successive_restarts: int              # REQUIRED - Num restarts
min_acceleration_threshold: float     # REQUIRED - Accel threshold
seconds_between_save: int             # REQUIRED - Autosave interval
```

**Config keywords**:
- `LocalOrientationGridRadius`, `MinLocalResolution`, `MaxLocalResolution`
- `MaxDiscreteCandidates`, `MaxAcceptedCost`, `MaxConvergenceCost`
- `MaxDeepeningHitRatio`, `MaxMCSteps`, `MCRadiusScaleFactor`
- `SuccessiveRestarts`, `MinAccelerationThreshold`, `SecondsBetweenSave`

#### Parameter Optimization (18 parameters)
```python
optimization_filename: str                   # REQUIRED
optimization_constrain_filename: str         # REQUIRED
detection_limit_filename: str                # REQUIRED
num_param_opt_steps: int                     # REQUIRED
num_element_per_pe: int                      # REQUIRED
param_mc_temperature: float                  # REQUIRED
orientation_search_method: SO3SearchMethod   # REQUIRED - ConstrainedEuler/UniformSO3
cooling_fraction: float                      # REQUIRED
thermalize_fraction: float                   # REQUIRED
parameter_refinements: int                   # REQUIRED
num_detectors: int                           # REQUIRED
detector_spacing: List[float]                # REQUIRED (num_detectors - 1 values)
detector_dist_deviation: float               # REQUIRED
detector_orient_deviation_euler: np.ndarray  # REQUIRED (deg → rad)
detector_orient_deviation_so3: float         # REQUIRED (deg → rad)
max_param_mc_local_restarts: int             # REQUIRED
max_param_mc_global_restarts: int            # REQUIRED
param_mc_global_search_elements: int         # REQUIRED
constrained_param_mc: bool                   # REQUIRED
search_vol_reduction_factor: float           # REQUIRED
```

**Config keywords**:
- `OptimizationFilename`, `OptimizationConstrainFilename`, `DetectionLimitFilename`
- `NumParameterOptimizationSteps`, `NumElementToOptimizePerPE`
- `ParameterMCInitTemperature`, `OrientationSearchMethod`
- `CoolingFraction`, `ThermalizeFraction`, `ParameterRefinements`
- `NumDetectors`, `DetectorSpacing` (multiple), `DetectorSpacingDeviation`
- `DetectorOrientationDeviationInEuler`, `DetectorOrientationDeviationInSO3`
- `ParamMCMaxLocalRestarts`, `ParamMCMaxGlobalRestarts`
- `ParamMCNumGlobalSearchElements`, `ConstrainedOptimization`
- `SearchVolumeReductionFactor`

#### Scattering Vector (2 parameters)
```python
consistency_error: float         # REQUIRED - Tolerance (deg → rad)
bragg_filter_tolerance: float    # REQUIRED - Omega tolerance (deg → rad)
```

**Config keywords**:
- `ConsistencyError`, `BraggFilterTolerance`

#### Pipeline Commands (8 flags)
```python
run_param_optimization: bool                        # OPTIONAL
run_reconstruction: bool                            # OPTIONAL
run_adp_reconstruction: bool                        # OPTIONAL
select_boundary_voxels: Optional[Tuple[...]]        # OPTIONAL - (cost, angle, radius)
intensity_decomposition: bool                       # OPTIONAL
lazy_bfs: bool                                      # OPTIONAL
lazy_strain: bool                                   # OPTIONAL
local_orientation_optimization: bool                # OPTIONAL
```

**Config keywords**:
- `RunParameterOptimization`, `RunReconstruction`, `RunAdpReconstruction`
- `SelectBoundaryVoxels` (3 values), `IntensityDecomposition`
- `LazyBFS`, `LazyStrain`, `LocalOrientationOptimization`

#### Strain Optimization (2 parameters)
```python
strain_enabled: bool             # OPTIONAL
strain_opt_config_filename: str  # OPTIONAL
```

**Config keywords**:
- `EnableStrain`, `StrainOptConfigFilename`

#### Partial Results (3 parameters)
```python
use_partial_result: bool         # OPTIONAL
partial_result_acceptance_conf: float  # OPTIONAL
partial_result_filename: str     # OPTIONAL
```

**Config keywords**:
- `PartialResultFilename`, `PartialResultAcceptanceConfidence`

#### Grid Type (1 parameter)
```python
mic_grid_type: GridType          # OPTIONAL - Default: TRIANGULAR
```

**Config keywords**:
- `GridType` (values: `Triangular`, `Square`)

## Implementation Strategy

### Phase 1: Core Data Structures (Day 1)

**Enumerations**:
```python
from enum import Enum, auto

class FileType(Enum):
    """Input/output file type."""
    BIN = 0
    TIF = 1
    ASCII = 2

class SO3SearchMethod(Enum):
    """Orientation search method."""
    CONSTRAINED_EULER = 0
    UNIFORM_QUATERNION = 1

class SymmetryType(Enum):
    """Crystal symmetry type."""
    CUBIC = 0
    HEXAGONAL = 1
    TETRAGONAL = 2
    NONE = 3

class GridType(Enum):
    """Microstructure grid type."""
    TRIANGULAR = 0
    SQUARE = 1
```

**Dataclass with all fields**:
```python
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import numpy as np

@dataclass
class ConfigFile:
    """Configuration file for IceNine experiments.

    C++ Reference: Src/ConfigFile.h, Src/ConfigFile.cpp
    """

    # Define all 80+ fields with defaults
    # Group by category for readability

    # File I/O
    file_start_num: int = 0
    file_end_num: int = 0
    # ... (all fields)

    def __post_init__(self):
        """Initialize numpy arrays that can't be in defaults."""
        if self.beam_direction is None:
            self.beam_direction = np.zeros(3)
        if self.sample_location is None:
            self.sample_location = np.zeros(3)
        # ... (all vector fields)
```

### Phase 2: Parsing Logic (Day 2)

**File reading and tokenization**:
```python
@classmethod
def from_file(cls, filename: str) -> 'ConfigFile':
    """Load configuration from .config file.

    Args:
        filename: Path to .config file

    Returns:
        ConfigFile instance with parsed parameters

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If parsing fails or required parameters missing
    """
    config = cls()

    # Read file
    with open(filename, 'r') as f:
        content = f.read()

    # Tokenize
    lines = config._tokenize(content)

    # Parse
    config._parse_lines(lines)

    # Validate
    config._validate()

    return config

def _tokenize(self, content: str) -> List[List[str]]:
    """Tokenize config file content.

    Rules:
    - Lines starting with '#' are comments (skip)
    - Empty lines are skipped
    - Split by whitespace (space, tab)
    - Returns list of token lists (one per line)
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
```

**Keyword dispatch**:
```python
def _parse_lines(self, lines: List[Tuple[int, List[str]]]):
    """Parse tokenized lines using keyword dispatch."""

    # Track which parameters were set
    self._initialized = set()

    # Keyword parser dictionary
    parsers = {
        'InfileBasename': self._parse_string('in_file_basename'),
        'BeamEnergy': self._parse_float('beam_energy'),
        'BeamDirection': self._parse_vector3('beam_direction'),
        'EtaLimit': self._parse_angle('eta_limit'),  # deg → rad
        'SampleSymmetry': self._parse_symmetry,
        'InFileType': self._parse_file_type,
        'DetectorSpacing': self._parse_detector_spacing,
        'SelectBoundaryVoxels': self._parse_boundary_voxels,
        # ... (80+ keyword parsers)
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
            raise ValueError(f"Error at line {line_num}: {e}")

def _parse_string(self, attr_name: str):
    """Create parser for string parameter."""
    def parser(tokens):
        if len(tokens) < 2:
            raise ValueError(f"Missing value for {tokens[0]}")
        setattr(self, attr_name, tokens[1])
    return parser

def _parse_float(self, attr_name: str):
    """Create parser for float parameter."""
    def parser(tokens):
        if len(tokens) < 2:
            raise ValueError(f"Missing value for {tokens[0]}")
        setattr(self, attr_name, float(tokens[1]))
    return parser

def _parse_int(self, attr_name: str):
    """Create parser for int parameter."""
    def parser(tokens):
        if len(tokens) < 2:
            raise ValueError(f"Missing value for {tokens[0]}")
        setattr(self, attr_name, int(tokens[1]))
    return parser

def _parse_bool(self, attr_name: str):
    """Create parser for bool parameter (0/1)."""
    def parser(tokens):
        if len(tokens) < 2:
            raise ValueError(f"Missing value for {tokens[0]}")
        setattr(self, attr_name, bool(int(tokens[1])))
    return parser

def _parse_angle(self, attr_name: str):
    """Create parser for angle (degrees → radians)."""
    def parser(tokens):
        if len(tokens) < 2:
            raise ValueError(f"Missing value for {tokens[0]}")
        angle_deg = float(tokens[1])
        setattr(self, attr_name, np.deg2rad(angle_deg))
    return parser

def _parse_vector3(self, attr_name: str, convert_angles: bool = False):
    """Create parser for 3D vector."""
    def parser(tokens):
        if len(tokens) < 4:
            raise ValueError(f"Vector3 requires 3 values: {tokens[0]}")
        vec = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])
        if convert_angles:
            vec = np.deg2rad(vec)
        setattr(self, attr_name, vec)
    return parser

def _parse_symmetry(self, tokens):
    """Parse symmetry enum."""
    if len(tokens) < 2:
        raise ValueError("Missing value for SampleSymmetry")

    symmetry_map = {
        'Cubic': SymmetryType.CUBIC,
        'Hexagonal': SymmetryType.HEXAGONAL,
        'Tetragonal': SymmetryType.TETRAGONAL,
        'None': SymmetryType.NONE
    }

    value = tokens[1]
    if value not in symmetry_map:
        raise ValueError(f"Unknown symmetry: {value}")

    self.sample_symmetry = symmetry_map[value]

def _parse_file_type(self, tokens):
    """Parse file type enum."""
    if len(tokens) < 2:
        raise ValueError("Missing value for InFileType")

    type_map = {
        'bin': FileType.BIN,
        'tif': FileType.TIF,
        'ascii': FileType.ASCII
    }

    value = tokens[1].lower()
    if value not in type_map:
        raise ValueError(f"Unknown file type: {value}")

    self.in_file_type = type_map[value]

def _parse_detector_spacing(self, tokens):
    """Parse detector spacing (index value pairs)."""
    if len(tokens) < 3:
        raise ValueError("DetectorSpacing requires index and value")

    if self.detector_spacing is None:
        self.detector_spacing = []

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
```

### Phase 3: Validation (Day 3)

**Required parameter checking**:
```python
def _validate(self):
    """Validate configuration.

    Checks:
    1. All required parameters are set
    2. Cross-parameter consistency
    3. Value ranges are valid

    Raises:
        ValueError: If validation fails
    """
    # Required keywords (from C++ ConfigFile.cpp lines 299-313)
    required = {
        # File I/O
        'InfileBasename', 'InfileExtension', 'InFileType',
        'InfileSerialLength', 'OutfileBasename', 'OutfileExtension',
        'OutfileSerialLength', 'OutStructureBasename',

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

        # Obsolete but still required
        'FileNumStart', 'FileNumEnd',
    }

    missing = required - self._initialized
    if missing:
        raise ValueError(f"Missing required parameters: {', '.join(sorted(missing))}")

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

    # ... more range checks

def _validate_consistency(self):
    """Validate cross-parameter consistency."""
    # Detector spacing must have num_detectors - 1 entries
    if self.detector_spacing is not None:
        expected = self.num_detectors - 1
        if len(self.detector_spacing) != expected:
            raise ValueError(
                f"DetectorSpacing requires {expected} entries, got {len(self.detector_spacing)}"
            )

    # If strain enabled, must have strain config file
    if self.strain_enabled and not self.strain_opt_config_filename:
        raise ValueError("StrainOptConfigFilename required when EnableStrain is 1")

    # ... more consistency checks
```

### Phase 4: Testing (Day 4)

**Test all example configs**:
```python
def test_parse_reconstruct_test_config():
    """Test parsing ReconstructTest.config."""
    config = ConfigFile.from_file('../../ConfigFiles/ReconstructTest.config')

    # Validate key parameters
    assert config.beam_energy == 50.02099
    assert np.allclose(config.beam_direction, [1, 0, 0])
    assert config.sample_symmetry == SymmetryType.CUBIC
    assert config.num_detectors == 2
    assert len(config.detector_spacing) == 1
    assert config.detector_spacing[0] == 2.0

def test_parse_all_parameter_types():
    """Test parsing all data types."""
    # String
    # Float
    # Int
    # Bool
    # Vector3
    # Vector3 with angle conversion
    # Enum
    # Multi-value (detector spacing)
    # Multi-value (boundary voxels)

def test_angle_conversion():
    """Test degree to radian conversion."""
    config = ConfigFile()
    config._parse_lines([(1, ['EtaLimit', '81'])])

    assert np.isclose(config.eta_limit, np.deg2rad(81))

def test_missing_required_parameter():
    """Test error on missing required parameter."""
    with pytest.raises(ValueError, match="Missing required parameters"):
        config = ConfigFile()
        config._validate()

def test_invalid_value_range():
    """Test error on invalid value."""
    config = ConfigFile()
    config.min_amplitude_fraction = 1.5  # > 1

    with pytest.raises(ValueError, match="must be in"):
        config._validate_ranges()

def test_consistency_check():
    """Test cross-parameter validation."""
    config = ConfigFile()
    config.num_detectors = 2
    config.detector_spacing = [1.0, 2.0, 3.0]  # Too many!

    with pytest.raises(ValueError, match="requires 1 entries"):
        config._validate_consistency()
```

## Implementation Checklist

### Phase 1: Core (Day 1)
- [ ] Create enum classes (FileType, SO3SearchMethod, SymmetryType, GridType)
- [ ] Create ConfigFile dataclass with all 80+ fields
- [ ] Add field defaults
- [ ] Add `__post_init__` for numpy array initialization
- [ ] Add docstrings with C++ references

### Phase 2: Parsing (Day 2)
- [ ] Implement `from_file()` class method
- [ ] Implement `_tokenize()` - file reading and tokenization
- [ ] Implement `_parse_lines()` - keyword dispatch
- [ ] Implement type-specific parsers:
  - [ ] `_parse_string()`
  - [ ] `_parse_float()`
  - [ ] `_parse_int()`
  - [ ] `_parse_bool()`
  - [ ] `_parse_angle()` (deg → rad)
  - [ ] `_parse_vector3()` (with optional angle conversion)
  - [ ] `_parse_symmetry()`
  - [ ] `_parse_file_type()`
  - [ ] `_parse_search_method()`
  - [ ] `_parse_grid_type()`
  - [ ] `_parse_detector_spacing()`
  - [ ] `_parse_boundary_voxels()`
- [ ] Create keyword → parser mapping (80+ entries)
- [ ] Add line number tracking for error messages

### Phase 3: Validation (Day 3)
- [ ] Implement `_validate()` - main validation
- [ ] Implement `_validate_ranges()` - value range checking
- [ ] Implement `_validate_consistency()` - cross-parameter checks
- [ ] Define required parameter set (60+ keywords)
- [ ] Add error messages matching C++ format
- [ ] Test with all example config files

### Phase 4: Testing (Day 4)
- [ ] Test parsing ReconstructTest.config
- [ ] Test parsing OneGrain.Config
- [ ] Test parsing StrainOpt.config (if available)
- [ ] Test all data type parsers individually
- [ ] Test angle conversion
- [ ] Test vector parsing
- [ ] Test enum parsing
- [ ] Test multi-value parameters
- [ ] Test error handling for:
  - [ ] Missing required parameters
  - [ ] Invalid keywords
  - [ ] Invalid values
  - [ ] Type conversion errors
  - [ ] Consistency violations
- [ ] Test default values
- [ ] Test optional parameters

### Documentation
- [ ] Add module docstring with overview
- [ ] Add class docstring with example usage
- [ ] Add method docstrings with C++ references
- [ ] Document all 80+ parameters with units and ranges
- [ ] Create example config file with comments
- [ ] Update icenine/__init__.py

## Test Files

### Test Config Examples

**Minimal config** (for testing):
```config
# Minimal valid configuration
InfileBasename Test_
InfileExtension d
InFileType ascii
InfileSerialLength 3
OutfileBasename Out_
OutfileExtension d
OutfileSerialLength 3
OutStructureBasename Struct_

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

FileNumStart 0
FileNumEnd 0
```

## Usage Examples

### Basic Usage
```python
from icenine.config_file import ConfigFile

# Load configuration
config = ConfigFile.from_file('ReconstructTest.config')

# Access parameters
print(f"Beam energy: {config.beam_energy} keV")
print(f"Sample file: {config.sample_filename}")
print(f"Num detectors: {config.num_detectors}")

# Use with ExperimentSetup
from icenine.experiment_setup import XDMExperimentSetup

setup = XDMExperimentSetup(config)
setup.initialize_experiment()
```

### Error Handling
```python
try:
    config = ConfigFile.from_file('invalid.config')
except FileNotFoundError:
    print("Config file not found")
except ValueError as e:
    print(f"Configuration error: {e}")
```

## Success Criteria

1. **Parse all example configs**: ReconstructTest.config, OneGrain.Config, StrainOpt.config
2. **Type conversion**: All data types parsed correctly (string, int, float, bool, vector, enum)
3. **Angle conversion**: All angle parameters converted from degrees to radians
4. **Validation**: Required parameters checked, ranges validated, consistency enforced
5. **Error messages**: Clear, helpful error messages with line numbers
6. **Compatibility**: Values match C++ parsing exactly (within float precision)

## Timeline

| Day | Tasks | Hours |
|-----|-------|-------|
| 1 | Enums, dataclass, field definitions | 6-8 |
| 2 | Parsing logic, keyword dispatch, type parsers | 6-8 |
| 3 | Validation, error handling | 4-6 |
| 4 | Testing, documentation | 4-6 |
| **Total** | **Full implementation** | **20-28 hours** |

## Risks and Mitigations

| Risk | Severity | Mitigation |
|------|----------|------------|
| Keyword name differences | Low | Cross-reference C++ array carefully |
| Angle conversion errors | Medium | Comprehensive tests for all angle parameters |
| Validation logic complexity | Medium | Break into small, testable functions |
| Missing required parameters | Low | Use C++ initialization check as reference |
| Type conversion edge cases | Low | Extensive type conversion tests |

## Next Steps After Completion

Once ConfigFile is complete:
1. **Port file I/O utilities** (InitFileIO namespace) - 3-4 days
2. **Port ExperimentSetup** - 4-5 days (now unblocked by both Sample and ConfigFile)

---

**References**:
- C++ Source: `Src/ConfigFile.h`, `Src/ConfigFile.cpp`
- Test Configs: `ConfigFiles/ReconstructTest.config`, `ConfigFiles/OneGrain.Config`
- C++ Usage: `Src/Driver.cpp`, `Src/ExperimentSetup.cpp`
