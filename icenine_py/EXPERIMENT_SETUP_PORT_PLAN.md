# ExperimentSetup Python Port Plan

**Status**: Planning Phase
**Created**: 2025-11-16
**Estimated Complexity**: HIGH (requires porting CSample first)

## Executive Summary

`CXDMExperimentSetup` is the central orchestration class that initializes all experimental parameters for IceNine simulations. It parses configuration files, reads external data files (detectors, omega ranges, crystal structures), and initializes sample geometry.

**Critical Blocker**: Requires porting `CSample` class first - this is the largest missing dependency.

## Architecture Overview

### Class Hierarchy

```
CExperimentSetup (Abstract Base)
├── Configuration file container (CConfigFile)
├── Beam parameters (energy, direction, limits)
└── Intensity thresholds

CXDMExperimentSetup (Concrete Implementation)
├── Detector list (vector<CDetector>)
├── Omega range system (CSimulationRange)
├── Optimization parameters (step sizes, constraints)
└── Sample initialization logic
```

### File Locations

**C++ Source**:
- `Src/ExperimentSetup.h` (298 lines) - Header with class declarations
- `Src/ExperimentSetup.cpp` (611 lines) - Implementation

**Python Target**:
- `icenine_py/icenine/experiment_setup.py` (new file)
- `icenine_py/tests/test_experiment_setup.py` (new file)

## Dependencies Analysis

### Already Ported ✅

| Component | Python Module | Status | Usage in ExperimentSetup |
|-----------|---------------|--------|--------------------------|
| CDetector | detector.py | Complete | vDetectorList storage, GetMaxQ() |
| CSimulationRange | simulation_range.py | Complete | oRangeToIndexMap for angle mapping |
| SRange | simulation_range.OmegaRange | Complete | vOmegaRangeList storage |
| SIntRange | simulation_range.FileRange | Complete | vFileRangeList storage |
| CUnitCell | crystal_structure.py | Complete | InitializeSample() crystal setup |
| SVector3 | geometry.py | Complete | Beam direction, positions |
| Symmetry | symmetry.py | Complete | GetSampleSymmetry() |
| DiffractionCore | diffraction_core.py | Complete | Reciprocal vector calculations |

### NOT Yet Ported ❌

| Component | File | Priority | Estimated LOC | Notes |
|-----------|------|----------|---------------|-------|
| **CSample** | Sample.h/cpp | **CRITICAL** | ~800 | Sample geometry, .mic file, crystal list |
| CConfigFile | ConfigFile.h | High | ~300 | Configuration file parser |
| InitFileIO | InitFilesIO.h/cpp | High | ~400 | File I/O utilities (omega, detector, structure) |
| CDetectorFile | DetectorFile.h | Medium | ~200 | Detector geometry file parser |
| SStepSizeInfo | - | Low | ~30 | Optimization parameter structure |

## Porting Strategy

### Phase 0: Prerequisites (MUST COMPLETE FIRST)

**Before starting ExperimentSetup, we must port:**

#### 1. CSample Class (CRITICAL BLOCKER)

**Location**: `Src/Sample.h`, `Src/Sample.cpp`

**Key Responsibilities**:
- Sample location and orientation in lab frame
- Microstructure voxel grid (via MicIO - already have mic_file.py)
- Crystal structure list management
- Sample-to-lab coordinate transformations
- Voxel neighbor queries

**Required for**:
- `InitializeSample()` - Cannot initialize without CSample
- All simulation workflows - Sample is central to reconstruction

**Estimated effort**: Large (separate task - recommend dedicated planning session)

#### 2. ConfigFile Parser

**Location**: `Src/ConfigFile.h`

**Key Responsibilities**:
- Parse .config files (key-value text format)
- Store ~50+ experimental parameters
- Type conversions (string → float, int, bool, vector)

**Implementation approach**:
```python
@dataclass
class ConfigFile:
    """Configuration file container for IceNine experiments."""

    # Beam parameters
    beam_energy: float
    beam_energy_width: float
    beam_direction: np.ndarray  # [x, y, z]
    beam_height: float

    # Sample parameters
    sample_filename: str
    sample_location: np.ndarray  # [x, y, z]
    sample_orientation: np.ndarray  # Euler angles [phi1, Phi, phi2]
    sample_radius: float

    # Detector parameters
    detector_filename: str
    num_detectors: int

    # Omega ranges
    omega_filename: str

    # Output paths
    output_prefix: str

    # ... 40+ more fields

    @classmethod
    def from_file(cls, filepath: str) -> 'ConfigFile':
        """Parse .config file and create ConfigFile instance."""
        pass
```

**Estimated effort**: Medium (2-3 days)

#### 3. File I/O Utilities (InitFileIO namespace)

**Location**: `Src/InitFilesIO.h`, `Src/InitFilesIO.cpp`

**Key Functions to Port**:
```python
# Already implemented in simulation_range.py:
# - read_rotation_interval_files() → ReadRotationIntervalFiles()

# Need to add:
def read_detector_file(filename: str) -> List[Detector]:
    """Read detector geometry from file."""
    pass

def read_crystal_structure_file(filename: str) -> CrystalStructure:
    """Read crystal structure from .dat file."""
    pass

def read_fundamental_zone_file(filename: str) -> np.ndarray:
    """Read FZ orientation list."""
    pass
```

**Estimated effort**: Medium (can be done incrementally)

### Phase 1: Core ExperimentSetup Class

Once prerequisites are complete, port ExperimentSetup in this order:

#### Step 1.1: Base Class (ExperimentSetup)

**File**: `icenine_py/icenine/experiment_setup.py`

```python
from dataclasses import dataclass, field
from typing import List, Optional
import numpy as np

from .detector import Detector
from .simulation_range import SimulationRange, OmegaRange, FileRange
from .crystal_structure import CrystalStructure
from .geometry import Vector3

@dataclass
class StepSizeInfo:
    """Optimization parameter structure (C++ SStepSizeInfo)."""
    euler_steps: np.ndarray  # [3] - Euler angle step sizes
    detector_pos: np.ndarray  # [3] - Detector position
    beam_center_j: float
    beam_center_k: float
    pixel_height: float
    pixel_width: float
    angular_radius: float


class ExperimentSetup:
    """Base class for X-ray diffraction experiment setup.

    Python port of CExperimentSetup from Src/ExperimentSetup.h.
    Defines minimum requirements for an X-ray experiment.
    """

    def __init__(self, config_file: Optional['ConfigFile'] = None):
        self.config_file = config_file
        self.initialized = False

        # Beam parameters
        self.beam_direction = np.array([0.0, 0.0, 1.0])
        self.beam_energy = 0.0  # keV
        self.beam_energy_width = 0.0
        self.beam_deflection_chi_laue = 0.0
        self.eta_limit = np.pi / 2.0  # 90 degrees

        # Thresholds
        self.min_accepted_intensity_fraction = 0.0

    def set_config_file(self, config: 'ConfigFile') -> None:
        """Set configuration and parse basic parameters.

        C++ Reference: ExperimentSetup.cpp:231-258
        """
        self.config_file = config

        # Parse basic beam parameters from config
        self.beam_energy = config.beam_energy
        self.beam_energy_width = config.beam_energy_width
        self.beam_direction = config.beam_direction.copy()

        # Normalize beam direction
        norm = np.linalg.norm(self.beam_direction)
        if norm > 0:
            self.beam_direction /= norm

        self.min_accepted_intensity_fraction = config.min_intensity_fraction
        self.eta_limit = config.eta_limit

    # Getter methods
    def get_beam_energy(self) -> float:
        return self.beam_energy

    def get_beam_direction(self) -> np.ndarray:
        return self.beam_direction.copy()

    # ... more getters


class XDMExperimentSetup(ExperimentSetup):
    """XDM/HEDM-specific experiment setup with full initialization.

    Python port of CXDMExperimentSetup from Src/ExperimentSetup.h.
    """

    def __init__(self, config_file: Optional['ConfigFile'] = None):
        super().__init__(config_file)

        # Detector setup
        self.detector_list: List[Detector] = []

        # Omega range system
        self.omega_range_list: List[OmegaRange] = []
        self.file_range_list: List[FileRange] = []
        self.range_to_index_map: Optional[SimulationRange] = None

        # Optimization parameters
        self.optimization_info: List[StepSizeInfo] = []
        self.optimization_constrains: List[StepSizeInfo] = []
        self.detection_sensitivity: List[StepSizeInfo] = []

    def initialize_experiment(self) -> None:
        """Initialize experiment by reading all external files.

        C++ Reference: ExperimentSetup.cpp:340-385

        Reads:
        - Detector geometry file
        - Omega range file
        - Optimization parameter files

        Sets up:
        - detector_list
        - omega_range_list, file_range_list
        - range_to_index_map
        - optimization parameters
        """
        if self.config_file is None:
            raise RuntimeError("ConfigFile not set")

        # Read detector file
        detector_filename = self.config_file.detector_filename
        self.detector_list = read_detector_file(detector_filename)

        # Read omega ranges
        omega_filename = self.config_file.omega_filename
        num_detectors = self.config_file.num_detectors

        from .simulation_range import read_rotation_interval_files
        self.omega_range_list, self.file_range_list = \
            read_rotation_interval_files(omega_filename, num_detectors)

        # Process omega ranges (replicate C++ ReadRotationInterval logic)
        self._process_omega_ranges()

        # Read optimization parameters (if specified)
        # ... optional

        self.initialized = True

    def _process_omega_ranges(self) -> None:
        """Process omega ranges and create SimulationRange mapper.

        C++ Reference: ExperimentSetup.cpp:170-224

        Handles:
        - Swapping reversed ranges (high < low)
        - Determining overall range bounds
        - Creating angle-to-file-number mapping
        """
        if not self.omega_range_list:
            raise RuntimeError("No omega ranges loaded")

        # Swap individual ranges if reversed
        for omega_range in self.omega_range_list:
            if omega_range.high < omega_range.low:
                omega_range.low, omega_range.high = omega_range.high, omega_range.low

        # Determine overall range (may be flipped for descending sequences)
        first_range = self.omega_range_list[0]
        last_range = self.omega_range_list[-1]

        if first_range.high > last_range.high:
            # Flipped sequence (descending)
            low = first_range.high
            high = last_range.low
        else:
            # Normal sequence (ascending)
            low = first_range.low
            high = last_range.high

        # Width from first range
        width = first_range.high - first_range.low

        # Create SimulationRange mapper
        self.range_to_index_map = SimulationRange(
            low=low,
            high=high,
            width=width,
            range_list=self.omega_range_list,
            start_file_num=0
        )

    def initialize_sample(self, sample: 'Sample', detector: Detector) -> None:
        """Initialize sample with crystal structure and detection limits.

        C++ Reference: ExperimentSetup.cpp:290-332

        Args:
            sample: Sample object to initialize
            detector: Detector for calculating detection limits
        """
        if self.config_file is None:
            raise RuntimeError("ConfigFile not set")

        # Load microstructure from .mic file
        sample.load_mic_file(self.config_file.sample_filename)

        # Set sample location and orientation from config
        sample.set_location(self.config_file.sample_location)
        sample.set_orientation(self.config_file.sample_orientation)

        # Read crystal structure file
        structure_filename = self.config_file.structure_filename
        crystal = read_crystal_structure_file(structure_filename)

        # Calculate detection limits based on detector geometry
        max_q = self.get_max_q(detector, sample)
        crystal.set_detection_limit(max_q)

        # Apply symmetry to reflection vectors
        symmetry = sample.get_sample_symmetry()
        crystal.apply_symmetry(symmetry)

        # Add crystal structure to sample
        sample.add_crystal_structure(crystal)

    def get_max_q(self, detector: Detector, sample: 'Sample') -> float:
        """Calculate maximum scattering vector magnitude.

        C++ Reference: ExperimentSetup.cpp:52-70

        Used to limit reciprocal lattice generation to observable reflections.

        Args:
            detector: Detector geometry
            sample: Sample with location in lab frame

        Returns:
            Maximum |Q| in inverse angstroms
        """
        from .diffraction_core import k_from_energy

        # Incident wave vector magnitude
        k_mag = k_from_energy(self.beam_energy)

        # Get detector corner position
        detector_origin = detector.get_detector_coordinate_origin()

        # Vector from sample to detector corner
        sample_to_detector = detector_origin - sample.get_location()

        # Maximum scattering angle occurs at detector edge
        # For elastic scattering: |Q| = 2k*sin(theta/2)
        distance = np.linalg.norm(sample_to_detector)
        detector_radius = detector.get_detector_radius()

        # Maximum scattering angle
        theta_max = np.arctan(detector_radius / distance)

        # Maximum Q
        max_q = 2.0 * k_mag * np.sin(theta_max / 2.0)

        return max_q

    def get_reciprocal_vector(self, scattered_direction: np.ndarray) -> np.ndarray:
        """Convert scattered beam direction to reciprocal lattice vector.

        C++ Reference: ExperimentSetup.cpp:479-488

        Physics: G = k_out - k_in (elastic scattering)

        Args:
            scattered_direction: Normalized scattered beam direction

        Returns:
            Reciprocal lattice vector G in inverse angstroms
        """
        from .diffraction_core import k_from_energy

        k_mag = k_from_energy(self.beam_energy)

        # k_in = k * beam_direction
        k_in = k_mag * self.beam_direction

        # k_out = k * scattered_direction
        k_out = k_mag * scattered_direction

        # G = k_out - k_in
        reciprocal_vector = k_out - k_in

        return reciprocal_vector

    # Accessor methods
    def get_detector_list(self) -> List[Detector]:
        return self.detector_list

    def get_omega_range_list(self) -> List[OmegaRange]:
        return self.omega_range_list

    def get_file_range_list(self) -> List[FileRange]:
        return self.file_range_list

    def get_range_to_index_map(self) -> SimulationRange:
        if self.range_to_index_map is None:
            raise RuntimeError("Range to index map not initialized")
        return self.range_to_index_map
```

**Estimated effort**: Medium (3-4 days once prerequisites done)

#### Step 1.2: File I/O Functions

Add to `icenine_py/icenine/file_io.py` (new module):

```python
"""File I/O utilities for IceNine experiment files.

Python port of InitFileIO namespace from Src/InitFilesIO.h/cpp.
"""

from typing import List, Tuple
import numpy as np

from .detector import Detector
from .crystal_structure import CrystalStructure


def read_detector_file(filename: str) -> List[Detector]:
    """Read detector geometry from file.

    C++ Reference: InitFilesIO.cpp:ReadDetectorFile()

    File format:
        NumDetectors <n>
        {
            BeamCenterJ <j>
            BeamCenterK <k>
            PixelHeight <h>
            PixelWidth <w>
            ...
        }

    Returns:
        List of Detector objects
    """
    # Parse custom detector file format
    # Create Detector objects from parameters
    pass


def read_crystal_structure_file(filename: str) -> CrystalStructure:
    """Read crystal structure from .dat file.

    C++ Reference: InitFilesIO.cpp:ReadStructureFile()

    File format (example gold.dat):
        4.0782  # Lattice parameter (Angstroms)
        FCC     # Lattice type
        ...

    Returns:
        CrystalStructure object
    """
    pass


def read_fundamental_zone_file(filename: str) -> np.ndarray:
    """Read fundamental zone orientation list.

    C++ Reference: InitFilesIO.cpp:ReadFundamentalZoneFile()

    Returns:
        Array of Euler angles (N x 3)
    """
    pass
```

**Estimated effort**: Medium (2-3 days)

### Phase 2: Testing and Validation

#### Test Strategy

**Approach**: Generate C++ reference data, validate Python output

**Test Harness** (C++):
```cpp
// icenine_py/cpp_harness/test_experiment_setup.cpp

void export_experiment_setup_state(const CXDMExperimentSetup& setup,
                                     const string& output_json) {
    ofstream os(output_json);
    os << "{\n";

    // Export beam parameters
    os << "  \"beam_energy\": " << setup.GetBeamEnergy() << ",\n";
    os << "  \"beam_direction\": ["
       << setup.GetBeamDirection().m_fVec[0] << ", "
       << setup.GetBeamDirection().m_fVec[1] << ", "
       << setup.GetBeamDirection().m_fVec[2] << "],\n";

    // Export detector list
    const vector<CDetector>& detectors = setup.GetDetectorList();
    os << "  \"num_detectors\": " << detectors.size() << ",\n";
    os << "  \"detectors\": [\n";
    for (size_t i = 0; i < detectors.size(); i++) {
        export_detector(detectors[i], os);
        if (i < detectors.size() - 1) os << ",\n";
    }
    os << "  ],\n";

    // Export omega ranges
    const vector<SRange>& omega_ranges = setup.GetOmegaRangeList();
    os << "  \"num_omega_ranges\": " << omega_ranges.size() << ",\n";

    // Export max Q calculation
    CSample sample;
    setup.InitializeSample(sample, detectors[0]);
    Float max_q = setup.GetMaxQ(detectors[0], sample);
    os << "  \"max_q\": " << max_q << "\n";

    os << "}\n";
}
```

**Python Test**:
```python
# icenine_py/tests/test_experiment_setup.py

import pytest
import json
import numpy as np
from icenine.experiment_setup import XDMExperimentSetup
from icenine.config_file import ConfigFile


@pytest.fixture
def cpp_reference_data():
    """Load C++ reference data."""
    with open('cpp_outputs/experiment_setup_test_data.json') as f:
        return json.load(f)


def test_initialize_experiment(cpp_reference_data):
    """Test experiment initialization against C++ reference."""

    # Load config file
    config = ConfigFile.from_file('../../ConfigFiles/ReconstructTest.config')

    # Create and initialize
    setup = XDMExperimentSetup(config)
    setup.initialize_experiment()

    # Validate beam parameters
    assert np.isclose(setup.get_beam_energy(),
                     cpp_reference_data['beam_energy'])

    assert np.allclose(setup.get_beam_direction(),
                      cpp_reference_data['beam_direction'])

    # Validate detector list
    detectors = setup.get_detector_list()
    assert len(detectors) == cpp_reference_data['num_detectors']

    # Validate omega ranges
    omega_ranges = setup.get_omega_range_list()
    assert len(omega_ranges) == cpp_reference_data['num_omega_ranges']


def test_max_q_calculation(cpp_reference_data):
    """Test max Q calculation against C++ reference."""

    config = ConfigFile.from_file('../../ConfigFiles/ReconstructTest.config')
    setup = XDMExperimentSetup(config)
    setup.initialize_experiment()

    # Need Sample class for this test
    # sample = Sample()
    # setup.initialize_sample(sample, setup.get_detector_list()[0])

    # max_q = setup.get_max_q(setup.get_detector_list()[0], sample)
    # assert np.isclose(max_q, cpp_reference_data['max_q'])
    pass  # Skip until Sample is ported
```

#### Test Files

Use existing test data:
- `ConfigFiles/ReconstructTest.config` - Configuration file
- `DataFiles/omega_180_2L.dat` - Omega ranges (180 wedges)
- `ConfigFiles/DetectorFile.txt` - Detector geometry
- `DataFiles/gold.dat` - Crystal structure
- `DataFiles/Au1007_small.mic` - Microstructure (requires CSample)

### Phase 3: Integration

Once ExperimentSetup is ported, integrate with:

1. **Forward Simulation** (future task)
   - Use ExperimentSetup to initialize simulation parameters
   - Feed detector list, omega ranges to simulator

2. **Reconstruction** (future task)
   - Use ExperimentSetup to load experimental data
   - Initialize sample for reconstruction

3. **Visualization** (future task)
   - Use ExperimentSetup to display experimental geometry

## Implementation Checklist

### Prerequisites
- [ ] Port CSample class (CRITICAL - separate task)
- [ ] Port ConfigFile parser
- [ ] Implement read_detector_file()
- [ ] Implement read_crystal_structure_file()
- [ ] Implement read_fundamental_zone_file()

### Core Implementation
- [ ] Create StepSizeInfo dataclass
- [ ] Implement ExperimentSetup base class
- [ ] Implement XDMExperimentSetup.__init__()
- [ ] Implement set_config_file()
- [ ] Implement initialize_experiment()
- [ ] Implement _process_omega_ranges()
- [ ] Implement initialize_sample()
- [ ] Implement get_max_q()
- [ ] Implement get_reciprocal_vector()
- [ ] Implement all getter methods

### Testing
- [ ] Create C++ test harness (test_experiment_setup.cpp)
- [ ] Generate reference data JSON
- [ ] Implement test_initialize_experiment()
- [ ] Implement test_omega_range_processing()
- [ ] Implement test_max_q_calculation()
- [ ] Implement test_reciprocal_vector()
- [ ] Validate against ReconstructTest.config

### Documentation
- [ ] Add docstrings to all methods (reference C++ line numbers)
- [ ] Create usage examples
- [ ] Update icenine/__init__.py
- [ ] Update README with ExperimentSetup capabilities

## Risk Assessment

| Risk | Severity | Mitigation |
|------|----------|------------|
| CSample not ported | **CRITICAL** | Port CSample first (separate task) |
| Config file format complexity | High | Incremental parsing, validate each field |
| Detector file format parsing | Medium | Test with multiple detector files |
| Max Q calculation errors | Medium | Extensive validation against C++ |
| Floating point precision | Low | Use high precision throughout |

## Success Criteria

1. **All tests pass**: Python output matches C++ reference data
2. **File compatibility**: Can read all C++ input files (.config, .dat, etc.)
3. **Numerical accuracy**: Max Q, reciprocal vectors within 1e-10 of C++
4. **Integration ready**: Can be used by forward simulation and reconstruction

## Estimated Timeline

| Phase | Duration | Dependencies |
|-------|----------|--------------|
| CSample port | 2-3 weeks | MicFile (done), CrystalStructure (done) |
| ConfigFile parser | 3-4 days | None |
| File I/O utilities | 3-4 days | None |
| ExperimentSetup core | 4-5 days | All prerequisites |
| Testing & validation | 3-4 days | C++ test harness |
| Documentation | 2 days | Implementation complete |
| **TOTAL** | **4-5 weeks** | - |

## Next Steps

1. **Create separate CSample porting plan** (highest priority)
2. **Start ConfigFile parser** (can be done in parallel)
3. **Design file I/O module structure**
4. **Set up C++ test harness infrastructure**

---

**References**:
- C++ Source: `Src/ExperimentSetup.h`, `Src/ExperimentSetup.cpp`
- Dependencies: `Src/Sample.h`, `Src/ConfigFile.h`, `Src/InitFilesIO.h`
- Already Ported: `detector.py`, `simulation_range.py`, `crystal_structure.py`
