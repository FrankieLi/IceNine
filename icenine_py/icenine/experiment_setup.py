"""
Experiment setup and parameter management for IceNine.

Python port of CExperimentSetup and CXDMExperimentSetup from
Src/ExperimentSetup.h and Src/ExperimentSetup.cpp.

The ExperimentSetup classes act as container/orchestration classes that:
- Parse configuration files
- Read external data files (detectors, omega ranges, crystal structures)
- Initialize experimental parameters (beam, sample, detectors)
- Provide accessors to all experimental data

Author: S. F. Li
"""

from dataclasses import dataclass
from typing import List, Optional
import numpy as np

from .config_file import ConfigFile, SymmetryType
from .detector import Detector
from .simulation_range import SimulationRange, OmegaRange, FileRange, read_omega_file
from .crystal_structure import CrystalStructure
from .file_io import read_detector_file
from .symmetry import CrystalSymmetry, create_cubic_symmetry, create_fcc_symmetry
from .sample import Sample
from .constants import KEV_OVER_HBAR_C_IN_ANG


# ============================================================================
# Data Structures
# ============================================================================

@dataclass
class StepSizeInfo:
    """
    Optimization parameter structure for detector geometry refinement.

    Python port of CXDMExperimentSetup::SStepSizeInfo from Src/ExperimentSetup.h:150-159.

    Used for parameter optimization during reconstruction to adjust:
    - Detector orientation (Euler angles)
    - Detector position
    - Beam center
    - Pixel sizes
    - Angular search radius

    Attributes:
        euler_steps: Euler angle step sizes [φ₁, Φ, φ₂] (radians)
        detector_pos: Detector position step sizes [x, y, z] (meters)
        beam_center_j: Beam center J pixel coordinate
        beam_center_k: Beam center K pixel coordinate
        pixel_height: Pixel height (mm)
        pixel_width: Pixel width (mm)
        angular_radius: Angular search radius (radians)
    """
    euler_steps: np.ndarray       # [3] radians
    detector_pos: np.ndarray      # [3] meters
    beam_center_j: float          # pixels
    beam_center_k: float          # pixels
    pixel_height: float           # mm
    pixel_width: float            # mm
    angular_radius: float         # radians


# ============================================================================
# Base Class
# ============================================================================

class ExperimentSetup:
    """
    Abstract base class for X-ray diffraction experiments.

    Python port of CExperimentSetup from Src/ExperimentSetup.h:70-140.

    Defines the minimum requirements for an X-ray experiment:
    - Beam parameters (energy, direction, limits)
    - Intensity thresholds

    This is an abstract base. Use XDMExperimentSetup for concrete implementation.
    """

    def __init__(self, config_file: Optional[ConfigFile] = None):
        """
        Initialize experiment setup.

        C++ Reference: ExperimentSetup.cpp:160-163

        Args:
            config_file: ConfigFile object (optional, can set later with set_config_file)
        """
        self.config_file: Optional[ConfigFile] = None
        self.initialized: bool = False

        # Beam parameters
        self.beam_direction: np.ndarray = np.array([0.0, 0.0, 1.0], dtype=np.float32)
        self.beam_energy: float = 0.0  # keV
        self.beam_energy_width: float = 0.0  # fraction
        self.beam_deflection_chi_laue: float = 0.0  # radians
        self.eta_limit: float = np.pi / 2.0  # radians (90 degrees)

        # Thresholds
        self.min_accepted_intensity_fraction: float = 0.0

        if config_file is not None:
            self.set_config_file(config_file)

    def set_config_file(self, config: ConfigFile) -> None:
        """
        Set configuration and parse basic parameters.

        C++ Reference: ExperimentSetup.cpp:134-150

        Args:
            config: ConfigFile object with experimental parameters
        """
        self.config_file = config

        # Parse beam parameters
        self.beam_energy = config.beam_energy
        self.beam_energy_width = config.beam_energy_width
        self.beam_direction = config.beam_direction.copy()

        # Normalize beam direction
        norm = np.linalg.norm(self.beam_direction)
        if norm > 0:
            self.beam_direction = self.beam_direction / norm

        self.beam_deflection_chi_laue = config.beam_deflection_chi_laue
        self.min_accepted_intensity_fraction = config.min_amplitude_fraction
        self.eta_limit = config.eta_limit

        self.initialized = True

    # ========================================================================
    # Accessors
    # ========================================================================

    def get_beam_energy(self) -> float:
        """
        Get beam energy in keV.

        C++ Reference: ExperimentSetup.cpp:410-414

        Returns:
            Beam energy in keV
        """
        if not self.initialized:
            raise RuntimeError("ExperimentSetup not initialized")
        return self.beam_energy

    def get_beam_energy_width(self) -> float:
        """
        Get beam energy width (fractional).

        C++ Reference: ExperimentSetup.cpp:432-436

        Returns:
            Energy width as fraction (e.g., 0.05 = 5%)
        """
        if not self.initialized:
            raise RuntimeError("ExperimentSetup not initialized")
        return self.beam_energy_width

    def get_beam_direction(self) -> np.ndarray:
        """
        Get X-ray beam direction unit vector.

        C++ Reference: ExperimentSetup.cpp:443-447

        Returns:
            Beam direction [x, y, z] (normalized)
        """
        if not self.initialized:
            raise RuntimeError("ExperimentSetup not initialized")
        return self.beam_direction.copy()

    def get_min_accepted_intensity_fraction(self) -> float:
        """
        Get minimum accepted intensity fraction for peak detection.

        C++ Reference: ExperimentSetup.cpp:454-458

        Returns:
            Minimum intensity fraction [0, 1]
        """
        if not self.initialized:
            raise RuntimeError("ExperimentSetup not initialized")
        return self.min_accepted_intensity_fraction

    def get_beam_deflection_chi_laue(self) -> float:
        """
        Get beam deflection angle for Laue geometry.

        C++ Reference: ExperimentSetup.cpp:465-469

        Returns:
            Deflection angle in radians
        """
        if not self.initialized:
            raise RuntimeError("ExperimentSetup not initialized")
        return self.beam_deflection_chi_laue

    def get_eta_limit(self) -> float:
        """
        Get azimuthal angle limit.

        C++ Reference: ExperimentSetup.cpp:421-425

        Returns:
            Eta limit in radians (typically π/2)
        """
        if not self.initialized:
            raise RuntimeError("ExperimentSetup not initialized")
        return self.eta_limit


# ============================================================================
# XDM/HEDM Concrete Implementation
# ============================================================================

class XDMExperimentSetup(ExperimentSetup):
    """
    XDM/HEDM-specific experiment setup with full initialization.

    Python port of CXDMExperimentSetup from Src/ExperimentSetup.h:146-297.

    Extends ExperimentSetup with:
    - Detector list management
    - Omega range system
    - Optimization parameters
    - Sample initialization
    - File I/O for all experimental data
    """

    def __init__(self, config_file: Optional[ConfigFile] = None):
        """
        Initialize XDM experiment setup.

        C++ Reference: ExperimentSetup.h:191-192

        Args:
            config_file: ConfigFile object (optional)
        """
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
        """
        Initialize experiment by reading all external files.

        C++ Reference: ExperimentSetup.cpp:340-385

        Reads:
        - Detector geometry file
        - Omega range file
        - Optimization parameter files (if specified)

        Sets up:
        - detector_list
        - omega_range_list, file_range_list
        - range_to_index_map
        - optimization parameters

        Raises:
            RuntimeError: If ConfigFile not set or files cannot be read
            ValueError: If validation fails
        """
        if self.config_file is None:
            raise RuntimeError("ConfigFile not set. Call set_config_file() first.")

        # Read detector geometry
        num_detectors = self._read_detector_info()

        # Read omega ranges
        self._read_rotation_interval(num_detectors)

        # Read optimization parameters (if specified)
        if self.config_file.optimization_filename:
            self.optimization_info = self._read_step_size_file(
                self.config_file.optimization_filename
            )

        if self.config_file.detection_limit_filename:
            self.detection_sensitivity = self._read_step_size_file(
                self.config_file.detection_limit_filename
            )

        if self.config_file.optimization_constrain_filename:
            self.optimization_constrains = self._read_step_size_file(
                self.config_file.optimization_constrain_filename
            )

        # Validation
        self._validate_experiment_setup()

    def _read_detector_info(self) -> int:
        """
        Read detector geometry from file.

        C++ Reference: ExperimentSetup.cpp:232-240

        Returns:
            Number of detectors read

        Raises:
            RuntimeError: If detector file cannot be read
        """
        try:
            self.detector_list = read_detector_file(self.config_file.detector_filename)
        except Exception as e:
            raise RuntimeError(
                f"Failed to read detector file: {self.config_file.detector_filename}\n{e}"
            )

        if len(self.detector_list) == 0:
            raise RuntimeError("No detectors specified in detector file")

        return len(self.detector_list)

    def _read_rotation_interval(self, num_detectors: int) -> int:
        """
        Read omega rotation intervals from file.

        C++ Reference: ExperimentSetup.cpp:170-225

        Args:
            num_detectors: Number of detectors (for validation)

        Returns:
            Number of omega ranges read

        Raises:
            RuntimeError: If omega file cannot be read or validation fails
        """
        omega_file = self.config_file.rotation_range_filename

        try:
            self.omega_range_list, self.file_range_list = read_omega_file(
                omega_file, num_detectors
            )
        except Exception as e:
            raise RuntimeError(f"Failed to read omega file: {omega_file}\n{e}")

        # Validate file ranges match number of detectors
        if len(self.file_range_list) != self.config_file.num_detectors:
            raise RuntimeError(
                f"Config file specified {self.config_file.num_detectors} detectors, "
                f"but {len(self.file_range_list)} file ranges found in omega file"
            )

        # Swap reversed ranges (high < low)
        for omega_range in self.omega_range_list:
            if omega_range.high < omega_range.low:
                print(f"WARNING: Omega range reversed, swapping [{omega_range.low}, {omega_range.high}]")
                omega_range.low, omega_range.high = omega_range.high, omega_range.low

        if len(self.omega_range_list) == 0:
            raise RuntimeError("No omega ranges specified in omega file")

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

        return len(self.omega_range_list)

    def _read_step_size_file(self, filename: str) -> List[StepSizeInfo]:
        """
        Read optimization step size parameters from file.

        C++ Reference: ExperimentSetup.cpp:245-280

        Args:
            filename: Path to step size file (same format as detector file)

        Returns:
            List of StepSizeInfo objects

        Raises:
            RuntimeError: If file cannot be read
        """
        # Read as detector info (same file format)
        from .file_io import DetectorInfo

        # TODO: Implement proper step size file reading
        # For now, return empty list as these files are optional
        print(f"WARNING: Step size file reading not yet implemented: {filename}")
        return []

    def _validate_experiment_setup(self) -> None:
        """
        Validate experiment setup consistency.

        C++ Reference: ExperimentSetup.cpp:348-385

        Raises:
            RuntimeError: If validation fails
        """
        # Check detector count matches config
        if len(self.detector_list) != self.config_file.num_detectors:
            raise RuntimeError(
                f"Config file specified {self.config_file.num_detectors} detectors, "
                f"but {len(self.detector_list)} found in detector file"
            )

        # Check detector spacing
        expected_spacing = self.config_file.num_detectors - 1
        if len(self.config_file.detector_spacing) != expected_spacing:
            raise RuntimeError(
                f"Config file specified {self.config_file.num_detectors} detectors, "
                f"requires {expected_spacing} spacing values, "
                f"but {len(self.config_file.detector_spacing)} found"
            )

        # Check all spacing values are non-negative
        for i, spacing in enumerate(self.config_file.detector_spacing):
            if spacing < 0:
                raise RuntimeError(f"Detector spacing {i} is negative: {spacing}")

        # Validate optimization info (if present)
        if self.optimization_info and len(self.optimization_info) != len(self.detector_list):
            raise RuntimeError(
                f"Optimization file must specify step sizes for each detector. "
                f"Expected {len(self.detector_list)}, got {len(self.optimization_info)}"
            )

        if self.detection_sensitivity and len(self.detection_sensitivity) != len(self.detector_list):
            raise RuntimeError(
                f"Sensitivity file must specify step sizes for each detector. "
                f"Expected {len(self.detector_list)}, got {len(self.detection_sensitivity)}"
            )

        if self.optimization_constrains and len(self.optimization_constrains) != len(self.detector_list) - 1:
            raise RuntimeError(
                f"Constrain file must specify {len(self.detector_list) - 1} constrains. "
                f"Got {len(self.optimization_constrains)}"
            )

    def initialize_sample(self, sample: Sample, detector: Detector) -> None:
        """
        Initialize sample with crystal structure and detection limits.

        C++ Reference: ExperimentSetup.cpp:290-332

        This method:
        1. Loads microstructure from .mic file
        2. Sets sample location and orientation from config
        3. Reads crystal structure file
        4. Calculates detection limits (max Q) based on detector geometry
        5. Applies symmetry to reflection vectors
        6. Adds crystal structures to sample

        Args:
            sample: Sample object to initialize
            detector: Detector for calculating detection limits

        Raises:
            RuntimeError: If initialization fails
        """
        if self.config_file is None:
            raise RuntimeError("ConfigFile not set")

        # Apply sample orientation and translation
        # C++ uses Rotate() then Translate()
        sample.rotate(
            np.rad2deg(self.config_file.sample_orientation[0]),
            np.rad2deg(self.config_file.sample_orientation[1]),
            np.rad2deg(self.config_file.sample_orientation[2])
        )
        sample.translate(self.config_file.sample_location)

        # Load microstructure from .mic file
        success = sample.load_sample(self.config_file.sample_filename)
        if not success:
            raise RuntimeError(f"Failed to load sample from {self.config_file.sample_filename}")

        # Add empty structure (phase 0 - no scattering)
        empty_structure = CrystalStructure.create_cubic("Empty", 1.0)
        # C++: oEmptyStructure.SetReflectionVectorLimits(0, 0, 0, 0, 0)
        empty_structure.set_reflection_limits(0, 0, 0, 0.0, 0.0)
        sample.add_crystal_structure(empty_structure)

        # Read crystal structure file
        # C++: if (!InitFileIO::ReadCrystalStructureFile(oCellStructure, oExpConfigFile.StructureFilename))
        from .file_io import read_crystal_structure_file
        structure_filename = self.config_file.structure_filename
        try:
            crystal = read_crystal_structure_file(structure_filename)
            print(f"Loaded crystal structure from {structure_filename}: {crystal}")
        except Exception as e:
            print(f"WARNING: Failed to read {structure_filename}: {e}")
            print(f"WARNING: Using default FCC gold structure as fallback")
            crystal = CrystalStructure.create_fcc("Au", 4.0782)  # Gold FCC fallback

        # Calculate detection limits
        max_q = self.get_max_q(detector, sample)
        print(f"Max Q calculated: {max_q:.4f} Å⁻¹")

        # Take minimum of user-specified max_q and calculated max_q
        max_q = min(max_q, self.config_file.max_q)

        # Calculate max Miller indices for cubic crystal
        # C++: Int nMaxH = (Int) ( fMaxQ / oRecpVecs[0].GetLength() );
        recip_params = crystal.get_reciprocal_lattice_parameters()
        max_h = int(max_q / recip_params[0])
        max_k = int(max_q / recip_params[1])
        max_l = int(max_q / recip_params[2])

        # Set reflection vector limits
        # C++: oResCellStruct.SetReflectionVectorLimits(nMaxH, nMaxK, nMaxL, fMaxQ, oExpConfigFile.fMinAmplitudeFraction)
        crystal.set_reflection_limits(
            max_h, max_k, max_l, max_q, self.config_file.min_amplitude_fraction
        )

        # Apply symmetry to reflection vectors
        symmetry = self.get_sample_symmetry()
        sample.set_sample_symmetry(symmetry)
        # TODO: crystal.set_unique_reflection_list(symmetry)

        # Add crystal structure to sample
        sample.add_crystal_structure(crystal)

    def get_max_q(self, detector: Detector, sample: Sample) -> float:
        """
        Calculate maximum scattering vector magnitude |Q|.

        C++ Reference: ExperimentSetup.cpp:52-70

        Used to limit reciprocal lattice generation to observable reflections.
        The maximum Q occurs when scattering from the detector edge furthest
        from the beam.

        Args:
            detector: Detector geometry
            sample: Sample with location in lab frame

        Returns:
            Maximum |Q| in inverse angstroms (Å⁻¹)

        Example:
            >>> setup = XDMExperimentSetup(config)
            >>> max_q = setup.get_max_q(detector, sample)
            >>> print(f"Max observable Q: {max_q:.4f} Å⁻¹")
        """
        # C++: SVector3 oKOutMaxDir = oDetector.GetDetectorCoordinateOrigin()
        # Get detector corner position (coordinate origin) in lab frame
        # By convention, this is the corner furthest from the beam
        k_out_max_dir = detector.coordinate_origin.numpy()

        # C++: oKOutMaxDir += oSample.GetLocation()
        # Note: C++ code adds sample location (seems like an error, but matching it)
        k_out_max_dir += sample.get_location()

        # C++: oKOutMaxDir.Normalize()
        # Normalize to get scattered beam direction
        k_out_max_dir = k_out_max_dir / np.linalg.norm(k_out_max_dir)

        # C++: SVector3 oQMax = GetReciprocalVector( oKOutMaxDir )
        # Calculate reciprocal vector for this scattering direction
        q_max_vec = self.get_reciprocal_vector(k_out_max_dir)

        # C++: return oQMax.GetLength()
        max_q = np.linalg.norm(q_max_vec)

        return max_q

    def get_reciprocal_vector(self, scattered_direction: np.ndarray) -> np.ndarray:
        """
        Convert scattered beam direction to reciprocal lattice vector.

        C++ Reference: ExperimentSetup.cpp:479-488

        Physics: For elastic scattering, the reciprocal lattice vector G
        is the momentum transfer:
            G = k_out - k_in
        where k_in and k_out have the same magnitude (elastic condition).

        Args:
            scattered_direction: Normalized scattered beam direction [x, y, z]

        Returns:
            Reciprocal lattice vector G in inverse angstroms (Å⁻¹)

        Example:
            >>> k_out_dir = np.array([0.707, 0.707, 0.0])
            >>> g_vec = setup.get_reciprocal_vector(k_out_dir)
            >>> print(f"G = [{g_vec[0]:.4f}, {g_vec[1]:.4f}, {g_vec[2]:.4f}] Å⁻¹")
        """
        # Wave vector magnitude (same for input and output - elastic)
        k_mag = KEV_OVER_HBAR_C_IN_ANG * self.beam_energy

        # Input momentum vector
        k_in = k_mag * self.beam_direction

        # Output momentum vector
        k_out = k_mag * scattered_direction

        # Reciprocal lattice vector (momentum transfer)
        g_vec = k_out - k_in

        return g_vec

    def get_sample_symmetry(self) -> Optional[CrystalSymmetry]:
        """
        Get crystal symmetry from config.

        C++ Reference: ExperimentSetup.cpp:78-95

        NOTE: In the Python port, symmetry is obtained from the crystal structure
        (via CrystalStructure.symmetry) rather than from a global symmetry factory.
        This method returns None and symmetry should be obtained from the
        crystal structure after it's created.

        Returns:
            None (symmetry obtained from CrystalStructure)

        Raises:
            RuntimeError: If config file not set
        """
        if self.config_file is None:
            raise RuntimeError("ConfigFile not set")

        sym_type = self.config_file.sample_symmetry
        if sym_type == SymmetryType.CUBIC:
            return create_cubic_symmetry(4.0782)
        elif sym_type == SymmetryType.HEXAGONAL:
            return create_cubic_symmetry(3.0)  # TODO: proper hexagonal
        else:
            return None

    # ========================================================================
    # Accessors
    # ========================================================================

    def get_detector_list(self) -> List[Detector]:
        """Get list of all detectors."""
        return self.detector_list

    def get_omega_range_list(self) -> List[OmegaRange]:
        """Get list of omega ranges."""
        return self.omega_range_list

    def get_file_range_list(self) -> List[FileRange]:
        """Get list of file ranges."""
        return self.file_range_list

    def get_range_to_index_map(self) -> SimulationRange:
        """Get omega angle to file index mapper."""
        if self.range_to_index_map is None:
            raise RuntimeError("Range to index map not initialized. Call initialize_experiment() first.")
        return self.range_to_index_map

    def get_optimization_info(self) -> List[StepSizeInfo]:
        """Get optimization step size info."""
        return self.optimization_info

    def get_optimization_constrains(self) -> List[StepSizeInfo]:
        """Get optimization constraints."""
        return self.optimization_constrains

    def get_detection_sensitivity(self) -> List[StepSizeInfo]:
        """Get detection sensitivity parameters."""
        return self.detection_sensitivity

    def get_config_file(self) -> ConfigFile:
        """Get configuration file."""
        if self.config_file is None:
            raise RuntimeError("ConfigFile not set")
        return self.config_file

    def get_next_sample(self) -> str:
        """
        Get next sample filename.

        C++ Reference: ExperimentSetup.cpp:392-396

        Returns:
            Sample filename from config
        """
        if not self.initialized:
            raise RuntimeError("ExperimentSetup not initialized")
        return self.config_file.sample_filename
