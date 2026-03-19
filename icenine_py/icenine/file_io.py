"""
File I/O utilities for IceNine experiment files.

Python port of InitFileIO namespace from Src/InitFilesIO.h/cpp and
DetectorFile parsing from Src/DetectorFile.h/cpp.

Handles reading:
- Detector geometry files (.txt format with { } blocks)
- Crystal structure files (via CrystalStructure factory methods)
- Omega range files (in simulation_range.py)

Author: S. F. Li
"""

from typing import List, Tuple, Dict, Optional
import numpy as np
import re

from .detector import Detector
from .geometry import euler_to_matrix
from .crystal_structure import CrystalStructure


# ============================================================================
# Detector File I/O
# ============================================================================

class DetectorInfo:
    """
    Detector geometry parameters from file.

    Python port of InitFileIO::CDetectorInfo from Src/DetectorFile.h.

    Attributes:
        j_unit_vector: J-axis direction in lab frame [x, y, z]
        k_unit_vector: K-axis direction in lab frame [x, y, z]
        beam_center_j: Beam center in J pixels
        beam_center_k: Beam center in K pixels
        lab_frame_location: Detector origin in lab frame (meters)
        lab_frame_orientation_euler: Detector orientation as Euler angles (radians)
        lab_frame_orientation_matrix: 3x3 rotation matrix
        num_j_pixels: Number of pixels in J direction
        num_k_pixels: Number of pixels in K direction
        pixel_width: Pixel width (J direction) in mm
        pixel_height: Pixel height (K direction) in mm
    """

    def __init__(self):
        """Initialize with default values."""
        self.j_unit_vector: np.ndarray = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        self.k_unit_vector: np.ndarray = np.array([0.0, -1.0, 0.0], dtype=np.float32)
        self.beam_center_j: float = 0.0
        self.beam_center_k: float = 0.0
        self.lab_frame_location: np.ndarray = np.zeros(3, dtype=np.float32)
        self.lab_frame_orientation_euler: np.ndarray = np.zeros(3, dtype=np.float32)
        self.lab_frame_orientation_matrix: np.ndarray = np.eye(3, dtype=np.float32)
        self.num_j_pixels: int = 0
        self.num_k_pixels: int = 0
        self.pixel_width: float = 0.0  # mm
        self.pixel_height: float = 0.0  # mm

    def parse_detector_block(self, block_text: str) -> bool:
        """
        Parse detector parameters from a { } block.

        C++ Reference: DetectorFile.cpp:72-200 (CDetectorInfo::ParseDetectorBlock)

        File format:
            {
            JUnitVector 1 0 0
            KUnitVector 0 -1 0
            BeamCenterJ 524.0225
            BeamCenterK 1009.933
            LabFrameLocation 4.60 0 0
            LabFrameOrientation 90 90 0
            NumJPixels 1024
            NumKPixels 1024
            PixelJLength 0.0040845890
            PixelKLength 0.0040845890
            }

        Args:
            block_text: Text content inside { } block

        Returns:
            True if parsing succeeded, False otherwise
        """
        lines = block_text.strip().split('\n')

        for line_num, line in enumerate(lines, 1):
            # Remove comments and whitespace
            line = line.split('#')[0].strip()
            if not line:
                continue

            tokens = line.split()
            if len(tokens) < 2:
                continue

            keyword = tokens[0]

            try:
                if keyword == "JUnitVector":
                    self.j_unit_vector = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])],
                                                  dtype=np.float32)
                elif keyword == "KUnitVector":
                    self.k_unit_vector = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])],
                                                  dtype=np.float32)
                elif keyword == "BeamCenterJ":
                    self.beam_center_j = float(tokens[1])
                elif keyword == "BeamCenterK":
                    self.beam_center_k = float(tokens[1])
                elif keyword == "LabFrameLocation":
                    self.lab_frame_location = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])],
                                                       dtype=np.float32)
                elif keyword == "LabFrameOrientation":
                    # Input is in degrees
                    euler_deg = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])
                    self.lab_frame_orientation_euler = np.deg2rad(euler_deg).astype(np.float32)
                    # Build rotation matrix using active Euler angles (ZXZ Bunge convention)
                    # euler_to_matrix takes DEGREES, so pass the original degree values
                    self.lab_frame_orientation_matrix = euler_to_matrix(
                        euler_deg[0], euler_deg[1], euler_deg[2]
                    )
                elif keyword == "NumJPixels":
                    self.num_j_pixels = int(tokens[1])
                elif keyword == "NumKPixels":
                    self.num_k_pixels = int(tokens[1])
                elif keyword == "PixelJLength":
                    self.pixel_width = float(tokens[1])  # mm
                elif keyword == "PixelKLength":
                    self.pixel_height = float(tokens[1])  # mm
                else:
                    print(f"Warning: Unknown keyword '{keyword}' on line {line_num}")

            except (ValueError, IndexError) as e:
                print(f"Error parsing line {line_num}: {line}")
                print(f"  {e}")
                return False

        return True

    def get_detector(self) -> Detector:
        """
        Create Detector object from parsed parameters.

        C++ Reference: DetectorFile.cpp:55-63 (CDetectorInfo::GetDetector)

        Returns:
            Detector object with geometry set from file
        """
        import torch

        # Create detector using factory method
        # C++ calls CXDMDetectorFactory::MakeDetector
        # Note: Detector uses (num_rows, num_cols) which corresponds to (K, J) in C++
        # J/K unit vectors come from the detector file and define the detector
        # coordinate system before orientation rotation is applied
        detector = Detector(
            num_rows=self.num_k_pixels,  # K direction = rows
            num_cols=self.num_j_pixels,  # J direction = columns
            beam_center_j=self.beam_center_j,
            beam_center_k=self.beam_center_k,
            pixel_width=self.pixel_width,
            pixel_height=self.pixel_height,
            position=torch.from_numpy(self.lab_frame_location),
            orientation=torch.from_numpy(self.lab_frame_orientation_matrix),
            j_unit_vector=torch.from_numpy(self.j_unit_vector),
            k_unit_vector=torch.from_numpy(self.k_unit_vector),
        )

        return detector


def read_detector_file(filename: str) -> List[Detector]:
    """
    Read detector geometry file.

    Python port of InitFileIO::CDetectorFile::Parse() from Src/DetectorFile.cpp:207-251.

    File format consists of one or more { } blocks, each defining a detector:
        # Comment lines start with #
        {
        JUnitVector 1 0 0
        KUnitVector 0 -1 0
        BeamCenterJ 524.0225
        BeamCenterK 1009.933
        LabFrameLocation 4.60 0 0
        LabFrameOrientation 90 90 0
        NumJPixels 1024
        NumKPixels 1024
        PixelJLength 0.0040845890
        PixelKLength 0.0040845890
        }

        {
        # Second detector...
        }

    Args:
        filename: Path to detector file

    Returns:
        List of Detector objects

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If parsing fails

    Example:
        >>> detectors = read_detector_file('ConfigFiles/DetectorFile.txt')
        >>> print(f"Loaded {len(detectors)} detectors")
        Loaded 2 detectors
    """
    try:
        with open(filename, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"Detector file not found: {filename}")

    # Find all { } blocks
    # C++ Reference: DetectorFile.cpp:262-303 (FindDetectorBlock)
    blocks = []
    brace_level = 0
    current_block = []

    for line in content.split('\n'):
        # Count braces
        open_braces = line.count('{')
        close_braces = line.count('}')

        if open_braces > 0:
            if brace_level > 0:
                raise ValueError(f"Nested {{ found in detector file: {filename}")
            brace_level += open_braces
            continue

        if close_braces > 0:
            if brace_level == 0:
                raise ValueError(f"Unmatched }} found in detector file: {filename}")
            brace_level -= close_braces
            # End of block - save it
            blocks.append('\n'.join(current_block))
            current_block = []
            continue

        if brace_level > 0:
            current_block.append(line)

    if brace_level != 0:
        raise ValueError(f"Unmatched {{ in detector file: {filename}")

    if not blocks:
        raise ValueError(f"No detector blocks found in {filename}")

    # Parse each block
    detectors = []
    for i, block_text in enumerate(blocks):
        detector_info = DetectorInfo()
        success = detector_info.parse_detector_block(block_text)
        if not success:
            raise ValueError(f"Failed to parse detector block {i+1} in {filename}")
        detectors.append(detector_info.get_detector())

    return detectors


# ============================================================================
# Crystal Structure File I/O
# ============================================================================

def read_crystal_structure_file(filename: str, element: str = "Au", lattice_type: str = "FCC") -> CrystalStructure:
    """
    Read crystal structure from .dat file.

    C++ Reference: InitFilesIO.cpp:275-336 ReadCrystalStructureFile()

    File format (comma/space/tab/# delimited):
        Line 1: a, b, c (unit cell lengths in Angstroms)
        Line 2: alpha, beta, gamma (unit cell angles in degrees)
        Line 3: number of atoms
        Lines 4+: Z, x, y, z (atomic number and fractional coordinates)

    Args:
        filename: Path to structure .dat file
        element: Element symbol (unused - read from file)
        lattice_type: Lattice type (unused - determined from file)

    Returns:
        CrystalStructure object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid

    Example:
        >>> copper = read_crystal_structure_file("DataFiles/copper.dat")
        >>> copper.lattice.a
        3.61
    """
    from pymatgen.core.structure import Structure
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.periodic_table import Element

    # Read file
    try:
        with open(filename, 'rb') as f:
            content = f.read().decode('utf-8', errors='ignore')
    except FileNotFoundError:
        raise FileNotFoundError(f"Crystal structure file not found: {filename}")

    # Tokenize by commas, #, spaces, tabs, newlines
    # C++: GeneralLib::Tokenize( vsTokens, sBuffStr, ",# \t\n");
    lines = []
    for line in content.split('\n'):
        # Split by delimiters and filter empty tokens
        tokens = []
        for token in line.replace(',', ' ').replace('#', ' ').split():
            token = token.strip()
            if token:
                tokens.append(token)
        if tokens:
            lines.append(tokens)

    if len(lines) < 4:
        raise ValueError(
            f"Crystal structure file too short: expected >= 4 lines, got {len(lines)}"
        )

    # Parse unit cell lengths (Line 1)
    # C++: oCell.SetUnitCellLength( atof(vsTokens[0][0]), ... )
    if len(lines[0]) < 3:
        raise ValueError(f"Failed to parse a, b, c from line 1: {lines[0]}")
    a = float(lines[0][0])
    b = float(lines[0][1])
    c = float(lines[0][2])

    # Parse unit cell angles (Line 2, in degrees)
    # C++: oCell.SetUnitCellBasisAngles( DEGREE_TO_RADIAN(...), ... )
    if len(lines[1]) < 3:
        raise ValueError(f"Failed to parse alpha, beta, gamma from line 2: {lines[1]}")
    alpha = float(lines[1][0])
    beta = float(lines[1][1])
    gamma = float(lines[1][2])

    # Parse number of atoms (Line 3)
    # C++: oCell.SetNumAtoms( atoi( vsTokens[2][0].c_str() ) )
    if len(lines[2]) < 1:
        raise ValueError(f"Failed to parse number of atoms from line 3: {lines[2]}")
    num_atoms = int(lines[2][0])

    # Parse atoms (Lines 4+)
    # C++: for each line: Z, x, y, z
    species = []
    coords = []
    for i in range(3, min(3 + num_atoms, len(lines))):
        if len(lines[i]) < 4:
            raise ValueError(
                f"Invalid atom specification at line {i+1}: expected Z,x,y,z, got {lines[i]}"
            )

        # C++: oCell.AddAtom( atoi( vsTokens[i][0].c_str() ), oPos )
        atomic_number = int(lines[i][0])
        x = float(lines[i][1])
        y = float(lines[i][2])
        z = float(lines[i][3])

        # Convert atomic number to element symbol
        element_symbol = Element.from_Z(atomic_number).symbol
        species.append(element_symbol)
        coords.append([x, y, z])

    # Create pymatgen Structure
    # C++: oCell.InitializeCoordinateSystem()
    lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    structure = Structure(lattice, species, coords)

    # Create CrystalStructure
    return CrystalStructure(structure)


# ============================================================================
# Helper Functions
# ============================================================================

def _extract_vector(tokens: List[str], line_num: int) -> np.ndarray:
    """
    Extract 3D vector from token list.

    Helper function matching C++ Parser::ExtractVector.

    Args:
        tokens: Token list [keyword, x, y, z, ...]
        line_num: Line number for error reporting

    Returns:
        3D vector as numpy array
    """
    if len(tokens) < 4:
        raise ValueError(f"Line {line_num}: Vector requires 3 values, got {len(tokens)-1}")
    return np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])], dtype=np.float32)


def _extract_real(tokens: List[str], line_num: int) -> float:
    """Extract float from token list."""
    if len(tokens) < 2:
        raise ValueError(f"Line {line_num}: Missing value")
    return float(tokens[1])


def _extract_int(tokens: List[str], line_num: int) -> int:
    """Extract integer from token list."""
    if len(tokens) < 2:
        raise ValueError(f"Line {line_num}: Missing value")
    return int(tokens[1])
