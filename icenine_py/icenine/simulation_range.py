"""
Simulation Range Module - Omega Range System for IceNine

This module implements the omega range system for handling discontinuous data
collection in synchrotron X-ray diffraction experiments. In physical experiments,
data cannot be collected continuously across all rotation angles due to:
- Beam shutter constraints
- Detector readout time
- Goniometer limitations

Data is collected in discrete angular "wedges" (e.g., [-90°, -85°], [-45°, -40°]).
This module provides:
- Data structures for angular ranges (OmegaRange) and file ranges (FileRange)
- Mapper (SimulationRange) to determine if a predicted omega angle is observable
- Omega file parser for loading experimental wedge configurations

C++ References:
- XDM++/libXDM/3dMath.h: SRange (line 474), SIntRange (line 459)
- Src/SimulationData.h: CSimulationRange (lines 288-565)
- Src/InitFilesIO.cpp: ReadRotationIntervalFiles (line 350)

Usage:
    >>> from icenine.simulation_range import OmegaRange, SimulationRange
    >>>
    >>> # Define experimental wedges (in radians)
    >>> wedges = [
    ...     OmegaRange(low=-1.5708, high=-1.4835),  # -90° to -85°
    ...     OmegaRange(low=0.6981, high=0.7854)     # 40° to 45°
    ... ]
    >>>
    >>> # Create mapper
    >>> mapper = SimulationRange(
    ...     low=-1.5708, high=1.5708, width=0.0175,
    ...     range_list=wedges
    ... )
    >>>
    >>> # Check if omega is observable
    >>> mapper.is_in_experimental_range(-1.5)  # True (in first wedge)
    >>> mapper.is_in_experimental_range(0.0)   # False (gap between wedges)
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np


@dataclass
class OmegaRange:
    """Angular range for data collection (radians).

    Represents a continuous angular wedge during which experimental data
    was collected. Samples rotate through omega (rotation about z-axis),
    but data collection is discontinuous due to physical constraints.

    C++ Reference: XDM++/libXDM/3dMath.h SRange (line 474)

    Attributes:
        low: Lower bound of angular range (radians)
        high: Upper bound of angular range (radians)
    """
    low: float
    high: float

    def contains(self, angle: float) -> bool:
        """Check if angle is within this range.

        Args:
            angle: Rotation angle in radians

        Returns:
            True if low <= angle <= high

        C++ Reference: 3dMath.h SRange::Contains (line 480)
        """
        return self.low <= angle <= self.high

    def width(self) -> float:
        """Calculate range width in radians.

        Returns:
            Width of angular range (high - low)
        """
        return self.high - self.low


@dataclass
class FileRange:
    """File number range for detector images.

    Associates detector image files with angular wedges. Each detector
    has a range of file numbers corresponding to collected images.

    C++ Reference: XDM++/libXDM/3dMath.h SIntRange (line 459)

    Attributes:
        low: Lower bound of file number range (inclusive)
        high: Upper bound of file number range (inclusive)
    """
    low: int
    high: int

    def contains(self, file_num: int) -> bool:
        """Check if file number is within this range.

        Args:
            file_num: Detector image file number

        Returns:
            True if low <= file_num <= high

        C++ Reference: 3dMath.h SIntRange::Contains (line 465)
        """
        return self.low <= file_num <= self.high


class SimulationRange:
    """Maps omega angles to wedge indices and file numbers.

    This class is CRITICAL for forward simulation and reconstruction.
    It determines whether a predicted omega angle (from Bragg condition)
    falls within experimental data collection wedges.

    **Why This Matters**:
    - Forward simulation must only generate peaks for omegas within wedges
    - Reconstruction must filter reflections whose omegas are unmeasured
    - Without this, simulations include spurious peaks at unmeasured angles

    **Algorithm**:
    - Discretizes overall angular range into uniform bins
    - Creates lookup table mapping bin index -> wedge index (or None for gaps)
    - Provides O(1) lookup for "is this omega observable?"

    C++ Reference: Src/SimulationData.h CSimulationRange (lines 288-565)

    Attributes:
        low: Overall minimum omega angle (radians)
        high: Overall maximum omega angle (radians)
        width: Width of each angular interval (radians)
        num_intervals: Total number of uniform bins
        start_file_num: Starting file number offset
        stop_file_num: Ending file number
        range_list: List of experimental omega wedges
        index_list: Lookup table (bin index -> wedge index or None)
    """

    def __init__(
        self,
        low: float,
        high: float,
        width: float,
        range_list: List[OmegaRange],
        start_file_num: int = 0
    ):
        """Initialize range mapper.

        Args:
            low: Overall minimum omega angle (radians)
            high: Overall maximum omega angle (radians)
            width: Width of each angular interval (radians)
            range_list: List of experimental omega wedges
            start_file_num: Starting file number offset (default: 0)

        Raises:
            ValueError: If range_list is empty or width <= 0
        """
        if not range_list:
            raise ValueError("range_list cannot be empty")
        if width <= 0:
            raise ValueError(f"width must be positive, got {width}")

        self.low = low
        self.high = high
        self.width = width

        # Handle backwards ranges (high < low)
        # C++ Reference: SimulationData.h Set() (lines 363-364)
        # C++ Code:
        #   if( h < l )
        #       fWidth = -fWidth;
        if high < low:
            self.width = -width

        self.num_intervals = round((high - low) / self.width)
        self.start_file_num = start_file_num
        self.stop_file_num = start_file_num + self.num_intervals - 1
        self.range_list = range_list

        # Build lookup table: bin index -> wedge index (or None)
        self.index_list: List[Optional[int]] = self._build_index_list(range_list)

    def angle_to_index(self, angle: float) -> int:
        """Convert omega angle to uniform bin index.

        Discretizes the overall angular range into uniform bins.

        Args:
            angle: Omega angle in radians

        Returns:
            Bin index (0 to num_intervals-1), or -1 if angle < low

        C++ Reference: SimulationData.h AngleToIndex (line 317)

        C++ Code:
            Float f = (fAngle - fLow) / fWidth;
            if (f < 0)
                return -1;
            else
                return Int(f);  // Truncation
        """
        f = (angle - self.low) / self.width
        if f < 0:
            return -1
        else:
            return int(f)  # Truncation (floor for positive)

    def to_file_number(self, angle: float) -> Optional[int]:
        """Map omega angle to file number.

        Converts an omega angle to the corresponding detector image
        file number. Returns None if angle is outside overall range.

        Args:
            angle: Omega angle in radians

        Returns:
            File number, or None if angle outside range

        C++ Reference: SimulationData.h ToFileNumber (line 406)

        C++ Code:
            Int n = AngleToIndex(fAngle);
            if (n > nStopFileNum || n < 0)
                return NoMatch;
            else
                return n + nStartFileNum;
        """
        n = self.angle_to_index(angle)
        if n > self.stop_file_num or n < 0:
            return None  # C++ NoMatch
        else:
            return n + self.start_file_num

    def angle_to_wedge_index(self, angle: float) -> Optional[int]:
        """Map omega angle to experimental wedge index.

        **CRITICAL METHOD**: Determines if an omega angle is within an
        experimental data collection wedge. Used to filter reflections
        during forward simulation and reconstruction.

        Args:
            angle: Omega angle in radians

        Returns:
            Wedge index (0 to len(range_list)-1), or None if:
            - Angle is outside overall range, OR
            - Angle is in a gap between wedges

        C++ Reference: SimulationData.h operator() (line 446)

        C++ Code:
            Int n = AngleToIndex(fAngle);
            if (n >= nNumIntervals || n < 0)
                return NoMatch;
            else
                return vIndexList[n];  // May be NoMatch if gap
        """
        n = self.angle_to_index(angle)
        if n < 0 or n >= self.num_intervals:
            return None  # Outside overall range
        else:
            return self.index_list[n]  # May be None if in gap

    def is_in_experimental_range(self, angle: float) -> bool:
        """Check if omega angle is within experimental wedges.

        Convenience method for forward simulation filtering.

        **Use this to filter reflections**:
            if exp_setup.is_omega_observable(omega):
                # Only simulate peaks for measured omegas
                simulate_peak(omega, ...)

        Args:
            angle: Omega angle in radians

        Returns:
            True if angle is within any experimental wedge
        """
        return self.angle_to_wedge_index(angle) is not None

    def get_wedge(self, wedge_idx: int) -> OmegaRange:
        """Get omega range for a specific wedge index.

        Args:
            wedge_idx: Wedge index (0 to len(range_list)-1)

        Returns:
            OmegaRange for this wedge

        Raises:
            IndexError: If wedge_idx is out of bounds
        """
        return self.range_list[wedge_idx]

    def index_to_interval(self, index: int) -> OmegaRange:
        """Get angular range for a specific bin index.

        Reconstructs the angular interval for a uniform bin.

        Args:
            index: Bin index (0 to num_intervals-1)

        Returns:
            OmegaRange corresponding to this bin

        C++ Reference: SimulationData.h IndexToInterval (line 420)
        """
        low_angle = self.low + index * self.width
        high_angle = low_angle + self.width
        return OmegaRange(low=low_angle, high=high_angle)

    def _build_index_list(self, range_list: List[OmegaRange]) -> List[Optional[int]]:
        """Build lookup table mapping bin indices to wedge indices.

        **IMPORTANT LIMITATION**: The C++ implementation only marks the CENTER
        bin of each wedge. If a wedge spans multiple bins, only the center
        bin is marked. This means `angle_to_wedge_index()` may return None
        for angles near wedge boundaries even though they are technically
        within the wedge.

        For proper forward simulation filtering, wedges should be comparable
        in width to bin width, or this method should be enhanced to mark
        ALL bins within each wedge (potential Phase 2 improvement).

        Args:
            range_list: List of experimental omega wedges

        Returns:
            List of length num_intervals, where each entry is:
            - wedge_idx (int) if bin is at the CENTER of a wedge
            - None if bin is not a wedge center

        C++ Reference: SimulationData.h Set() (line 362)

        C++ Code:
            for (Size_Type i = 0; i < vRange.size(); i++) {
                Float fMid = (vRange[i].fHigh + vRange[i].fLow) / 2.0;
                Int nIndex = AngleToIndex(fMid);
                if (nIndex >= 0 && nIndex < nNumIntervals) {
                    vIndexList[nIndex] = i;
                }
            }
        """
        # Initialize all bins to None (no wedge)
        index_list: List[Optional[int]] = [None] * self.num_intervals

        # For each wedge, mark ONLY its center bin (matches C++ behavior)
        for wedge_idx, omega_range in enumerate(range_list):
            # Calculate center of wedge
            f_mid = (omega_range.high + omega_range.low) / 2.0
            n_index = self.angle_to_index(f_mid)

            # Mark this bin if valid
            if n_index >= 0 and n_index < self.num_intervals:
                index_list[n_index] = wedge_idx

        return index_list


def read_omega_file(
    filename: str,
    num_detectors: int
) -> Tuple[List[OmegaRange], List[FileRange]]:
    """Parse C++ omega range file.

    Omega files specify:
    1. File number ranges for each detector
    2. Angular wedges (in degrees) for data collection

    File format (binary .dat, parsed as text):
        Line 1: Header (ignored)
        Lines 2 to (1+num_detectors): File ranges as int pairs (nLow nHigh)
        Remaining lines: Omega ranges in DEGREES (fLow fHigh)

    **Important**: Omega angles are stored in DEGREES in file,
    but converted to RADIANS for internal use.

    C++ Reference: Src/InitFilesIO.cpp ReadRotationIntervalFiles (line 350)

    C++ Code:
        // Skip first line
        getline(iss, sBuf);

        // Read file ranges
        for (Size_Type i = 0; i < nNumFileRange; i++) {
            getline(iss, sBuf);
            vector<string> oTokens = Parser::Tokenize(sBuf);
            SIntRange s;
            s.nLow = atoi(oTokens[0].c_str());
            s.nHigh = atoi(oTokens[1].c_str());
            oFileRange.push_back(s);
        }

        // Read omega ranges (degrees -> radians)
        while (getline(iss, sBuf)) {
            vector<string> oTokens = Parser::Tokenize(sBuf);
            SRange s;
            s.fLow = DEGREE_TO_RADIAN(atof(oTokens[0].c_str()));
            s.fHigh = DEGREE_TO_RADIAN(atof(oTokens[1].c_str()));
            oRotationRange.push_back(s);
        }

    Args:
        filename: Path to omega .dat file
        num_detectors: Expected number of detectors

    Returns:
        Tuple of (omega_ranges, file_ranges) where:
        - omega_ranges: List of angular wedges in RADIANS
        - file_ranges: List of file number ranges per detector

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    try:
        # Read binary file and decode as UTF-8
        with open(filename, 'rb') as f:
            content = f.read().decode('utf-8', errors='ignore')
    except FileNotFoundError:
        raise FileNotFoundError(f"Omega file not found: {filename}")

    # Split into lines and remove empty lines
    lines = [line.strip() for line in content.split('\n') if line.strip()]

    if len(lines) < 1 + num_detectors:
        raise ValueError(
            f"Omega file too short: expected at least {1 + num_detectors} lines, "
            f"got {len(lines)}"
        )

    # Skip header (first line)
    lines = lines[1:]

    # Parse file ranges (next num_detectors lines)
    file_ranges = []
    for i in range(num_detectors):
        tokens = lines[i].split()
        if len(tokens) < 2:
            raise ValueError(
                f"Invalid file range at line {i+2}: expected 2 tokens, "
                f"got {len(tokens)}"
            )

        low = int(tokens[0])
        high = int(tokens[1])
        file_ranges.append(FileRange(low, high))

    # Parse omega ranges (remaining lines, degrees -> radians)
    omega_ranges = []
    for i, line in enumerate(lines[num_detectors:], start=num_detectors+2):
        tokens = line.split()
        if len(tokens) < 2:
            raise ValueError(
                f"Invalid omega range at line {i}: expected 2 tokens, "
                f"got {len(tokens)}"
            )

        low_deg = float(tokens[0])
        high_deg = float(tokens[1])

        # Convert degrees to radians
        omega_ranges.append(OmegaRange(
            low=np.deg2rad(low_deg),
            high=np.deg2rad(high_deg)
        ))

    return omega_ranges, file_ranges
