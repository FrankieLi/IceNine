"""
Tests for simulation_range module (omega range system).

Validates Python implementation against C++ CSimulationRange using
test data generated from C++ harness.
"""

import pytest
import numpy as np
import json
from pathlib import Path
from typing import Dict, Any, List

from icenine.simulation_range import (
    OmegaRange,
    FileRange,
    SimulationRange,
    read_omega_file
)


# Path to C++ generated test data
TEST_DATA_DIR = Path(__file__).parent.parent / "cpp_outputs"
CPP_TEST_DATA_FILE = TEST_DATA_DIR / "omega_range_test_data.json"


@pytest.fixture
def cpp_test_data() -> Dict[str, Any]:
    """Load C++ generated test data."""
    if not CPP_TEST_DATA_FILE.exists():
        pytest.skip(
            f"C++ test data not found: {CPP_TEST_DATA_FILE}\n"
            "Run: cd icenine_py/cpp_harness && make test_omega_ranges"
        )

    with open(CPP_TEST_DATA_FILE, 'r') as f:
        return json.load(f)


# ============================================================================
# Test OmegaRange
# ============================================================================

class TestOmegaRange:
    """Test OmegaRange dataclass and methods."""

    def test_contains_inside_range(self):
        """Test contains() for angle inside range."""
        omega_range = OmegaRange(low=-1.0, high=1.0)
        assert omega_range.contains(0.0)
        assert omega_range.contains(0.5)
        assert omega_range.contains(-0.5)

    def test_contains_outside_range(self):
        """Test contains() for angle outside range."""
        omega_range = OmegaRange(low=-1.0, high=1.0)
        assert not omega_range.contains(-1.5)
        assert not omega_range.contains(1.5)

    def test_contains_boundary(self):
        """Test contains() at exact boundaries."""
        omega_range = OmegaRange(low=-1.0, high=1.0)
        assert omega_range.contains(-1.0)  # Low boundary
        assert omega_range.contains(1.0)   # High boundary

    def test_width_calculation(self):
        """Test width() method."""
        omega_range = OmegaRange(low=-1.0, high=1.0)
        assert np.isclose(omega_range.width(), 2.0)

        omega_range2 = OmegaRange(low=0.0, high=np.pi)
        assert np.isclose(omega_range2.width(), np.pi)

    def test_negative_range(self):
        """Test range with negative angles."""
        omega_range = OmegaRange(low=-np.pi, high=-np.pi/2)
        assert omega_range.contains(-3.0)
        assert not omega_range.contains(0.0)

    def test_wraparound_range(self):
        """Test range that could wrap around (though not handled specially)."""
        # Note: Current implementation doesn't handle wraparound
        omega_range = OmegaRange(low=np.pi/2, high=3*np.pi/2)
        assert omega_range.contains(np.pi)
        assert not omega_range.contains(0.0)


# ============================================================================
# Test FileRange
# ============================================================================

class TestFileRange:
    """Test FileRange dataclass and methods."""

    def test_contains_inside(self):
        """Test contains() for file number inside range."""
        file_range = FileRange(low=0, high=100)
        assert file_range.contains(0)
        assert file_range.contains(50)
        assert file_range.contains(100)

    def test_contains_outside(self):
        """Test contains() for file number outside range."""
        file_range = FileRange(low=0, high=100)
        assert not file_range.contains(-1)
        assert not file_range.contains(101)

    def test_contains_boundary(self):
        """Test contains() at exact boundaries."""
        file_range = FileRange(low=10, high=20)
        assert file_range.contains(10)  # Low boundary
        assert file_range.contains(20)  # High boundary

    def test_zero_based_range(self):
        """Test zero-based file range."""
        file_range = FileRange(low=0, high=10)
        assert file_range.contains(0)
        assert file_range.contains(10)


# ============================================================================
# Test SimulationRange Basics
# ============================================================================

class TestSimulationRangeBasics:
    """Test basic SimulationRange functionality."""

    def test_initialization(self):
        """Test SimulationRange initialization."""
        omega_ranges = [OmegaRange(low=-1.0, high=1.0)]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        assert sim_range.low == -1.0
        assert sim_range.high == 1.0
        assert sim_range.width == 0.1
        assert sim_range.num_intervals == 20
        assert len(sim_range.range_list) == 1

    def test_angle_to_index_positive(self):
        """Test angle_to_index() for positive angles."""
        omega_ranges = [OmegaRange(low=0.0, high=1.0)]
        sim_range = SimulationRange(
            low=0.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        assert sim_range.angle_to_index(0.0) == 0
        assert sim_range.angle_to_index(0.5) == 5
        assert sim_range.angle_to_index(0.95) == 9

    def test_angle_to_index_negative(self):
        """Test angle_to_index() for negative angles."""
        omega_ranges = [OmegaRange(low=-1.0, high=1.0)]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        assert sim_range.angle_to_index(-1.0) == 0
        assert sim_range.angle_to_index(0.0) == 10
        # int(0.9 / 0.1) = int(9.0) but (0.9 - (-1.0)) / 0.1 = 19.0 → int(19) but actually 18 due to rounding
        assert sim_range.angle_to_index(0.9) == 18  # (0.9 - (-1.0)) / 0.1 = 1.9 / 0.1 = 19.0, but truncation gives 18

    def test_angle_to_index_out_of_range(self):
        """Test angle_to_index() for angles outside overall range."""
        omega_ranges = [OmegaRange(low=0.0, high=1.0)]
        sim_range = SimulationRange(
            low=0.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        assert sim_range.angle_to_index(-0.5) == -1
        assert sim_range.angle_to_index(-1.0) == -1

    def test_num_intervals_calculation(self):
        """Test that num_intervals is calculated correctly."""
        omega_ranges = [OmegaRange(low=-np.pi, high=np.pi)]
        sim_range = SimulationRange(
            low=-np.pi, high=np.pi, width=np.pi/180,  # 1 degree bins
            range_list=omega_ranges
        )

        # 360 degrees total
        assert sim_range.num_intervals == 360

    def test_file_number_mapping(self):
        """Test to_file_number() basic functionality."""
        omega_ranges = [OmegaRange(low=0.0, high=1.0)]
        sim_range = SimulationRange(
            low=0.0, high=1.0, width=0.1,
            range_list=omega_ranges,
            start_file_num=100
        )

        # File number = index + start_file_num
        assert sim_range.to_file_number(0.0) == 100  # index 0 + 100
        assert sim_range.to_file_number(0.5) == 105  # index 5 + 100

        # Out of range returns None
        assert sim_range.to_file_number(-0.5) is None


# ============================================================================
# Test SimulationRange Wedge Lookup (CRITICAL)
# ============================================================================

class TestSimulationRangeWedgeLookup:
    """Test wedge lookup functionality (critical for forward simulation)."""

    def test_single_wedge_all_angles_observable(self):
        """Test single wedge covering entire range - all angles observable."""
        omega_ranges = [OmegaRange(low=-1.0, high=1.0)]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        # Center of range should be observable (wedge 0)
        assert sim_range.angle_to_wedge_index(0.0) == 0
        assert sim_range.is_in_experimental_range(0.0)

    def test_multiple_wedges_with_gaps(self):
        """Test multiple wedges with gaps - some angles not observable."""
        # Wedges: [-1.0, -0.5] and [0.5, 1.0]
        # Gap: (-0.5, 0.5)
        omega_ranges = [
            OmegaRange(low=-1.0, high=-0.5),
            OmegaRange(low=0.5, high=1.0)
        ]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        # Angles in wedges should be observable
        assert sim_range.angle_to_wedge_index(-0.75) == 0
        assert sim_range.angle_to_wedge_index(0.75) == 1

        # Angles in gap should NOT be observable
        assert sim_range.angle_to_wedge_index(0.0) is None

    def test_angle_in_gap_returns_none(self):
        """Test that angles in gaps return None for wedge index."""
        omega_ranges = [
            OmegaRange(low=-1.0, high=-0.5),
            OmegaRange(low=0.5, high=1.0)
        ]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        # Gap between wedges
        assert sim_range.angle_to_wedge_index(0.0) is None
        assert sim_range.angle_to_wedge_index(-0.25) is None
        assert sim_range.angle_to_wedge_index(0.25) is None

    def test_angle_outside_overall_range(self):
        """Test angles outside overall range return None."""
        omega_ranges = [OmegaRange(low=-1.0, high=1.0)]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        assert sim_range.angle_to_wedge_index(-1.5) is None
        assert sim_range.angle_to_wedge_index(1.5) is None

    def test_is_in_experimental_range_true(self):
        """Test is_in_experimental_range() for observable angles.

        NOTE: Due to C++ limitation, only the CENTER bin of each wedge is marked.
        So only the center angle (0.0) is observable, not -0.5 or 0.5.
        """
        omega_ranges = [OmegaRange(low=-1.0, high=1.0)]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        # Only center of wedge is marked (bin 10 = angle 0.0)
        assert sim_range.is_in_experimental_range(0.0)
        # Other angles not marked even though technically in wedge
        assert not sim_range.is_in_experimental_range(-0.5)
        assert not sim_range.is_in_experimental_range(0.5)

    def test_is_in_experimental_range_false(self):
        """Test is_in_experimental_range() for non-observable angles."""
        omega_ranges = [OmegaRange(low=-1.0, high=-0.5)]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        # Outside wedge
        assert not sim_range.is_in_experimental_range(0.0)
        assert not sim_range.is_in_experimental_range(0.5)

        # Outside overall range
        assert not sim_range.is_in_experimental_range(-1.5)

    def test_boundary_angles(self):
        """Test wedge lookup at boundaries."""
        omega_ranges = [OmegaRange(low=-0.5, high=0.5)]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        # Center should be in wedge
        assert sim_range.angle_to_wedge_index(0.0) == 0

    def test_get_wedge_by_index(self):
        """Test get_wedge() accessor."""
        omega_ranges = [
            OmegaRange(low=-1.0, high=-0.5),
            OmegaRange(low=0.5, high=1.0)
        ]
        sim_range = SimulationRange(
            low=-1.0, high=1.0, width=0.1,
            range_list=omega_ranges
        )

        wedge0 = sim_range.get_wedge(0)
        assert wedge0.low == -1.0
        assert wedge0.high == -0.5

        wedge1 = sim_range.get_wedge(1)
        assert wedge1.low == 0.5
        assert wedge1.high == 1.0


# ============================================================================
# Test C++ Validation
# ============================================================================

class TestSimulationRangeCppValidation:
    """Validate Python implementation against C++ test data."""

    def test_single_wedge_cpp_match(self, cpp_test_data):
        """Test single wedge case matches C++ exactly.

        NOTE: We compute angles from degrees directly to avoid JSON precision loss.
        """
        test_case = cpp_test_data["test_cases"][0]
        assert test_case["name"] == "single_wedge_full_range"

        # Create Python SimulationRange with same config
        config = test_case["config"]
        omega_ranges = [
            OmegaRange(low=r["low"], high=r["high"])
            for r in config["range_list"]
        ]
        sim_range = SimulationRange(
            low=config["low"],
            high=config["high"],
            width=config["width"],
            range_list=omega_ranges
        )

        # Test each test case
        for test in test_case["tests"]:
            # Compute angle from degrees to match C++ precision
            angle = np.deg2rad(test["angle_deg"])
            expected_file = test["expected_file"]
            expected_wedge = test["expected_wedge"]

            # Check file number
            py_file = sim_range.to_file_number(angle)
            assert py_file == expected_file, \
                f"File number mismatch at {test['angle_deg']}°: " \
                f"Python={py_file}, C++={expected_file}"

            # Check wedge index
            py_wedge = sim_range.angle_to_wedge_index(angle)
            assert py_wedge == expected_wedge, \
                f"Wedge index mismatch at {test['angle_deg']}°: " \
                f"Python={py_wedge}, C++={expected_wedge}"

    def test_multiple_wedges_cpp_match(self, cpp_test_data):
        """Test multiple wedges case matches C++ exactly.

        NOTE: We compute angles from degrees directly to avoid JSON precision loss.
        """
        test_case = cpp_test_data["test_cases"][1]
        assert test_case["name"] == "multiple_wedges_with_gaps"

        # Create Python SimulationRange with same config
        config = test_case["config"]
        omega_ranges = [
            OmegaRange(low=r["low"], high=r["high"])
            for r in config["range_list"]
        ]
        sim_range = SimulationRange(
            low=config["low"],
            high=config["high"],
            width=config["width"],
            range_list=omega_ranges
        )

        # Test each test case
        for test in test_case["tests"]:
            # Compute angle from degrees to match C++ precision
            angle = np.deg2rad(test["angle_deg"])
            expected_file = test["expected_file"]
            expected_wedge = test["expected_wedge"]

            # Check file number
            py_file = sim_range.to_file_number(angle)
            assert py_file == expected_file, \
                f"File number mismatch at {test['angle_deg']}°: " \
                f"Python={py_file}, C++={expected_file}"

            # Check wedge index (CRITICAL for forward simulation)
            py_wedge = sim_range.angle_to_wedge_index(angle)
            assert py_wedge == expected_wedge, \
                f"Wedge index mismatch at {test['angle_deg']}°: " \
                f"Python={py_wedge}, C++={expected_wedge}"

    def test_fine_resolution_cpp_match(self, cpp_test_data):
        """Test fine resolution case matches C++ exactly.

        NOTE: We compute angles from degrees directly to avoid JSON precision loss.
        """
        test_case = cpp_test_data["test_cases"][2]
        assert test_case["name"] == "fine_angular_resolution"

        # Create Python SimulationRange with same config
        config = test_case["config"]
        omega_ranges = [
            OmegaRange(low=r["low"], high=r["high"])
            for r in config["range_list"]
        ]
        sim_range = SimulationRange(
            low=config["low"],
            high=config["high"],
            width=config["width"],
            range_list=omega_ranges
        )

        # Test each test case
        for test in test_case["tests"]:
            # Compute angle from degrees to match C++ precision
            angle = np.deg2rad(test["angle_deg"])
            expected_wedge = test["expected_wedge"]

            # Check wedge index
            py_wedge = sim_range.angle_to_wedge_index(angle)
            assert py_wedge == expected_wedge, \
                f"Wedge index mismatch at {test['angle_deg']}°: " \
                f"Python={py_wedge}, C++={expected_wedge}"

    def test_experimental_omega_5000_random(self, cpp_test_data):
        """Test experimental omega file with 5000 random angles.

        This test validates against actual experimental data from
        omega_180_2L.dat with 180 omega ranges and 5000 random
        test angles at 1° resolution.
        """
        test_case = cpp_test_data["test_cases"][3]
        assert test_case["name"] == "experimental_omega_5000_random"

        # Create Python SimulationRange with same config
        config = test_case["config"]
        omega_ranges = [
            OmegaRange(low=r["low"], high=r["high"])
            for r in config["range_list"]
        ]
        sim_range = SimulationRange(
            low=config["low"],
            high=config["high"],
            width=config["width"],
            range_list=omega_ranges
        )

        # Validate configuration
        assert len(omega_ranges) == 180, "Expected 180 omega ranges from experimental file"
        assert len(test_case["tests"]) == 5000, "Expected 5000 random test angles"

        # Test all 5000 random angles
        num_matches = 0
        num_gaps = 0

        for i, test in enumerate(test_case["tests"]):
            # Compute angle from degrees to match C++ precision
            # With setprecision(17), we have sufficient precision
            angle = np.deg2rad(test["angle_deg"])
            expected_file = test["expected_file"]
            expected_wedge = test["expected_wedge"]

            # Check file number
            py_file = sim_range.to_file_number(angle)
            assert py_file == expected_file, \
                f"File number mismatch at test {i} ({test['angle_deg']}°): " \
                f"Python={py_file}, C++={expected_file}"

            # Check wedge index (CRITICAL for forward simulation)
            py_wedge = sim_range.angle_to_wedge_index(angle)
            assert py_wedge == expected_wedge, \
                f"Wedge index mismatch at test {i} ({test['angle_deg']}°): " \
                f"Python={py_wedge}, C++={expected_wedge}"

            # Track statistics
            if expected_wedge is not None:
                num_matches += 1
            else:
                num_gaps += 1

            # Progress indicator for long test
            if (i + 1) % 1000 == 0:
                print(f"    Validated {i + 1}/5000 test angles...")

        # Report statistics
        # Note: With 101 omega ranges covering most of 80-180°, it's possible
        # all random angles fall in wedges (no gaps)
        print(f"  Statistics: {num_matches} in wedges, {num_gaps} in gaps")
        assert num_matches + num_gaps == 5000, "All angles should be accounted for"

    def test_all_cpp_cases_covered(self, cpp_test_data):
        """Verify we tested all C++ test cases."""
        assert len(cpp_test_data["test_cases"]) == 4
        for i, test_case in enumerate(cpp_test_data["test_cases"]):
            assert "name" in test_case
            assert "config" in test_case
            assert "tests" in test_case


# ============================================================================
# Test Omega File Parser (Mock Data)
# ============================================================================

class TestOmegaFileParser:
    """Test omega file parsing functionality."""

    def test_parse_error_handling(self):
        """Test error handling for invalid files."""
        with pytest.raises(FileNotFoundError):
            read_omega_file("nonexistent_file.dat", num_detectors=1)

    # Note: Full file parsing tests would require creating mock omega files
    # This is deferred for now since we have C++ validation of the core logic
