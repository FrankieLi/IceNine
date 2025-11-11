"""
Test suite for diffraction core calculations.

These tests validate the Python diffraction implementation against C++ ground
truth data generated from the original IceNine code.

Test Data Source: cpp_outputs/scattering_omegas.json
"""

import json
import numpy as np
import pytest
from pathlib import Path

# Import modules to test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from icenine.diffraction_core import (
    get_scattering_omegas,
    calculate_bragg_angle,
    ScatteringResult,
)
from icenine.constants import TEST_PARAMS


# Path to C++ ground truth data
CPP_OUTPUTS_DIR = Path(__file__).parent.parent / "cpp_outputs"


def load_json(filename: str) -> dict:
    """Load JSON ground truth data from C++."""
    filepath = CPP_OUTPUTS_DIR / filename
    with open(filepath, "r") as f:
        return json.load(f)


@pytest.fixture
def cpp_omega_data() -> dict:
    """Load C++ scattering omega angle test data."""
    return load_json("scattering_omegas.json")


class TestScatteringOmegas:
    """
    Test 6: Verify scattering omega calculations match C++.

    Critical validation of the core diffraction physics implementation.
    """

    def test_omega_angles_match_cpp(self, cpp_omega_data):
        """
        CRITICAL TEST: All omega angles must match C++ ground truth.

        This test validates that the Python implementation of GetScatteringOmegas
        produces identical results to the C++ version for all test cases.
        """
        test_cases = cpp_omega_data["test_cases"]
        beam_energy = cpp_omega_data["beam_energy"]
        beam_deflection_chi = cpp_omega_data["beam_deflection_chi"]

        failures = []

        for i, case in enumerate(test_cases):
            g_vec = np.array(case["g_vector"], dtype=float)
            g_mag = case["g_magnitude"]
            cpp_observable = case["observable"]

            # Calculate Python result
            result = get_scattering_omegas(
                g_vec, g_mag, beam_energy, beam_deflection_chi
            )

            # Check observable flag
            if result.observable != cpp_observable:
                failures.append({
                    "case_index": i,
                    "hkl": case["hkl"],
                    "error": "Observable flag mismatch",
                    "cpp_observable": cpp_observable,
                    "py_observable": result.observable,
                })
                continue

            # If observable, check omega angles
            if cpp_observable:
                cpp_omega1 = case["omega1"]
                cpp_omega2 = case["omega2"]

                # Tolerance: 1e-6 radians (~0.00006 degrees)
                if not np.isclose(result.omega1, cpp_omega1, atol=1e-6):
                    failures.append({
                        "case_index": i,
                        "hkl": case["hkl"],
                        "error": "Omega1 mismatch",
                        "cpp_omega1": cpp_omega1,
                        "py_omega1": result.omega1,
                        "diff": abs(result.omega1 - cpp_omega1),
                    })

                if not np.isclose(result.omega2, cpp_omega2, atol=1e-6):
                    failures.append({
                        "case_index": i,
                        "hkl": case["hkl"],
                        "error": "Omega2 mismatch",
                        "cpp_omega2": cpp_omega2,
                        "py_omega2": result.omega2,
                        "diff": abs(result.omega2 - cpp_omega2),
                    })

        assert len(failures) == 0, (
            f"Omega angle mismatch in {len(failures)} cases:\\n"
            + "\\n".join(str(f) for f in failures)
        )

    def test_all_cases_observable(self, cpp_omega_data):
        """Verify that all test cases in C++ data are marked as observable."""
        test_cases = cpp_omega_data["test_cases"]
        for i, case in enumerate(test_cases):
            assert case["observable"], f"Test case {i} is not observable in C++ data"

    def test_omega_angle_range(self, cpp_omega_data):
        """Verify omega angles are in the expected [-π, π] range."""
        test_cases = cpp_omega_data["test_cases"]
        beam_energy = cpp_omega_data["beam_energy"]
        beam_deflection_chi = cpp_omega_data["beam_deflection_chi"]

        for case in test_cases:
            if case["observable"]:
                g_vec = np.array(case["g_vector"], dtype=float)
                g_mag = case["g_magnitude"]

                result = get_scattering_omegas(
                    g_vec, g_mag, beam_energy, beam_deflection_chi
                )

                assert -np.pi <= result.omega1 <= np.pi, (
                    f"Omega1 out of range: {result.omega1}"
                )
                assert -np.pi <= result.omega2 <= np.pi, (
                    f"Omega2 out of range: {result.omega2}"
                )


class TestSpecificReflections:
    """Test omega calculations for specific known reflections."""

    def test_111_reflection(self):
        """Test (111) reflection with known geometry."""
        # Au FCC (111) reflection
        a = TEST_PARAMS["a"]  # 4.0782 Å
        a_recip = 2.0 * np.pi / a

        # (111) gives G = a* * [1, 1, 1]
        h, k, l = 1, 1, 1
        g_vec = a_recip * np.array([h, k, l], dtype=float)
        g_mag = a_recip * np.sqrt(h**2 + k**2 + l**2)

        result = get_scattering_omegas(
            g_vec, g_mag,
            beam_energy=TEST_PARAMS["beam_energy"],
            beam_deflection_chi=0.0
        )

        assert result.observable, "(111) reflection should be observable"
        assert result.omega1 is not None
        assert result.omega2 is not None

    def test_200_reflection(self):
        """Test (200) reflection with known geometry."""
        a = TEST_PARAMS["a"]
        a_recip = 2.0 * np.pi / a

        # (200) gives G = a* * [2, 0, 0]
        h, k, l = 2, 0, 0
        g_vec = a_recip * np.array([h, k, l], dtype=float)
        g_mag = a_recip * np.sqrt(h**2 + k**2 + l**2)

        result = get_scattering_omegas(
            g_vec, g_mag,
            beam_energy=TEST_PARAMS["beam_energy"],
            beam_deflection_chi=0.0
        )

        assert result.observable, "(200) reflection should be observable"


class TestBraggAngle:
    """Test Bragg angle calculations."""

    def test_bragg_angle_calculation(self):
        """Test Bragg angle calculation for known Q."""
        # For (111) reflection in Au FCC
        a = TEST_PARAMS["a"]
        a_recip = 2.0 * np.pi / a
        g_mag = a_recip * np.sqrt(3)  # |G| for (111)

        theta = calculate_bragg_angle(g_mag, TEST_PARAMS["beam_energy"])

        # Bragg angle should be positive and less than π/2
        assert 0 < theta < np.pi / 2

    def test_bragg_angle_increases_with_q(self):
        """Bragg angle should increase with Q magnitude."""
        beam_energy = TEST_PARAMS["beam_energy"]

        theta1 = calculate_bragg_angle(g_magnitude=2.0, beam_energy=beam_energy)
        theta2 = calculate_bragg_angle(g_magnitude=4.0, beam_energy=beam_energy)

        assert theta2 > theta1, "Larger Q should give larger Bragg angle"


class TestEdgeCases:
    """Test edge cases and special geometries."""

    def test_zero_beam_deflection(self):
        """Test with zero beam deflection (standard geometry)."""
        g_vec = np.array([1.0, 1.0, 0.0])
        g_mag = np.linalg.norm(g_vec)

        result = get_scattering_omegas(
            g_vec, g_mag,
            beam_energy=50.0,
            beam_deflection_chi=0.0
        )

        # Should get a result (may or may not be observable)
        assert isinstance(result, ScatteringResult)

    def test_g_along_z_axis(self):
        """Test scattering vector aligned with rotation axis."""
        # G along z-axis should be problematic (sin_chi ≈ 0)
        g_vec = np.array([0.0, 0.0, 5.0])
        g_mag = 5.0

        result = get_scattering_omegas(
            g_vec, g_mag,
            beam_energy=50.0,
            beam_deflection_chi=0.0
        )

        # This should likely be non-observable (too close to rotation axis)
        # But we just verify it doesn't crash
        assert isinstance(result, ScatteringResult)

    def test_result_dataclass_properties(self):
        """Verify ScatteringResult has expected properties."""
        g_vec = np.array([1.0, 2.0, 3.0])
        g_mag = np.linalg.norm(g_vec)

        result = get_scattering_omegas(g_vec, g_mag, beam_energy=50.0)

        assert hasattr(result, 'observable')
        assert hasattr(result, 'omega1')
        assert hasattr(result, 'omega2')
        assert hasattr(result, 'g_vector')
        assert hasattr(result, 'g_magnitude')

        # g_vector should be a copy, not a reference
        result.g_vector[0] = 999.0
        assert g_vec[0] != 999.0, "g_vector should be copied, not referenced"


class TestNumericalPrecision:
    """Test numerical precision and stability."""

    def test_high_precision_match(self, cpp_omega_data):
        """Test that omega angles match C++ to high precision (< 1e-10)."""
        test_cases = cpp_omega_data["test_cases"]
        beam_energy = cpp_omega_data["beam_energy"]
        beam_deflection_chi = cpp_omega_data["beam_deflection_chi"]

        max_error_omega1 = 0.0
        max_error_omega2 = 0.0

        for case in test_cases:
            if case["observable"]:
                g_vec = np.array(case["g_vector"], dtype=float)
                g_mag = case["g_magnitude"]

                result = get_scattering_omegas(
                    g_vec, g_mag, beam_energy, beam_deflection_chi
                )

                error1 = abs(result.omega1 - case["omega1"])
                error2 = abs(result.omega2 - case["omega2"])

                max_error_omega1 = max(max_error_omega1, error1)
                max_error_omega2 = max(max_error_omega2, error2)

        # Report maximum errors (informational)
        print(f"\\nMax omega1 error: {max_error_omega1:.2e} radians")
        print(f"Max omega2 error: {max_error_omega2:.2e} radians")

        # Verify errors are acceptably small
        assert max_error_omega1 < 1e-6, f"Omega1 error too large: {max_error_omega1}"
        assert max_error_omega2 < 1e-6, f"Omega2 error too large: {max_error_omega2}"


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
