"""
Comprehensive validation tests for MIC file I/O against all available test files.

This test suite validates the Python implementation against all .mic files in
the repository to ensure broad compatibility with the C++ implementation.
"""

import pytest
import numpy as np
import tempfile
import shutil
from pathlib import Path
from typing import List, Dict, Tuple

from icenine.mic_file import MicFile, euler_to_matrix, matrix_to_euler


# ==========================================================================================
# Test Discovery
# ==========================================================================================


def find_all_mic_files(base_dir: Path) -> List[Path]:
    """Find all .mic files in the repository."""
    mic_files = []

    # Check DataFiles.back/
    datafiles_back = base_dir / "DataFiles.back"
    if datafiles_back.exists():
        mic_files.extend(datafiles_back.glob("*.mic"))

    # Check DataFiles/
    datafiles = base_dir / "DataFiles"
    if datafiles.exists():
        mic_files.extend(datafiles.glob("*.mic"))

    return sorted(mic_files)


@pytest.fixture(scope="module")
def repo_root():
    """Get repository root directory."""
    # Assume tests are in icenine_py/tests/
    return Path(__file__).parent.parent.parent


@pytest.fixture(scope="module")
def all_mic_files(repo_root):
    """Find all .mic files in the repository."""
    return find_all_mic_files(repo_root)


@pytest.fixture
def temp_dir():
    """Create temporary directory for test outputs."""
    temp = tempfile.mkdtemp()
    yield Path(temp)
    shutil.rmtree(temp)


# ==========================================================================================
# Validation Helpers
# ==========================================================================================


def validate_rotation_matrix(R: np.ndarray, tolerance: float = 1e-4) -> Tuple[bool, str]:
    """
    Validate that a matrix is a proper rotation matrix.

    Returns:
        (is_valid, error_message)
    """
    # Check shape
    if R.shape != (3, 3):
        return False, f"Wrong shape: {R.shape}, expected (3, 3)"

    # Check orthogonality: R^T @ R = I
    identity_check = R.T @ R
    if not np.allclose(identity_check, np.eye(3), atol=tolerance):
        max_error = np.max(np.abs(identity_check - np.eye(3)))
        return False, f"Not orthogonal (max error: {max_error:.2e})"

    # Check determinant = 1 (proper rotation, not reflection)
    det = np.linalg.det(R)
    if not np.isclose(det, 1.0, atol=tolerance):
        return False, f"Determinant {det:.6f} != 1.0"

    return True, ""


def validate_voxel_data(mic: MicFile) -> Dict[str, any]:
    """
    Validate voxel data integrity.

    Returns:
        Dictionary with validation results
    """
    results = {
        "valid": True,
        "errors": [],
        "warnings": [],
        "n_voxels": len(mic.voxels),
        "n_invalid_rotations": 0,
        "n_invalid_positions": 0,
    }

    for i, voxel in enumerate(mic.voxels):
        # Check rotation matrix
        is_valid, error_msg = validate_rotation_matrix(voxel.orientation)
        if not is_valid:
            results["valid"] = False
            results["errors"].append(f"Voxel {i}: Invalid rotation - {error_msg}")
            results["n_invalid_rotations"] += 1

        # Check position is finite
        if not np.all(np.isfinite(voxel.position)):
            results["valid"] = False
            results["errors"].append(f"Voxel {i}: Non-finite position {voxel.position}")
            results["n_invalid_positions"] += 1

        # Check generation is reasonable
        if voxel.generation < 0 or voxel.generation > 20:
            results["warnings"].append(
                f"Voxel {i}: Unusual generation {voxel.generation}"
            )

        # Check phase is reasonable
        if voxel.phase < 0 or voxel.phase > 10:
            results["warnings"].append(f"Voxel {i}: Unusual phase {voxel.phase}")

        # Check confidence is in [0, 1]
        if voxel.confidence < 0.0 or voxel.confidence > 1.0:
            results["warnings"].append(
                f"Voxel {i}: Confidence {voxel.confidence:.3f} out of [0, 1] range"
            )

    return results


# ==========================================================================================
# File Reading Tests
# ==========================================================================================


class TestFileReading:
    """Test reading all available .mic files."""

    def test_read_all_files(self, all_mic_files):
        """Test that all .mic files can be read without errors."""
        if len(all_mic_files) == 0:
            pytest.skip("No .mic files found in repository")

        results = []
        for mic_file in all_mic_files:
            try:
                mic = MicFile.read(str(mic_file))
                results.append({
                    "file": mic_file.name,
                    "success": True,
                    "n_voxels": len(mic.voxels),
                    "side_length": mic.initial_side_length,
                    "error": None,
                })
            except Exception as e:
                results.append({
                    "file": mic_file.name,
                    "success": False,
                    "n_voxels": 0,
                    "side_length": 0.0,
                    "error": str(e),
                })

        # Print summary
        successful = sum(1 for r in results if r["success"])
        total = len(results)
        print(f"\n{'='*70}")
        print(f"MIC File Reading Summary: {successful}/{total} successful")
        print(f"{'='*70}")

        for r in results:
            status = "✓" if r["success"] else "✗"
            if r["success"]:
                print(f"{status} {r['file']:40s} - {r['n_voxels']:5d} voxels")
            else:
                print(f"{status} {r['file']:40s} - FAILED: {r['error']}")

        # All files should be readable
        failures = [r for r in results if not r["success"]]
        if failures:
            pytest.fail(
                f"{len(failures)}/{total} files failed to read:\n"
                + "\n".join(f"  - {r['file']}: {r['error']}" for r in failures)
            )

    def test_validate_all_files(self, all_mic_files):
        """Test that all voxel data in .mic files is valid."""
        if len(all_mic_files) == 0:
            pytest.skip("No .mic files found in repository")

        validation_results = []

        for mic_file in all_mic_files:
            try:
                mic = MicFile.read(str(mic_file))
                validation = validate_voxel_data(mic)
                validation_results.append({
                    "file": mic_file.name,
                    "validation": validation,
                })
            except Exception as e:
                # Skip files that can't be read (handled by previous test)
                continue

        # Print validation summary
        print(f"\n{'='*70}")
        print("MIC File Validation Summary")
        print(f"{'='*70}")

        for r in validation_results:
            v = r["validation"]
            if v["valid"]:
                status = "✓"
                msg = f"{v['n_voxels']} voxels, all valid"
            else:
                status = "✗"
                msg = f"{len(v['errors'])} errors, {len(v['warnings'])} warnings"

            print(f"{status} {r['file']:40s} - {msg}")

            # Print errors if any
            if v["errors"]:
                for error in v["errors"][:3]:  # Show first 3 errors
                    print(f"    ERROR: {error}")
                if len(v["errors"]) > 3:
                    print(f"    ... and {len(v['errors']) - 3} more errors")

        # Check for invalid files
        invalid_files = [r for r in validation_results if not r["validation"]["valid"]]
        if invalid_files:
            pytest.fail(
                f"{len(invalid_files)}/{len(validation_results)} files have invalid data"
            )


# ==========================================================================================
# Round-Trip Tests
# ==========================================================================================


class TestRoundTrip:
    """Test round-trip (read → write → read) for all files."""

    @pytest.mark.parametrize(
        "test_files",
        [
            ["Au1007_small.mic"],  # Known good file
            ["empty.mic"],  # Edge case
            ["oneTriangle2.mic"],  # Small file
        ],
    )
    def test_roundtrip_specific_files(self, repo_root, test_files, temp_dir):
        """Test round-trip on specific known files."""
        for filename in test_files:
            # Try both DataFiles and DataFiles.back
            mic_file = repo_root / "DataFiles" / filename
            if not mic_file.exists():
                mic_file = repo_root / "DataFiles.back" / filename

            if not mic_file.exists():
                pytest.skip(f"Test file not found: {filename}")
                continue

            # Read original
            mic_original = MicFile.read(str(mic_file))

            # Write to temp
            temp_file = temp_dir / f"roundtrip_{filename}"
            mic_original.write(str(temp_file))

            # Read back
            mic_roundtrip = MicFile.read(str(temp_file))

            # Compare
            assert len(mic_roundtrip.voxels) == len(mic_original.voxels), (
                f"Voxel count mismatch: {len(mic_roundtrip.voxels)} != "
                f"{len(mic_original.voxels)}"
            )

            assert mic_roundtrip.initial_side_length == pytest.approx(
                mic_original.initial_side_length, abs=1e-6
            )

            # Compare voxels
            for i, (v_orig, v_rt) in enumerate(
                zip(mic_original.voxels, mic_roundtrip.voxels)
            ):
                # Position should match closely
                np.testing.assert_allclose(
                    v_rt.position, v_orig.position, atol=1e-5,
                    err_msg=f"{filename} voxel {i}: position mismatch"
                )

                # Orientation matrix should match (tolerance for Euler conversion)
                np.testing.assert_allclose(
                    v_rt.orientation, v_orig.orientation, atol=1e-3,
                    err_msg=f"{filename} voxel {i}: orientation mismatch"
                )

                # Metadata should match exactly
                assert v_rt.generation == v_orig.generation
                assert v_rt.phase == v_orig.phase
                assert v_rt.points_up == v_orig.points_up

    def test_roundtrip_sample_of_all_files(self, all_mic_files, temp_dir):
        """Test round-trip on a sample of all available files."""
        if len(all_mic_files) == 0:
            pytest.skip("No .mic files found in repository")

        # Test a sample (e.g., every 3rd file, or all if fewer than 20)
        if len(all_mic_files) <= 20:
            sample_files = all_mic_files
        else:
            sample_files = all_mic_files[::3]  # Every 3rd file

        results = []

        for mic_file in sample_files:
            try:
                # Read original
                mic_original = MicFile.read(str(mic_file))

                # Write to temp
                temp_file = temp_dir / f"roundtrip_{mic_file.name}"
                mic_original.write(str(temp_file))

                # Read back
                mic_roundtrip = MicFile.read(str(temp_file))

                # Check voxel count preserved
                voxel_match = len(mic_roundtrip.voxels) == len(mic_original.voxels)

                # Check positions preserved
                if len(mic_roundtrip.voxels) > 0:
                    position_errors = []
                    for v_orig, v_rt in zip(mic_original.voxels, mic_roundtrip.voxels):
                        pos_diff = np.max(np.abs(v_rt.position - v_orig.position))
                        position_errors.append(pos_diff)
                    max_position_error = max(position_errors)
                    position_match = max_position_error < 1e-5
                else:
                    max_position_error = 0.0
                    position_match = True

                results.append({
                    "file": mic_file.name,
                    "success": True,
                    "voxel_match": voxel_match,
                    "position_match": position_match,
                    "max_position_error": max_position_error,
                    "error": None,
                })

            except Exception as e:
                results.append({
                    "file": mic_file.name,
                    "success": False,
                    "voxel_match": False,
                    "position_match": False,
                    "max_position_error": float("inf"),
                    "error": str(e),
                })

        # Print summary
        print(f"\n{'='*70}")
        print(f"Round-Trip Test Summary: {len(sample_files)} files tested")
        print(f"{'='*70}")

        for r in results:
            if r["success"] and r["voxel_match"] and r["position_match"]:
                status = "✓"
                msg = f"max position error: {r['max_position_error']:.2e}"
            else:
                status = "✗"
                if not r["success"]:
                    msg = f"FAILED: {r['error']}"
                elif not r["voxel_match"]:
                    msg = "Voxel count mismatch"
                else:
                    msg = f"Position mismatch (error: {r['max_position_error']:.2e})"

            print(f"{status} {r['file']:40s} - {msg}")

        # Check for failures
        failures = [r for r in results if not (r["success"] and r["voxel_match"] and r["position_match"])]
        if failures:
            pytest.fail(
                f"{len(failures)}/{len(results)} files failed round-trip test"
            )


# ==========================================================================================
# File Statistics
# ==========================================================================================


class TestFileStatistics:
    """Gather statistics about .mic files in the repository."""

    def test_file_size_distribution(self, all_mic_files):
        """Report statistics on file sizes."""
        if len(all_mic_files) == 0:
            pytest.skip("No .mic files found in repository")

        sizes = []
        for mic_file in all_mic_files:
            try:
                mic = MicFile.read(str(mic_file))
                sizes.append(len(mic.voxels))
            except:
                pass

        if len(sizes) == 0:
            pytest.skip("No files could be read")

        print(f"\n{'='*70}")
        print("MIC File Size Distribution")
        print(f"{'='*70}")
        print(f"Total files analyzed: {len(sizes)}")
        print(f"Total voxels: {sum(sizes)}")
        print(f"Min voxels: {min(sizes)}")
        print(f"Max voxels: {max(sizes)}")
        print(f"Mean voxels: {np.mean(sizes):.1f}")
        print(f"Median voxels: {np.median(sizes):.1f}")

        # Histogram
        bins = [0, 1, 10, 100, 1000, 10000, max(sizes) + 1]
        hist, _ = np.histogram(sizes, bins=bins)

        print("\nSize distribution:")
        for i in range(len(hist)):
            if hist[i] > 0:
                print(f"  {bins[i]:6d} - {bins[i+1]:6d} voxels: {hist[i]:3d} files")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
