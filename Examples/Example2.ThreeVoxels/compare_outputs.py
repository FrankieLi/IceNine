#!/usr/bin/env python3
"""
Compare C++ and Python IceNine detector image outputs.

Generates statistical comparison and identifies discrepancies.

Usage:
    cd Examples/Example2.ThreeVoxels
    python compare_outputs.py
"""

import sys
from pathlib import Path
import numpy as np
import torch

# Add icenine_py to Python path
REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "icenine_py"))

from icenine.image_data import ImageData


def load_detector_image(filepath: Path) -> torch.Tensor:
    """Load detector image from ASCII file."""
    try:
        # 2048x2048 from StandardGeometry.2Det
        img = ImageData(num_rows=2048, num_cols=2048, mode='dense')
        img.load_ascii(str(filepath))
        return img._pixels_dense
    except Exception as e:
        print(f"  Error loading {filepath.name}: {e}")
        return None


def compare_images(cpp_data: torch.Tensor, py_data: torch.Tensor, filename: str) -> dict:
    """Compare two detector images and return metrics."""

    # Compute differences
    abs_diff = torch.abs(cpp_data - py_data)
    rel_diff = abs_diff / (torch.abs(cpp_data) + 1e-10)

    # Compute statistics
    metrics = {
        'filename': filename,
        'max_abs_diff': abs_diff.max().item(),
        'mean_abs_diff': abs_diff.mean().item(),
        'max_rel_diff': rel_diff.max().item(),
        'mean_rel_diff': rel_diff.mean().item(),
        'cpp_max': cpp_data.max().item(),
        'py_max': py_data.max().item(),
        'cpp_sum': cpp_data.sum().item(),
        'py_sum': py_data.sum().item(),
        'cpp_nonzero': (cpp_data > 0).sum().item(),
        'py_nonzero': (py_data > 0).sum().item(),
    }

    # Correlation (only if both have non-zero pixels)
    if metrics['cpp_nonzero'] > 0 and metrics['py_nonzero'] > 0:
        # Flatten and compute correlation
        cpp_flat = cpp_data.flatten()
        py_flat = py_data.flatten()
        correlation = torch.corrcoef(torch.stack([cpp_flat, py_flat]))[0, 1].item()
        metrics['correlation'] = correlation
    else:
        metrics['correlation'] = 0.0

    return metrics


def main():
    """Run comparison between C++ and Python outputs."""

    example_dir = Path(__file__).parent
    cpp_dir = example_dir / "ScatteringData"
    python_dir = example_dir / "ScatteringData_Python"

    print("=" * 80)
    print("IceNine C++ vs Python Output Comparison - Example2.ThreeVoxels")
    print("=" * 80)
    print(f"\nC++ output: {cpp_dir}")
    print(f"Python output: {python_dir}\n")

    # Check directories exist
    if not cpp_dir.exists():
        print(f"✗ Error: C++ output directory not found: {cpp_dir}")
        return 1

    if not python_dir.exists():
        print(f"✗ Error: Python output directory not found: {python_dir}")
        print("  Run run_python_simulation.py first!")
        return 1

    # Get file lists
    cpp_files = sorted(cpp_dir.glob("*.d*"))
    python_files = sorted(python_dir.glob("**/*.d*"))

    print(f"Found {len(cpp_files)} C++ output files")
    print(f"Found {len(python_files)} Python output files\n")

    if len(cpp_files) == 0:
        print("✗ No C++ files found!")
        return 1

    if len(python_files) == 0:
        print("✗ No Python files found!")
        return 1

    # Compare each file
    print("Comparing detector images...")
    print("-" * 80)

    results = []
    missing_python = []
    missing_cpp = []

    for cpp_file in cpp_files:
        # Find corresponding Python file
        python_file = python_dir / cpp_file.name

        if not python_file.exists():
            missing_python.append(cpp_file.name)
            continue

        # Load images
        cpp_data = load_detector_image(cpp_file)
        py_data = load_detector_image(python_file)

        if cpp_data is None or py_data is None:
            continue

        # Check dimensions match
        if cpp_data.shape != py_data.shape:
            print(f"⚠ {cpp_file.name}: Shape mismatch C++{cpp_data.shape} vs Py{py_data.shape}")
            continue

        # Compare
        metrics = compare_images(cpp_data, py_data, cpp_file.name)
        results.append(metrics)

    # Check for extra Python files
    cpp_names = {f.name for f in cpp_files}
    for py_file in python_files:
        if py_file.name not in cpp_names:
            missing_cpp.append(py_file.name)

    print(f"Compared {len(results)} file pairs\n")

    # Report missing files
    if missing_python:
        print(f"⚠ {len(missing_python)} files in C++ but not in Python:")
        for name in missing_python[:5]:
            print(f"    {name}")
        if len(missing_python) > 5:
            print(f"    ... and {len(missing_python) - 5} more")
        print()

    if missing_cpp:
        print(f"⚠ {len(missing_cpp)} files in Python but not in C++:")
        for name in missing_cpp[:5]:
            print(f"    {name}")
        if len(missing_cpp) > 5:
            print(f"    ... and {len(missing_cpp) - 5} more")
        print()

    if len(results) == 0:
        print("✗ No files were successfully compared!")
        return 1

    # Summary statistics
    print("=" * 80)
    print("Summary Statistics")
    print("=" * 80)

    max_abs_diffs = [r['max_abs_diff'] for r in results]
    mean_abs_diffs = [r['mean_abs_diff'] for r in results]
    max_rel_diffs = [r['max_rel_diff'] for r in results]
    correlations = [r['correlation'] for r in results]

    print(f"\nAbsolute Differences:")
    print(f"  Max across all files:  {max(max_abs_diffs):.6e}")
    print(f"  Mean max difference:   {np.mean(max_abs_diffs):.6e}")
    print(f"  Median max difference: {np.median(max_abs_diffs):.6e}")
    print(f"  Mean avg difference:   {np.mean(mean_abs_diffs):.6e}")

    print(f"\nRelative Differences:")
    print(f"  Max relative error:    {max(max_rel_diffs):.6e}")
    print(f"  Mean relative error:   {np.mean(max_rel_diffs):.6e}")

    print(f"\nCorrelation:")
    print(f"  Mean correlation:      {np.mean(correlations):.6f}")
    print(f"  Min correlation:       {min(correlations):.6f}")

    # Pass/fail criteria
    print("\n" + "=" * 80)
    print("Pass/Fail Criteria")
    print("=" * 80)

    THRESHOLDS = {
        'max_abs_diff': 1e-3,
        'mean_abs_diff': 1e-5,
        'max_rel_diff': 0.01,  # 1%
        'min_correlation': 0.99
    }

    passed = True

    if max(max_abs_diffs) < THRESHOLDS['max_abs_diff']:
        print(f"✓ Max absolute difference < {THRESHOLDS['max_abs_diff']}")
    else:
        print(f"✗ Max absolute difference >= {THRESHOLDS['max_abs_diff']}")
        passed = False

    if np.mean(mean_abs_diffs) < THRESHOLDS['mean_abs_diff']:
        print(f"✓ Mean absolute difference < {THRESHOLDS['mean_abs_diff']}")
    else:
        print(f"✗ Mean absolute difference >= {THRESHOLDS['mean_abs_diff']}")
        passed = False

    if max(max_rel_diffs) < THRESHOLDS['max_rel_diff']:
        print(f"✓ Max relative error < {THRESHOLDS['max_rel_diff']*100}%")
    else:
        print(f"✗ Max relative error >= {THRESHOLDS['max_rel_diff']*100}%")
        passed = False

    if min(correlations) > THRESHOLDS['min_correlation']:
        print(f"✓ Minimum correlation > {THRESHOLDS['min_correlation']}")
    else:
        print(f"✗ Minimum correlation <= {THRESHOLDS['min_correlation']}")
        passed = False

    print()

    if passed:
        print("=" * 80)
        print("✓ ALL TESTS PASSED - Python and C++ outputs match!")
        print("=" * 80)
        return 0
    else:
        print("=" * 80)
        print("✗ TESTS FAILED - Outputs differ beyond threshold")
        print("=" * 80)

        # Show worst offenders
        print("\nWorst 5 files by absolute difference:")
        worst = sorted(results, key=lambda r: r['max_abs_diff'], reverse=True)[:5]
        for r in worst:
            print(f"  {r['filename']:30s}  max_diff={r['max_abs_diff']:.6e}")

        return 1


if __name__ == "__main__":
    sys.exit(main())
