#!/usr/bin/env python3
"""
Compare C++ and Python IceNine detector image outputs for Example2.ManyGrains.

Generates statistical comparison and identifies discrepancies.

Usage:
    cd Examples/Example2.ManyGrains
    uv run python compare_outputs.py
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
        img = ImageData(num_rows=2048, num_cols=2048, mode='dense')
        img.load_ascii(str(filepath))
        return img._pixels_dense
    except Exception as e:
        print(f"  Error loading {filepath.name}: {e}")
        return None


def compare_images(cpp_data: torch.Tensor, py_data: torch.Tensor, filename: str) -> dict:
    """Compare two detector images and return metrics."""

    abs_diff = torch.abs(cpp_data - py_data)

    cpp_nonzero = (cpp_data > 0).sum().item()
    py_nonzero = (py_data > 0).sum().item()

    # Pixel overlap analysis
    cpp_mask = cpp_data > 0
    py_mask = py_data > 0
    both_nonzero = (cpp_mask & py_mask).sum().item()
    cpp_only = (cpp_mask & ~py_mask).sum().item()
    py_only = (~cpp_mask & py_mask).sum().item()

    # Relative diff only where both are nonzero
    if both_nonzero > 0:
        overlap_mask = cpp_mask & py_mask
        rel_diff_at_overlap = (abs_diff[overlap_mask] / (torch.abs(cpp_data[overlap_mask]) + 1e-10))
        max_rel_diff = rel_diff_at_overlap.max().item()
        mean_rel_diff = rel_diff_at_overlap.mean().item()
    else:
        max_rel_diff = 0.0
        mean_rel_diff = 0.0

    metrics = {
        'filename': filename,
        'max_abs_diff': abs_diff.max().item(),
        'mean_abs_diff': abs_diff.mean().item(),
        'max_rel_diff': max_rel_diff,
        'mean_rel_diff': mean_rel_diff,
        'cpp_max': cpp_data.max().item(),
        'py_max': py_data.max().item(),
        'cpp_sum': cpp_data.sum().item(),
        'py_sum': py_data.sum().item(),
        'cpp_nonzero': cpp_nonzero,
        'py_nonzero': py_nonzero,
        'both_nonzero': both_nonzero,
        'cpp_only': cpp_only,
        'py_only': py_only,
    }

    # Correlation (only if both have non-zero pixels)
    if cpp_nonzero > 0 and py_nonzero > 0:
        cpp_flat = cpp_data.flatten()
        py_flat = py_data.flatten()
        correlation = torch.corrcoef(torch.stack([cpp_flat, py_flat]))[0, 1].item()
        metrics['correlation'] = correlation
    else:
        metrics['correlation'] = float('nan') if (cpp_nonzero == 0 and py_nonzero == 0) else 0.0

    return metrics


def main():
    """Run comparison between C++ and Python outputs."""

    example_dir = Path(__file__).parent
    cpp_dir = example_dir / "ScatteringData"
    python_dir = example_dir / "ScatteringData_Python"

    print("=" * 80)
    print("IceNine C++ vs Python Output Comparison - Example2.ManyGrains")
    print("=" * 80)
    print(f"\nC++ output: {cpp_dir}")
    print(f"Python output: {python_dir}\n")

    if not cpp_dir.exists():
        print(f"Error: C++ output directory not found: {cpp_dir}")
        return 1

    if not python_dir.exists():
        print(f"Error: Python output directory not found: {python_dir}")
        print("  Run run_python_simulation.py first!")
        return 1

    # Get file lists
    cpp_files = sorted(cpp_dir.glob("*.d*"))
    python_files = sorted(python_dir.glob("*.d*"))

    print(f"Found {len(cpp_files)} C++ output files")
    print(f"Found {len(python_files)} Python output files\n")

    if len(cpp_files) == 0:
        print("No C++ files found!")
        return 1

    if len(python_files) == 0:
        print("No Python files found!")
        return 1

    # Compare each file
    print("Comparing detector images...")
    print("-" * 80)

    results = []
    missing_python = []
    files_with_diffs = []
    empty_both = 0
    progress_interval = max(1, len(cpp_files) // 20)

    for i, cpp_file in enumerate(cpp_files):
        if (i + 1) % progress_interval == 0:
            print(f"  Progress: {i+1}/{len(cpp_files)} files compared...")

        python_file = python_dir / cpp_file.name

        if not python_file.exists():
            missing_python.append(cpp_file.name)
            continue

        cpp_data = load_detector_image(cpp_file)
        py_data = load_detector_image(python_file)

        if cpp_data is None or py_data is None:
            continue

        if cpp_data.shape != py_data.shape:
            print(f"  {cpp_file.name}: Shape mismatch C++{cpp_data.shape} vs Py{py_data.shape}")
            continue

        metrics = compare_images(cpp_data, py_data, cpp_file.name)
        results.append(metrics)

        # Track empty images
        if metrics['cpp_nonzero'] == 0 and metrics['py_nonzero'] == 0:
            empty_both += 1

        # Track images with differences
        if metrics['cpp_only'] > 0 or metrics['py_only'] > 0 or metrics['max_abs_diff'] > 1e-6:
            files_with_diffs.append(metrics)

    # Check for extra Python files
    cpp_names = {f.name for f in cpp_files}
    missing_cpp = [f.name for f in python_files if f.name not in cpp_names]

    print(f"\nCompared {len(results)} file pairs\n")

    if missing_python:
        print(f"  {len(missing_python)} files in C++ but not in Python:")
        for name in missing_python[:5]:
            print(f"    {name}")
        if len(missing_python) > 5:
            print(f"    ... and {len(missing_python) - 5} more")
        print()

    if missing_cpp:
        print(f"  {len(missing_cpp)} files in Python but not in C++:")
        for name in missing_cpp[:5]:
            print(f"    {name}")
        if len(missing_cpp) > 5:
            print(f"    ... and {len(missing_cpp) - 5} more")
        print()

    if len(results) == 0:
        print("No files were successfully compared!")
        return 1

    # Filter to non-empty results for statistics
    nonempty_results = [r for r in results if r['cpp_nonzero'] > 0 or r['py_nonzero'] > 0]

    print("=" * 80)
    print("Summary Statistics")
    print("=" * 80)

    print(f"\nFile counts:")
    print(f"  Total compared:       {len(results)}")
    print(f"  Both empty:           {empty_both}")
    print(f"  With content:         {len(nonempty_results)}")
    print(f"  With differences:     {len(files_with_diffs)}")

    if nonempty_results:
        # Pixel overlap
        total_cpp_px = sum(r['cpp_nonzero'] for r in nonempty_results)
        total_py_px = sum(r['py_nonzero'] for r in nonempty_results)
        total_both = sum(r['both_nonzero'] for r in nonempty_results)
        total_cpp_only = sum(r['cpp_only'] for r in nonempty_results)
        total_py_only = sum(r['py_only'] for r in nonempty_results)

        print(f"\nPixel Overlap (across all images):")
        print(f"  C++ nonzero pixels:   {total_cpp_px}")
        print(f"  Python nonzero pixels: {total_py_px}")
        print(f"  Both nonzero:         {total_both}")
        print(f"  C++ only:             {total_cpp_only}")
        print(f"  Python only:          {total_py_only}")
        if total_cpp_px > 0:
            print(f"  Match rate (vs C++):  {total_both / total_cpp_px * 100:.4f}%")
        if total_py_px > 0:
            print(f"  Pixel ratio (Py/C++): {total_py_px / total_cpp_px:.4f}")

        max_abs_diffs = [r['max_abs_diff'] for r in nonempty_results]
        mean_abs_diffs = [r['mean_abs_diff'] for r in nonempty_results]
        max_rel_diffs = [r['max_rel_diff'] for r in nonempty_results if r['both_nonzero'] > 0]
        correlations = [r['correlation'] for r in nonempty_results
                        if not (isinstance(r['correlation'], float) and np.isnan(r['correlation']))]

        print(f"\nAbsolute Differences:")
        print(f"  Max across all files:  {max(max_abs_diffs):.6e}")
        print(f"  Mean max difference:   {np.mean(max_abs_diffs):.6e}")
        print(f"  Median max difference: {np.median(max_abs_diffs):.6e}")
        print(f"  Mean avg difference:   {np.mean(mean_abs_diffs):.6e}")

        if max_rel_diffs:
            print(f"\nRelative Differences (at overlapping pixels):")
            print(f"  Max relative error:    {max(max_rel_diffs):.6e}")
            print(f"  Mean max rel error:    {np.mean(max_rel_diffs):.6e}")
            print(f"  Median max rel error:  {np.median(max_rel_diffs):.6e}")

        if correlations:
            valid_corr = [c for c in correlations if not np.isnan(c)]
            if valid_corr:
                print(f"\nCorrelation:")
                print(f"  Mean correlation:      {np.mean(valid_corr):.6f}")
                print(f"  Min correlation:       {min(valid_corr):.6f}")

        # Intensity comparison
        total_cpp_sum = sum(r['cpp_sum'] for r in nonempty_results)
        total_py_sum = sum(r['py_sum'] for r in nonempty_results)
        print(f"\nTotal Intensity:")
        print(f"  C++ total:             {total_cpp_sum:.6f}")
        print(f"  Python total:          {total_py_sum:.6f}")
        if total_cpp_sum > 0:
            print(f"  Ratio (Py/C++):        {total_py_sum / total_cpp_sum:.6f}")

    # Show worst offenders
    if files_with_diffs:
        print(f"\nWorst 10 files by pixel mismatch (C++-only + Py-only):")
        worst_px = sorted(files_with_diffs,
                          key=lambda r: r['cpp_only'] + r['py_only'], reverse=True)[:10]
        for r in worst_px:
            print(f"  {r['filename']:30s}  cpp_only={r['cpp_only']:4d}  "
                  f"py_only={r['py_only']:4d}  max_rel={r['max_rel_diff']:.2e}")

        print(f"\nWorst 10 files by max relative diff:")
        worst_rel = sorted(files_with_diffs,
                           key=lambda r: r['max_rel_diff'], reverse=True)[:10]
        for r in worst_rel:
            print(f"  {r['filename']:30s}  max_rel={r['max_rel_diff']:.6e}  "
                  f"cpp_px={r['cpp_nonzero']:5d}  py_px={r['py_nonzero']:5d}")

    # Pass/fail
    print("\n" + "=" * 80)
    print("Pass/Fail Assessment")
    print("=" * 80)

    passed = True

    if total_cpp_px > 0:
        match_rate = total_both / total_cpp_px
        if match_rate > 0.99:
            print(f"  PASS: Pixel match rate {match_rate*100:.2f}% > 99%")
        else:
            print(f"  FAIL: Pixel match rate {match_rate*100:.2f}% <= 99%")
            passed = False

    if max_rel_diffs:
        max_rel = max(max_rel_diffs)
        if max_rel < 0.01:
            print(f"  PASS: Max relative error {max_rel:.2e} < 1%")
        else:
            print(f"  FAIL: Max relative error {max_rel:.2e} >= 1%")
            passed = False

    if total_cpp_px > 0:
        px_ratio = total_py_px / total_cpp_px
        if 0.99 < px_ratio < 1.01:
            print(f"  PASS: Pixel count ratio {px_ratio:.4f} within 1%")
        else:
            print(f"  FAIL: Pixel count ratio {px_ratio:.4f} outside 1%")
            passed = False

    print()
    if passed:
        print("=" * 80)
        print("ALL TESTS PASSED - Python and C++ outputs match!")
        print("=" * 80)
        return 0
    else:
        print("=" * 80)
        print("TESTS FAILED - Outputs differ beyond threshold")
        print("=" * 80)
        return 1


if __name__ == "__main__":
    sys.exit(main())
