#!/usr/bin/env python3
"""
Compare C++ and Python forward simulation output by pixel overlap.

Counts how many non-zero C++ pixels are also non-zero in Python output,
and vice versa. This is the key metric — intensity values will differ
due to triangle size differences, but the same pixels should be hit.

Usage:
    cd icenine_py
    uv run python ../Examples/Example2.ThreeVoxels/compare_pixel_overlap.py
"""

import sys
from pathlib import Path
import numpy as np
import torch

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "icenine_py"))

from icenine.image_data import ImageData


def load_nonzero_pixels(filepath: Path) -> set:
    """Load detector image and return set of (row, col) tuples for non-zero pixels."""
    try:
        img = ImageData(num_rows=2048, num_cols=2048, mode='dense')
        img.load_ascii(str(filepath))
        data = img._pixels_dense
        nonzero = torch.nonzero(data, as_tuple=False)
        return {(r.item(), c.item()) for r, c in nonzero}
    except Exception as e:
        return None


def main():
    example_dir = Path(__file__).parent
    cpp_dir = example_dir / "ScatteringData"
    python_dir = example_dir / "ScatteringData_Python"

    print("=" * 70)
    print("Pixel Overlap Comparison: C++ vs Python")
    print("=" * 70)

    cpp_files = sorted(cpp_dir.glob("*.d*"))
    python_files = sorted(python_dir.glob("*.d*"))

    print(f"C++ files: {len(cpp_files)}")
    print(f"Python files: {len(python_files)}\n")

    # Aggregate stats
    total_cpp_pixels = 0
    total_py_pixels = 0
    total_matched = 0
    total_cpp_only = 0
    total_py_only = 0
    files_with_data = 0

    # Per-detector stats
    d0_cpp = 0; d0_py = 0; d0_match = 0
    d1_cpp = 0; d1_py = 0; d1_match = 0

    worst_files = []

    for cpp_file in cpp_files:
        py_file = python_dir / cpp_file.name
        if not py_file.exists():
            continue

        cpp_pixels = load_nonzero_pixels(cpp_file)
        py_pixels = load_nonzero_pixels(py_file)

        if cpp_pixels is None or py_pixels is None:
            continue

        if len(cpp_pixels) == 0 and len(py_pixels) == 0:
            continue

        files_with_data += 1
        matched = cpp_pixels & py_pixels
        cpp_only = cpp_pixels - py_pixels
        py_only = py_pixels - cpp_pixels

        total_cpp_pixels += len(cpp_pixels)
        total_py_pixels += len(py_pixels)
        total_matched += len(matched)
        total_cpp_only += len(cpp_only)
        total_py_only += len(py_only)

        # Per-detector
        if cpp_file.name.endswith('.d0'):
            d0_cpp += len(cpp_pixels)
            d0_py += len(py_pixels)
            d0_match += len(matched)
        elif cpp_file.name.endswith('.d1'):
            d1_cpp += len(cpp_pixels)
            d1_py += len(py_pixels)
            d1_match += len(matched)

        # Track worst files
        if len(cpp_pixels) > 0:
            recall = len(matched) / len(cpp_pixels)
            if recall < 1.0 and len(cpp_only) > 0:
                worst_files.append((cpp_file.name, len(cpp_pixels), len(py_pixels),
                                    len(matched), len(cpp_only), len(py_only)))

    # Summary
    print(f"Files with non-zero pixels: {files_with_data}")
    print()

    print("=" * 70)
    print("Overall Pixel Overlap")
    print("=" * 70)
    print(f"  C++ non-zero pixels:    {total_cpp_pixels}")
    print(f"  Python non-zero pixels: {total_py_pixels}")
    print(f"  Matched (both non-zero):{total_matched}")
    print(f"  C++ only (missed by Py):{total_cpp_only}")
    print(f"  Python only (extra):    {total_py_only}")
    print()

    if total_cpp_pixels > 0:
        recall = total_matched / total_cpp_pixels * 100
        print(f"  Recall (C++ pixels found in Py): {recall:.1f}% ({total_matched}/{total_cpp_pixels})")
    if total_py_pixels > 0:
        precision = total_matched / total_py_pixels * 100
        print(f"  Precision (Py pixels in C++):     {precision:.1f}% ({total_matched}/{total_py_pixels})")
    print(f"  Py/C++ pixel ratio:               {total_py_pixels/max(total_cpp_pixels,1):.2f}x")

    # Per-detector
    print()
    print("=" * 70)
    print("Per-Detector Breakdown")
    print("=" * 70)
    for name, cpp, py, match in [("Detector 0 (.d0)", d0_cpp, d0_py, d0_match),
                                   ("Detector 1 (.d1)", d1_cpp, d1_py, d1_match)]:
        if cpp > 0:
            print(f"  {name}: C++={cpp}, Py={py}, Match={match}, "
                  f"Recall={match/cpp*100:.1f}%, Ratio={py/cpp:.2f}x")

    # Worst files
    if worst_files:
        worst_files.sort(key=lambda x: x[3] / max(x[1], 1))  # Sort by recall ascending
        print()
        print("=" * 70)
        print("Worst 10 Files by Recall")
        print("=" * 70)
        for name, cpp, py, match, cpp_only, py_only in worst_files[:10]:
            recall = match / max(cpp, 1) * 100
            print(f"  {name:30s} C++={cpp:4d} Py={py:4d} Match={match:4d} "
                  f"Recall={recall:.0f}% CppOnly={cpp_only} PyOnly={py_only}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
