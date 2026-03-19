#!/usr/bin/env python3
"""
Detector geometry diagnostic: print all detector geometric properties
for comparison against C++ values.

Usage:
    cd Examples/Example2.ThreeVoxels
    uv run python debug_detector_geometry.py
"""

import sys
from pathlib import Path
import numpy as np
import torch

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "icenine_py"))

from icenine.file_io import read_detector_file


def main():
    example_dir = Path(__file__).parent
    detector_file = example_dir / "ConfigFiles" / "StandardGeometry.2Det"

    print("=" * 70)
    print("Detector Geometry Diagnostic")
    print("=" * 70)
    print(f"File: {detector_file}\n")

    detectors = read_detector_file(str(detector_file))
    print(f"Number of detectors: {len(detectors)}\n")

    for i, det in enumerate(detectors):
        print(f"{'=' * 70}")
        print(f"Detector {i}")
        print(f"{'=' * 70}")

        # Basic parameters
        print(f"  Dimensions: {det.num_cols} x {det.num_rows} (cols x rows)")
        print(f"  Pixel size: {det.pixel_width:.10f} x {det.pixel_height:.10f} mm (w x h)")
        print(f"  Beam center: J={det.beam_center_j:.4f}, K={det.beam_center_k:.4f} pixels")
        print(f"  Beam center: J={det.beam_center_j * det.pixel_width:.6f}, "
              f"K={det.beam_center_k * det.pixel_height:.6f} mm")

        # Position
        print(f"\n  Position (lab frame): {det._position.numpy()}")

        # Orientation matrix
        print(f"  Orientation matrix:")
        for row in range(3):
            vals = det._orientation[row].numpy()
            print(f"    [{vals[0]:10.6f} {vals[1]:10.6f} {vals[2]:10.6f}]")

        # Detector-frame basis vectors
        print(f"\n  Det-frame J basis: {det._det_frame_basis_j.numpy()}")
        print(f"  Det-frame K basis: {det._det_frame_basis_k.numpy()}")
        print(f"  Det-frame coord origin: {det._det_frame_coord_origin.numpy()}")

        # Lab-frame basis vectors (after orientation rotation)
        print(f"\n  Lab-frame J basis: {det._lab_frame_basis_j.numpy()}")
        print(f"  Lab-frame K basis: {det._lab_frame_basis_k.numpy()}")
        print(f"  Lab-frame coord origin: {det._lab_frame_coord_origin.numpy()}")

        # Detector plane
        plane = det._detector_plane
        print(f"\n  Plane normal: {plane.normal.numpy()}")
        print(f"  Plane D: {plane.d.item():.10f}")
        print(f"  Plane coeffs [A,B,C,D]: {plane.coeffs.numpy()}")

        # Verify: beam center should map to detector position
        # i.e., if we shoot a ray from origin along the detector position vector,
        # it should hit near beam center pixels
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
