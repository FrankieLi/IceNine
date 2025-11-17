"""
Python test harness to validate against C++ Detector implementation

Run:
    python validate_py_detector.py > py_detector_output.txt

Compare with C++:
    diff cpp_detector_output.txt py_detector_output.txt
"""

import torch
import sys
sys.path.insert(0, '/Users/sfli/Research/IceNine/icenine_py')

from icenine.detector import Detector
from icenine.geometry import Ray


def print_vector(name, v):
    """Print vector in same format as C++"""
    print(f"{name}: [{v[0]:.10f}, {v[1]:.10f}, {v[2]:.10f}]")


def print_matrix(name, m):
    """Print matrix in same format as C++"""
    print(f"{name}:")
    for i in range(3):
        row = m[i]
        print(f"  [{row[0]:.10f}, {row[1]:.10f}, {row[2]:.10f}]")


def main():
    print("# Python Detector Validation Output")
    print("# Format: key: value")
    print()

    # Test 1: Basic detector at origin
    print("## Test 1: Basic detector at origin")

    detector = Detector(
        num_rows=2048,
        num_cols=2048,
        pixel_height=0.2,
        pixel_width=0.2,
        beam_center_j=1024.0,  # pixels
        beam_center_k=1024.0,  # pixels
        position=torch.tensor([100.0, 0.0, 0.0]),
        orientation=torch.eye(3),
        dtype=torch.float64  # Use float64 for precision matching
    )

    print(f"num_rows: {detector.num_rows}")
    print(f"num_cols: {detector.num_cols}")
    print(f"pixel_width: {detector.pixel_width}")
    print(f"pixel_height: {detector.pixel_height}")
    print(f"beam_center_j: {detector.beam_center_j}")
    print(f"beam_center_k: {detector.beam_center_k}")

    print_vector("position", detector.position)
    print_matrix("orientation", detector.orientation)

    # Test coordinate transformations
    print()
    print("### Coordinate Transformations")

    # Lab to detector coordinate at detector center
    test_point1 = torch.tensor([100.0, 0.0, 0.0], dtype=torch.float64)
    j1, k1 = detector.lab_to_detector_coordinate(test_point1)
    print(f"lab_to_detector([100, 0, 0]): j={j1.item():.10f}, k={k1.item():.10f}")

    # Lab to pixel at detector center
    row1, col1 = detector.lab_to_pixel(test_point1)
    print(f"lab_to_pixel([100, 0, 0]): row={row1.item():.10f}, col={col1.item():.10f}")

    # Lab to detector coordinate at offset point
    test_point2 = torch.tensor([100.0, 10.0, 5.0], dtype=torch.float64)
    j2, k2 = detector.lab_to_detector_coordinate(test_point2)
    print(f"lab_to_detector([100, 10, 5]): j={j2.item():.10f}, k={k2.item():.10f}")

    # Detector to lab coordinate
    lab_point = detector.detector_to_lab_coordinate(
        torch.tensor(10.0, dtype=torch.float64),
        torch.tensor(5.0, dtype=torch.float64)
    )
    print_vector("detector_to_lab(j=10, k=5)", lab_point)

    # Pixel to lab coordinate
    lab_from_pixel = detector.pixel_to_lab_coordinate(
        torch.tensor(1024.0, dtype=torch.float64),
        torch.tensor(1024.0, dtype=torch.float64)
    )
    print_vector("pixel_to_lab(col=1024, row=1024)", lab_from_pixel)

    # Get basis vectors
    j_basis, k_basis = detector.basis_vectors
    print_vector("j_basis", j_basis)
    print_vector("k_basis", k_basis)

    print()

    # Test 2: Detector with Euler angle rotation
    print("## Test 2: Detector with Euler rotation (10°, 5°, 0°)")

    detector2 = Detector(
        num_rows=1024,
        num_cols=1024,
        pixel_height=0.2,
        pixel_width=0.2,
        beam_center_j=512.0,
        beam_center_k=512.0,
        position=torch.tensor([150.0, 10.0, -5.0], dtype=torch.float64),
        dtype=torch.float64
    )

    # Set orientation using Euler angles
    detector2.set_orientation_euler(phi=10.0, theta=5.0, psi=0.0)

    print_vector("position", detector2.position)
    print_matrix("orientation", detector2.orientation)

    # Test transformations with rotated detector
    test_point = torch.tensor([150.0, 0.0, 0.0], dtype=torch.float64)
    j, k = detector2.lab_to_detector_coordinate(test_point)
    print(f"lab_to_detector([150, 0, 0]): j={j.item():.10f}, k={k.item():.10f}")

    # Get rotated basis vectors
    j_basis_rot, k_basis_rot = detector2.basis_vectors
    print_vector("j_basis_rotated", j_basis_rot)
    print_vector("k_basis_rotated", k_basis_rot)

    print()

    # Test 3: Ray intersection
    print("## Test 3: Ray-detector intersection")

    detector3 = Detector(
        num_rows=2048,
        num_cols=2048,
        pixel_height=0.2,
        pixel_width=0.2,
        beam_center_j=1024.0,
        beam_center_k=1024.0,
        position=torch.tensor([100.0, 0.0, 0.0], dtype=torch.float64),
        orientation=torch.eye(3, dtype=torch.float64),
        dtype=torch.float64
    )

    # Ray along X-axis
    ray = Ray(
        origin=torch.tensor([0.0, 0.0, 0.0], dtype=torch.float64),
        direction=torch.tensor([1.0, 0.0, 0.0], dtype=torch.float64)
    )

    intersects, t = detector3.intersect_ray(ray)

    print(f"ray_intersects: {'true' if intersects.item() else 'false'}")
    print(f"ray_t: {t.item():.10f}")

    if intersects.item():
        hit_point = ray.at(t)
        print_vector("hit_point", hit_point)

        row, col = detector3.lab_to_pixel(hit_point)
        print(f"hit_pixel: row={row.item():.10f}, col={col.item():.10f}")

    print()


if __name__ == '__main__':
    main()
