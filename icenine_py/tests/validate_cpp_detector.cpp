/**
 * C++ test harness to validate Python Detector implementation
 *
 * Compile with:
 *   g++ -std=c++11 -I../XDM++/libXDM -I../Src validate_cpp_detector.cpp \
 *       ../XDM++/libXDM/3dMath.cpp ../Src/Detector.cpp -o validate_cpp_detector
 *
 * Run:
 *   ./validate_cpp_detector > cpp_detector_output.txt
 */

#include <iostream>
#include <iomanip>
#include "../Src/Detector.h"
#include "../XDM++/libXDM/3dMath.h"

using namespace std;
using namespace GeneralLib;

void print_vector(const char* name, const SVector3& v) {
    cout << name << ": [" << setprecision(10)
         << v.m_fX << ", " << v.m_fY << ", " << v.m_fZ << "]" << endl;
}

void print_matrix(const char* name, const SMatrix3x3& m) {
    cout << name << ":" << endl;
    for (int i = 0; i < 3; i++) {
        cout << "  [" << setprecision(10);
        for (int j = 0; j < 3; j++) {
            cout << m.m[i][j];
            if (j < 2) cout << ", ";
        }
        cout << "]" << endl;
    }
}

int main() {
    cout << "# C++ Detector Validation Output" << endl;
    cout << "# Format: key: value" << endl << endl;

    // Test 1: Basic detector at origin
    {
        cout << "## Test 1: Basic detector at origin" << endl;

        BBox2D range(Point(0, 0), Point(2048, 2048));
        SVector3 j_unit(0, 1, 0);
        SVector3 k_unit(0, 0, 1);
        SVector3 position(100.0, 0.0, 0.0);
        SMatrix3x3 orientation;
        orientation.SetIdentity();

        CDetector detector = CXDMDetectorFactory::MakeDetector(
            2048, 2048,           // num_j_pixels, num_k_pixels
            position,             // position
            1024.0, 1024.0,       // beam_center_j, beam_center_k (in pixels)
            j_unit, k_unit,       // basis vectors
            0.2, 0.2,             // pixel_width, pixel_height
            orientation           // orientation matrix
        );

        cout << "num_rows: " << detector.GetNumRows() << endl;
        cout << "num_cols: " << detector.GetNumCols() << endl;
        cout << "pixel_width: " << detector.GetPixelWidth() << endl;
        cout << "pixel_height: " << detector.GetPixelHeight() << endl;
        cout << "beam_center_j: " << detector.GetBeamCenterJ() << endl;
        cout << "beam_center_k: " << detector.GetBeamCenterK() << endl;

        print_vector("position", detector.GetLocation());
        print_matrix("orientation", detector.GetOrientationMatrix());

        // Test coordinate transformations
        cout << endl << "### Coordinate Transformations" << endl;

        // Lab to detector coordinate at detector center
        SVector3 test_point1(100.0, 0.0, 0.0);
        Float j1, k1;
        detector.LabToDetectorCoordinate(j1, k1, test_point1);
        cout << "lab_to_detector([100, 0, 0]): j=" << j1 << ", k=" << k1 << endl;

        // Lab to pixel at detector center
        Float row1, col1;
        detector.LabToPixel(row1, col1, test_point1);
        cout << "lab_to_pixel([100, 0, 0]): row=" << row1 << ", col=" << col1 << endl;

        // Lab to detector coordinate at offset point
        SVector3 test_point2(100.0, 10.0, 5.0);
        Float j2, k2;
        detector.LabToDetectorCoordinate(j2, k2, test_point2);
        cout << "lab_to_detector([100, 10, 5]): j=" << j2 << ", k=" << k2 << endl;

        // Detector to lab coordinate
        SVector3 lab_point = detector.DetectorToLabCoordinate(10.0, 5.0);
        print_vector("detector_to_lab(j=10, k=5)", lab_point);

        // Pixel to lab coordinate
        SVector3 lab_from_pixel = detector.PixelToLabCoordinate(1024.0, 1024.0);
        print_vector("pixel_to_lab(col=1024, row=1024)", lab_from_pixel);

        // Get basis vectors
        SVector3 j_basis, k_basis;
        detector.GetCoordinateBasis(j_basis, k_basis);
        print_vector("j_basis", j_basis);
        print_vector("k_basis", k_basis);

        cout << endl;
    }

    // Test 2: Detector with Euler angle rotation
    {
        cout << "## Test 2: Detector with Euler rotation (10°, 5°, 0°)" << endl;

        BBox2D range(Point(0, 0), Point(1024, 1024));
        SVector3 j_unit(0, 1, 0);
        SVector3 k_unit(0, 0, 1);
        SVector3 position(150.0, 10.0, -5.0);
        SMatrix3x3 orientation;

        // Convert degrees to radians (C++ expects radians)
        Float phi = 10.0 * PI / 180.0;
        Float theta = 5.0 * PI / 180.0;
        Float psi = 0.0;
        orientation.BuildActiveEulerMatrix(phi, theta, psi);

        CDetector detector = CXDMDetectorFactory::MakeDetector(
            1024, 1024,
            position,
            512.0, 512.0,
            j_unit, k_unit,
            0.2, 0.2,
            orientation
        );

        print_vector("position", detector.GetLocation());
        print_matrix("orientation", detector.GetOrientationMatrix());

        // Test transformations with rotated detector
        SVector3 test_point(150.0, 0.0, 0.0);
        Float j, k;
        detector.LabToDetectorCoordinate(j, k, test_point);
        cout << "lab_to_detector([150, 0, 0]): j=" << j << ", k=" << k << endl;

        // Get rotated basis vectors
        SVector3 j_basis, k_basis;
        detector.GetCoordinateBasis(j_basis, k_basis);
        print_vector("j_basis_rotated", j_basis);
        print_vector("k_basis_rotated", k_basis);

        cout << endl;
    }

    // Test 3: Ray intersection
    {
        cout << "## Test 3: Ray-detector intersection" << endl;

        BBox2D range(Point(0, 0), Point(2048, 2048));
        SVector3 j_unit(0, 1, 0);
        SVector3 k_unit(0, 0, 1);
        SVector3 position(100.0, 0.0, 0.0);
        SMatrix3x3 orientation;
        orientation.SetIdentity();

        CDetector detector = CXDMDetectorFactory::MakeDetector(
            2048, 2048,
            position,
            1024.0, 1024.0,
            j_unit, k_unit,
            0.2, 0.2,
            orientation
        );

        // Ray along X-axis
        CRay ray(SVector3(0, 0, 0), SVector3(1, 0, 0));
        Float t;
        bool intersects = detector.Intersects(ray, t);

        cout << "ray_intersects: " << (intersects ? "true" : "false") << endl;
        cout << "ray_t: " << t << endl;

        if (intersects) {
            SVector3 hit_point = ray.GetStart() + t * ray.GetDirection();
            print_vector("hit_point", hit_point);

            Float row, col;
            detector.LabToPixel(row, col, hit_point);
            cout << "hit_pixel: row=" << row << ", col=" << col << endl;
        }

        cout << endl;
    }

    return 0;
}
