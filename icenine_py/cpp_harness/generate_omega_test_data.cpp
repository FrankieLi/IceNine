///////////////////////////////////////////////////////////////////////////////
//
//  generate_omega_test_data.cpp
//
//  Purpose: Generate test data for validating Python SimulationRange
//           implementation against C++ CSimulationRange.
//
//  Outputs: JSON file with test cases containing:
//           - Range configurations
//           - Expected results for angle_to_index(), to_file_number(), etc.
//
//  Usage:   ./generate_omega_test_data <output_json_path>
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "Src/SimulationData.h"
#include "Src/InitFilesIO.h"
#include "XDM++/libXDM/3dMath.h"

using namespace std;
using namespace XDMSimulation;

///////////////////////////////////////////////////////////////////////////////
//  Helper: Convert radians to degrees
///////////////////////////////////////////////////////////////////////////////
double rad_to_deg(double rad) {
    return rad * 180.0 / M_PI;
}

///////////////////////////////////////////////////////////////////////////////
//  Helper: Convert degrees to radians
///////////////////////////////////////////////////////////////////////////////
double deg_to_rad(double deg) {
    return deg * M_PI / 180.0;
}

///////////////////////////////////////////////////////////////////////////////
//  Helper: Write JSON escaped string
///////////////////////////////////////////////////////////////////////////////
void write_json_string(ostream& os, const string& s) {
    os << "\"" << s << "\"";
}

///////////////////////////////////////////////////////////////////////////////
//  Test Case 1: Single Wedge Covering Entire Range
///////////////////////////////////////////////////////////////////////////////
void generate_single_wedge_test(ostream& os) {
    cout << "Generating Test Case 1: Single Wedge..." << endl;

    // Configuration: -90° to +90°, 1° width, single wedge covering all
    Float fLow = deg_to_rad(-90.0);
    Float fHigh = deg_to_rad(90.0);
    Float fWidth = deg_to_rad(1.0);

    vector<SRange> vRangeList;
    vRangeList.push_back(SRange(fLow, fHigh));

    // Create CSimulationRange
    CSimulationRange oRange;
    oRange.Set(fLow, fHigh, fWidth, vRangeList);

    os << "    {\n";
    os << "      \"name\": \"single_wedge_full_range\",\n";
    os << "      \"description\": \"Single wedge covering entire range -90 to +90 degrees\",\n";
    os << "      \"config\": {\n";
    os << "        \"low\": " << fLow << ",\n";
    os << "        \"high\": " << fHigh << ",\n";
    os << "        \"width\": " << fWidth << ",\n";
    os << "        \"range_list\": [\n";
    os << "          {\"low\": " << fLow << ", \"high\": " << fHigh << "}\n";
    os << "        ]\n";
    os << "      },\n";
    os << "      \"tests\": [\n";

    // Test various angles
    double test_angles_deg[] = {-90.0, -45.0, 0.0, 45.0, 89.0};
    int num_tests = sizeof(test_angles_deg) / sizeof(test_angles_deg[0]);

    for (int i = 0; i < num_tests; i++) {
        Float fAngle = deg_to_rad(test_angles_deg[i]);
        Size_Type nFileNum = oRange.ToFileNumber(fAngle);
        Size_Type nWedgeIdx = oRange(fAngle);

        os << "        {\n";
        os << "          \"angle_deg\": " << test_angles_deg[i] << ",\n";
        os << "          \"angle_rad\": " << fAngle << ",\n";
        os << "          \"expected_file\": " << nFileNum << ",\n";

        if (nWedgeIdx == XDMSimulation::NoMatch) {
            os << "          \"expected_wedge\": null\n";
        } else {
            os << "          \"expected_wedge\": " << nWedgeIdx << "\n";
        }

        os << "        }";
        if (i < num_tests - 1) os << ",";
        os << "\n";
    }

    os << "      ]\n";
    os << "    }";
}

///////////////////////////////////////////////////////////////////////////////
//  Test Case 2: Multiple Wedges with Gaps
///////////////////////////////////////////////////////////////////////////////
void generate_multiple_wedges_test(ostream& os) {
    cout << "Generating Test Case 2: Multiple Wedges with Gaps..." << endl;

    // Configuration: -90° to +90°, 1° width
    // Wedges: [-90, -85], [-45, -40], [40, 45], [85, 90]
    Float fLow = deg_to_rad(-90.0);
    Float fHigh = deg_to_rad(90.0);
    Float fWidth = deg_to_rad(1.0);

    vector<SRange> vRangeList;
    vRangeList.push_back(SRange(deg_to_rad(-90.0), deg_to_rad(-85.0)));
    vRangeList.push_back(SRange(deg_to_rad(-45.0), deg_to_rad(-40.0)));
    vRangeList.push_back(SRange(deg_to_rad(40.0), deg_to_rad(45.0)));
    vRangeList.push_back(SRange(deg_to_rad(85.0), deg_to_rad(90.0)));

    // Create CSimulationRange
    CSimulationRange oRange;
    oRange.Set(fLow, fHigh, fWidth, vRangeList);

    os << ",\n";
    os << "    {\n";
    os << "      \"name\": \"multiple_wedges_with_gaps\",\n";
    os << "      \"description\": \"Four wedges with gaps between them\",\n";
    os << "      \"config\": {\n";
    os << "        \"low\": " << fLow << ",\n";
    os << "        \"high\": " << fHigh << ",\n";
    os << "        \"width\": " << fWidth << ",\n";
    os << "        \"range_list\": [\n";
    os << "          {\"low\": " << deg_to_rad(-90.0) << ", \"high\": " << deg_to_rad(-85.0) << "},\n";
    os << "          {\"low\": " << deg_to_rad(-45.0) << ", \"high\": " << deg_to_rad(-40.0) << "},\n";
    os << "          {\"low\": " << deg_to_rad(40.0) << ", \"high\": " << deg_to_rad(45.0) << "},\n";
    os << "          {\"low\": " << deg_to_rad(85.0) << ", \"high\": " << deg_to_rad(90.0) << "}\n";
    os << "        ]\n";
    os << "      },\n";
    os << "      \"tests\": [\n";

    // Test angles: in wedges and in gaps
    double test_angles_deg[] = {
        -87.5,  // Wedge 0
        -42.5,  // Wedge 1
        0.0,    // Gap!
        42.5,   // Wedge 2
        87.5,   // Wedge 3
        -60.0,  // Gap between wedge 0 and 1
        20.0    // Gap between wedge 1 and 2
    };
    int num_tests = sizeof(test_angles_deg) / sizeof(test_angles_deg[0]);

    for (int i = 0; i < num_tests; i++) {
        Float fAngle = deg_to_rad(test_angles_deg[i]);
        Size_Type nFileNum = oRange.ToFileNumber(fAngle);
        Size_Type nWedgeIdx = oRange(fAngle);

        os << "        {\n";
        os << "          \"angle_deg\": " << test_angles_deg[i] << ",\n";
        os << "          \"angle_rad\": " << fAngle << ",\n";

        if (nFileNum == XDMSimulation::NoMatch) {
            os << "          \"expected_file\": null,\n";
        } else {
            os << "          \"expected_file\": " << nFileNum << ",\n";
        }

        if (nWedgeIdx == XDMSimulation::NoMatch) {
            os << "          \"expected_wedge\": null\n";
        } else {
            os << "          \"expected_wedge\": " << nWedgeIdx << "\n";
        }

        os << "        }";
        if (i < num_tests - 1) os << ",";
        os << "\n";
    }

    os << "      ]\n";
    os << "    }";
}

///////////////////////////////////////////////////////////////////////////////
//  Test Case 3: Fine Angular Resolution
///////////////////////////////////////////////////////////////////////////////
void generate_fine_resolution_test(ostream& os) {
    cout << "Generating Test Case 3: Fine Angular Resolution..." << endl;

    // Configuration: -10° to +10°, 0.1° width (200 bins)
    // Single wedge: [-5, +5]
    Float fLow = deg_to_rad(-10.0);
    Float fHigh = deg_to_rad(10.0);
    Float fWidth = deg_to_rad(0.1);

    vector<SRange> vRangeList;
    vRangeList.push_back(SRange(deg_to_rad(-5.0), deg_to_rad(5.0)));

    // Create CSimulationRange
    CSimulationRange oRange;
    oRange.Set(fLow, fHigh, fWidth, vRangeList);

    os << ",\n";
    os << "    {\n";
    os << "      \"name\": \"fine_angular_resolution\",\n";
    os << "      \"description\": \"Fine 0.1 degree bins, single wedge -5 to +5\",\n";
    os << "      \"config\": {\n";
    os << "        \"low\": " << fLow << ",\n";
    os << "        \"high\": " << fHigh << ",\n";
    os << "        \"width\": " << fWidth << ",\n";
    os << "        \"range_list\": [\n";
    os << "          {\"low\": " << deg_to_rad(-5.0) << ", \"high\": " << deg_to_rad(5.0) << "}\n";
    os << "        ]\n";
    os << "      },\n";
    os << "      \"tests\": [\n";

    // Test angles at boundaries and outside
    double test_angles_deg[] = {-9.9, -5.0, -4.9, 0.0, 4.9, 5.0, 9.9};
    int num_tests = sizeof(test_angles_deg) / sizeof(test_angles_deg[0]);

    for (int i = 0; i < num_tests; i++) {
        Float fAngle = deg_to_rad(test_angles_deg[i]);
        Size_Type nWedgeIdx = oRange(fAngle);

        os << "        {\n";
        os << "          \"angle_deg\": " << test_angles_deg[i] << ",\n";
        os << "          \"angle_rad\": " << fAngle << ",\n";

        if (nWedgeIdx == XDMSimulation::NoMatch) {
            os << "          \"expected_wedge\": null\n";
        } else {
            os << "          \"expected_wedge\": " << nWedgeIdx << "\n";
        }

        os << "        }";
        if (i < num_tests - 1) os << ",";
        os << "\n";
    }

    os << "      ]\n";
    os << "    }";
}

///////////////////////////////////////////////////////////////////////////////
//  Test Case 4: Experimental Data - 5000 Random Angles
///////////////////////////////////////////////////////////////////////////////
void generate_experimental_omega_test(ostream& os) {
    cout << "Generating Test Case 4: Experimental Omega File (5000 random angles)..." << endl;

    // Read actual experimental omega file
    string omega_file = "../../DataFiles/omega_180_2L.dat";
    vector<SRange> vOmegaRangeList;
    vector<SIntRange> vFileRangeList;
    Size_Type nNumDetectors = 2;

    bool success = InitFileIO::ReadRotationIntervalFiles(
        vOmegaRangeList,
        vFileRangeList,
        nNumDetectors,
        omega_file
    );

    if (!success) {
        cerr << "ERROR: Failed to read omega file: " << omega_file << endl;
        cerr << "Skipping experimental test case." << endl;
        return;
    }

    cout << "  Read " << vOmegaRangeList.size() << " omega ranges from " << omega_file << endl;
    cout << "  File ranges: " << vFileRangeList.size() << " detectors" << endl;

    // Replicate CXDMExperimentSetup::ReadRotationInterval logic
    // Lines 189-193: Swap individual ranges if reversed
    for (size_t i = 0; i < vOmegaRangeList.size(); i++) {
        if (vOmegaRangeList[i].fHigh < vOmegaRangeList[i].fLow) {
            std::swap(vOmegaRangeList[i].fHigh, vOmegaRangeList[i].fLow);
        }
    }

    // Lines 210-219: Determine overall range (may be flipped)
    Size_Type nLastIndex = vOmegaRangeList.size() - 1;
    Float fHigh, fLow;

    if (vOmegaRangeList[0].fHigh > vOmegaRangeList[nLastIndex].fHigh) {
        // Flipped sequence (descending)
        fHigh = vOmegaRangeList[nLastIndex].fLow;   // flipped on purpose
        fLow = vOmegaRangeList[0].fHigh;
    } else {
        // Normal sequence (ascending)
        fHigh = vOmegaRangeList[nLastIndex].fHigh;
        fLow = vOmegaRangeList[0].fLow;
    }

    // Line 221: Width from first range
    Float fWidth = vOmegaRangeList[0].fHigh - vOmegaRangeList[0].fLow;

    // Create CSimulationRange
    CSimulationRange oRange;
    oRange.Set(fLow, fHigh, fWidth, vOmegaRangeList);

    cout << "  Overall range: " << rad_to_deg(fLow) << "° to " << rad_to_deg(fHigh) << "°" << endl;
    cout << "  Angular width: " << rad_to_deg(fWidth) << "°" << endl;

    os << ",\n";
    os << "    {\n";
    os << "      \"name\": \"experimental_omega_5000_random\",\n";
    os << "      \"description\": \"Real experimental omega file with 5000 random test angles at 1° resolution\",\n";
    os << "      \"omega_file\": \"omega_180_2L.dat\",\n";
    os << "      \"config\": {\n";
    os << "        \"low\": " << std::setprecision(17) << fLow << ",\n";
    os << "        \"high\": " << std::setprecision(17) << fHigh << ",\n";
    os << "        \"width\": " << std::setprecision(17) << fWidth << ",\n";
    os << "        \"range_list\": [\n";

    // Write all omega ranges
    for (size_t i = 0; i < vOmegaRangeList.size(); i++) {
        os << "          {\"low\": " << std::setprecision(17) << vOmegaRangeList[i].fLow
           << ", \"high\": " << std::setprecision(17) << vOmegaRangeList[i].fHigh << "}";
        if (i < vOmegaRangeList.size() - 1) os << ",";
        os << "\n";
    }

    os << "        ]\n";
    os << "      },\n";
    os << "      \"tests\": [\n";

    // Generate 5000 random test angles uniformly distributed across the range
    srand(42);  // Fixed seed for reproducibility
    const int num_tests = 5000;

    cout << "  Generating " << num_tests << " random test angles..." << endl;

    for (int i = 0; i < num_tests; i++) {
        // Random angle in the overall range
        double random_fraction = (double)rand() / RAND_MAX;
        double angle_deg = rad_to_deg(fLow) + random_fraction * (rad_to_deg(fHigh) - rad_to_deg(fLow));
        Float fAngle = deg_to_rad(angle_deg);

        Size_Type nFileNum = oRange.ToFileNumber(fAngle);
        Size_Type nWedgeIdx = oRange(fAngle);

        if (i > 0) {
            os << ",\n";
        }

        os << "        {\n";
        os << "          \"angle_deg\": " << std::setprecision(17) << angle_deg << ",\n";
        os << "          \"angle_rad\": " << std::setprecision(17) << fAngle << ",\n";

        if (nFileNum == XDMSimulation::NoMatch) {
            os << "          \"expected_file\": null,\n";
        } else {
            os << "          \"expected_file\": " << nFileNum << ",\n";
        }

        if (nWedgeIdx == XDMSimulation::NoMatch) {
            os << "          \"expected_wedge\": null\n";
        } else {
            os << "          \"expected_wedge\": " << nWedgeIdx << "\n";
        }

        os << "        }";

        // Progress indicator
        if ((i + 1) % 1000 == 0) {
            cout << "    Generated " << (i + 1) << " / " << num_tests << " test points..." << endl;
        }
    }

    os << "\n      ]\n";
    os << "    }";

    cout << "  Completed experimental omega test case." << endl;
}

///////////////////////////////////////////////////////////////////////////////
//  Main
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <output_json_path>" << endl;
        return 1;
    }

    string output_file = argv[1];

    cout << "IceNine Omega Range Test Data Generator" << endl;
    cout << "========================================" << endl;
    cout << "Output: " << output_file << endl;
    cout << endl;

    // Open output file
    ofstream os(output_file.c_str());
    if (!os.is_open()) {
        cerr << "ERROR: Failed to open output file: " << output_file << endl;
        return 1;
    }

    // Write JSON header
    os << "{\n";
    os << "  \"description\": \"C++ CSimulationRange validation test data\",\n";
    os << "  \"generator\": \"generate_omega_test_data.cpp\",\n";
    os << "  \"test_cases\": [\n";

    // Generate test cases
    generate_single_wedge_test(os);
    generate_multiple_wedges_test(os);
    generate_fine_resolution_test(os);
    generate_experimental_omega_test(os);

    // Write JSON footer
    os << "\n  ]\n";
    os << "}\n";

    os.close();

    cout << endl;
    cout << "SUCCESS: Generated test data in " << output_file << endl;

    return 0;
}
