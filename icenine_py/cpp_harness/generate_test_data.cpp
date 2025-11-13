/////////////////////////////////////////////////////////////////
//
//  File:    generate_test_data.cpp
//  Author:  IceNine Migration Project
//
//  Purpose: Generate ground truth test data from C++ IceNine
//           for validating Python/PyTorch port
//
//  Outputs JSON files:
//    - cubic_symmetry_24ops.json: 24 cubic symmetry matrices
//    - vector_equivalence_tests.json: Test cases for vector equivalence
//    - au_fcc_all_hkl_qmax8.json: All Miller indices before symmetry reduction
//    - au_fcc_unique_hkl_qmax8.json: Unique reflections after reduction
//    - q_magnitudes.json: Q magnitudes for test reflections
//    - scattering_omegas.json: ~5000 scattering omega angles with uniform SO(3) rotations
//
/////////////////////////////////////////////////////////////////

#include "XDM++/libXDM/Symmetry.h"
#include "XDM++/libXDM/CrystalStructure.h"
#include "XDM++/libXDM/3dMath.h"
#include "XDM++/libXDM/Quaternion.h"
#include "XDM++/libXDM/PhysicalConstants.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace GeneralLib;
using namespace LatticeSymmetry;

// Test parameters from ConfigFiles/ReconstructTest.config
const Float AU_LATTICE_A = 4.0782;  // Angstroms
const Float BEAM_ENERGY = 50.02099; // keV
const Float WAVELENGTH = 0.2478;    // Angstroms
const Float MAX_Q = 8.0;            // Å^-1
const Float TWO_PI = 6.28318530718;
const Float BEAM_DEFLECTION_CHI = 0.0;  // radians (no beam deflection)

//------------------------------------------------------------------------------
// GetScatteringOmegas
//
// Standalone implementation extracted from Simulation.cpp (lines 65-119)
// for testing purposes.
//
// This function calculates the omega angles at which the Bragg condition
// is satisfied for a specific scattering vector.
//
// Returns true if observable, false if too close to rotation axis.
//------------------------------------------------------------------------------
bool GetScatteringOmegas(Float &fOmegaRes1, Float &fOmegaRes2,
                         const SVector3 &oScatteringVec,
                         const Float &fScatteringVecMag,
                         Float fBeamEnergy,
                         Float fBeamDeflectionChiLaue)
{
    // NOTE: reciprocal vectors are measured in angstrom
    // k = 2pi/lambda = E / (h-bar c)

    Float fWavenumber = PhysicalConstants::keV_over_hbar_c_in_ang * fBeamEnergy;

    Float fSinTheta = fScatteringVecMag / ((Float)2.0 * fWavenumber);   // Bragg angle
    Float fCosChi = oScatteringVec.m_fZ / fScatteringVecMag;             // Tilt angle of G relative to z-axis
    Float fSinChi = sqrt((Float)1.0 - fCosChi * fCosChi);

    Float fSinChiLaue = sin(fBeamDeflectionChiLaue);     // Tilt angle of k_i (+ means up)
    Float fCosChiLaue = cos(fBeamDeflectionChiLaue);

    Float fNumerator = fSinTheta + fCosChi * fSinChiLaue;
    Float fDenom     = fSinChi * fCosChiLaue;

    if (fabs(fNumerator) <= fabs(fDenom))
    {
        // [-pi:pi]: angle to bring G to nominal position along +y-axis
        Float fDeltaOmega0 = atan2(oScatteringVec.m_fX, oScatteringVec.m_fY);

        //  [0:pi/2] since arg >0: phi goes from above to Bragg angle
        Float fDeltaOmega_b1 = asin(fNumerator / fDenom);

        Float fDeltaOmega_b2 = PI - fDeltaOmega_b1;

        fOmegaRes1 = fDeltaOmega_b1 + fDeltaOmega0;  // oScatteringVec.m_fY > 0
        fOmegaRes2 = fDeltaOmega_b2 + fDeltaOmega0;  // oScatteringVec.m_fY < 0

        if (fOmegaRes1 > PI)          // range really doesn't matter
            fOmegaRes1 -= 2.f * PI;

        if (fOmegaRes1 < -PI)
            fOmegaRes1 += 2.f * PI;

        if (fOmegaRes2 > PI)
            fOmegaRes2 -= 2.f * PI;

        if (fOmegaRes2 < -PI)
            fOmegaRes2 += 2.f * PI;

        return true;
    }
    else
    {
        fOmegaRes1 = fOmegaRes2 = 0;     // too close to rotation axis to be illuminated
        return false;
    }
}

// JSON helpers
void write_matrix3x3_json(ofstream& out, const SMatrix3x3& m, bool last = false) {
    out << "    [";
    for (int i = 0; i < 3; i++) {
        out << "[";
        for (int j = 0; j < 3; j++) {
            out << fixed << setprecision(15) << m.m[i][j];
            if (j < 2) out << ", ";
        }
        out << "]";
        if (i < 2) out << ", ";
    }
    out << "]";
    if (!last) out << ",";
    out << "\n";
}

void write_vector3_json(ofstream& out, const SVector3& v) {
    out << "[" << fixed << setprecision(15)
        << v.m_fX << ", " << v.m_fY << ", " << v.m_fZ << "]";
}

// Generate all Miller indices up to max_q
struct HKL {
    int h, k, l;
    Float q_mag;
};

vector<HKL> generate_all_hkl(Float lattice_a, Float max_q) {
    vector<HKL> hkls;

    // Reciprocal lattice parameter for cubic: a* = 2π/a
    Float a_recip = TWO_PI / lattice_a;

    // Maximum h,k,l to check
    int max_index = (int)ceil(max_q / a_recip) + 1;

    for (int h = -max_index; h <= max_index; h++) {
        for (int k = -max_index; k <= max_index; k++) {
            for (int l = -max_index; l <= max_index; l++) {
                if (h == 0 && k == 0 && l == 0) continue;

                // For FCC: h,k,l all even or all odd (systematic absences)
                bool all_even = (h % 2 == 0) && (k % 2 == 0) && (l % 2 == 0);
                bool all_odd = (h % 2 != 0) && (k % 2 != 0) && (l % 2 != 0);
                if (!(all_even || all_odd)) continue;

                // Calculate |Q| = a* * sqrt(h^2 + k^2 + l^2)
                Float q_mag = a_recip * sqrt(Float(h*h + k*k + l*l));

                if (q_mag <= max_q) {
                    HKL hkl = {h, k, l, q_mag};
                    hkls.push_back(hkl);
                }
            }
        }
    }

    return hkls;
}

// Check if two vectors are equivalent under symmetry
bool are_equivalent(const CCubicSymmetry& sym, const SVector3& v1, const SVector3& v2, Float tol = 1e-5) {
    const vector<SMatrix3x3>& ops = sym.GetOperatorList();

    for (size_t i = 0; i < ops.size(); i++) {
        SVector3 rotated = ops[i] * v1;
        SVector3 diff = rotated - v2;
        if (diff.GetLength() < tol) {
            return true;
        }
    }
    return false;
}

// Get unique reflections after symmetry reduction
vector<HKL> get_unique_reflections(const vector<HKL>& all_hkls, const CCubicSymmetry& sym) {
    vector<HKL> unique;

    for (size_t i = 0; i < all_hkls.size(); i++) {
        SVector3 v_current(all_hkls[i].h, all_hkls[i].k, all_hkls[i].l);
        bool is_unique = true;

        for (size_t j = 0; j < unique.size(); j++) {
            SVector3 v_unique(unique[j].h, unique[j].k, unique[j].l);
            if (are_equivalent(sym, v_current, v_unique)) {
                is_unique = false;
                break;
            }
        }

        if (is_unique) {
            unique.push_back(all_hkls[i]);
        }
    }

    return unique;
}

int main() {
    cout << "IceNine C++ Test Data Generator\n";
    cout << "================================\n\n";

    // 1. Generate cubic symmetry operators
    cout << "Generating cubic symmetry operators..." << endl;
    CCubicSymmetry& cubic_sym = CCubicSymmetry::Get();
    const vector<SMatrix3x3>& sym_ops = cubic_sym.GetOperatorList();

    cout << "  Found " << sym_ops.size() << " symmetry operators\n";

    // Write to JSON
    ofstream out_sym("../cpp_outputs/cubic_symmetry_24ops.json");
    out_sym << "{\n";
    out_sym << "  \"description\": \"24 cubic symmetry operators (point group m-3m)\",\n";
    out_sym << "  \"count\": " << sym_ops.size() << ",\n";
    out_sym << "  \"matrices\": [\n";
    for (size_t i = 0; i < sym_ops.size(); i++) {
        write_matrix3x3_json(out_sym, sym_ops[i], i == sym_ops.size() - 1);
    }
    out_sym << "  ]\n";
    out_sym << "}\n";
    out_sym.close();
    cout << "  Written: cubic_symmetry_24ops.json\n\n";

    // 2. Generate vector equivalence test cases
    cout << "Generating vector equivalence test cases..." << endl;
    ofstream out_equiv("../cpp_outputs/vector_equivalence_tests.json");
    out_equiv << "{\n";
    out_equiv << "  \"description\": \"Test cases for vector equivalence under cubic symmetry\",\n";
    out_equiv << "  \"test_cases\": [\n";

    struct TestCase { int h1, k1, l1, h2, k2, l2; bool expected; };
    vector<TestCase> tests = {
        {1, 1, 1,  1, -1, -1, true},   // Equivalent for cubic
        {1, 1, 1, -1,  1, -1, true},   // Equivalent for cubic
        {1, 0, 0,  0,  1,  0, true},   // Equivalent for cubic
        {2, 0, 0,  0,  0,  2, true},   // Equivalent for cubic
        {1, 0, 0,  1,  1,  0, false},  // Not equivalent
        {1, 1, 1,  2,  0,  0, false},  // Not equivalent
    };

    for (size_t i = 0; i < tests.size(); i++) {
        SVector3 v1(tests[i].h1, tests[i].k1, tests[i].l1);
        SVector3 v2(tests[i].h2, tests[i].k2, tests[i].l2);
        bool result = are_equivalent(cubic_sym, v1, v2);

        out_equiv << "    {\n";
        out_equiv << "      \"v1\": [" << tests[i].h1 << ", " << tests[i].k1 << ", " << tests[i].l1 << "],\n";
        out_equiv << "      \"v2\": [" << tests[i].h2 << ", " << tests[i].k2 << ", " << tests[i].l2 << "],\n";
        out_equiv << "      \"expected\": " << (tests[i].expected ? "true" : "false") << ",\n";
        out_equiv << "      \"computed\": " << (result ? "true" : "false") << "\n";
        out_equiv << "    }";
        if (i < tests.size() - 1) out_equiv << ",";
        out_equiv << "\n";
    }
    out_equiv << "  ]\n";
    out_equiv << "}\n";
    out_equiv.close();
    cout << "  Written: vector_equivalence_tests.json\n\n";

    // 3. Generate all Miller indices for Au FCC
    cout << "Generating Miller indices for Au FCC (Q_max = " << MAX_Q << " Å⁻¹)..." << endl;
    vector<HKL> all_hkls = generate_all_hkl(AU_LATTICE_A, MAX_Q);
    cout << "  Generated " << all_hkls.size() << " reflections (before symmetry reduction)\n";

    // Write all HKLs to JSON
    ofstream out_all_hkl("../cpp_outputs/au_fcc_all_hkl_qmax8.json");
    out_all_hkl << "{\n";
    out_all_hkl << "  \"description\": \"All Miller indices for Au FCC before symmetry reduction\",\n";
    out_all_hkl << "  \"lattice_a\": " << AU_LATTICE_A << ",\n";
    out_all_hkl << "  \"max_q\": " << MAX_Q << ",\n";
    out_all_hkl << "  \"count\": " << all_hkls.size() << ",\n";
    out_all_hkl << "  \"reflections\": [\n";
    for (size_t i = 0; i < all_hkls.size(); i++) {
        out_all_hkl << "    {\"h\": " << all_hkls[i].h
                    << ", \"k\": " << all_hkls[i].k
                    << ", \"l\": " << all_hkls[i].l
                    << ", \"q_mag\": " << fixed << setprecision(15) << all_hkls[i].q_mag << "}";
        if (i < all_hkls.size() - 1) out_all_hkl << ",";
        out_all_hkl << "\n";
    }
    out_all_hkl << "  ]\n";
    out_all_hkl << "}\n";
    out_all_hkl.close();
    cout << "  Written: au_fcc_all_hkl_qmax8.json\n\n";

    // 4. Get unique reflections after symmetry reduction
    cout << "Applying symmetry reduction..." << endl;
    vector<HKL> unique_hkls = get_unique_reflections(all_hkls, cubic_sym);
    cout << "  Unique reflections: " << unique_hkls.size() << "\n";

    // Write unique HKLs to JSON
    ofstream out_unique_hkl("../cpp_outputs/au_fcc_unique_hkl_qmax8.json");
    out_unique_hkl << "{\n";
    out_unique_hkl << "  \"description\": \"Unique Miller indices for Au FCC after symmetry reduction\",\n";
    out_unique_hkl << "  \"lattice_a\": " << AU_LATTICE_A << ",\n";
    out_unique_hkl << "  \"max_q\": " << MAX_Q << ",\n";
    out_unique_hkl << "  \"count\": " << unique_hkls.size() << ",\n";
    out_unique_hkl << "  \"reflections\": [\n";
    for (size_t i = 0; i < unique_hkls.size(); i++) {
        out_unique_hkl << "    {\"h\": " << unique_hkls[i].h
                       << ", \"k\": " << unique_hkls[i].k
                       << ", \"l\": " << unique_hkls[i].l
                       << ", \"q_mag\": " << fixed << setprecision(15) << unique_hkls[i].q_mag << "}";
        if (i < unique_hkls.size() - 1) out_unique_hkl << ",";
        out_unique_hkl << "\n";
    }
    out_unique_hkl << "  ]\n";
    out_unique_hkl << "}\n";
    out_unique_hkl.close();
    cout << "  Written: au_fcc_unique_hkl_qmax8.json\n\n";

    // 5. Generate Q magnitudes for specific test reflections
    cout << "Generating Q magnitudes for test reflections..." << endl;
    vector<HKL> test_reflections = {
        {1, 1, 1, 0}, {2, 0, 0, 0}, {2, 2, 0, 0}, {3, 1, 1, 0}, {2, 2, 2, 0}
    };

    Float a_recip = TWO_PI / AU_LATTICE_A;
    ofstream out_q("../cpp_outputs/q_magnitudes.json");
    out_q << "{\n";
    out_q << "  \"description\": \"Q magnitudes for test reflections\",\n";
    out_q << "  \"lattice_a\": " << AU_LATTICE_A << ",\n";
    out_q << "  \"a_reciprocal\": " << fixed << setprecision(15) << a_recip << ",\n";
    out_q << "  \"reflections\": [\n";
    for (size_t i = 0; i < test_reflections.size(); i++) {
        int h = test_reflections[i].h;
        int k = test_reflections[i].k;
        int l = test_reflections[i].l;
        Float q_mag = a_recip * sqrt(Float(h*h + k*k + l*l));

        out_q << "    {\"h\": " << h << ", \"k\": " << k << ", \"l\": " << l
              << ", \"q_mag\": " << fixed << setprecision(15) << q_mag << "}";
        if (i < test_reflections.size() - 1) out_q << ",";
        out_q << "\n";
    }
    out_q << "  ]\n";
    out_q << "}\n";
    out_q.close();
    cout << "  Written: q_magnitudes.json\n\n";

    // 6. Generate scattering omega angles for test reflections
    cout << "Generating scattering omega angles with random rotations..." << endl;

    // Test reflections with various orientations to get comprehensive test coverage
    struct OmegaTestCase {
        int h, k, l;
        Float gx, gy, gz;  // scattering vector in sample frame
        Float qw, qx, qy, qz;  // quaternion that produced this rotation
    };

    vector<OmegaTestCase> omega_tests;

    // Test on unique reflections with identity orientation (no rotation)
    cout << "  Adding " << unique_hkls.size() << " test cases with identity orientation..." << endl;
    for (size_t i = 0; i < unique_hkls.size(); i++) {
        OmegaTestCase test;
        test.h = unique_hkls[i].h;
        test.k = unique_hkls[i].k;
        test.l = unique_hkls[i].l;

        // G vector = a* * (h, k, l) for cubic
        test.gx = a_recip * test.h;
        test.gy = a_recip * test.k;
        test.gz = a_recip * test.l;

        // Identity quaternion
        test.qw = 1.0;
        test.qx = 0.0;
        test.qy = 0.0;
        test.qz = 0.0;

        omega_tests.push_back(test);
    }

    // Add ~5000 random rotations applied to a subset of reflections
    // We'll use a few representative reflections and apply many random rotations to each
    const int NUM_RANDOM_TESTS = 5000;
    cout << "  Generating " << NUM_RANDOM_TESTS << " test cases with random SO(3) rotations..." << endl;

    // Select a few representative reflections to rotate
    vector<HKL> base_reflections;
    if (unique_hkls.size() > 0) base_reflections.push_back(unique_hkls[0]);  // First unique
    if (unique_hkls.size() > 1) base_reflections.push_back(unique_hkls[1]);  // Second unique
    if (unique_hkls.size() > unique_hkls.size()/2) {
        base_reflections.push_back(unique_hkls[unique_hkls.size()/2]);  // Middle
    }
    if (unique_hkls.size() > 2) {
        base_reflections.push_back(unique_hkls[unique_hkls.size()-1]);  // Last unique
    }

    // Add some specific interesting reflections if available
    for (size_t i = 0; i < unique_hkls.size(); i++) {
        const HKL& hkl = unique_hkls[i];
        if ((hkl.h == 1 && hkl.k == 1 && hkl.l == 1) ||
            (hkl.h == 2 && hkl.k == 0 && hkl.l == 0) ||
            (hkl.h == 2 && hkl.k == 2 && hkl.l == 0)) {
            // Check if not already added
            bool already_added = false;
            for (const auto& br : base_reflections) {
                if (br.h == hkl.h && br.k == hkl.k && br.l == hkl.l) {
                    already_added = true;
                    break;
                }
            }
            if (!already_added) {
                base_reflections.push_back(hkl);
            }
        }
    }

    // Initialize random rotation generator
    CRandomRotationGenerator random_rotation_gen;

    // Generate random rotations
    int tests_per_reflection = NUM_RANDOM_TESTS / base_reflections.size();
    if (tests_per_reflection < 1) tests_per_reflection = 1;

    for (size_t i = 0; i < base_reflections.size(); i++) {
        int num_tests = (i == base_reflections.size() - 1)
            ? (NUM_RANDOM_TESTS - omega_tests.size() + unique_hkls.size())  // Adjust last batch to hit target
            : tests_per_reflection;

        for (int j = 0; j < num_tests; j++) {
            // Generate random rotation
            SQuaternion random_quat = random_rotation_gen.GetRandomQuaternion();
            SMatrix3x3 rotation_matrix = random_quat.GetRotationMatrix3x3();

            // Base g-vector in crystal frame
            SVector3 g_crystal(a_recip * base_reflections[i].h,
                              a_recip * base_reflections[i].k,
                              a_recip * base_reflections[i].l);

            // Rotate g-vector to sample frame
            SVector3 g_sample = rotation_matrix * g_crystal;

            OmegaTestCase test;
            test.h = base_reflections[i].h;
            test.k = base_reflections[i].k;
            test.l = base_reflections[i].l;
            test.gx = g_sample.m_fX;
            test.gy = g_sample.m_fY;
            test.gz = g_sample.m_fZ;
            test.qw = random_quat.m_fW;
            test.qx = random_quat.m_fX;
            test.qy = random_quat.m_fY;
            test.qz = random_quat.m_fZ;

            omega_tests.push_back(test);
        }
    }

    cout << "  Total test cases: " << omega_tests.size() << endl;

    ofstream out_omega("../cpp_outputs/scattering_omegas.json");
    out_omega << "{\n";
    out_omega << "  \"description\": \"Scattering omega angles for test reflections with random SO(3) rotations\",\n";
    out_omega << "  \"beam_energy\": " << BEAM_ENERGY << ",\n";
    out_omega << "  \"beam_deflection_chi\": " << BEAM_DEFLECTION_CHI << ",\n";
    out_omega << "  \"beam_direction\": [1, 0, 0],\n";
    out_omega << "  \"num_test_cases\": " << omega_tests.size() << ",\n";
    out_omega << "  \"test_cases\": [\n";

    int observable_count = 0;
    for (size_t i = 0; i < omega_tests.size(); i++) {
        SVector3 g_vec(omega_tests[i].gx, omega_tests[i].gy, omega_tests[i].gz);
        Float g_mag = g_vec.GetLength();

        Float omega1, omega2;
        bool observable = GetScatteringOmegas(omega1, omega2, g_vec, g_mag,
                                               BEAM_ENERGY, BEAM_DEFLECTION_CHI);

        if (observable) observable_count++;

        out_omega << "    {\n";
        out_omega << "      \"hkl\": [" << omega_tests[i].h << ", "
                  << omega_tests[i].k << ", " << omega_tests[i].l << "],\n";
        out_omega << "      \"quaternion\": [" << fixed << setprecision(15)
                  << omega_tests[i].qw << ", " << omega_tests[i].qx << ", "
                  << omega_tests[i].qy << ", " << omega_tests[i].qz << "],\n";
        out_omega << "      \"g_vector\": [" << fixed << setprecision(15)
                  << g_vec.m_fX << ", " << g_vec.m_fY << ", " << g_vec.m_fZ << "],\n";
        out_omega << "      \"g_magnitude\": " << fixed << setprecision(15) << g_mag << ",\n";
        out_omega << "      \"observable\": " << (observable ? "true" : "false");

        if (observable) {
            out_omega << ",\n";
            out_omega << "      \"omega1\": " << fixed << setprecision(15) << omega1 << ",\n";
            out_omega << "      \"omega2\": " << fixed << setprecision(15) << omega2 << "\n";
        } else {
            out_omega << "\n";
        }

        out_omega << "    }";
        if (i < omega_tests.size() - 1) out_omega << ",";
        out_omega << "\n";
    }

    out_omega << "  ]\n";
    out_omega << "}\n";
    out_omega.close();

    cout << "  Generated " << omega_tests.size() << " omega angle test cases\n";
    cout << "  Observable peaks: " << observable_count << "\n";
    cout << "  Written: scattering_omegas.json\n\n";

    cout << "Test data generation complete!\n";
    cout << "All files written to cpp_outputs/\n";

    return 0;
}
