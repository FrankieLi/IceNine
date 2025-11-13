///////////////////////////////////////////////////////////////////////////////
//
//  test_neighbors.cpp
//
//  Test harness to export neighbor information from C++ MicFile for
//  validation against Python implementation.
//
//  Author: S. F. Li
//  Date: 2025-11-12
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "MicIO.h"
#include "XDMVoxel.h"

using namespace std;
using namespace XDMUtility;

///////////////////////////////////////////////////////////////////////////////
//  Calculate distance between two voxels
///////////////////////////////////////////////////////////////////////////////
double VoxelDistance(const SVoxel& v1, const SVoxel& v2) {
    double dx = v1.GetCenter().m_fX - v2.GetCenter().m_fX;
    double dy = v1.GetCenter().m_fY - v2.GetCenter().m_fY;
    double dz = v1.GetCenter().m_fZ - v2.GetCenter().m_fZ;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

///////////////////////////////////////////////////////////////////////////////
//  Find neighbors within radius for a voxel
///////////////////////////////////////////////////////////////////////////////
vector<int> FindNeighborsWithinRadius(
    const vector<SVoxel>& voxels,
    size_t query_idx,
    double radius
) {
    vector<int> neighbors;
    const SVoxel& query_voxel = voxels[query_idx];

    for (size_t i = 0; i < voxels.size(); ++i) {
        if (i == query_idx) continue;  // Skip self

        double dist = VoxelDistance(query_voxel, voxels[i]);
        if (dist <= radius) {
            neighbors.push_back(static_cast<int>(i));
        }
    }

    return neighbors;
}

///////////////////////////////////////////////////////////////////////////////
//  Write neighbors for each voxel in a MIC file
///////////////////////////////////////////////////////////////////////////////
void WriteNeighbors(const string& mic_file, const string& output_file, double radius) {
    // Read MIC file
    MicFile<SVoxel> mic;
    if (!mic.Read(mic_file)) {
        cerr << "ERROR: Failed to read MIC file: " << mic_file << endl;
        return;
    }

    const vector<SVoxel>& voxels = mic.GetVoxels();
    cout << "Read " << voxels.size() << " voxels from " << mic_file << endl;
    cout << "Search radius: " << radius << " m" << endl;

    // Open output file
    ofstream out(output_file);
    if (!out.is_open()) {
        cerr << "ERROR: Failed to open output file: " << output_file << endl;
        return;
    }

    // Write header
    out << "# Neighbor lists from C++ distance-based search" << endl;
    out << "# Format: voxel_index num_neighbors neighbor_0 neighbor_1 ..." << endl;
    out << "# Radius: " << radius << " m" << endl;
    out << "# Total voxels: " << voxels.size() << endl;

    // For each voxel, find neighbors and write to file
    int total_neighbors = 0;
    for (size_t i = 0; i < voxels.size(); ++i) {
        vector<int> neighbors = FindNeighborsWithinRadius(voxels, i, radius);

        // Write voxel index and neighbor count
        out << i << " " << neighbors.size();

        // Write neighbor indices
        for (int neighbor_idx : neighbors) {
            out << " " << neighbor_idx;
        }
        out << endl;

        total_neighbors += neighbors.size();
    }

    cout << "Wrote neighbor information to " << output_file << endl;
    cout << "Average neighbors per voxel: "
         << (double)total_neighbors / voxels.size() << endl;

    out.close();
}

///////////////////////////////////////////////////////////////////////////////
//  Main
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        cerr << "Usage: " << argv[0] << " <input.mic> <output_neighbors.txt> [radius]" << endl;
        cerr << "Example: " << argv[0] << " DataFiles/Au1007_small.mic neighbors_cpp.txt 0.05" << endl;
        cerr << "  radius: search radius in meters (default: 0.05)" << endl;
        return 1;
    }

    string mic_file = argv[1];
    string output_file = argv[2];
    double radius = 0.05;  // Default 50 mm

    if (argc == 4) {
        radius = atof(argv[3]);
    }

    cout << "C++ Neighbor Validation Test" << endl;
    cout << "============================" << endl;
    cout << "Input:  " << mic_file << endl;
    cout << "Output: " << output_file << endl;
    cout << "Radius: " << radius << " m" << endl << endl;

    WriteNeighbors(mic_file, output_file, radius);

    return 0;
}
