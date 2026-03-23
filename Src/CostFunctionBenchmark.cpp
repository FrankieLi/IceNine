//==============================================================================
//  CostFunctionBenchmark.cpp
//
//  Standalone benchmark for profiling the IceNine cost function at per-operation
//  granularity. Measures:
//    1. Single GetScatteringOmegas call
//    2. GetObservablePeaks (all reciprocal vectors)
//    3. Single peak overlap: 0 detectors hit
//    4. Single peak overlap: 1 detector hit
//    5. Single peak overlap: 2 detectors hit
//    6. Full VoxelCostFunction::operator() (single voxel)
//
//  Usage: ./CostFunctionBenchmark <ConfigFile>
//  Example: ./CostFunctionBenchmark Examples/Example2.ThreeVoxels/ConfigFiles/ReconstructBenchmark.config
//==============================================================================

#include "ConfigFile.h"
#include "ExperimentSetup.h"
#include "ReconstructionSetup.h"
#include "Simulation.h"
#include "SearchTraits.h"
#include "CostFunctions.h"
#include "OverlapInfo.h"
#include "Voxel.h"
#include "Sample.h"
#include "PeakFilters.h"
#include "DiffractionCore.h"
#include "MicIO.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace CostFunctions;

//------------------------------------------------------------------------------
//  Timer helper
//------------------------------------------------------------------------------
struct BenchResult {
  double mean_us;    // microseconds
  double stddev_us;
  double min_us;
  int n_iters;
};

template<typename Fn>
BenchResult benchmark(const string& name, int n_iters, Fn fn)
{
  // Warmup
  for (int i = 0; i < std::min(3, n_iters); i++)
    fn();

  vector<double> times;
  times.reserve(n_iters);

  for (int i = 0; i < n_iters; i++) {
    auto t0 = chrono::high_resolution_clock::now();
    fn();
    auto t1 = chrono::high_resolution_clock::now();
    double us = chrono::duration<double, micro>(t1 - t0).count();
    times.push_back(us);
  }

  double sum = accumulate(times.begin(), times.end(), 0.0);
  double mean = sum / n_iters;

  double sq_sum = 0;
  for (auto t : times) sq_sum += (t - mean) * (t - mean);
  double stddev = sqrt(sq_sum / n_iters);

  double min_val = *min_element(times.begin(), times.end());

  cout << fixed << setprecision(2);
  cout << "  " << setw(45) << left << name
       << "  mean=" << setw(10) << right << mean << " us"
       << "  stddev=" << setw(8) << stddev << " us"
       << "  min=" << setw(10) << min_val << " us"
       << "  (n=" << n_iters << ")" << endl;

  return {mean, stddev, min_val, n_iters};
}

//------------------------------------------------------------------------------
//  Main
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <ConfigFile>" << endl;
    return 1;
  }

  string sConfigFile = argv[1];
  cout << "=== IceNine Cost Function Benchmark ===" << endl;
  cout << "Config: " << sConfigFile << endl;
  cout << endl;

  //--------------------------------------------
  // Step 1: Initialize everything
  //--------------------------------------------
  cout << "--- Initialization ---" << endl;

  CConfigFile oConfigFile;
  if (!oConfigFile.InputConfigParameters(sConfigFile)) {
    cerr << "ERROR: Failed to read config file: " << sConfigFile << endl;
    return 1;
  }

  Reconstruction::ReconstructionSetup oSetup;
  oSetup.InitializeWithDataFiles(oConfigFile);
  cout << "  Config loaded, data files read." << endl;

  const CXDMExperimentSetup& oExpSetup = oSetup.ExperimentalSetup();
  CSample& oSample = oSetup.SampleGeometry();
  const CSimulationData& oExpData = oSetup.Data();
  const vector<CDetector>& oDetList = oExpSetup.GetDetectorList();
  const XDMSimulation::CSimulationRange& oRangeMap = oExpSetup.GetRangeToIndexMap();

  CSimulation oSimulator(oExpSetup);
  cout << "  Simulator initialized." << endl;
  cout << "  Detectors: " << oDetList.size() << endl;

  //--------------------------------------------
  // Step 2: Get voxel and reciprocal vectors
  //--------------------------------------------
  auto pMicBase = oSample.GetMic();
  auto pMic = std::dynamic_pointer_cast<CMic>(pMicBase);
  if (!pMic) {
    cerr << "ERROR: Could not cast MicIOBase to CMic" << endl;
    return 1;
  }
  const vector<SVoxel>& voxels = pMic->GetVoxels();
  cout << "  Voxels in sample: " << voxels.size() << endl;

  if (voxels.size() < 3) {
    cerr << "ERROR: Need at least 3 voxels (using voxel index 2)" << endl;
    return 1;
  }

  // Use voxel 2 (the one that reconstructs correctly in validation)
  SVoxel oTestVoxel = voxels[2];
  cout << "  Using voxel 2, phase=" << oTestVoxel.nPhase << endl;
  cout << "  Voxel center: ("
       << oTestVoxel.GetCenter().m_fX << ", "
       << oTestVoxel.GetCenter().m_fY << ", "
       << oTestVoxel.GetCenter().m_fZ << ")" << endl;
  cout << "  Orientation matrix:" << endl;
  for (int r = 0; r < 3; r++) {
    cout << "    [";
    for (int c = 0; c < 3; c++)
      cout << setw(10) << setprecision(5) << oTestVoxel.oOrientMatrix.m[r][c];
    cout << " ]" << endl;
  }

  const vector<CRecpVector>& oRecipVectors =
    oSample.GetStructureList()[oTestVoxel.nPhase].GetReflectionVectorList();
  cout << "  Reciprocal vectors: " << oRecipVectors.size() << endl;

  //--------------------------------------------
  // Step 3: Build cost function (same type chain as Reconstructor.h)
  //--------------------------------------------
  typedef Reconstruction::ReconstructionSetup::EtaAngularFilter PeakFilterT;
  typedef OrientationSearch::GeneralSearchCostFn<
    PeakFilterT, CSimulation,
    CostFunctions::TrianglePixelOverlapCounter,
    CostFunctions::DetectorOverlapCounter,
    CostFunctions::XDMSampleVertexOverlapCounter,
    CostFunctions::VoxelToVertices,
    CostFunctions::VoxelCostFunction, SVoxel
  > LocalSearchCostFunctions;

  PeakFilterT oFilter = oSetup.EtaThresholdFilter();
  LocalSearchCostFunctions oSearchCostFn(oFilter, oSimulator);
  auto& oVoxelCostFn = oSearchCostFn.VoxelCostFn();

  cout << endl;

  //--------------------------------------------
  // Benchmark 1: Single GetScatteringOmegas call
  //--------------------------------------------
  cout << "--- Benchmark 1: Single GetScatteringOmegas ---" << endl;
  {
    SVector3 g = oTestVoxel.oOrientMatrix * oRecipVectors[0].v;
    Float fMag = oRecipVectors[0].fMag;
    Float omega1, omega2;
    int N = 100000;

    benchmark("GetScatteringOmegas (single)", N, [&]() {
      oSimulator.GetScatteringOmegas(omega1, omega2, g, fMag);
    });

    // Also benchmark all reciprocal vectors sequentially
    benchmark("GetScatteringOmegas (all " + to_string(oRecipVectors.size()) + " recip vecs)", 1000, [&]() {
      Float o1, o2;
      for (size_t i = 0; i < oRecipVectors.size(); i++) {
        SVector3 gi = oTestVoxel.oOrientMatrix * oRecipVectors[i].v;
        oSimulator.GetScatteringOmegas(o1, o2, gi, oRecipVectors[i].fMag);
      }
    });
  }
  cout << endl;

  //--------------------------------------------
  // Benchmark 2: GetObservablePeaks (full)
  //--------------------------------------------
  cout << "--- Benchmark 2: GetObservablePeaks ---" << endl;
  vector<SPeakInfo> oPeakInfoList;
  {
    int N = 1000;
    benchmark("GetObservablePeaks (all recip vecs)", N, [&]() {
      oPeakInfoList.clear();
      oPeakInfoList.reserve(oRecipVectors.size() * 2);
      oSimulator.GetObservablePeaks(oPeakInfoList, oTestVoxel.oOrientMatrix,
                                     oRecipVectors, oDetList);
    });
    cout << "  Observable peaks generated: " << oPeakInfoList.size() << endl;
  }
  cout << endl;

  //--------------------------------------------
  // Benchmark 3-5: Per-peak overlap by detector hit count
  //
  // We run GenerateProjectedPixels + DetOverlapCounter for individual peaks
  // to measure the cost of overlap computation by detector hit category.
  //--------------------------------------------
  cout << "--- Benchmark 3-5: Per-peak overlap (categorized by detector hits) ---" << endl;
  {
    // First, classify peaks by how many detectors they actually hit
    vector<SVector3> oVertexList = oSearchCostFn.VertexExtractorFn()(oTestVoxel);

    int n_hit_0 = 0, n_hit_1 = 0, n_hit_2 = 0;
    int peak_0_det = -1, peak_1_det = -1, peak_2_det = -1;

    // Use the full cost function internals to classify peaks
    const SMatrix3x3 oCurOrientation = oSample.GetOrientationMatrix();
    for (size_t p = 0; p < oPeakInfoList.size(); p++) {
      if (!oPeakInfoList[p].bObservable) continue;
      Size_Type nOmegaIndex = oRangeMap(oPeakInfoList[p].fOmega);
      if (nOmegaIndex == XDMSimulation::NoMatch) continue;

      // Save/restore sample orientation
      oSample.RotateZ(oPeakInfoList[p].fOmega);
      const SVector3& oNormal = oPeakInfoList[p].oScatteringDir;

      int dets_hit = 0;
      for (size_t d = 0; d < oDetList.size(); d++) {
        bool all_hit = true;
        for (int vi = 0; vi < 3; vi++) {
          CRay oReflectedRay = DiffractionCore::BuildReflectedRay(
            oSample, oVertexList[vi],
            DiffractionCore::GetReflectedRayDir(oSample, oNormal, oExpSetup.GetXrayBeamDirection()));
          Point pixel;
          bool hit = DiffractionCore::GetIlluminatedPixel(pixel, oDetList[d], oReflectedRay);
          if (!hit) { all_hit = false; break; }
        }
        if (all_hit) dets_hit++;
      }

      if (dets_hit == 0 && peak_0_det < 0) peak_0_det = p;
      if (dets_hit == 1 && peak_1_det < 0) peak_1_det = p;
      if (dets_hit == 2 && peak_2_det < 0) peak_2_det = p;
      if (dets_hit == 0) n_hit_0++;
      else if (dets_hit == 1) n_hit_1++;
      else n_hit_2++;

      oSample.SetOrientation(oCurOrientation);
    }

    cout << "  Peak classification: 0-det=" << n_hit_0
         << " 1-det=" << n_hit_1 << " 2-det=" << n_hit_2 << endl;

    // Now benchmark individual peak evaluations for each category
    // We create single-peak PeakInfoLists and run CalculateDiffractionOverlap
    auto& oVertexOverlapCounter = oSearchCostFn.VertexCostFn();

    auto benchSinglePeak = [&](const string& label, int peakIdx, int N) {
      if (peakIdx < 0) {
        cout << "  " << setw(45) << left << label << "  SKIPPED (no such peak)" << endl;
        return;
      }
      vector<SPeakInfo> singlePeak = { oPeakInfoList[peakIdx] };
      benchmark(label, N, [&]() {
        oVertexOverlapCounter.CalculateDiffractionOverlap(
          oVertexList, singlePeak, oSample, oDetList, oRangeMap, oExpData);
      });
    };

    benchSinglePeak("Single peak: 0 detectors hit", peak_0_det, 10000);
    benchSinglePeak("Single peak: 1 detector hit", peak_1_det, 10000);
    benchSinglePeak("Single peak: 2 detectors hit", peak_2_det, 10000);

    // Benchmark: all peaks via CalculateDiffractionOverlap (without GetObservablePeaks)
    benchmark("CalculateDiffractionOverlap (all " + to_string(oPeakInfoList.size()) + " peaks)", 100, [&]() {
      oVertexOverlapCounter.CalculateDiffractionOverlap(
        oVertexList, oPeakInfoList, oSample, oDetList, oRangeMap, oExpData);
    });
  }
  cout << endl;

  //--------------------------------------------
  // Benchmark 6: Full VoxelCostFunction::operator()
  //--------------------------------------------
  cout << "--- Benchmark 6: Full VoxelCostFunction (single voxel) ---" << endl;
  {
    SOverlapInfo oResult;
    int N = 100;

    benchmark("VoxelCostFunction::operator()", N, [&]() {
      oResult = oVoxelCostFn(oTestVoxel, oSample, oDetList, oRangeMap, oExpData);
    });

    Float fCost = 1.0f - oResult.fQuality;
    Float fHitRatio = (oResult.nPixelOnDetector > 0) ?
      Float(oResult.nPixelOverlap) / Float(oResult.nPixelOnDetector) : 0.0f;

    cout << "  Result: quality=" << oResult.fQuality
         << "  cost=" << fCost
         << "  hit_ratio=" << fHitRatio << endl;
    cout << "  Pixels: overlap=" << oResult.nPixelOverlap
         << "  on_det=" << oResult.nPixelOnDetector << endl;
    cout << "  Peaks: overlap=" << oResult.nPeakOverlap
         << "  on_det=" << oResult.nPeakOnDetector << endl;
  }
  cout << endl;

  //--------------------------------------------
  // Per-peak diagnostic output (for C++/Python validation)
  //--------------------------------------------
  cout << "--- Per-Peak Diagnostics ---" << endl;
  {
    vector<SVector3> oVertexList = oSearchCostFn.VertexExtractorFn()(oTestVoxel);
    auto& oVertexOverlapCounter = oSearchCostFn.VertexCostFn();

    // Count observable peaks that pass omega filter
    int nValidPeaks = 0;
    for (size_t p = 0; p < oPeakInfoList.size(); p++) {
      if (!oPeakInfoList[p].bObservable) continue;
      Size_Type nOmegaIndex = oRangeMap(oPeakInfoList[p].fOmega);
      if (nOmegaIndex == XDMSimulation::NoMatch) continue;
      nValidPeaks++;
    }
    cout << "  Observable peaks (after eta filter): " << oPeakInfoList.size()
         << " (valid omega: " << nValidPeaks << ")" << endl;
    cout << endl;

    cout << "  " << setw(4) << right << "Peak" << "  "
         << setw(10) << right << "omega_deg" << "  "
         << setw(8) << right << "pix_ovlp" << "  "
         << setw(8) << right << "pix_det" << "  "
         << setw(7) << right << "pk_ovlp" << "  "
         << setw(6) << right << "pk_det" << "  "
         << setw(10) << right << "n_det_ovlp" << "  "
         << setw(10) << right << "quality_i" << endl;

    cout << "  " << setw(4) << right << "----" << "  "
         << setw(10) << right << "----------" << "  "
         << setw(8) << right << "--------" << "  "
         << setw(8) << right << "--------" << "  "
         << setw(7) << right << "-------" << "  "
         << setw(6) << right << "------" << "  "
         << setw(10) << right << "----------" << "  "
         << setw(10) << right << "----------" << endl;

    int peakCount = 0;
    for (size_t p = 0; p < oPeakInfoList.size(); p++) {
      if (!oPeakInfoList[p].bObservable) continue;

      vector<SPeakInfo> singlePeak = { oPeakInfoList[p] };
      SOverlapInfo oPeakResult = oVertexOverlapCounter.CalculateDiffractionOverlap(
        oVertexList, singlePeak, oSample, oDetList, oRangeMap, oExpData);

      Float omega_deg = oPeakInfoList[p].fOmega * 180.0 / M_PI;
      Float quality_i = oPeakResult.fQuality;
      if (oPeakResult.nPixelOnDetector == 0)
        quality_i = -1.0;

      cout << "  " << setw(4) << right << peakCount << "  "
           << setw(10) << setprecision(4) << fixed << right << omega_deg << "  "
           << setw(8) << right << oPeakResult.nPixelOverlap << "  "
           << setw(8) << right << oPeakResult.nPixelOnDetector << "  "
           << setw(7) << right << oPeakResult.nPeakOverlap << "  "
           << setw(6) << right << oPeakResult.nPeakOnDetector << "  "
           << setw(10) << right << oPeakResult.nDetectorsOverlap << "  "
           << setw(10) << setprecision(6) << right << quality_i << endl;

      peakCount++;
    }

    cout << endl;

    // Recompute aggregate quality from per-peak data
    cout << "--- Aggregate Quality Recomputation ---" << endl;
    Float runningQuality = 0.0;
    int nPoints = 0;
    peakCount = 0;
    for (size_t p = 0; p < oPeakInfoList.size(); p++) {
      if (!oPeakInfoList[p].bObservable) continue;

      vector<SPeakInfo> singlePeak = { oPeakInfoList[p] };
      SOverlapInfo oPeakResult = oVertexOverlapCounter.CalculateDiffractionOverlap(
        oVertexList, singlePeak, oSample, oDetList, oRangeMap, oExpData);

      if (oPeakResult.nPixelOnDetector > 0) {
        Float curQuality = oPeakResult.fQuality;
        runningQuality += (curQuality - runningQuality) / Float(nPoints + 1);
        nPoints++;
      }
      peakCount++;
    }

    cout << "  Recomputed quality: " << setprecision(6) << runningQuality
         << "  (n_points=" << nPoints << ")" << endl;

    // Re-run full eval to show comparison
    SOverlapInfo oFinalResult = oVoxelCostFn(oTestVoxel, oSample, oDetList, oRangeMap, oExpData);
    cout << "  VoxelCostFunction quality: " << setprecision(6) << oFinalResult.fQuality
         << "  (n_points implied from peaks)" << endl;
  }
  cout << endl;

  //--------------------------------------------
  // Summary
  //--------------------------------------------
  cout << "=== Benchmark Complete ===" << endl;

  return 0;
}
