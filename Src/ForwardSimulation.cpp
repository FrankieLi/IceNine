//==============================================================================
// Copyright (c) 2014, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by S. F. Li (li31@llnl.gov)
// LLNL-CODE-657639
// All rights reserved.
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the Lawrence Livermore National Lab nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL LAB BE LIABLE FOR
// ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==============================================================================

//------------------------------------------------------------------------------------
//  Author:  S. F. Li (Frankie)
//  e-mail:  li31@llnl.gov; sfli@cmu.edu
//------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//  File:   ForwardSimulation.cpp
//
//
//--------------------------------------------------------------------------------------------------------

#include "ForwardSimulation.h"

//--------------------------------------------------------------------------------------------------------
//  Public:
//  CXDMForwardSimulation::CXDMForwardSimulation(const string & sConfigFilename)
//  Constructor
//
//  input:  sConfigFilename - contains file name of the config file with
//  experimental setup and simulation
//                           geometry
//
//--------------------------------------------------------------------------------------------------------
CXDMForwardSimulation::CXDMForwardSimulation(const CConfigFile &oConfigFile) {
  oSetupFile = oConfigFile;  // to be removed
  oExpSetup.SetConfigFile(oConfigFile);
}

//--------------------------------------------------------------------------------------------------------
//
//  CXDMForwardSimulation:: SimulateDetectorImagesOptimized
//
//--------------------------------------------------------------------------------------------------------
bool CXDMForwardSimulation::SimulateDetectorImagesOptimized(
    bool bStrainEnabled) {
  CSample oCurrentLayer;
  oExpSetup.InitializeExperiment();  // read in data files

  const vector<SRange> &vOmegaRangeList = oExpSetup.GetOmegaRangeList();
  const vector<SIntRange> &vFileRangeList = oExpSetup.GetFileRangeList();
  const vector<CDetector> &oDetectorList = oExpSetup.GetDetectorList();
  const CSimulationRange &oRangeToIndexMap = oExpSetup.GetRangeToIndexMap();

  oSimulator.Initialize(oExpSetup);

  oExpSetup.InitializeSample(
      oCurrentLayer,
      oDetectorList[0]);  // binding CurrentLayer to the first detector's limit

  //----------------------------
  //  DEBUG
  //----------------------------
  //   std::cout << "Recip Metric Tensor" << std::endl;
  //   const vector<CUnitCell> & oCryStructList =
  //   oCurrentLayer.GetStructureList();

  //   for( Size_Type i = 0; i < oCryStructList.size(); i ++ )
  //   {
  //     std::cout << i << "--------------------- " << std::endl;
  //     std::cout << oCryStructList[i].GetRecipocalMetric() << std::endl;
  //   }

  //----------------------------

  ImageMap oSimData;
  oSimData.resize(boost::extents[vOmegaRangeList.size()][oDetectorList.size()]);
  for (Size_Type i = 0; i < vOmegaRangeList.size(); i++)
    for (Size_Type j = 0; j < oDetectorList.size(); j++) {
      oSimData[i][j].Resize(oDetectorList[j].GetNumCols(),
                            oDetectorList[j].GetNumRows());
      oSimData[i][j].Fill(0);
    }
  std::cout << "Begin Simulation " << std::endl;
  if (bStrainEnabled)
    SimulatePeaksWithStrain_Debug(oSimData, oDetectorList, oCurrentLayer,
                                  oRangeToIndexMap);
  else
    SimulatePeaks(oSimData, oDetectorList, oCurrentLayer, oRangeToIndexMap);
  std::cout << "Finished Simulation " << std::endl;

  //
  //  Output
  //
  for (Size_Type nDetNum = 0; nDetNum < oDetectorList.size(); nDetNum++) {
#pragma omp parallel for  // parallize through open mp
    for (Size_Type i = 0; i < vOmegaRangeList.size(); i++) {
      Int nCurrentFileIndex = vFileRangeList[nDetNum].nLow + i;
      stringstream tmpSS;
      tmpSS << oSetupFile.OutFileBasename
            << InitFileIO::NumToSuffix(nCurrentFileIndex,
                                       oSetupFile.OutFileSerialLength)
            << "." << oSetupFile.OutFileExt << nDetNum;

      oSimData[i][nDetNum].PrintRaster(tmpSS.str());
      std::cout << tmpSS.str() << std::endl;
    }
  }

  return true;
}

//--------------------------------------------------------------------------------------------------------
//
//  CXDMForwardSimulation:: SimulateDetectorImagesOptimized
//
//--------------------------------------------------------------------------------------------------------
bool CXDMForwardSimulation::SimulateDetectorImagesOptimized(
    ImageMap &oSimData, const std::vector<SVoxel> &TestVoxelList,
    const bool bStrainEnabled) {
  CSample oCurrentLayer;
  oExpSetup.InitializeExperiment();  // read in data files

  const vector<SRange> &vOmegaRangeList = oExpSetup.GetOmegaRangeList();
  const vector<SIntRange> &vFileRangeList = oExpSetup.GetFileRangeList();
  const vector<CDetector> &oDetectorList = oExpSetup.GetDetectorList();
  const CSimulationRange &oRangeToIndexMap = oExpSetup.GetRangeToIndexMap();

  oSimulator.Initialize(oExpSetup);
  oExpSetup.InitializeSample(
      oCurrentLayer,
      oDetectorList[0]);  // binding CurrentLayer to the first detector's limit

  boost::shared_ptr<CMic> pMic =
      boost::dynamic_pointer_cast<CMic>(oCurrentLayer.GetMic());

  pMic->SetVoxelList(TestVoxelList);

  oSimData.resize(boost::extents[vOmegaRangeList.size()]
                                [oDetectorList.size()]);  // This is a feature!

  std::cout << "Before Zero-ing " << std::endl;
  for (Size_Type i = 0; i < vOmegaRangeList.size(); i++)
    for (Size_Type j = 0; j < oDetectorList.size(); j++) {
      oSimData[i][j].Resize(oDetectorList[j].GetNumCols(),
                            oDetectorList[j].GetNumRows());
      oSimData[i][j].Fill(0);
    }
  std::cout << "Begin Simulation " << std::endl;
  if (bStrainEnabled)
    SimulatePeaksWithStrain_Debug(oSimData, oDetectorList, oCurrentLayer,
                                  oRangeToIndexMap);
  else
    SimulatePeaks(oSimData, oDetectorList, oCurrentLayer, oRangeToIndexMap);
  std::cout << "Finished Simulation " << std::endl;

  return true;
}

//--------------------------------------------------------------------------------------------------------
//  Random accept function for peaks
//  Example Usage:
//             bool bHit = ProjectVoxel( oCurImageList, vDetectorList,
//             oCurrentLayer,
//                                       *voxelIterator, oScatteringDir,
//                                       bind( RandomAccept, _1, fIntensity ) );
//--------------------------------------------------------------------------------------------------------
std::pair<Bool, Float> RandomAccept(const SVector3 &oScatteringDir,
                                    Float fIntensity) {
  return std::pair<Bool, Float>(true, fIntensity);
}

//--------------------------------------------------------------------------------------------------------
//
//  CXDMForwardSimulation:: SimulatePeaks   -- Fast version
//
//  Note for reader:  GetReflectionVectorList is "lazy."  In another words, it
//  doesn't recalculate unless
//  it needs to.  Therefore its complexity is really O(1).
//
//--------------------------------------------------------------------------------------------------------
void CXDMForwardSimulation::SimulatePeaks(
    ImageMap &oSimData, const DetectorListT &vDetectorList,
    CSample oCurrentLayer, const CSimulationRange &oRangeToIndexMap) {
  HEDM::XDMEtaAcceptFn FAcceptFn(-oExpSetup.GetEtaLimit(),
                                 oExpSetup.GetEtaLimit());

  // i.e., E = 64 => 64 keV/(hbar c) in final unit of ang
  Float fWavenumber =
      PhysicalConstants::keV_over_hbar_c_in_ang *
      oExpSetup.GetBeamEnergy();  // (Energy (keV)/( hbar c ) - in ang.
  const vector<CUnitCell> &oCryStructList = oCurrentLayer.GetStructureList();

  //
  //  Note:  This loop cannot be parallelized due to the sparse matrix used.
  //  (Memory access and all)
  //
  boost::shared_ptr<CMic> pMic =
      boost::dynamic_pointer_cast<CMic>(oCurrentLayer.GetMic());

  int nVoxelCount = 0;
  for (vector<SVoxel>::const_iterator pCurVoxel = pMic->VoxelListBegin();
       pCurVoxel != pMic->VoxelListEnd(); pCurVoxel++) {
    nVoxelCount++;

    if (nVoxelCount % 10000 == 0) std::cout << nVoxelCount << std::endl;
    //------------------------------------------
    Int nCryStructIndex = pCurVoxel->nPhase;
    const vector<CRecpVector> &oRecipVectors =
        oCryStructList[nCryStructIndex].GetReflectionVectorList();
    //-------------------------------------------

    for (Size_Type nRecipIndex = 0; nRecipIndex < oRecipVectors.size();
         nRecipIndex++) {
      SVector3 oScatteringVec = oRecipVectors[nRecipIndex].v;
      oScatteringVec.Transform(pCurVoxel->oOrientMatrix);  // g_hkl' = O * g_hkl

      // ----------- CHECK THIS
      //
      //  This is starin in crystal frame
      //
      //  -- TO TEST THIS
      //  --  Check principal axis strains
      //  --  Check components
      //  --  Check trivial rotations
      //
      //
      // Apply strain here   g_hkl'' =  (I + S) * g_hkl -- note the *LEFT*
      // multiply -- we've taken the inverse of (I+S)
      ///

      Float fOmegaRes[2];
      //------------------------------------
      //  Magnitude is precomputed.  Rotation does not change magintude
      //------------------------------------
      bool bPeakObservable = oSimulator.GetScatteringOmegas(
          fOmegaRes[0], fOmegaRes[1], oScatteringVec,
          oRecipVectors[nRecipIndex].fMag);

      Float fSinTheta =
          oRecipVectors[nRecipIndex].fMag / (Float(2.0) * fWavenumber);
      FAcceptFn.fSin2Theta =
          sin(Float(2) * asin(fSinTheta));  // change this function object to
                                            // save x,y points instead

      if (bPeakObservable) {
        oScatteringVec.Normalize();
        const SVector3 &oScatteringDir = oScatteringVec;
        for (int i = 0; i < 2;
             i++)  // using this for loop to enforce uniformity
        {
          Size_Type nOmegaIndex = oRangeToIndexMap(fOmegaRes[i]);
          if (nOmegaIndex != XDMSimulation::NoMatch) {
            const SVector3 oCurOrientation = oCurrentLayer.GetOrientation();
            oCurrentLayer.RotateZ(fOmegaRes[i]);
            FAcceptFn.fFormIntensity = oRecipVectors[nRecipIndex].fIntensity;
            typedef boost::multi_array_types::index_range range;
            range rDetRange = range(0, vDetectorList.size());
            ImageMap::array_view<1>::type oCurImageList =
                oSimData[boost::indices[nOmegaIndex][rDetRange]];
            oSimulator.ProjectVoxel(oCurImageList, vDetectorList, oCurrentLayer,
                                    *pCurVoxel, oScatteringDir, FAcceptFn);
            oCurrentLayer.SetOrientation(
                oCurOrientation.m_fX, oCurOrientation.m_fY,
                oCurOrientation.m_fZ);  // use this to reduce numerical errors
          }
        }
      }

    }  // end for each Reciprocal Vector
  }
}

//--------------------------------------------------------------------------------------------------------
//
//  CXDMForwardSimulation:: SimulatePeaks   -- Fast version
//
//  Note for reader:  GetReflectionVectorList is "lazy."  In another words, it
//  doesn't recalculate unless
//  it needs to.  Therefore its complexity is really O(1).
//
//--------------------------------------------------------------------------------------------------------
void CXDMForwardSimulation::SimulatePeaksWithStrain_Debug(
    ImageMap &oSimData, const DetectorListT &vDetectorList,
    CSample oCurrentLayer, const CSimulationRange &oRangeToIndexMap) {
  HEDM::XDMEtaAcceptFn FAcceptFn(-oExpSetup.GetEtaLimit(),
                                 oExpSetup.GetEtaLimit());

  // i.e., E = 64 => 64 keV/(hbar c) in final unit of ang
  Float fWavenumber =
      PhysicalConstants::keV_over_hbar_c_in_ang *
      oExpSetup.GetBeamEnergy();  // (Energy (keV)/( hbar c ) - in ang.

  const vector<CUnitCell> &oCryStructList = oCurrentLayer.GetStructureList();

  //
  //  Note:  This loop cannot be parallelized due to the sparse matrix used.
  //  (Memory access and all)
  //
  boost::shared_ptr<CMic> pMic =
      boost::dynamic_pointer_cast<CMic>(oCurrentLayer.GetMic());
  for (vector<SVoxel>::const_iterator pCurVoxel = pMic->VoxelListBegin();
       pCurVoxel != pMic->VoxelListEnd(); pCurVoxel++) {
    //------------------------------------------
    Int nCryStructIndex = pCurVoxel->nPhase;
    vector<CRecpVector> oRecipVectors =
        oCryStructList[nCryStructIndex].GetReflectionVectorList();

    //-------------------------------------------

    for (Size_Type nRecipIndex = 0; nRecipIndex < oRecipVectors.size();
         nRecipIndex++) {
      SVector3 oScatteringVec = oRecipVectors[nRecipIndex].v;

      //------------------------------------------------------------------------------------
      //  Accomodate for strain debug
      //
      //  ** In this case  oDeformation is M, where M * A_0 = A, A is the "new"
      //  B matrix
      //
      //------------------------------------------------------------------------------------
      oScatteringVec =
          pCurVoxel->oDeformation * oScatteringVec;  // accomodate for strain
      oScatteringVec.Transform(pCurVoxel->oOrientMatrix);  // g_hkl' = O * g_hkl

      // ----------- CHECK THIS
      //
      //  This is starin in crystal frame
      //
      //  -- TO TEST THIS
      //  --  Check principal axis strains
      //  --  Check components
      //  --  Check trivial rotations
      //
      //
      // Apply strain here   g_hkl'' =  (I + S) * g_hkl -- note the *LEFT*
      // multiply -- we've taken the inverse of (I+S)
      ///

      Float fOmegaRes[2];
      //------------------------------------
      //  Magnitude is precomputed.  Rotation does not change magintude
      //------------------------------------
      bool bPeakObservable = oSimulator.GetScatteringOmegas(
          fOmegaRes[0], fOmegaRes[1], oScatteringVec,
          oScatteringVec.GetLength());

      Float fSinTheta =
          oRecipVectors[nRecipIndex].fMag / (Float(2.0) * fWavenumber);
      FAcceptFn.fSin2Theta =
          sin(Float(2) * asin(fSinTheta));  // change this function object to
                                            // save x,y points instead

      if (bPeakObservable) {
        oScatteringVec.Normalize();
        const SVector3 &oScatteringDir = oScatteringVec;
        for (int i = 0; i < 2;
             i++)  // using this for loop to enforce uniformity
        {
          Size_Type nOmegaIndex = oRangeToIndexMap(fOmegaRes[i]);
          if (nOmegaIndex != XDMSimulation::NoMatch) {
            const SVector3 oCurOrientation = oCurrentLayer.GetOrientation();
            oCurrentLayer.RotateZ(fOmegaRes[i]);
            FAcceptFn.fFormIntensity = oRecipVectors[nRecipIndex].fIntensity;
            typedef boost::multi_array_types::index_range range;
            range rDetRange = range(0, vDetectorList.size());
            ImageMap::array_view<1>::type oCurImageList =
                oSimData[boost::indices[nOmegaIndex][rDetRange]];
            oSimulator.ProjectVoxel(oCurImageList, vDetectorList, oCurrentLayer,
                                    *pCurVoxel, oScatteringDir, FAcceptFn);
            oCurrentLayer.SetOrientation(
                oCurOrientation.m_fX, oCurOrientation.m_fY,
                oCurOrientation.m_fZ);  // use this to reduce numerical errors
          }
        }
      }

    }  // end for each Reciprocal Vector
  }
}
