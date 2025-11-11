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
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL LAB BE LIABLE FOR ANY
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
////////////////////////////////////////////////////////////
//
//  ReconstructionAnalysis.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//   Purpose:  Analysis tool for HEDM recontructions. 
////////////////////////////////////////////////////////////



#ifndef _RECONSTRUCTION_ANALYSIS_
#define _RECONSTRUCTION_ANALYSIS_

#include "SampleReconstruction.h"
#include "XDMVoxel.h"
// For annoyomous functions


#include <tuple>

namespace ReconstructionAnalysis
{
  class CMicCostCalculator : public CXDMReconstruction
  {

  private:
    
  protected:

    //------------------------------------------------------------------------
    //  Default constructor -- use with caution
    //------------------------------------------------------------------------
    CMicCostCalculator( ) {}


  public:

    //------------------------------------------------------------------------
    //  ctor
    //------------------------------------------------------------------------
    CMicCostCalculator( const CConfigFile & oConfigFile )
      : CXDMReconstruction( oConfigFile ) {}

    //------------------------------------------------------------------------
    //  InitializeWithDataFiles
    //------------------------------------------------------------------------    
    virtual void InitializeWithDataFiles( )
    {
      CXDMReconstruction::InitializeWithDataFiles();
    }

    //------------------------------------------------------------------------
    //  EvaluateCost
    //------------------------------------------------------------------------    
    void EvaluateCost()
    {
      std::cerr << "Evaluating Cost of Mic File " << std::endl;

      const vector<CDetector> &  vDetectorList = oExpSetup.GetDetectorList();
      const CSimulationRange & oRangeToIndexMap  = oExpSetup.GetRangeToIndexMap();
      
      Float fInitialSideLength = GetMicSideLength();
      
      CXDMOverlapFunction oOverlapFn( this, oReconstructionLayer, vDetectorList,
                                      oRangeToIndexMap, oExpData );

      vector<SVoxel>::iterator pVoxelIter = oReconstructionLayer.VoxelListBegin();

      typedef std::pair< SVoxel, SOverlapInfo > VoxelCostPair;
      vector< VoxelCostPair > oVoxelCostList;
      
      for( ; pVoxelIter != oReconstructionLayer.VoxelListEnd(); pVoxelIter ++ )
      {
        SOverlapInfo oOverlapInfo = GetOverlapProperty( *pVoxelIter,
                                                        oReconstructionLayer,
                                                        vDetectorList,
                                                        oRangeToIndexMap,
                                                        oExpData );
        oVoxelCostList.push_back( std::make_pair( *pVoxelIter, oOverlapInfo ) );
      }

      stringstream ss;
      ss << oSetupFile.OutStructureBasename << ".mic.cost";

      ofstream os( ss.str().c_str() );

      // Writing a generalized Mic file -- need to be incorpoerated into the
      // MicIO class

      os << fInitialSideLength << std::endl;
      for( Size_Type i = 0; i < oVoxelCostList.size(); i ++ )
      {
        os << XDMUtility::CXDMVoxel::LeftVertex( oVoxelCostList[i].first ) << " "
           << XDMUtility::CXDMVoxel::Label( oVoxelCostList[i].first ) << " "
           <<  oVoxelCostList[i].first.nGeneration << " "
           <<  oVoxelCostList[i].first.nPhase << " " 
           <<  RadianToDegree( oVoxelCostList[i].first.oOrientMatrix.GetEulerAngles() ) << " "
          // -------- cost related items;

           << oVoxelCostList[i].second.fQuality << " "
           << oVoxelCostList[i].second.nPixelOverlap << " "
           << oVoxelCostList[i].second.nPixelOnDetector << " "
           << oVoxelCostList[i].second.nPeakOverlap << " "
           << oVoxelCostList[i].second.nPeakOnDetector << ""
           << std::endl;
      }
      
    }
    
    //------------------------------------------------------------------------
    //  Serializer for MPI (required by CDXDMReconstruction, but may not be used )
    //------------------------------------------------------------------------
    virtual Bool Save   ( CSerializer   & oSerialBuf ) const
    {
      return true;
    }

    //------------------------------------------------------------------------
    //  Restore
    //------------------------------------------------------------------------
    virtual Bool Restore( CDeserializer & oSerialBuf )
    {
      return true;
    }
    
  };
  
} //  namespace ReconstructionAnalysis

#endif
