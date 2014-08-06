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
//  Reconstructor.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//   Purpose:  Interface and implementations of simple reconstructors
//
////////////////////////////////////////////////////////////

#ifndef _RECONSTRUCTOR_H
#define _RECONSTRUCTOR_H

#include "DiscreteSearch.h"
#include "ConfigFile.h"
#include "ExperimentSetup.h"
#include "Voxel.h"
#include "ReconstructionSetup.h"
#include "Simulation.h"
#include "OverlapInfo.h"
#include "ReconstructionStrategies.h"
#include "ReconstructionFilters.h"
#include "ContinuousSearch.h"
#include "SearchTraits.h"
#include <ctime>
#include "boost/tuple/tuple.hpp"


namespace Reconstruction
{

  namespace Details
  {
    
    
  }
  
  //-----------------------------
  //  The actual driver that runs the reconstruction  (i.e., will have multiple reconstructor objects)
  //-----------------------------
  class BasicVoxelReconstructor
  {
    typedef CostFunctions::SOverlapInfo           SOverlapInfo;
    typedef ReconstructionSetup::EtaAngularFilter PeakFilterT;
    typedef OrientationSearch
    ::GeneralSearchCostFn< PeakFilterT, CSimulation,
                           CostFunctions::TrianglePixelOverlapCounter,
                           CostFunctions::DetectorOverlapCounter,
                           CostFunctions::XDMSampleVertexOverlapCounter,
                           CostFunctions::VoxelToVertices,
                           CostFunctions::VoxelCostFunction, SVoxel >  LocalSearchCostFunctions;
    
    typedef OrientationSearch
    ::GeneralSearchCostFn< PeakFilterT, CSimulation,
                           CostFunctions::PixelBasedPeakOverlapCounter,
                           CostFunctions::DetectorOverlapCounter,
                           CostFunctions::XDMSampleVertexOverlapCounter,
                           CostFunctions::VoxelToCenter,
                           CostFunctions::InterpolatedVoxelCostFunction, SVoxel > GlobalSearchCostFunctions;
  public:
    typedef SearchDetails::SCandidate       SCandidate;
    typedef SearchDetails::SSearchParameter SSearchParameter;
    
  private:  
    const CSimulation    & oSimulator;  // holding my own copy
    ReconstructionSetup  & oSetup;      // Definitely no copy - this holds the entire data
                                        // But the only reason why this isn't const is because
                                        // of sample geometry
    const CConfigFile         & InputParam;    
    
    vector< vector<SMatrix3x3> > AllLevelLocalGrids;

    //---------------------------------------------------------
    //  RunDiscreteSearch
    //---------------------------------------------------------
    vector<SCandidate> RunDiscreteSearch( const SVoxel & oVoxel,
                                          const vector<SMatrix3x3> & oLocalGrid );

    //----------------------------------------
    //  Disable default constructor
    //----------------------------------------
    BasicVoxelReconstructor(); 

  public:

    BasicVoxelReconstructor( const CSimulation & oSim,  ReconstructionSetup & ReconSetup )
      : oSimulator( oSim ),
        oSetup( ReconSetup ),
        InputParam( ReconSetup.InputParameters() )
    {

      for ( int i = 0; i <= InputParam.nMaxLocalResolution; ++ i )
      {
        vector<SMatrix3x3> oLocalGrid;
        OrientationSearch
          ::Utilities
          ::GenerateLocalGrid( oLocalGrid, InputParam.fLocalOrientationGridRadius, i );
        AllLevelLocalGrids.push_back( oLocalGrid );
      }
    }
    
    //----------------------------------------
    //  ReconstructVoxel
    //----------------------------------------
    std::pair<SVoxel, Int> ReconstructVoxel( const SVoxel & oInputVoxel );
    
    //----------------------------------------
    //  LocalOptimization
    //----------------------------------------
    SOverlapInfo LocalOptimization( SVoxel & oVoxel );

    
    //----------------------------------------
    //  EvaluateOverlapInfo
    //----------------------------------------
    SOverlapInfo EvaluateOverlapInfo( const SVoxel & oVoxel );

    //----------------------------------------
    //  DEBUG_OutputVoxelCost
    //----------------------------------------
    void DEBUG_OutputVoxelCost( const SVoxel & TestVoxel );
  };

}

#endif
