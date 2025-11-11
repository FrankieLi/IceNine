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
//  DiscreteAdapative.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//   Purpose:  Interface and implementations of simple reconstructors
//
////////////////////////////////////////////////////////////

#ifndef _DISCRETE_ADAPTIVE_H
#define _DISCRETE_ADAPTIVE_H


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
#include <tuple>
namespace Reconstruction
{
    
  //---------------------------------------------
  //  Adapator to allow continuous search cost function
  //  to be used in discrete search
  //---------------------------------------------
  template< class OverlapCostFn >
  class ContinuousToDiscreteAdapater
  {

  public:
    typedef typename OverlapCostFn::SamplePointT SamplePointT;
  private:
    OverlapCostFn  * CostFn;
    CSample        & oSample;
    const XDMSimulation::CSimulationRange & oRangeToIndexMap;
    const DetectorListT   & oDetList;
    const CSimulationData & oExpData;
    
    typedef CostFunctions::SOverlapInfo SOverlapInfo;

  public:
    
    ContinuousToDiscreteAdapater( OverlapCostFn *CostFn_, CSample & s_,
                                  const XDMSimulation::CSimulationRange & oRangeMap_,
                                  const DetectorListT & DetList_,
                                  const CSimulationData & oData_ )
      : CostFn( CostFn_ ), oSample( s_ ),
        oRangeToIndexMap( oRangeMap_ ),
        oDetList( DetList_ ),
        oExpData( oData_ ) {}
    
    //------------------------------------------
    //  Actual wrapper
    //------------------------------------------
    SOverlapInfo operator()( const SamplePointT & v )
    { 
      SOverlapInfo oOverlapInfo
        = this->CostFn->operator()( v, oSample, oDetList, oRangeToIndexMap, oExpData );

      return oOverlapInfo;
    }
  };

  
  //-----------------------------
  //  The actual driver that runs the reconstruction  (i.e., will have multiple reconstructor objects)
  //
  //  TODO:  Make this take two parameters as templates (LocalSearch vs. Global Search)
  // 
  //-----------------------------
  template<class SamplePointT>
  class DiscreteRefinement
  {
    typedef CostFunctions::SOverlapInfo           SOverlapInfo;
    typedef ReconstructionSetup::EtaAngularFilter PeakFilterT;
    typedef typename CostFunctions::ShapePixelOverlapCounter<SamplePointT> ShapePixelCounterT;
    typedef typename OrientationSearch
    ::GeneralSearchCostFn< PeakFilterT, CSimulation,
                           ShapePixelCounterT,     // <--- need to change to template-template parameter
                           CostFunctions::DetectorOverlapCounter,
                           CostFunctions::XDMSampleVertexOverlapCounter,
                           CostFunctions::VoxelToVertices,
                           CostFunctions::VoxelCostFunction,
                           SamplePointT >  LocalSearchCostFunctions;
    
    typedef typename OrientationSearch
    ::GeneralSearchCostFn< PeakFilterT, CSimulation,
                           CostFunctions::PixelBasedPeakOverlapCounter,
                           CostFunctions::DetectorOverlapCounter,
                           CostFunctions::XDMSampleVertexOverlapCounter,
                           CostFunctions::VoxelToCenter,
                           CostFunctions::InterpolatedVoxelCostFunction,
                           SamplePointT> GlobalSearchCostFunctions;
    typedef typename LocalSearchCostFunctions::VoxelOverlapFn LocalVoxelCostFn;
    typedef ContinuousToDiscreteAdapater< LocalVoxelCostFn > AdaptedLocalVoxelCostFn;    // A shim so that discrete search
                                                                                         // can use continuous search's cost function
    
  public:
    typedef SearchDetails::SCandidate       SCandidate;
    typedef SearchDetails::SSearchParameter SSearchParameter;

  protected:  
    const CSimulation    & oSimulator;  // holding my own copy
    ReconstructionSetup  & oSetup;      // Definitely no copy - this holds the entire data
                                        // But the only reason why this isn't const is because
                                        // of sample geometry
    const CConfigFile         & InputParam;    
    
    vector<SCandidate> RunDiscreteSearch( const SamplePointT & oVoxel,
                                          vector<SMatrix3x3>::const_iterator pFirst,
                                          vector<SMatrix3x3>::const_iterator pEnd,
                                          const vector<SMatrix3x3> & oLocalGrid,
                                          int nQMax, LocalVoxelCostFn & CompleteVoxelCostFn,
                                          Float fAngularRadius );
    //----------------------------------------
    //  Disable default constructor
    //----------------------------------------
    DiscreteRefinement(); 
    
  public:
    
    DiscreteRefinement( const CSimulation & oSim,  ReconstructionSetup & ReconSetup )
      : oSimulator( oSim ),
        oSetup( ReconSetup ),
        InputParam( ReconSetup.InputParameters() )
    {  }

    std::pair<SamplePointT, Int> ReconstructVoxel( const SamplePointT & oInputVoxel );

    //-------------------------
    //  To be removed
    //-------------------------
    SOverlapInfo EvaluateOverlapInfo( const SamplePointT & oVoxel );
    SOverlapInfo LocalOptimization( SamplePointT & oInputVoxel );
  };


}


#include "DiscreteAdaptive.tmpl.cpp"

#endif
