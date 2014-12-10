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
//  GrainReconstruction.h
//   Author:   Frankie Li (li31@llnl.gov)
//
//   Purpose:  Simple serial grain reconstruction
//
////////////////////////////////////////////////////////////

#ifndef _GRAIN_COST_FUNCTION_H
#define _GRAIN_COST_FUNCTION_H

#include "Reconstructor.h"
#include "DiscreteAdaptive.h"
#include <ctime>
#include <boost/shared_ptr.hpp>

namespace Reconstruction
{
  //----------------------------------------
  //
  //  GrainReconstruction
  //
  //  This is the root of all serial reconstruction calls, at least
  //  for the simple version.
  //
  //----------------------------------------
  class GrainCostFunction
  {

  public:  
    // parameterization of GrainCostFunction later
    typedef SVoxel SamplePointT;
    typedef MicAnalysis::CMicGrid SamplePointGrid; 
    
    typedef MicFile<SamplePointT>                                 Mic;
    
    typedef ReconstructionStrategies::BreadthFirstStrategy<SamplePointT, SamplePointGrid> ReconstructionStrategy;
    
    typedef ReconstructionStrategy::SamplePointPtr SamplePointPtr;
    typedef boost::shared_ptr< ReconstructionSetup >   ReconstructionSetupPtr;
    typedef boost::shared_ptr< CSimulation >           SimulatorPtr;

  private:    
    ReconstructionSetupPtr  pSetup;
    SimulatorPtr            pSimulator;
    
    GrainCostFunction();

      
  public:

    GrainCostFunction( ReconstructionSetupPtr _pSetup,
		       SimulatorPtr          _pSimulator )
      : pSetup( _pSetup ), pSimulator( _pSimulator )
      {  }
    
    //----------------------------------------
    //  GetSetup
    //----------------------------------------
    ReconstructionSetupPtr GetSetup() 
    {
      return pSetup;
    }

    //----------------------------------------
    //  GetSimulator
    //----------------------------------------
    SimulatorPtr GetSimulator() 
    {
      return pSimulator;
    }
    
    //----------------------------------------
    //  GetAssociatedPixelList
    //  - Return the set of pixels associated with
    //    the grain specified by the list of voxels.
    //----------------------------------------

    //----------------------------------------
    //  GetAverageVoxelOverlapRatio
    //----------------------------------------
    Float GetAverageVoxelOverlapRatio( std::vector<SamplePointT> & VoxelList )
    {
      if( VoxelList.size() == 0 )
	return 0;

      DiscreteRefinement<SamplePointT>  AdpReconstructor( *pSimulator, *pSetup );
      
      Float fAveragedOverlapRatio = 0;
      for( int i = 0; i < VoxelList.size(); i ++ )
      {
	CostFunctions::SOverlapInfo oInfo = AdpReconstructor.EvaluateOverlapInfo( VoxelList[i] );
	fAveragedOverlapRatio += CostFunctions::Utilities::GetHitRatio  ( oInfo );
      }
      
      return fAveragedOverlapRatio / static_cast<Float>( VoxelList.size() );
    }

    
  };

  
}


#endif
