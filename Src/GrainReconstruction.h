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

#ifndef _Grain_RECONSTRUCTION_H
#define _Grain_RECONSTRUCTION_H

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
  class GrainReconstruction
  {

  public:  
    // parameterization of GrainReconstruction later
    typedef SVoxel SamplePointT;
    typedef MicAnalysis::CMicGrid SamplePointGrid; 
    
    typedef MicFile<SamplePointT>                                 Mic;
    
    typedef ReconstructionStrategies::BreadthFirstStrategy<SamplePointT, SamplePointGrid> ReconstructionStrategy;
    
    typedef ReconstructionStrategy::SamplePointPtr SamplePointPtr;
    typedef boost::shared_ptr< ReconstructionSetup >   ReconstructionSetupPtr;
  private:    
    ReconstructionSetupPtr  pSetup;
    ReconstructionStrategy  VoxelQueue;   // acts as a queue
    CSimulation oSimulator;
    
    GrainReconstruction();
  public:

    GrainReconstruction( const CConfigFile & oConfigFile )
    {
      pSetup = ReconstructionSetupPtr( new ReconstructionSetup() );
      RUNTIME_ASSERT( oConfigFile.nMicGridType == eTriangular, 
		      "Non-triangular grid not implemented for serial reconstruction\n" );
      pSetup->InitializeWithDataFiles( oConfigFile );
      boost::shared_ptr<Mic> pMic= boost::dynamic_pointer_cast<Mic>( pSetup->ReconstructionRegion());
      VoxelQueue.Initialize( *pMic,
                             pSetup->MinSideLength() );
      oSimulator.Initialize( pSetup->ExperimentalSetup() );
    }
    
    //----------------------------------------
    //  ReconstructGrain
    //     Given a voxel marking the central seed point, reconstruct outward 
    //     in a breadth first sense.
    //----------------------------------------
    vector<SamplePointPtr> ReconstructGrain( SamplePointT Center )
    {
      using namespace ReconstructionStrategies::MultiStagedDetails;
      DiscreteRefinement<SamplePointT>  AdpReconstructor( oSimulator, *pSetup );

      VoxelQueue.Reset();
      
      int nCenterCode;
      
      boost::tie( Center, nCenterCode ) = AdpReconstructor.ReconstructVoxel( Center  );    // Fit center
      CostFunctions::SOverlapInfo oInfo = AdpReconstructor.EvaluateOverlapInfo( Center );
      Center.fConfidence        = CostFunctions::Utilities::GetConfidence( oInfo );
      Center.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio  ( oInfo );
      
      Center.nID = REFIT;
      VoxelQueue.Push( Center );
      
      if( nCenterCode != SearchDetails::CONVERGED
	  && Center.fPixelOverlapRatio < pSetup->InputParameters().fMinAccelerationThreshold )
	return VoxelQueue.SolutionVector();
      
      VoxelQueue.InsertSeed( Center );
      Float fBestConf = Center.fPixelOverlapRatio;  // begin fitting outward, breadth first

      while( VoxelQueue.Size() > 0 )
      {
	SamplePointT CurVoxel = VoxelQueue.First();
	VoxelQueue.Pop();
	oInfo = AdpReconstructor.LocalOptimization( CurVoxel );
	fBestConf = std::max( CurVoxel.fPixelOverlapRatio, fBestConf );

	if( ( CurVoxel.fPixelOverlapRatio / fBestConf ) > 0.9  )  // need more sophisticated acceptance
	{
	  CurVoxel.nID = FITTED;
	  VoxelQueue.Push( CurVoxel );
	  VoxelQueue.InsertSeed( CurVoxel );
	}
	else
	{
	  CurVoxel.nID = REFIT; // reset to unfitted
	  VoxelQueue.Push( CurVoxel );
	}
      }    
      return VoxelQueue.SolutionVector();
    }
    
    //----------------------------------------
    //  ReconstructRandomGrain
    //   Choose a random point and apply reconstruct grain
    //----------------------------------------
    vector<SamplePointPtr> ReconstructRandomGrain( )
    {
      typedef ReconstructionStrategies::RestrictedStratifiedGrid< SamplePointT,SamplePointGrid > 
	ReconstructionStrategy;
      
      ReconstructionStrategy InitialPoints;
      
      boost::shared_ptr< Mic > pMic 
	= boost::dynamic_pointer_cast<Mic>( pSetup->ReconstructionRegion() );
      
      InitialPoints.Initialize( *pMic,
				pSetup->MinSideLength(),
				pSetup->InputParameters().SampleCenter,
				pSetup->InputParameters().SampleRadius );
      
      vector<SamplePointPtr> FittedGrain = ReconstructGrain( InitialPoints.First() );
      
      return FittedGrain;
    }
    
    //----------------------------------------
    //  GetDetectorPtr
    //----------------------------------------
    ReconstructionSetupPtr GetSetup() const
    {
      return pSetup;
    }

    //----------------------------------------
    //  GetAssociatedPixelList
    //----------------------------------------

  };
}


#endif
