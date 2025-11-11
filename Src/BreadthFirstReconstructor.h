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
//  BreadthFirstReconstruction.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of general breadth first reconstruction
//
//  NOTE:     Currently this is very similar to IntensityDecomposition, but
//            this will change as intensity decomposition gets modified.
////////////////////////////////////////////////////////////


#ifndef BREADTH_FIRST_RECON_H
#define BREADTH_FIRST_RECON_H

#include "ParallelReconstructor.h"
#include "DiscreteAdaptive.h"
#include <memory>
#include <queue>
#include "XDMCommCore.h"
#include <limits>
#include "SimulationData.h"


namespace ParallelReconstructor
{
  //-----------------------------------
  //  IntensityFitCmd
  //
  //  Common functions for both intensity fits
  //-----------------------------------
  class BFSFitCmd
  {
  protected:
    enum Commands
      {
        BEGIN_BFS_FIT = 1,
        FIT           = 2,
        REPORT        = 3,
        RESTART_FIT   = 5,
        PROCESS_DONE  = -2
      };
  };

  //-----------------------------------
  // LazyBFSClient
  //
  // Purpose:  Threshold Decomposition
  //-----------------------------------
  template< class SamplePointT = SVoxel,
            class Reconstructor = Reconstruction::DiscreteRefinement<SamplePointT>,
            class SamplePointGrid = MicAnalysis::CMicGrid  >
  class LazyBFSClient
    : public ParallelReconstructorClient,
      public BFSFitCmd
  {
  public:
    typedef std::shared_ptr< Reconstructor >  ReconstructorPtr;
    typedef ReconstructionStrategies::BreadthFirstStrategy<SamplePointT, SamplePointGrid> ReconstructionStrategy;
    typedef PBRMath::BBox2D                                       BBox2D;
    typedef typename ReconstructionStrategy::SamplePointPtr       SamplePointPtr;
    typedef MicFile<SamplePointT>                                 Mic;

  private:
    ReconstructionSetup   & LocalSetup;
    
    CSimulation            Simulator;
    ReconstructionStrategy VoxelQueue;
    ReconstructorPtr       pReconstructor;
    
    LogStream & osLogFile;
    XDMParallel::XDMCommunicator Comm;
    Int nMyID;

    vector<SamplePointPtr> Refit( const SamplePointT & oCenter );


    //--------
    //  The second parameter is a hack to not write a new communication
    //  section for restart from existing files.
    //--------
    vector<SamplePointPtr> Fit( const SamplePointT & oCenter,
                                bool bSkipDiscreteSearch = false );
    

  public:
    
    LazyBFSClient( ReconstructionSetup & oReconstructionSetup,
                   int nID, LogStream & os )
      : LocalSetup( oReconstructionSetup ),
        osLogFile( os ), nMyID( nID ) {}
      
    void Initialize( );
    void Process();
    
  };

  //-----------------------------------
  //  LazyBFSServer
  //
  //-----------------------------------
  template< class SamplePointT = SVoxel, class SamplePointGrid = MicAnalysis::CMicGrid >
  class LazyBFSServer
    : public ParallelReconstructorServer,
      public BFSFitCmd
  {
  public:
    typedef ReconstructionStrategies::RestrictedStratifiedGrid< SamplePointT,SamplePointGrid > ReconstructionStrategy;
    typedef XDMParallel::SingleElementDistribution                                             SingleElementDistributor;
    typedef MicFile<SamplePointT>                                                              Mic;
    typedef std::shared_ptr<MicIOBase>                                                       MicPtr;
  protected:
    
    int nMyID;
    int nProcessingElements;
    const ReconstructionSetup & LocalSetup;  // Server doesn't do anything to this other than to send
    LogStream & osLogFile;

    XDMParallel::XDMCommunicator Comm;
    ReconstructionStrategy VoxelQueue;
    
    void WritePartialFitResult( time_t & oTime, const string & posfix );
    
  public:
    
    LazyBFSServer( const ReconstructionSetup & oReconstructionSetup,
                   int nID, int nPE,
                   LogStream & os )
      : nMyID( nID ), nProcessingElements( nPE ),
        LocalSetup( oReconstructionSetup ),
        osLogFile( os ) {}

    //--------------------------------------------------
    //  This is the interface layer
    //--------------------------------------------------
    void Initialize( );
    
    //------------------------------------------
    //  Initialize With Partial Fit
    //   TODO:  Change CMic -> iterator
    //    -- note:  Every point in the start point will
    //      still be refitted
    //------------------------------------------
    void Initialize( MicPtr StartPoint );
    
    //------------------------------------------
    //  To be changed
    //   TODO:  Change CMic -> iterator
    //------------------------------------------
    void InitializeRestart( MicPtr PartialRestartPtr );
    
    void Process();

    //------------------------------------------
    //  TODO:  return iterator pair
    //------------------------------------------
    Mic Solution() { return VoxelQueue.Solution(); }

  };  
  


}


#include "BreadthFirstReconstructor.tmpl.cpp"

#endif
