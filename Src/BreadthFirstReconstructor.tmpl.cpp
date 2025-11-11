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
//  BreadthFirstReconstruction.cpp
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of general breadth first reconstruction
//
//  NOTE:     Currently this is very similar to IntensityDecomposition, but
//            this will change as intensity decomposition gets modified.
////////////////////////////////////////////////////////////



namespace ParallelReconstructor
{


  //--------------------------------------------------------------------------------------
  //                  C L I E N T
  //--------------------------------------------------------------------------------------


  
  //-----------------------------------
  //  Initialize
  //-----------------------------------
  template< class SamplePointT, class R, class G >
  void LazyBFSClient<SamplePointT, R, G>::Initialize( )
  {
    typedef typename LazyBFSClient<SamplePointT, R, G>::Mic Mic;
    std::shared_ptr< Mic > pMic= std::dynamic_pointer_cast<Mic>( LocalSetup.ReconstructionRegion() );
    VoxelQueue.Initialize( *pMic,
                           LocalSetup.MinSideLength() );

    std::cout << "Client:  Finished with initialization" << std::endl;
    // broadcast
    int nCommand;
    Comm.BcastRecvCommand( 0, & nCommand );
    RUNTIME_ASSERT( nCommand ==  BEGIN_BFS_FIT,
                    "Unexpected command recv'd in LazyBFSClient, Initialize\n");
    Simulator.Initialize( LocalSetup.ExperimentalSetup() );
    pReconstructor = ReconstructorPtr ( new R( Simulator, LocalSetup ) );
  }

  
  //-----------------------------------
  // Refit
  //-----------------------------------
  template< class SamplePointT, class R, class G >
  vector< typename  LazyBFSClient< SamplePointT, R, G >::SamplePointPtr >
  LazyBFSClient<SamplePointT, R, G>::Refit( const SamplePointT & oCenter )
  {
    using namespace ReconstructionStrategies::MultiStagedDetails;
    SamplePointT v = oCenter;
    
    //--------------------------
    //  Might want to give a command to choose LocalOptimization vs. EvaluateOverlapInfo
    //--------------------------
    CostFunctions::SOverlapInfo oInfo = pReconstructor->LocalOptimization( v );

    if( CostFunctions::Utilities::GetConfidence( oInfo )
        >= LocalSetup.InputParameters().fPartialResultAcceptanceConf )
    {
      return Fit( v, true );  // skip discrete search
    }
    else
    {
      return Fit( v, false ); // do not skip discrete search
    }
  }
  
  //-----------------------------------
  //  Fit
  //-----------------------------------
  template< class SamplePointT, class R, class G >
  vector< typename  LazyBFSClient< SamplePointT, R, G >::SamplePointPtr >
  LazyBFSClient<SamplePointT, R, G>::Fit( const SamplePointT & oCenter, bool bSkipDiscreteSearch )
  {
    using namespace ReconstructionStrategies::MultiStagedDetails;
    VoxelQueue.Reset();

    int nCenterCode;
    SamplePointT oResult;
    
    if( bSkipDiscreteSearch )
    {
      oResult = oCenter;
      nCenterCode = SearchDetails::PARTIAL;   // restart case
    }
    else
    {
      time_t oStartTime, oStopTime;
      time( &oStartTime );
      std::tie( oResult, nCenterCode ) = pReconstructor->ReconstructVoxel( oCenter );
      time( &oStopTime );
      double oTimeDiff = difftime( oStopTime, oStartTime );
      GET_LOG( osLogFile ) << "  First Voxel Required: " << oTimeDiff << " Sec " << std::endl;
      GET_LOG( osLogFile ) << "    Orientation: " << RadianToDegree( oResult.oOrientMatrix.GetEulerAngles() ) << std::endl;
    }
    
    CostFunctions::SOverlapInfo oInfo = pReconstructor->EvaluateOverlapInfo( oResult );  
    
    oResult.fConfidence        = CostFunctions::Utilities::GetConfidence( oInfo );
    oResult.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio  ( oInfo );

    GET_LOG( osLogFile ) << " [Conf, HitRatio ] " << oResult.fConfidence << " " << oResult.fPixelOverlapRatio << std::endl;
    if( nCenterCode != SearchDetails::CONVERGED
      	&& oResult.fPixelOverlapRatio < LocalSetup.InputParameters().fMinAccelerationThreshold )
    {
      //      std::cout << oResult.fPixelOverlapRatio << " " << LocalSetup.InputParameters().fMinAccelerationThreshold << std::endl;
      oResult.nID                = REFIT;
      VoxelQueue.Push( oResult );
      return VoxelQueue.SolutionVector();
    }
    
    GET_LOG( osLogFile ) << "Fitting neighbors via Breadth First Search" << std::endl;
    
    oResult.nID = FITTED;
    VoxelQueue.Push( oResult );   // only unique elements will be saved in the solution
    VoxelQueue.InsertSeed( oResult );
    Float fBestConf = oResult.fPixelOverlapRatio;
    while( VoxelQueue.Size() > 0 )
    {
      SamplePointT CurVoxel = VoxelQueue.First();
      VoxelQueue.Pop();
      oInfo = pReconstructor->LocalOptimization( CurVoxel );
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
      GET_LOG( osLogFile ) << "QueueSize " << VoxelQueue.Size() << " Fitted "
                           << VoxelQueue.NumFitted() << " continuous voxels " << std::endl 
                           << "  [ Cost ] " << CurVoxel.fCost 
                           << "  " <<oInfo.nPixelOverlap << "/" << oInfo.nPixelOnDetector << " "
                           << oInfo.nDetectorsOverlap << " "
                           << oInfo.nPeakOverlap << "/" << oInfo.nPeakOnDetector << " " 
                           << oInfo.fQuality << " |=| " 
                           << fBestConf << " " << CurVoxel.fPixelOverlapRatio 
                           << "    Original : " << RadianToDegree( oResult.oOrientMatrix.GetEulerAngles() )
                           << "    Optimized: " << RadianToDegree( CurVoxel.oOrientMatrix.GetEulerAngles() )
                           << " VoxelLocation " << CurVoxel.GetCenter() << std::endl;
    }
    
    return VoxelQueue.SolutionVector();
  }
  
  //-----------------------------------
  //  Process
  //-----------------------------------
  template< class SamplePointT, class R, class SamplePointGrid >
  void LazyBFSClient< SamplePointT, R, SamplePointGrid>::Process()
  {
    GET_LOG( osLogFile ) << "[Client] Recving Work Unit:" << std::endl;
    SamplePointT oCenter;
    const int nServerPE = 0;
    
    int nCommand = PROCESS_DONE;
    do
    {
      time_t oStartTime, oStopTime;
      int nSourcePE;
      time( &oStartTime );
      Comm.RecvCommand( &nSourcePE, &nCommand );

      if( nCommand != PROCESS_DONE )
      {
        time( &oStopTime );
        double oTimeDiff = difftime( oStopTime, oStartTime );
        GET_LOG( osLogFile ) << "Time spent waiting for command (sec): " << oTimeDiff << std::endl;
        
        time( &oStartTime );
        Comm.RecvWorkUnit( nServerPE, &oCenter );
        RUNTIME_ASSERT( nServerPE == 0,
                        "Driver::FitElements: Recving messages from units other than ROOT!  STOP \n" );
        
        vector<SamplePointT> oFittedVoxels;
        vector<SamplePointPtr> oTmp;
        switch( nCommand )
        {
          case FIT:
            oTmp = Fit( oCenter );
            break;
          case RESTART_FIT:
            oTmp = Refit( oCenter );
            break;
          case PROCESS_DONE:
            break;
          default:
            RUNTIME_ASSERT( 0, "[LazyBFSClient] UNEXPECTED COMMAND\n" );
        } 
        for( Size_Type i = 0; i < oTmp.size(); i ++ )
        {
          RUNTIME_ASSERT(  SamplePointGrid::IsValid( oTmp[i] ), "Solution is not valid\n" );
          oFittedVoxels.push_back( * oTmp[i] );
        }
        int nFitted = oTmp.size();
        time( &oStopTime );
        oTimeDiff = difftime( oStopTime, oStartTime );
        GET_LOG( osLogFile ) << "  Time elapsed for this set of voxels (sec): " << oTimeDiff << std::endl;
        GET_LOG( osLogFile ) << "  Time per voxel (sec): " << oTimeDiff / nFitted << std::endl;
        
        Comm.SendCommand( nMyID, nServerPE, REPORT );
        Comm.SendWorkUnitList( nServerPE, oFittedVoxels );
      }
      else
      {
        std::cout << "Procss done command recieved " << std::endl;
      }
    } while( nCommand != PROCESS_DONE );
    
    RUNTIME_ASSERT( nCommand == PROCESS_DONE,
                    "Error!  [LazyBFSClient]: Unexpected command recieved\n");
  }
}



//--------------------------------------------------------------------------------------
//                  S E R V E R
//--------------------------------------------------------------------------------------


namespace ParallelReconstructor
{
  //-----------------------------------
  //  WritePartialFitResult
  //-----------------------------------
  template< class S, class G>
  void LazyBFSServer<S, G>
  ::WritePartialFitResult( time_t & oRefTime,
                           const string & postfix )
  {
    time_t oCurrentTime;
    time( & oCurrentTime );
    if( difftime( oCurrentTime, oRefTime )
        > LocalSetup.InputParameters().nSecondsBetweenSave )
    {
      time( &oRefTime );
      
      // no reduction to fundamental zone at this moment
      stringstream ss;
      ss << LocalSetup.InputParameters().OutStructureBasename <<  ".mic" << postfix;
      VoxelQueue.Solution().Write( ss.str() );
      GET_LOG( osLogFile ) << "Wrote to " << ss.str() << std::endl;
      
    } 
  }

  //------------------------------------------
  //  InitializeRestart
  //------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void LazyBFSServer<SamplePointT, SamplePointGrid>
  ::InitializeRestart( typename LazyBFSServer<SamplePointT, SamplePointGrid>::MicPtr PartialRestartPtr_ )
  {
    RUNTIME_ASSERT( 0, "InitializeRestart from LazyBFS is deprecated" );
  }
  
  //-----------------------------------
  //  Initialize
  //-----------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void LazyBFSServer<SamplePointT, SamplePointGrid>::Initialize( )
  {
    typedef typename LazyBFSServer<SamplePointT, SamplePointGrid>::Mic Mic;
    std::shared_ptr< Mic > pMic     = std::dynamic_pointer_cast<Mic>( LocalSetup.ReconstructionRegion() );
    std::shared_ptr< Mic > pPartial = std::dynamic_pointer_cast<Mic>( LocalSetup.PartialResult() );
    if ( LocalSetup.InputParameters().bUsePartialResult )
    {
      GET_LOG( osLogFile ) << "Restarting with partial result " << std::endl;
      VoxelQueue.InitializeRestart( *pMic,
                                    *pPartial,
                                    LocalSetup.InputParameters().fPartialResultAcceptanceConf,
                                    LocalSetup.MinSideLength(),
                                    LocalSetup.InputParameters().SampleCenter,
                                    LocalSetup.InputParameters().SampleRadius );
      std::cout << " Number of voxels in queue: " << VoxelQueue.Size()
                << " " << LocalSetup.InputParameters().SampleRadius << std::endl;
    }
    else
    {
      GET_LOG( osLogFile ) << "Start fitting from empty mic file " << std::endl;
      VoxelQueue.Initialize( *pMic,
                             LocalSetup.MinSideLength(),
                             LocalSetup.InputParameters().SampleCenter,
                             LocalSetup.InputParameters().SampleRadius );
      std::cout << " Number of voxels in queue: " << VoxelQueue.Size() << std::endl;
    }
    // broadcast
    Comm.BcastSendCommand( nMyID, BEGIN_BFS_FIT );
  }
  
  //-----------------------------------
  //  Initialize
  //   Note:  All points will be refitted.  This is
  //          designed for intensity fitting
  //-----------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void LazyBFSServer<SamplePointT, SamplePointGrid>
  ::Initialize( typename LazyBFSServer<SamplePointT, SamplePointGrid>::MicPtr StartPoint )
  {
    typedef typename LazyBFSServer<SamplePointT, SamplePointGrid>::Mic Mic;
    std::shared_ptr< Mic > pStartPoint = std::dynamic_pointer_cast<Mic>( StartPoint );
    VoxelQueue.Initialize( *pStartPoint,
                           pStartPoint->GetMinSideLength(),
                           LocalSetup.InputParameters().SampleCenter,
                           LocalSetup.InputParameters().SampleRadius );
    // broadcast
    Comm.BcastSendCommand( nMyID, BEGIN_BFS_FIT );
  }
  
  //-----------------------------------
  //  Process
  //-----------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void LazyBFSServer<SamplePointT, SamplePointGrid>::Process()
  {
    int FitCommand = FIT; 
    if( LocalSetup.InputParameters().bUsePartialResult )
    {
      GET_LOG( osLogFile ) << "RESTART_FIT in progress: " << std::endl;
      FitCommand = RESTART_FIT;
    }
    std::cout << "VoxelQueue Size " << VoxelQueue.Size() << " With NumProcessingElements: " << nProcessingElements << std::endl;
    if( VoxelQueue.Empty() )
    {
      Comm.SendCommand( nMyID, 1, nProcessingElements - 1, PROCESS_DONE  );
      GET_LOG( osLogFile ) <<  "WARNING!  Unexpected exit!  Reconstruction Region is empty"  << std::endl;
    }
    
    const Int nTotalClients = nProcessingElements - 1;
    std::queue<int> WaitQueue;
    for( int i = 1; i <= nTotalClients; i ++ )  // start from 1, since 0 is Server
      WaitQueue.push( i );

    std::cout << " WaitQueue Size " << WaitQueue.size() << std::endl;
    SingleElementDistributor VoxelDistributor;
    std::cout << "Before distribution " << std::endl;
    Utilities::WorkUnitDistribution( VoxelQueue, VoxelDistributor,
                                     WaitQueue, Comm, osLogFile,
                                     nMyID, FitCommand );
    
    GET_LOG( osLogFile ) << "nProcessingElements " << nProcessingElements << std::endl;
    
    DEBUG_ASSERT( static_cast<int>( WaitQueue.size() ) < nTotalClients,
                  "Haven't sent one voxel out.  Check number of voxels!" );
    GET_LOG( osLogFile ) << "\n\n Done distributing.   Num clients working: "
                         << nTotalClients - WaitQueue.size()  << std::endl;
    
    time_t oStartTime;
    time( &oStartTime );

    int nFileExt = 0;
    int nCollected = 0;

    
    while( static_cast<int>( WaitQueue.size() ) < nTotalClients
           || ! VoxelQueue.Empty() )         //  While someone's still working or voxel queue is not empty
    {
      Int nClientPE;    
      vector<SamplePointT> vResultVoxelList;
      vector<Int>    vResultCodeList;
      Int nCommand;
            
      Comm.RecvCommand( &nClientPE, &nCommand );
      WaitQueue.push( nClientPE );
      RUNTIME_ASSERT( nCommand == REPORT, "\nDriver::Server: Unknown command from client recv \n" );
      Comm.RecvWorkUnitList( nClientPE, vResultVoxelList );
      for( Size_Type n = 0; n < vResultVoxelList.size(); n ++ )
      {
        VoxelQueue.Push( vResultVoxelList[n] );
        nCollected ++;
      }
      if( vResultVoxelList.size() > 20 )   // probably shouldn't be hard coded
        VoxelQueue.RefreshSampleCenter();
      Utilities::WorkUnitDistribution( VoxelQueue, VoxelDistributor, WaitQueue,
                                       Comm, osLogFile, nMyID,  FitCommand );
      stringstream ss;
      ss << ".part" << nFileExt;
      WritePartialFitResult( oStartTime, ss.str() );
      nFileExt ++;  
    }
    
    GET_LOG( osLogFile ) << "Exit Wait Queue Size " << WaitQueue.size() << std::endl;
    GET_LOG( osLogFile ) << "VoxelQueue is  " << VoxelQueue.Empty() << std::endl;

    // empty WaitQueue
    while( WaitQueue.size() > 0 )
    {
      Comm.SendCommand( nMyID, WaitQueue.front(), PROCESS_DONE );
      WaitQueue.pop();      
    }
  }
  
} // end namespace
