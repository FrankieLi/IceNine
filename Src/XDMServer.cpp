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
//  XDMClient.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of the client side of the paralle IceNine
// 
//
//
////////////////////////////////////////////////////////////

#include "XDMServer.h"

namespace XDMParallel
{
  //---------------------------------------------------------------------------
  //  BroadcastExpData
  //
  //  Purpose:  Distribute initial experimental data along
  //            with parameters
  //---------------------------------------------------------------------------
  void CParallelServer::BroadcastExpData()
  {
    CSerializer oSimulationCheckpoint;
    Bool bSuccess = Save( oSimulationCheckpoint );
    RUNTIME_ASSERT( bSuccess, "\nCatastrophic Error!  Serialization failed!  Unable to send via MPI\n" );
    GET_LOG( osLogFile ) << "Server Started, Prepared to send:" << std::endl;
    char * pSendBuf = const_cast< char * >( oSimulationCheckpoint.GetBuffer() );
    Size_Type nBufSize = oSimulationCheckpoint.GetSize();
    GET_LOG( osLogFile ) << "nBufSize " << nBufSize  << std::endl;
   

    MPI_Bcast( &nBufSize, sizeof( nBufSize), MPI_CHAR, 0, MPI_COMM_WORLD); 
    MPI_Bcast( pSendBuf, nBufSize, MPI_CHAR, 0, MPI_COMM_WORLD);
    GET_LOG( osLogFile ) << "Initial data Sent Sent: [bytes sent] "
                         << nBufSize <<  std::endl;

    //---------------------------------------
    oSetup.Data().BCast_Send();
    //---------------------------------------
    GET_LOG( osLogFile ) << " Detector Data sent " <<  std::endl;
  }

  //---------------------------------------------------------------------------
  //  WriteFitResult
  //
  //  Parameter:  Postfix is used to for intermediate files
  //---------------------------------------------------------------------------
  void CParallelServer::WriteFitResult( boost::shared_ptr<MicIOBase> pMic, const string & oPostFix,
                                        bool bReduceToFZ )

  {
    CSymmetry * pSym = NULL;
    
    switch( oLocalConfigFile.SampleSymmetry )
    {
      case LatticeSymmetry::eCubic:
        pSym =  & ( LatticeSymmetry::CCubicSymmetry::Get() );
        break;
      case LatticeSymmetry::eHexagonal:
        pSym =  & ( LatticeSymmetry::CHexagonalSymmetry::Get() );
        break;
      case LatticeSymmetry::eTetragonal:
        pSym =  & ( LatticeSymmetry::CTetragonalSymmetry::Get() );
        break;
      default:
        cerr << "WARNING!!!!  Not using any symmetry!!" << std::endl;
        break;
    }
    
    if( pSym != NULL && bReduceToFZ)    // Reduce to fundamental zone
    {
      typedef MicFile< SVoxel >      TriangleMic;
      typedef MicFile< SquareVoxel > SquareMic;
      switch( oSetup.InputParameters().nMicGridType  )
      {
        case eTriangular:
          {
            XDMUtility
              ::ReduceToFundamentalZone( *( boost::dynamic_pointer_cast<TriangleMic>( pMic ) ), *pSym );
            break;
          }
        case eSquare:
          {
            XDMUtility
            ::ReduceToFundamentalZone( *( boost::dynamic_pointer_cast<SquareMic>  ( pMic ) ), *pSym );
          break;
          }
        default:
          {
            RUNTIME_ASSERT(0, "Reduction to FZ:  Unexpected grid option used for LazyBFS\n");
          } 
      }
    }
    
    stringstream ss;
    ss << oLocalConfigFile.OutStructureBasename <<  ".mic" << oPostFix;
    pMic->Write( ss.str() );
    GET_LOG( osLogFile ) << "Wrote to " << ss.str() << std::endl;
  }
  
  //---------------------------------------------------------------------------
  //  WriteDetectorFile
  //
  //  Parameter:  Postfix is used to for intermediate files
  //---------------------------------------------------------------------------
  void CParallelServer::WriteDetectorFile( const string & oPostFix )
  {
    stringstream ss;
    ss << oLocalConfigFile.OutStructureBasename <<  ".opt_det" << oPostFix;
    ofstream os;
    os.open( ss.str().c_str() );
    oSetup.ExperimentalSetup().WriteDetectorFile( os );
    os.close();
    GET_LOG( osLogFile ) << "Wrote to " << ss.str() << std::endl;
  }

  //------------------------------------------------------------------------
  //  RunIntensityDecomposition
  //------------------------------------------------------------------------
  void CParallelServer::RunIntensityDecomposition()
  {
    std::cout << "DEBUG RUN " << std::endl;
    RUNTIME_ASSERT( 0, "RunIntensityDecomposition disabled\n");
//     Float fMinThreshold = 0.1;
//     Float fMaxThreshold = 0.4;
//     for( int i = 0; i  < 1; i ++ )
//     {
//       std::cout << "Run Intensity Decomposition" << std::endl;
//       Comm.BcastSendCommand( nMyID, FIT_INT_DECOMP );
      
//       typedef ParallelReconstructor::IntensityDecompositionServer IDS;
//       IDS IDSServer( oSetup, nMyID, nProcessingElements, osLogFile );
//       IDSServer.Initialize( fMinThreshold, fMaxThreshold );  
//       IDSServer.Process();
      
//       CMic Result = IDSServer.Solution();
//       GET_LOG( osLogFile ) << "Wrote Result from Intensity Decomposition to " << std::endl;
//       stringstream ss;
//       ss << ".IDS" << i;
//       WriteFitResult( Result, ss.str(), true );
//       fMinThreshold = fMinThreshold + 0.2;
//       fMaxThreshold = fMaxThreshold + 0.2;
//     }
  }

  //------------------------------------------------------------------------
  //  RunLazyBFS
  //------------------------------------------------------------------------
  void CParallelServer::RunLazyBFS()
  {
    typedef CSample::MicIOPtrT     MicIOPtrT;
    std::cout << "Run LazyBFS" << std::endl;
    Comm.BcastSendCommand( nMyID, FIT_LBFS );
    switch( oSetup.InputParameters().nMicGridType  )
    {
      case eTriangular:
        {
          typedef ParallelReconstructor::LazyBFSServer<> LBFS;
          LBFS Server( oSetup, nMyID, nProcessingElements, osLogFile );
          Server.Initialize( );
          std::cout << "begin processing, triangle " << std::endl;
          Server.Process();
          MicIOPtrT pResult = oSetup.PartialResult();
          *( boost::dynamic_pointer_cast<LBFS::Mic>( pResult )) = Server.Solution();
        }
        break;
      case eSquare:
        {
          typedef ParallelReconstructor::LazyBFSServer<SquareVoxel, MicAnalysis::RectMicGrid<SquareVoxel> > LBFS;
          LBFS Server( oSetup, nMyID, nProcessingElements, osLogFile );
          Server.Initialize( );
          std::cout << "begin processing, square " << std::endl;
          Server.Process();
          MicIOPtrT pResult = oSetup.PartialResult();
          *( boost::dynamic_pointer_cast<LBFS::Mic>(pResult) ) = Server.Solution();
          break;
        }
      default:
        {
          RUNTIME_ASSERT(0, "Unexpected grid option used for LazyBFS\n");
        }
    }
    GET_LOG( osLogFile ) << "Wrote Result from LazyBFS to " << std::endl;
    stringstream ss;
    ss << ".LBFS";
    WriteFitResult( oSetup.PartialResult(), ss.str(), true );
  }


    //------------------------------------------------------------------------
  //  RunLazyBFS
  //------------------------------------------------------------------------
  void CParallelServer::RunParameterOptimization()
  {
    GET_LOG( osLogFile ) << " Begin Parameter Monte Carlo " << std::endl;
    
    RUNTIME_ASSERT( oSetup.InputParameters().nMicGridType == eTriangular,
                    "Error:  Only triangular voxels accepted for parameter MC right now\n" );
    typedef ParallelReconstructor:: ParameterOptimization::ParameterOptimizationServer<> PMC;
    PMC Server( oSetup, nMyID, nProcessingElements, osLogFile );
    Server.Initialize( );
    
    Server.Process();
    

    WriteFitResult( Server.Solution(), ".opt");
  }
  
  //------------------------------------------------------------------------
  //  RunLocalOrientationOptimization
  //------------------------------------------------------------------------
  void CParallelServer::RunLocalOrientationOptimization()
  {
    std::cout << "Run Local Optimization" << std::endl;
    Comm.BcastSendCommand( nMyID, FIT_LOCAL_OPTIMIZE );
    
    typedef ParallelReconstructor::LazyBFSServer<> LBFS;    // select server type - perhaps more than one voxel at a time
    LBFS Server( oSetup, nMyID, nProcessingElements, osLogFile );
    if ( oSetup.InputParameters().bUsePartialResult )
    {
//       typedef LBFS::Mic TriangleMic;
//       boost::shared_ptr< TriangleMic > pPartialResult= boost::dynamic_pointer_cast<TriangleMic>( oSetup.PartialResult() );
//       Server.InitializeRestart( *pPartialResult );
      Server.InitializeRestart( oSetup.PartialResult() );
      std::cout << " Restarting fit from partial result " << std::endl;
    }
    else
    {
      Server.Initialize( );
    }
    Server.Process();
    *( oSetup.PartialResult() ) = Server.Solution();
    GET_LOG( osLogFile ) << "Wrote Result from Local Optimization to " << std::endl;
    stringstream ss;
    ss << ".LBFS";
    WriteFitResult( oSetup.PartialResult(), ss.str(), true );
  }
  
  //------------------------------------------------------------------------
  //  RunIntensityDecomposition
  //------------------------------------------------------------------------
  void CParallelServer::RunStrainLazyBFS()
  {
    std::cout << "Run Strain Reconstruction" << std::endl;
    RUNTIME_ASSERT( 0, "Strain LazyBFS disabled " );
    // DEBUG
    std::cout << " -- DEBUG -- StrainLazyBFS: Replacing Reconstruction Region with Partial Result! " << std::endl;

    //oSetup.ReconstructionRegion() = PartialResult;
    
//     Comm.BcastSendCommand( nMyID, FIT_STRAIN_LBFS );
    
//     typedef ParallelReconstructor::LazyBFSServer<> LBFS;
//     LBFS Server( oSetup, nMyID, nProcessingElements, osLogFile );
//     Server.Initialize( PartialResult );  
//     Server.Process();
//     PartialResult = Server.Solution();
//     GET_LOG( osLogFile ) << "Wrote Result from Strain based LazyBFS to " << std::endl;
//     stringstream ss;
//     ss << ".strain";
//     WriteFitResult( PartialResult, ss.str(), true );
  }
  //------------------------------------------------------------------------
  //
  //  ProcessCommands  ( ParallelVersion )
  //  - Process all commands ( reconstruct, save file, reload, parameter
  //                           optimization, and so on... )
  //
  //------------------------------------------------------------------------
  void CParallelServer::ProcessCommands(  )
  {

//     int nDebug = 1;
//     while ( nDebug == 1 );
      
    
    GET_LOG( osLogFile ) << "-- Start Server --" << std::endl;
    GET_LOG( osLogFile ) << "Number of commands in the pipeline: " << oCommandList.size() << std::endl;

    for( Size_Type i = 0; i < oCommandList.size(); i ++ )
    {
      switch( oCommandList[i] )
      {
        case FIT:
          {
            RUNTIME_ASSERT( 0, "Single Voxel Reconstruction is deprecated.  Please choose amongst other options.\n");
            break;
          }
          //--------------------------------------------------
          //  Stuff below this point is using a backward compatible version
          //  of IceNine-2.  Need to be rewritten, but I just don't have the
          //  time.
          //--------------------------------------------------
        case OPT_PARAM:
          {
            RunParameterOptimization();
            stringstream ss;
            ss << i;
            GET_LOG( osLogFile ) << "Wrote New Detector file to " << ss.str() << std::endl;
            WriteDetectorFile( ss.str() );
            break;
          }
        case USE_BND:
          {
//             vVoxelList = XDMUtility::GetBoundaryVoxels( * LatticeSymmetry
//                                                         ::Utilities
//                                                         ::GetSymmetryObj( oSetup.SampleGeometry().GetSampleSymmetry() ),
//                                                         VoxelQueue.Solution().GetInitialSideLength(),
//                                                         vVoxelList,
//                                                         oLocalConfigFile.fMaxMisorientation,
//                                                         oLocalConfigFile.fMaxCost,
//                                                         oLocalConfigFile.fMaxRadius );
//             GET_LOG( osLogFile ) << "Ready to use boundary voxels"
//                                  << " " << oLocalConfigFile.fMaxMisorientation
//                                  << " " << oLocalConfigFile.fMaxCost
//                                  << " " << oLocalConfigFile.fMaxRadius
//                                  << std::endl;
            
//             RUNTIME_ASSERT( vVoxelList.size() > 0 ,
//                             " ERROR:  USE_BND - Boundary criteria does not yield any voxels ");

            RUNTIME_ASSERT( 0, "Boundary Optimization Removed for refactoring\n");
            break;
          }
        case FIT_INT_DECOMP:
          {
            RunIntensityDecomposition();
            break;
          }
        case FIT_LBFS:
          {
            RunLazyBFS();
            break;
          }
        case FIT_LOCAL_OPTIMIZE:
          {
            RunLocalOrientationOptimization();
            break;
          }
        case FIT_STRAIN_LBFS:
          {
            RunStrainLazyBFS();
            break;
          }
        default:
          RUNTIME_ASSERT(0, "\nStart Server: Unknown Command\n");
          exit(0);
      }
      GET_LOG( osLogFile ) << "Done with command: " << i << std::endl;
    }
    Comm.BcastSendCommand( 0, ALLDONE );
  }

  
  //------------------------------------------------------------------------
  //
  //  ProcessCommands  ( ParallelVersion )
  //  - Process all commands ( reconstruct, save file, reload, parameter
  //                           optimization, and so on... )
  //
  //------------------------------------------------------------------------
  bool CParallelServer::Run( )
  {
    stringstream tmpSS;
    tmpSS <<  "ParallelReconLog_0.log";
    osLogFile.SetWritable();        
    osLogFile.open( tmpSS.str().c_str() );

    time_t oStartTime, oStopTime;
    struct tm * pTimeObj;
    
    time ( &oStartTime );
    pTimeObj = localtime ( &oStartTime );
    GET_LOG( osLogFile ) << "StartTime: " <<  asctime( pTimeObj )  << std::endl;
    
    BroadcastExpData();
        
    //---------------------------------------------
    //  Process overall commands
    //---------------------------------------------
    ProcessCommands();
    //---------------------------------------------
    //  Clean up log
    //---------------------------------------------
    time ( &oStopTime );
    pTimeObj = localtime ( &oStopTime );
    GET_LOG( osLogFile ) << "Stop: " <<  asctime( pTimeObj )  << std::endl;
    GET_LOG( osLogFile ) << "Total Time elapsed (sec): "
                         << difftime( oStopTime, oStartTime ) << std::endl;
    osLogFile.close();
    return true;
  }
  
}
