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
//  XDMParallel.cpp
//  Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of basic communication wrapper via MPI.
//
////////////////////////////////////////////////////////////
#include "XDMClient.h"


namespace XDMParallel
{
 
  //----------------------------------
  //  RecvExpParameters
  //----------------------------------
  void CParallelClient::RecvExpParameters( Int nServerPE,
                                           SEnergyOpt & oEnergyParam,
                                           vector<SDetParamMsg> &vDetectorParams )
  {
    Comm.RecvWorkUnit    ( nServerPE, &oEnergyParam  );
    Comm.RecvWorkUnitList( nServerPE, vDetectorParams );
  }

  //---------------------------------------------------------------------------
  //  RunLazyBFS
  //
  //---------------------------------------------------------------------------
  int CParallelClient::RunLazyBFS( Int nMyID )
  {

    //--------------------------------
    //  TODO:  Make a virtual superclass that
    //         encompasses all LBFS clients.
    //         Use dynamic dispatch to clean this up.
    //--------------------------------
    //    typedef Reconstruction::DiscreteRefinement<SVoxel> DR;
    time_t oStartTime, oStopTime;
    switch( oSetup.InputParameters().nMicGridType  )
    {
      case eTriangular:
        {
          typedef ParallelReconstructor::LazyBFSClient<SVoxel> LBFS;
          LBFS Client( oSetup, nMyID, osLogFile );
          Client.Initialize();
          GET_LOG( osLogFile ) << "Beginning Lazy BFS: Triangle" << std::endl;
          time( &oStartTime );
          Client.Process();
          time( &oStopTime );
          break;
        }
      case eSquare:
        {
          typedef Reconstruction::DiscreteRefinement<SquareVoxel> DR;
          typedef MicAnalysis::RectMicGrid<SquareVoxel> Grid;
          typedef ParallelReconstructor::LazyBFSClient<SquareVoxel, DR, Grid> LBFS;
          LBFS Client( oSetup, nMyID, osLogFile );
          Client.Initialize();
          GET_LOG( osLogFile ) << "Beginning Lazy BFS: Square" << std::endl;
          time( &oStartTime );
          Client.Process();
          time( &oStopTime );
          break;
        }
    }
    double oTimeDiff = difftime( oStopTime, oStartTime );
    GET_LOG( osLogFile ) << "Total Processing Time: (sec) " << oTimeDiff << "... Client Exiting " << std::endl;
    return PROCESS_DONE;
  }

  //---------------------------------------------------------------------------
  //  RunParameterOptimization
  //
  //---------------------------------------------------------------------------
  int CParallelClient::RunParameterOptimization( Int nMyID )
  {
    RUNTIME_ASSERT( oSetup.InputParameters().nMicGridType == eTriangular,
                    "Error:  Only triangular voxels accepted for parameter MC right now\n" );
    typedef ParallelReconstructor::ParameterOptimization::ParameterOptimizationClient<> PMC;
    PMC  Client( oSetup, nMyID, osLogFile );
    //  Client.Initialize();
    GET_LOG( osLogFile ) << "Running Parameter Optimization" << std::endl;
    time_t oStartTime, oStopTime;
    time( &oStartTime );
    
    Client.Process();
    time( &oStopTime );
    double oTimeDiff = difftime( oStopTime, oStartTime );
    GET_LOG( osLogFile ) << "Total Processing Time: (sec) " << oTimeDiff << "... Client Exiting " << std::endl;
    return PROCESS_DONE;
  }
  
  //---------------------------------------------------------------------------
  //  RunLocalOrientationOptimization
  //
  //---------------------------------------------------------------------------
  int CParallelClient::RunLocalOrientationOptimization( Int nMyID )
  {
    typedef Reconstruction::DiscreteRefinement<SVoxel> DR;
    typedef Reconstruction::LocalReconstructionAdaptor< DR > LDR;
    typedef ParallelReconstructor::LazyBFSClient<SVoxel, LDR> LBFS;
    LBFS Client( oSetup, nMyID, osLogFile );
    Client.Initialize();
    
    GET_LOG( osLogFile ) << "Begin local optimization" << std::endl;
    
    time_t oStartTime, oStopTime;
    time( &oStartTime );
    Client.Process();
    time( &oStopTime );
    
    double oTimeDiff = difftime( oStopTime, oStartTime );
          
    GET_LOG( osLogFile ) << "Total Processing Time: (sec) " << oTimeDiff << "... Client Exiting " << std::endl;
    
    return PROCESS_DONE;
  }
  
  //---------------------------------------------------------------------------
  //  InitializeClient
  //---------------------------------------------------------------------------
  void CParallelClient::InitializeClient()
  {
    Size_Type nBufSize;
    GET_LOG( osLogFile ) << "Client Started, Prepared to recv:" << std::endl;
    MPI_Bcast( &nBufSize, sizeof(nBufSize), MPI_CHAR, 0,  MPI_COMM_WORLD); 
    CDeserializer oSimulationCheckpoint( nBufSize );
    char * pBuf = oSimulationCheckpoint.GetBuffer();  
    MPI_Bcast( pBuf, nBufSize, MPI_CHAR, 0,  MPI_COMM_WORLD );
    GET_LOG( osLogFile ) << "Recv'd " << nBufSize << "bytes, begin processing commands "  << std::endl;

    //---------------------------------------
    GET_LOG( osLogFile ) << "Recv-ing detector data: [bytes sent] "
                         <<  std::endl;
    oSetup.Data().BCast_Recv();
    //---------------------------------------
    
    Restore( oSimulationCheckpoint );
    GET_LOG( osLogFile ) << "... Checkpoint Restored "<< std::endl;
    SpecifyOptions( oSetup.InputParameters() );
    GET_LOG( osLogFile ) << "... Parameter Specified "<< std::endl;
    GET_LOG( osLogFile ) << "... Max Local Resolution  " << oSetup.InputParameters().nMaxLocalResolution << std::endl;
    GET_LOG( osLogFile ) << "... Simulator Initialized "<< std::endl;
  
  }

  //-----------------------------------------------------------------
  //  StartClient
  //-----------------------------------------------------------------
  void CParallelClient::Run( int nMyID )
  {
    //--------------
    //  TODO change this to something custom
    //--------------
    stringstream tmpSS;
    tmpSS <<  "ParallelReconLog_" << nMyID << ".log";
    if( nMyID < 3 )                // really could be changed
      osLogFile.SetWritable();        
    osLogFile.open( tmpSS.str().c_str() );

    struct tm * pTimeObj;
    time_t oStartTime, oStopTime;
    time ( &oStartTime );
    pTimeObj  = localtime ( &oStartTime );
    GET_LOG( osLogFile ) << "StartTime: " <<  asctime( pTimeObj )  << std::endl;

    InitializeClient();

    //---------------------------------------------
    //  Process overall commands
    //---------------------------------------------
    Int nCommand;
    do
    {
      Comm.BcastRecvCommand( 0, &nCommand );
      switch( nCommand ) 
      {
        GET_LOG( osLogFile ) << "Processing command: " << nCommand << std::endl;
        case FIT:
        case FIT_ADP:
          RUNTIME_ASSERT(0 , "ReconstructElements is no longer available\n");
          break;
        case OPT_PARAM:
          nCommand = RunParameterOptimization( nMyID );
          break;
        case FIT_LBFS:
          nCommand = RunLazyBFS( nMyID );
          break;
        case FIT_LOCAL_OPTIMIZE:
          nCommand = RunLocalOrientationOptimization( nMyID );
          break;
        case ALLDONE:
          break;
        case PROCESS_DONE:
          break;
        case WAIT:
          break;
        default:
          RUNTIME_ASSERT( 0, "\n StartClient: Recvn-ed unknown command\n" );
          break;
      }
    } while ( nCommand != ALLDONE );

    //---------------------------------------------
    //  Clean up log
    //---------------------------------------------
    time ( &oStopTime );
    pTimeObj = localtime ( &oStopTime );
    GET_LOG( osLogFile ) << "Stop: " <<  asctime( pTimeObj )  << std::endl;
    GET_LOG( osLogFile ) << "Total Time elapsed (sec): "
                         << difftime( oStopTime, oStartTime ) << std::endl;
 
    osLogFile.close();
  }
  
}

