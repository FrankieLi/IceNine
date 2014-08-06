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
/////////////////////////////////////////////////////
//
//  ParallelDriver.cpp
//
//  Author:  Frankie Li (sfli@cmu.edu)
//   
//  Purpose:  ParallelDriver
//
/////////////////////////////////////////////////////



#ifndef _P_DRIVER_
#define _P_DRIVER_

#include "XDMClient.h"
#include "XDMServer.h"
//#include "LazyBFS.h"

#include <string>

namespace XDMParallel
{
  class CParallelProcessing
  {
  public:
    void Process( int argc, char* args[],
                  const string & sConfigFilename )
    {
      int nMyRank, nProcessingElements; 
      MPI_Init( &argc, &args ); 
      MPI_Comm_rank( MPI_COMM_WORLD, & nMyRank );
      MPI_Comm_size( MPI_COMM_WORLD, &nProcessingElements );
      
      if( nMyRank == 0 )
      {
        CConfigFile oConfigFile;
        RUNTIME_ASSERT( oConfigFile.InputConfigParameters( sConfigFilename ),
                        "CParallelProcessing:  Failed to reading in Config file\n  ");
        CParallelServer oServer( oConfigFile, nMyRank, nProcessingElements );
        //LazyBFSServer oServer( oConfigFile, nMyRank, nProcessingElements );
        oServer.Run();
      }
      else
      {
        CParallelClient oClient;
        //LazyBFSClient oClient;
        oClient.Run( nMyRank );
      }
      
      MPI_Barrier( MPI_COMM_WORLD );
      MPI_Finalize();
    }
  };
  
  
}

#endif
