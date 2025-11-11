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
//  XDMServer.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of the server side of the paralle IceNine
// 
//
//
////////////////////////////////////////////////////////////


#ifndef _XDM_SERVER_H_
#define _XDM_SERVER_H_
#include "XDMCommCore.h"
#include "BreadthFirstReconstructor.h"
#include "MicIO.h"
#include <memory>
#include "ParameterOptimization.h"
namespace XDMParallel
{
  using namespace core;
  //------------------------------------------
  //  Implementation of the server of IceNine
  //
  //------------------------------------------
  class CParallelServer : public core::CParallelDetails
  {
  private:
    typedef std::shared_ptr<MicIOBase>    MicFilePtrT;
    
  protected:
    CConfigFile oLocalConfigFile;
    Int         nMyID;
    Int         nProcessingElements;

    //----------------------------------
    //  BroadcastExpData
    //
    //  Purpose:  Distribute initial experimental data along
    //            with parameters
    //----------------------------------
    void BroadcastExpData( );
    
    //----------------------------------
    //
    //  WriteFitResult
    //
    //  Parameter:  Postfix is used to for intermediate files
    //
    //----------------------------------
    void WriteFitResult( MicFilePtrT pMic , const string & oPostFix = string(""),
                         bool bReduceToFZ = false );
    
    //----------------------------------
    //  WriteDetectorFiles
    //
    //  Output optimized detector parameters
    //----------------------------------
    void WriteDetectorFile( const string & oPostFix = string("") );
    void ProcessCommands();
        
    //---------------------------------------------------------------//
    //                   Available functions
    //---------------------------------------------------------------//
    void RunIntensityDecomposition();
    void RunLazyBFS();
    void RunStrainLazyBFS();
    void RunLocalOrientationOptimization();
    void RunParameterOptimization();
    
  public:

    CParallelServer( const CConfigFile & oConfigFile,
                     int nCommRank, int nPE ) : CParallelDetails( oConfigFile ),
                                                oLocalConfigFile( oConfigFile ),
                                                nMyID( nCommRank ),
                                                nProcessingElements( nPE )
    {
      std::cout << "Running parallel initialization" << std::endl;
      oSetup.InitializeWithDataFiles( oConfigFile );
      //oSetup.InitializeWithDataFiles( oConfigFile, true );   // run server mode if intensity and sparse matrix are needed
      oSampleCenter = oSetup.InputParameters().SampleCenter;
      fSampleRadius = oSetup.InputParameters().SampleRadius;

      SpecifyOptions( oSetup.InputParameters() );
      RUNTIME_ASSERT( nCommRank == 0, "Only rank 0 should be running the server code\n");
    }
    
    //------------------------------------------------------------------------
    //
    //  ProcessCommands  ( ParallelVersion )
    //  - Process all commands ( reconstruct, save file, reload, parameter
    //                           optimization, and so on... )
    //
    //------------------------------------------------------------------------
    bool Run(  );

  };

}


#endif
