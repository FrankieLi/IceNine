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

#include "XDMParallel.h"

namespace XDMParallel
{
  
  const char*  pXDMCommands[] =
    {
      "FIT",
      "REPORT",
      "FIT_VECTOR",
      "FIT_REGRID",
      "REPORT_VECTOR",
      "FIT_MC_COMPLETE",
      "FIT_MC",
      "FIT_MC_LIST",
      "OPT_PARAM",
      "EVAL_OVERLAP",
      "REPORT_MC",
      "REPORT_MC_LIST",
      "SET_EXP_PARAM",
      "GENERAL_TAG",
      "ALLDONE",
      "WAIT",
      "PROCESS_DONE"   // gotta formalize this -- maybe wait has to be used more carefully.
    };
  
  //----------------------------------
  //  SendCommand
  //----------------------------------
  void XDMCommunicator::SendCommand( Int nSrcPE, Int nDestPE, Int nCommand )
  {
    SSimpleCommand oCommand;
    oCommand.nCommand = nCommand;
    oCommand.nSource = nSrcPE;
    SendPODBuffered( nDestPE, oCommand );
    ClearSendEvents();   // This is necessary because we don't need to check sent items
  }

  //----------------------------------
  //----------------------------------
  void XDMCommunicator::DEBUG_SendCommand( Int nSrcPE, Int nDestPE, Int nCommand, ofstream & os )
  {
    os << " Command sent: " << pXDMCommands[ nCommand ]  << std::endl;
    SendCommand( nSrcPE, nDestPE, nCommand );
  }
  
  //----------------------------------
  //  BcastCommand
  //----------------------------------
  Bool XDMCommunicator::BcastSendCommand( Int nSrcPE, Int nCommand )
  {
    SSimpleCommand oCommand;
    oCommand.nCommand = nCommand;
    oCommand.nSource = nSrcPE;
    Int nErrCode = MPI_Bcast( &oCommand, sizeof( oCommand ), MPI_CHAR, nSrcPE, MPI_COMM_WORLD );
    RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "\nDriver::SendCommand: ERROR:  unable to send signal\n" );
    return true;
  }
  
  //----------------------------------
  // RecvCommand
  //----------------------------------
  void XDMCommunicator::RecvCommand( Int *nSrcPE, Int *nCommand )
  {
    SSimpleCommand oCommand;    
    RecvPODBuffered( MPI_ANY_SOURCE, oCommand );
    WaitRecv(); // blocking recv still
    *nCommand = oCommand.nCommand;
    *nSrcPE   = oCommand.nSource;
  }
  
  //----------------------------------
  //  BcastCommand
  //----------------------------------
  Bool XDMCommunicator::BcastRecvCommand( Int nSrcPE, Int *nCommand )
  {
    SSimpleCommand oCommand;
    Int nErrCode = MPI_Bcast( &oCommand, sizeof( oCommand ), MPI_CHAR, nSrcPE, MPI_COMM_WORLD );
    RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "\nDriver::SendCommand: ERROR:  unable to send signal\n" );
    *nCommand = oCommand.nCommand;
    return true;
  }
  
  //----------------------------------
  //
  //  Send signal to the processing element
  //  nFirstPE, nLastPE inclusively.  The command
  //  sent will be nCommand
  //----------------------------------
  void XDMCommunicator::SendCommand( Int nSrcPE, Int nFirstPE, Int nLastPE, Int nCommand )
  {
    for( Int nCurPE = nFirstPE; nCurPE <= nLastPE; nCurPE ++ )
      SendCommand( nSrcPE, nCurPE, nCommand );
  }
  
 
}
