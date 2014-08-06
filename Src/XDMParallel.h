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
//  XDMParallel.h
//  Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of basic communication wrapper via MPI.
//
////////////////////////////////////////////////////////////

#ifndef _XDMPARALLEL_H_
#define _XDMPARALLEL_H_

#include <sstream>
#include "Serializer.h"
#include "MicIO.h"
#include <mpi.h>
#include <cstdio>
#include <ctime>
#include <map>
#include "AsynchronousMPI.h"


// debug
#include <ostream>
using std::ofstream;
namespace XDMParallel
{

  //--------------------------
  //
  //  Enumerated type for commands  (Change to map of function)
  //
  //--------------------------
  enum EXDMCommands
    {
      FIT,
      LBFS_REFIT,
      FIT_ADP,
      FIT_INT_DECOMP,
      FIT_LBFS,
      FIT_STRAIN_LBFS,
      FIT_LOCAL_OPTIMIZE,
      REPORT,
      FIT_VECTOR,
      FIT_REGRID,
      REPORT_VECTOR,
      FIT_MC_COMPLETE,
      FIT_MC,
      FIT_MC_LIST,
      OPT_PARAM,
      EVAL_OVERLAP,
      USE_BND,
      REPORT_MC,
      REPORT_MC_LIST,
      SET_EXP_PARAM,
      GENERAL_TAG,
      ALLDONE,
      WAIT,
      PROCESS_DONE   // gotta formalize this -- maybe wait has to be used more carefully.
    };


  
  //--------------------------
  //  SSimpleCommand
  //
  //  A simple, compact type that
  //  contains a command and the origin
  //--------------------------
  struct SSimpleCommand
  {
    int nCommand;
    int nSource;
  };

  //----------------------------------
  //  XDMCommunicator
  //
  //----------------------------------
  class XDMCommunicator : XDMMPI::AMpi
  {
    
  public:
    //----------------------------------
    //
    // SignalPE
    //
    // Precondition:  nFirstPE -> nLastPe exist
    //
    // Purpose:  Send signal to the processing element
    //           nFirstPE, nLastPE inclusively.  The command
    //           sent will be nCommand
    //
    // TODO:     Write MPI_ERRR_HANDLER!!!
    //
    //----------------------------------
    void SendCommand( Int nSrcPE, Int nFirstPE, Int nLastPE, Int nCommand );
    void SendCommand( Int nSrcPE, Int nDestPE, Int nCommand );

    void DEBUG_SendCommand( Int nSrcPE, Int nDestPE, Int nCommand, ofstream & os );
  
    //----------------------------------
    // RecvCommand
    //----------------------------------
    void RecvCommand( Int *nSrcPE, Int *nCommand );
  
    //----------------------------------
    //  BcastCommand
    //----------------------------------
    Bool BcastSendCommand( Int SrcPE, Int nCommand);
    Bool BcastRecvCommand( Int nSrcPE, Int *nCommand );
    
    //----------------------------------
    //  SendWorkUnit
    //  Send workunits over MPI
    //  (Designed for POD objects only)
    //----------------------------------
    template< typename WorkUnitT >
    void SendWorkUnit( Int nDestPE, const WorkUnitT & oWorkUnit );

    //----------------------------------
    //  SendWorkUnitList
    //  Send vectors of workunit over MPI
    //  (Designed for vector<POD> objects only)
    //----------------------------------
    template< typename WorkUnitT >
    void SendWorkUnitList( Int nDestPE, const vector<WorkUnitT> & oWorkUnit );

    //----------------------------------
    //  RecvWorkUnit
    //
    //  Return the command acompanied by the work unit
    //
    //  (everything's returned by reference parameters -- may change to
    //  tuples in the future)
    //  (Designed for POD objects only)
    //----------------------------------
    template< typename WorkUnitT >
    Bool RecvWorkUnit( Int nSourcePE, WorkUnitT  *oWorkUnit );
  
    //----------------------------------
    //  Save as RecvWorkUnit, but the vector form
    //----------------------------------
    template< typename WorkUnitT >
    Bool RecvWorkUnitList( Int nSourcePE, vector<WorkUnitT> &oWorkUnit );
  
    //----------------------------------
    //  BcastSend
    //  Purpose:  Sending general messages via Bcast
    //----------------------------------
    template< typename WorkUnitT >
    Bool BcastSend( Int nSrcPE, const WorkUnitT & oWorkUnit );
  
    //----------------------------------
    //  BcastSendList
    //  Purpose:  Sending general messages via Bcast
    //----------------------------------
    template< typename WorkUnitT >
    Bool BcastSendList( Int nSrcPE, const vector<WorkUnitT> & oWorkUnit );

    //----------------------------------
    //  BcastSend
    //  Purpose:  Sending general messages via Bcast
    //----------------------------------
    template< typename WorkUnitT >
    Bool BcastRecv( Int nSrcPE, WorkUnitT *oWorkUnit );
  
    //----------------------------------
    //  BcastSendList
    //  Purpose:  Sending general messages via Bcast
    //----------------------------------
    template< typename WorkUnitT >
    Bool BcastRecvList( Int nSrcPE, vector<WorkUnitT> & oWorkUnit );
  
  };
}

#include "XDMParallel.tmpl.cpp"

#endif
