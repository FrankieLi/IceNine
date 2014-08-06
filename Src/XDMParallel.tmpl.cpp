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
//  XDMParallel.tmpl.cpp
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Interface for templated send and recv functions
//
//
////////////////////////////////////////////////////////////




namespace XDMParallel
{
  
  //----------------------------------
  //  SendWorkUnit
  //
  //  Work unit is sent
  //
  //  NEED ERROR HANDLING AND RECOVERY
  //----------------------------------
  template< typename WorkUnitT >
  void XDMCommunicator::SendWorkUnit( Int nDestPE, const WorkUnitT & oWorkUnit )
  {
    XDMMPI::SendBuffered * pEvent = new XDMMPI::SendBuffered;
    pEvent->Serializer().InsertCompactObj( oWorkUnit );
    pEvent->Process( nDestPE );
    oSendEventQueue.push_back( pEvent );
    ClearSendEvents();   // This is necessary because we don't need to check sent items
  }

  //----------------------------------
  //  The list version
  //----------------------------------
  template< typename WorkUnitT >
  void XDMCommunicator::SendWorkUnitList( Int nDestPE, const vector<WorkUnitT> & oWorkUnitList )
  {
    XDMMPI::SendBuffered * pEvent = new XDMMPI::SendBuffered;
    pEvent->Serializer().InsertCompactVector( oWorkUnitList );
    pEvent->Process( nDestPE );
    oSendEventQueue.push_back( pEvent );
    ClearSendEvents();   // This is necessary because we don't need to check sent items
  }
 
  //----------------------------------
  //  BcastSendWorkUnit
  //  Using broadcast to send work units out
  //----------------------------------
  template< typename WorkUnitT >
  Bool XDMCommunicator::BcastSend( Int nSrcPE, const WorkUnitT & oWorkUnit )
  {
    Int nErrCode;
    CSerializer oSendBuf;
    oSendBuf.InsertCompactObj( oWorkUnit );
    char * pWorkUnitBuf = const_cast< char * >( oSendBuf.GetBuffer() );
    Size_Type nBufSize = oSendBuf.GetSize();
    nErrCode = MPI_Bcast( &nBufSize, sizeof( nBufSize), MPI_CHAR, nSrcPE, MPI_COMM_WORLD ); 
    nErrCode = MPI_Bcast( pWorkUnitBuf, nBufSize, MPI_CHAR, nSrcPE, MPI_COMM_WORLD);
    RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "\nDriver::SendWorkUnit: Error!  Cannot send work unit!\n" );
  
    return true;
  }
  
  //----------------------------------
  //  BcastSendList
  //  Purpose:  Sending general messages via Bcast
  //----------------------------------
  template< typename WorkUnitT >
  Bool XDMCommunicator::BcastSendList( Int nSrcPE, const vector<WorkUnitT> & oWorkUnitList )
  {
    Int nErrCode;
    CSerializer oSendBuf;
    oSendBuf.InsertCompactVector( oWorkUnitList );
    char * pWorkUnitBuf = const_cast< char * >( oSendBuf.GetBuffer() );
    Size_Type nBufSize = oSendBuf.GetSize();
    nErrCode = MPI_Bcast( &nBufSize, sizeof( nBufSize), MPI_CHAR, nSrcPE, MPI_COMM_WORLD ); 
    nErrCode = MPI_Bcast( pWorkUnitBuf, nBufSize, MPI_CHAR, nSrcPE, MPI_COMM_WORLD);
    RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "\nDriver::SendWorkUnit: Error!  Cannot send work unit!\n" );
  
    return true;
  }

  //----------------------------------
  //  BcastRecv
  //  Purpose:  Sending general messages via Bcast
  //----------------------------------
  template< typename WorkUnitT >
  Bool XDMCommunicator::BcastRecv( Int nSrcPE, WorkUnitT *oWorkUnit )
  {
    Int nErrCode;
    Size_Type nBufSize;

    nErrCode = MPI_Bcast( &nBufSize, sizeof( nBufSize), MPI_CHAR, nSrcPE, MPI_COMM_WORLD ); 
    CDeserializer oClientUnit( nBufSize );
    char *pClientBuf = oClientUnit.GetBuffer();
    nErrCode = MPI_Bcast( pClientBuf, nBufSize, MPI_CHAR, nSrcPE, MPI_COMM_WORLD);
    RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "\nDriver::RecvWorkUnit: Error!  Cannot recv work unit!\n" );
  
    return oClientUnit.GetCompactObj( oWorkUnit );
  }
  
  //----------------------------------
  //  BcastRecvList
  //  Purpose:  Sending general messages via Bcast
  //----------------------------------
  template< typename WorkUnitT >
  Bool XDMCommunicator:: BcastRecvList( Int nSrcPE, vector<WorkUnitT> & oWorkUnit )
  {
    Int nErrCode;
    Size_Type nBufSize;

    nErrCode = MPI_Bcast( &nBufSize, sizeof( nBufSize), MPI_CHAR, nSrcPE, MPI_COMM_WORLD ); 
    CDeserializer oClientUnit( nBufSize );
    char *pClientBuf = oClientUnit.GetBuffer();
    nErrCode = MPI_Bcast( pClientBuf, nBufSize, MPI_CHAR, nSrcPE, MPI_COMM_WORLD);
    RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "\nDriver::RecvWorkUnit: Error!  Cannot recv work unit!\n" );
  
    return oClientUnit.GetCompactVector( oWorkUnit );
  }
 
  //----------------------------------
  //  RecvWorkUnit
  //
  //  Return the command acompanied by the work unit
  //
  //  (everything's returned by reference parameters -- may change to
  //  tuples in the future)
  //
  //  TODO:  handle the error and the status
  //
  //----------------------------------
  template< typename WorkUnitT >
  Bool XDMCommunicator::RecvWorkUnit( Int nSourcePE, WorkUnitT *oWorkUnit )
  {
    XDMMPI::RecvBuffered Event;
    Event.Process( nSourcePE );
    Event.Wait();  // blocking on recv of work unit now - will alter in the future
    bool bDeserialized = Event.Deserializer().GetCompactObj( oWorkUnit );
    return bDeserialized;
  }

  //----------------------------------
  //  The list version
  //----------------------------------
  template< typename WorkUnitT >
  Bool XDMCommunicator::RecvWorkUnitList( Int nSourcePE, vector<WorkUnitT> &oWorkUnitList )
  {
    XDMMPI::RecvBuffered Event;
    Event.Process( nSourcePE );
    Event.Wait();  // blocking on recv of work unit now - will alter in the future
    bool bDeserialized = Event.Deserializer().GetCompactVector( oWorkUnitList );
    return bDeserialized;
  }
  
}
