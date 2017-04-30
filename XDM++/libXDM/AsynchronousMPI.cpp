//------------------------------------------------------
//  File:     AsychronousMPI.cpp
//  Author:   Frankie Li
//  e-mail:   sfli@cmu.edu
//  Purpose:  A simple wrapper for Asychronous MPI calls.
//
//------------------------------------------------------

#include "AsynchronousMPI.h"

namespace XDMMPI
{

  //----------------------------------------------------------
  // Send
  //----------------------------------------------------------
  //  void AMpi::Send( int nDest, SerializerPtr pSerializer )
  //  {
//     oSendInfoQueue.push( SCommInfo() );
//     SCommInfo & oCommInfo = oSendInfoQueue.back();
//     oCommInfo.nBufferSize = pSerializer->GetSize();
//     oCommInfo.pData       = pSerializer;

//     MPI_Request  oSizeSentReq;   // to be thrown away
//     int nErrCode = MPI_Isend( &( oCommInfo.nBufferSize ), sizeof( oCommInfo.nBufferSize ),
//                               MPI_CHAR, nDest, DEFAULT_TAG, MPI_COMM_WORLD,
//                               &oSizeSentReq );
    
//     RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "AMpi::Send: Error!  Cannot set serialized buffer !\n" );
//     char * pBuf = const_cast< char * >( pSerializer->GetBuffer() );
    
//     nErrCode = MPI_Isend( pBuf, oCommInfo.nBufferSize, MPI_CHAR,
//                           nDest, DEFAULT_TAG, MPI_COMM_WORLD, &( oCommInfo.oCommRequest) );
//     RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "AMpi::Send: Error!  Cannot set serialized buffer !\n" );
//  }
  
  //----------------------------------------------------------
  // Recv
  //----------------------------------------------------------
//   void AMpi::Recv( int nSrc, GeneralLib::CDeserializer & oSerializer )
//   {
  //   int nErrCode;

//     Size_Type nBufSize;
//     MPI_Request  oSizeRecvReq;   // to be thrown away
//     nErrCode = MPI_Irecv( &nBufSize, sizeof( nBufSize ), MPI_CHAR,
//                           nSrc, DEFAULT_TAG, MPI_COMM_WORLD, &oSizeRecvReq );
//     MPI_Status oStatus;
//     RUNTIME_ASSERT( MPI_SUCCESS != MPI_Wait( &oSizeRecvReq, &oStatus ),
//                     "AMpi::Recv:  communication problems with recv.\n" );  // have to block on buffer size
    
//     char *pClientBuf = oSerializer.GetBuffer();
//     nErrCode = MPI_Irecv( pClientBuf, nBufSize, MPI_CHAR,
//                           nSourcePE, DEFAULT_TAG, MPI_COMM_WORLD, &oRecvStatus );
//     RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "\nDriver::RecvWorkUnit: Error!  Cannot recv work unit!\n" );
//  }
}
