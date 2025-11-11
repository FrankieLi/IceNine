//------------------------------------------------------
//  File:     AsychronousMPI.h
//  Author:   Frankie Li
//  e-mail:   sfli@cmu.edu
//  Purpose:  A simple wrapper for Asychronous MPI calls.
//
//------------------------------------------------------
#ifndef _ASYCHRONOUS_MPI_H
#define _ASYCHRONOUS_MPI_H

#include <mpi.h>
#include <vector>
#include "Serializer.h"
#include <boost/circular_buffer.hpp>
#include <memory>
#include <deque>
#include "Types.h"
#include "Error.h"

namespace XDMMPI
{

  //------------------------------------
  //  General SendEvents
  //------------------------------------
  class CommEvent
  {
  protected:
    const static int DEFAULT_TAG = 0; 
    MPI_Status  oStatus;
    MPI_Request oRequest;
    Size_Type   nSize;
  public:

    virtual ~CommEvent() {}
    virtual bool Completed()
    {
      int nFlag;
      int nErrCode = MPI_Test( &oRequest, &nFlag, &oStatus );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Check complete failed !\n" );
      return nFlag;
    }

    virtual void Wait()
    {
      int nErrCode = MPI_Wait( &oRequest, &oStatus );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Check complete failed !\n" );
    }

    virtual void Advance() { }
    virtual void Process() { }
  };

  //------------------------------------
  //  SendSingleBuffered
  //
  //  Send buffer pBuf, buffered without
  //  first sending size information.
  //  Use this if the buffer size is known
  //  ahead of time.
  //------------------------------------
  template <typename T >
  class SendSingleBuffered : public CommEvent
  {
  private:
    T *pBuf;
  public:
    ~SendSingleBuffered()
    {
      delete pBuf;
    }
        
    void Process ( int nDest, const T & POD )   // object must be a POD
    {
      pBuf  =  new T;
      *pBuf = POD;
      int nErrCode = MPI_Isend( pBuf, sizeof( *pBuf ),  MPI_CHAR,
                                nDest, DEFAULT_TAG,
                                MPI_COMM_WORLD, &oRequest );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "SendSingleBuffered Failed !\n" );
    } 
  };

  class CmplxSend : public CommEvent
  {
  protected:
    MPI_Status  _ByteStatus;
    MPI_Request _ByteRequest;
  public:
    ~CmplxSend()
    {  }
    
    bool Completed()
    {
      int nFlag;
      int nErrCode = MPI_Test( &_ByteRequest, &nFlag, &_ByteStatus );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Check complete failed !\n" );
      if( nFlag )
        return CommEvent::Completed( );
      else
        return false;
    }
    void Wait()
    {
      int nErrCode = MPI_Wait( &_ByteRequest, &_ByteStatus );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Wait failed !\n" );
      CommEvent::Wait();
    }
  };
  
  //------------------------------------
  //  SendBuffered
  //------------------------------------
  class SendBuffered : public CmplxSend
  {
  private:
    GeneralLib::CSerializer _SerialBuffer;
  public:   

    ~SendBuffered()
    {
      //      std::cerr << "Send Buffered being destroyed " << std::endl;
    }

    GeneralLib::CSerializer & Serializer() { return _SerialBuffer; }    
    void Process ( int nDest )
    {
      char * pBuf = const_cast< char * >( _SerialBuffer.GetBuffer() );
      nSize = _SerialBuffer.GetSize();
      int nErrCode = MPI_Isend( &nSize, sizeof( nSize ),
                                MPI_CHAR, nDest, DEFAULT_TAG,
                                MPI_COMM_WORLD, &_ByteRequest );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Send Size Failed !\n" );
      nErrCode = MPI_Isend( pBuf, nSize,  MPI_CHAR,
                            nDest, DEFAULT_TAG,
                            MPI_COMM_WORLD, &oRequest );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Recv Buffer failed !\n" );
    }
  };

  //------------------------------------
  //  SendUnbuffered
  //
  //  Send the buffer pBuf directly through
  //  two commands.  Use this if the buffer
  //  size is unknown to the reciever
  //
  //------------------------------------
  class SendUnbuffered : public CmplxSend
  {
  public:
    ~SendUnbuffered()
    {
      //      std::cerr << "SendUnbuffered being destroeyd " << std::endl;
    }
    void Process ( int nDest, char * pBuf, Size_Type n )
    {
      nSize = n;
      int nErrCode = MPI_Isend( &nSize, sizeof( nSize ),
                                MPI_CHAR, nDest, DEFAULT_TAG,
                                MPI_COMM_WORLD, &_ByteRequest );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Send Size Failed !\n" );
      nErrCode = MPI_Isend( pBuf, nSize,  MPI_CHAR,
                            nDest, DEFAULT_TAG,
                            MPI_COMM_WORLD, &oRequest );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Recv Buffer failed !\n" );
    }
  };

  class CmplxRecv : public CommEvent
  {
  protected:
    int       nSrc;
    int       nState;   // state of this communication
                        // 0 - before send
                        // 1 - sending size information
                        // 2 - size sent
                        // 3 - sending buffer
                        // 4 - completed    

    virtual void RecvAction( int nState ) { RUNTIME_ASSERT(0, "Not implemented" );}
    virtual bool TestAction( )   // Default Test action
    {
      int nFlag;
      int nErrCode = MPI_Test( &oRequest, &nFlag, &oStatus );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Recv Buffer failed !\n" );
      return nFlag;
    }
    
    virtual void Process_( )
    {
      switch ( nState )
      {
        case 0:
          {
            RecvAction( nState );
            nState ++;
            return;
          }
        case 1:
        case 3:
          {
            if( TestAction() )
              nState ++;
            return;
          }
        case 2:
          {
            RecvAction( nState );
            nState ++;
            return;
          }
        case 4:  // do nothing
          return;
      }
    }
  public:
    CmplxRecv() : nState( 0 ) {}
    bool Completed()
    {
      if( nState == 4 )
        return true;
      
      Process_();
      return false;
    }
    void Wait()
    {
      while ( nState != 4 )
        Process_();
    }
    void Advance()  { Process_(); }
  };
  
  //------------------------------------
  //  RecvUnbuffered
  //
  //  Send the buffer pBuf directly through
  //  two commands.  Use this if the buffer
  //  size is unknown to the reciever
  //  NOTE:  Because pBuf is being sent directly,
  //         pBuf is assumed to be held by someone else.
  //         There will be no clean up by RecvUnbuffered
  //------------------------------------
  class RecvUnbuffered : public CmplxRecv
  {
  private:
    char    * pBuf;
    void RecvAction( int nState )
    {
      if( nState == 0 )
      {
        int nErrCode = MPI_Irecv( &nSize, sizeof( nSize ), MPI_CHAR, nSrc,
                                  DEFAULT_TAG, MPI_COMM_WORLD, &oRequest );
        RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Recv Size Failed !\n" );
        return;
      }
      if( nState == 2 )
      {
        int nErrCode = MPI_Irecv( pBuf, nSize,  MPI_CHAR, nSrc,
                                  DEFAULT_TAG, MPI_COMM_WORLD, &oRequest );
        RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Recv Buffer failed !\n" );
        return;
      }
      RUNTIME_ASSERT(0, "Unexpected state in RecvAction\n" );
    }

  public:
    void Process ( int src, char * buf, Size_Type n )
    {
      nSrc  = src;
      pBuf  = buf;
      nSize = n;
      Process_();
    }
  };
  
  //------------------------------------
  //  RecvBuffered
  //  - Item is first buffered then sent
  //
  //------------------------------------
  class RecvBuffered : public CmplxRecv
  {
  private:
    GeneralLib::CDeserializer _SerialBuf;
    void RecvAction( int nState )
    {
      if( nState == 0 )
      {
        int nErrCode = MPI_Irecv( &nSize, sizeof( nSize ), MPI_CHAR, nSrc,
                                  DEFAULT_TAG, MPI_COMM_WORLD, &oRequest );
        RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Recv Size Failed !\n" );

        return;
      }
      if( nState == 2 )
      {
        _SerialBuf.Resize( nSize );
        char * pBuf  = _SerialBuf.GetBuffer();
        int nErrCode = MPI_Irecv( pBuf, nSize,  MPI_CHAR, nSrc,
                                  DEFAULT_TAG, MPI_COMM_WORLD, &oRequest );
        RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "Recv Buffer failed !\n" );
        return;
      }
      RUNTIME_ASSERT(0, "Unexpected state in RecvAction\n" );
    }    
  public:
    GeneralLib::CDeserializer & Deserializer() { return _SerialBuf; }
    void Process( int src )
    {
      nSrc = src;
      Process_();
    }
  };

  
  //------------------------------------
  //  RecvSingleBuffered  - note:  The reciving
  //  buffer is empty until "Completed" or "Wait" is called
  //------------------------------------
  template< typename T >
  class RecvSingleBuffered : public CommEvent
  {
  private:
    T  Buf;
    T* ReturnPtr;
  public:
    void Process( int nSrc, T * pPOD )   // object must be a POD
    {
      ReturnPtr = pPOD;
      int nErrCode = MPI_Irecv( &Buf, sizeof( Buf ),  MPI_CHAR,
                                nSrc, DEFAULT_TAG,
                                MPI_COMM_WORLD, &oRequest );
      RUNTIME_ASSERT( nErrCode == MPI_SUCCESS, "SendSingleBuffered Failed !\n" );
    }
    
    bool Completed()
    {
      if ( CommEvent::Completed() )
      {
        *ReturnPtr = Buf;
        return true;
      }
      return false;
    }

    void Wait()
    {
      CommEvent::Wait();
      *ReturnPtr = Buf;
    }
  };


  //--------------------------------------------------------------------------------------------------

  
  //------------------------------------------------------
  //  Unbuffered AsycnhronousMPI wrapper class
  //  - Use with data that'd stay around during the
  //    entire duration of the communication. (I.e., won't go
  //    out of scope or get destroyed.)
  //------------------------------------------------------
  class AMpi
  {
  public:
    typedef std::shared_ptr< GeneralLib::CSerializer > SerializerPtr;

  protected:

    std::deque< CommEvent * > oSendEventQueue;
    std::deque< CommEvent * > oRecvEventQueue;

    void ClearEventDeque( std::deque< CommEvent * > & Deque )
    {
      typedef std::deque< CommEvent * >::iterator DequeIter;
      DequeIter pCur  = Deque.begin();
      DequeIter pLast = Deque.begin();
      for(; pCur != Deque.end(); ++ pCur )
        if( ! (*pCur)->Completed() )
        {
          std::swap( *pLast, *pCur );
          ++ pLast;
        }
        else
        {
          delete *pCur;
        }
      Deque.erase( pLast, Deque.end() );
    }
    
  public:
    
    ~AMpi()
    {
      typedef std::deque< CommEvent * >::iterator DequeIter;
      DequeIter pCur = oSendEventQueue.begin();
      for(int i = 0; pCur != oSendEventQueue.end(); ++ pCur, ++i )
        delete *pCur ;
      
      pCur = oRecvEventQueue.begin();
      for(; pCur != oRecvEventQueue.end(); ++ pCur )
        delete *pCur ;
    }
    
    //-----------------------------
    //  GetNumActiveRequests
    //   - Return the number of active requests.
    //   - Also automatically clear requests
    //     that are completed (removed from vector)
    //-----------------------------
    int  GetNumActiveRequests();

    bool WaitAllSendRequest( );
    bool WaitAnySendRequest( );

    template <typename T > 
    void SendPODBuffered( int nDest, const T & POD )
    {
      SendSingleBuffered<T> * pEvent = new SendSingleBuffered<T>;
      pEvent->Process( nDest, POD );
      oSendEventQueue.push_back( pEvent );
    }

    template <typename T > 
    void RecvPODBuffered( int nSrc, T & POD )
    {
      RecvSingleBuffered<T> * pEvent = new RecvSingleBuffered<T>;
      pEvent->Process( nSrc, &POD );
      oRecvEventQueue.push_back( pEvent );
    }

    void SendNUnbuffered( int nDest, char * pBuf, Size_Type nSize )
    {
      SendUnbuffered * pEvent = new SendUnbuffered;
      pEvent->Process( nDest, pBuf, nSize );
      oSendEventQueue.push_back( pEvent );
    }
    
    void RecvNUnbuffered( int nSrc, char * pBuf, Size_Type nSize )
    {
      RecvUnbuffered * pEvent = new RecvUnbuffered;
      pEvent->Process( nSrc, pBuf, nSize );
      oRecvEventQueue.push_back( pEvent );
    }
    
    void WaitSent()
    {
      RUNTIME_ASSERT( oSendEventQueue.size() > 0, " There is NO event\n");
      oSendEventQueue.front()->Wait();

      delete oSendEventQueue.front();
      oSendEventQueue.pop_front();
    }
    
    void WaitRecv()
    {
      RUNTIME_ASSERT( oRecvEventQueue.size() > 0, " There is NO event\n");
      oRecvEventQueue.front()->Wait();
      delete oRecvEventQueue.front();
      oRecvEventQueue.pop_front();
    }

    void WaitAndAdvanceRecv()
    {
      RUNTIME_ASSERT( oRecvEventQueue.size() > 0, " There is NO event\n");
      
      typedef std::deque< CommEvent * >::iterator DequeIter;
      for( DequeIter pCur = oRecvEventQueue.begin();
           pCur != oRecvEventQueue.end(); ++ pCur )
        (*pCur)->Advance();
      
      oRecvEventQueue.front()->Wait();
      delete oRecvEventQueue.front();
      oRecvEventQueue.pop_front();
    }

    //--------------------------
    //  Check to see if the earliest event
    //  is completed
    //--------------------------
    bool SendCompleted()
    {
      if( oSendEventQueue.size() == 0 )
        return true;
      if( oSendEventQueue.front()->Completed() )
      {
        delete oSendEventQueue.front();
        oSendEventQueue.pop_front();
        return true;
      }
      return false;
    }

    bool RecvCompleted()
    {
      if( oRecvEventQueue.size() == 0 )
        return true;
      if( oRecvEventQueue.front()->Completed() )
      {
        delete oRecvEventQueue.front();
        oRecvEventQueue.pop_front();
        return true;
      }
      return false;
    }
        
    //-----------------------------
    // Recv
    // Desc:  A wrapper around recv that
    //        takes care of data size and
    //        serialization
    //-----------------------------
//     void Recv( int nSrc, GeneralLib::CDeserializer & oSerializer )
//     {
      
//     }
    
    //-----------------------------
    //  Wait for recv requests
    //-----------------------------
    bool WaitAllRecvRequest( );
    bool WaitAnyRecvRequest( );

    void ClearSendEvents( )
    {
      ClearEventDeque( oSendEventQueue );
    }
  };



  
}


#endif
