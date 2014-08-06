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
//  SimulationData.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//   Purpose:  This organizes the simulation data and associates them properly.
//             In another words, a set of angular range (omega) is associated
//             with each detector.  There's a one-to-one corrospondance of angular
//             range with image.
//
//   Note:     Currently this is quite limited.
//
////////////////////////////////////////////////////////////


#ifndef SIMULATION_DATA_H_
#define SIMULATION_DATA_H_

#include "Detector.h"
#include "ImageData.h"
#include "Serializer.h"
#include "boost/multi_array.hpp"
#include "IteratorAdapter.h"

namespace XDMSimulation
{
  const Size_Type NoMatch = std::string::npos;
  
  //---------------------------------------------------
  //
  //  Class Simulation Data 
  //
  //  TODO:  Make this into a template so that a switch
  //  between CImageData and CSearchableImageData can
  //  be interchangable
  //---------------------------------------------------
  template< class RasterType>
  class SimulationData
  {
  public:
    typedef SimulationData< RasterType >                    Self;
    typedef CSearchableImageData                            RasterT;
    typedef boost::multi_array<CSearchableImageData, 2>     ImageMapT;
    typedef ImageMapT::index                                ImageMapIndexT;
    //  typedef ImageMapT::iterator                             ImageMapIteratorT;
    
    typedef boost::multi_array_types::index_range range;


    typedef typename XDMUtility::Const_SerializingIterator< const ImageMapT >   Const_ImageIteratorT;
    typedef typename XDMUtility::SerializingIterator< ImageMapT >               ImageIteratorT;
    
  public:
    ImageMapIndexT nNumIntervals;
    ImageMapIndexT nNumDetectors;
    Int nImageWidth;
    Int nImageHeight;
    ImageMapT mImageMap;

    SimulationData(){};
    
    SimulationData( ImageMapIndexT nIntervals, ImageMapIndexT nDetectors,
                    Int nImWidth, Int nImHeight ): nNumIntervals( nIntervals ),
                                                   nNumDetectors( nDetectors ),
                                                   nImageWidth( nImWidth ),
                                                   nImageHeight(nImHeight )
    {
      Initialize( nIntervals, nDetectors, nImWidth, nImHeight );
    }

    //---------------------------------------------
    //  Initialize
    //---------------------------------------------
    void Initialize( ImageMapIndexT nIntervals, ImageMapIndexT nDetectors,
                     Int nImWidth, Int nImHeight )
    {
      nNumIntervals = nIntervals;
      nNumDetectors = nDetectors;
      nImageWidth = nImWidth;
      nImageHeight = nImHeight;
      mImageMap.resize( boost::extents[ nIntervals ][ nDetectors ] );
      
      for( ImageMapIndexT i = 0; i < nNumIntervals; i ++)
      {
        for ( ImageMapIndexT j = 0; j < nNumDetectors; j ++ )
        {
          mImageMap[i][j].Resize( nImWidth, nImHeight );
        }
      }
    };
    
    //---------------------------------------------
    //  Clear
    //---------------------------------------------
    void Clear( )
    { 
      for( ImageMapIndexT i = 0; i < nNumIntervals; i ++)
      {
        for ( ImageMapIndexT j = 0; j < nNumDetectors; j ++ )
        {
          mImageMap[i][j].ClearImage();
        }
      }
      mImageMap.resize( boost::extents[ 0 ][ 0 ] );
    };
    
    //---------------------------------------------
    //  Save
    //---------------------------------------------
    Bool Save( CSerializer & oSerialBuf ) const
    {

      Bool bSuccess = true;
      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nNumIntervals );
      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nNumDetectors );
     
      for( ImageMapIndexT i = 0; i < nNumIntervals; i ++)
      {
        for ( ImageMapIndexT j = 0; j < nNumDetectors; j ++ )
        {
          bSuccess = bSuccess && mImageMap[i][j].Save( oSerialBuf );
        }
      }
      return bSuccess;
    }

    //---------------------------------------------
    //  Restore
    //
    //  Note:  bBuildSearchTree can be turned off in cases where a search
    //         tree is not useful.
    //---------------------------------------------
    Bool Restore( CDeserializer & oSerialBuf, bool bBuildSearchTree = true )
    {
      Bool bSuccess = true;
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nNumIntervals );
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nNumDetectors );
      mImageMap.resize( boost::extents[ nNumIntervals ][ nNumDetectors ] );

      for( ImageMapIndexT i = 0; i < nNumIntervals; i ++)
      {
        for ( ImageMapIndexT j = 0; j < nNumDetectors; j ++ )
        {
          bSuccess = bSuccess && mImageMap[i][j].Restore( oSerialBuf, bBuildSearchTree );
        }
      }
      return bSuccess; 
    }


    //---------------------------------------------
    //  BCast_Send
    //
    //  Note:  This is breaking the design rule - the only reason
    //         this is allowed is because the data structure
    //         is so large.
    //---------------------------------------------
    Bool BCast_Send() const
    {
      // send "header" information
      CSerializer Header;
      bool bSuccess =             Header.InsertCompactObj( nNumIntervals );
      bSuccess      = bSuccess && Header.InsertCompactObj( nNumDetectors );
      bSuccess      = bSuccess && Header.InsertCompactObj( nImageWidth   );
      bSuccess      = bSuccess && Header.InsertCompactObj( nImageHeight  );
      
      RUNTIME_ASSERT( bSuccess, "\nCatastrophic Error!  Serialization failed!  Unable to send via MPI\n" );
      char * pSendBuf = const_cast< char * >( Header.GetBuffer() );
      Size_Type nBufSize = Header.GetSize();
      MPI_Bcast( &nBufSize, sizeof( nBufSize ), MPI_CHAR, 0, MPI_COMM_WORLD ); 
      MPI_Bcast( pSendBuf, nBufSize, MPI_CHAR, 0, MPI_COMM_WORLD);
  
      // send actual data
      std::cerr << "Sending raw data ";
      
      int nImages = 0;
      for( Const_ImageIteratorT pCur = begin();
           pCur != end(); ++ pCur )
      {
        nImages ++;
        CSerializer DataSlice;
        Bool SaveSuccess = DataSlice.InsertCompactObj( pCur.Index2() );
        SaveSuccess      = SaveSuccess && DataSlice.InsertCompactObj( pCur.Index1() );
        SaveSuccess      = SaveSuccess && pCur->Save( DataSlice );
        
        RUNTIME_ASSERT( SaveSuccess, "Failed ot save detector image");
        pSendBuf = const_cast< char * >( DataSlice.GetBuffer() );
        nBufSize = DataSlice.GetSize();
        MPI_Bcast( &nBufSize, sizeof( nBufSize ), MPI_CHAR, 0, MPI_COMM_WORLD); 
        MPI_Bcast( pSendBuf, nBufSize, MPI_CHAR, 0, MPI_COMM_WORLD);
        std::cerr << " . ";
      }
      std::cerr << std::endl << " Number of Images Sent "  << nImages << std::endl;
      return true;
    }
    
    //---------------------------------------------
    //  BCast_Recv
    //
    //---------------------------------------------
    void BCast_Recv( bool bBuildSearchTree = true )
    {
      // recv "header" information
      //---------------------------------------
      Size_Type nBufSize;
      MPI_Bcast( &nBufSize, sizeof(nBufSize), MPI_CHAR, 0,  MPI_COMM_WORLD); 
      CDeserializer Header( nBufSize );
      char * pBuf = Header.GetBuffer();  
      MPI_Bcast( pBuf, nBufSize, MPI_CHAR, 0,  MPI_COMM_WORLD );

      bool bSuccess =             Header.GetCompactObj( &nNumIntervals );
      bSuccess      = bSuccess && Header.GetCompactObj( &nNumDetectors );
      bSuccess      = bSuccess && Header.GetCompactObj( &nImageWidth   );
      bSuccess      = bSuccess && Header.GetCompactObj( &nImageHeight  );
      mImageMap.resize( boost::extents[ nNumIntervals ][ nNumDetectors ] );

      // recv actual data
      for( int i = 0; i < nNumIntervals * nNumDetectors; i ++  )
      {
        MPI_Bcast( &nBufSize, sizeof(nBufSize), MPI_CHAR, 0,  MPI_COMM_WORLD); 
        CDeserializer DataSlice( nBufSize );
        pBuf = DataSlice.GetBuffer();
        MPI_Bcast( pBuf, nBufSize, MPI_CHAR, 0,  MPI_COMM_WORLD );

        ImageMapIndexT nDet, nInterval;
        Bool SaveSuccess =                DataSlice.GetCompactObj( &nDet );
        SaveSuccess      = SaveSuccess && DataSlice.GetCompactObj( &nInterval );
        SaveSuccess      = SaveSuccess && mImageMap[nInterval][nDet].Restore( DataSlice, bBuildSearchTree );
      }
    }
    
    
    Const_ImageIteratorT begin() const { return Const_ImageIteratorT( mImageMap, nNumIntervals, nNumDetectors, 0, 0 ); }
    Const_ImageIteratorT end  () const { return Const_ImageIteratorT( mImageMap, nNumIntervals, nNumDetectors, nNumIntervals, nNumDetectors ); }

    ImageIteratorT begin() { return ImageIteratorT( mImageMap, nNumIntervals, nNumDetectors, 0, 0 ); }
    ImageIteratorT end  () { return ImageIteratorT( mImageMap, nNumIntervals, nNumDetectors, nNumIntervals, nNumDetectors ); }
    
    
  };


  //---------------------------------------------------
  //
  //  Class Simulation Range
  //
  //  Provide an abstraction to 1) store the simulation angle
  //  information and 2) an interface to different kinds of look up.
  //  For example, non uniform look up will require interval search.
  //
  //  (Note that currently only support constant interval ranges)
  //---------------------------------------------------
  class CSimulationRange
  {
  public:
    typedef vector<Size_Type>::const_iterator IndexIter;
    
  private:
    
    //-------------------------------------
    //
    //  A range is generated for the interval ( fHigh, fLow ]
    //  with width being the width of individual intervals.
    //  In another words, there are a total of (fHigh - fLow) / fWidth
    //  intervals.
    //
    //  Note that one is to specify number of intervals, high, and low.
    //
    //-------------------------------------
    Float fLow;
    Float fHigh;          
    Float fWidth;
    Int   nNumIntervals;
    
    Int   nStartFileNum;     // starting number of the file
    Int   nStopFileNum;
  
    vector<SRange>    vRangeList; 
    vector<Size_Type> vIndexList;   // a mapping of [ l, h) -> unique wedge structures
                                    // This gives the index for a specific angle range on vRangeList
    
    inline Int AngleToIndex( Float fAngle ) const
    {
      Float f = ( fAngle - fLow ) / fWidth;   // this is needed to emulate floor()

      // below is a fast version of floor( f )
      if (  f < 0 )
        return -1;
      else
        return Int( f );                      // note that type demotion is a truncation,
                                              // and truncation of negative is the ceil() instead of floor()
    }
    
  public:
    
    //---------------------------------------------------------
    // Constructors
    //---------------------------------------------------------
    CSimulationRange()
      :fLow(12059), fHigh(-1940), fWidth(0),
       nNumIntervals(1020501), nStartFileNum(-2),
       nStopFileNum(-1093), vRangeList(), vIndexList() // deliberate gibberish
    {
      
    }

    //---------------------------------------------------------
    // General Constructors that maps the wedge structure
    //---------------------------------------------------------
    CSimulationRange( Float l, Float h,
                      Float _fWidth,
                      const vector<SRange> &vRange,
                      Int s = 0)
    {
      Set( l, h, _fWidth, vRange, s );
    }
    
    //---------------------------------------------------------
    //  Set
    //---------------------------------------------------------
    void Set( Float l, Float h,  Float _fWidth,
              const vector<SRange> &vRange, Int s = 0 )
    {
      fLow  = l;
      fHigh = h;
      fWidth = _fWidth;

      if( h < l )
        fWidth = -fWidth;
      
      nNumIntervals = round( ( h - l ) / fWidth );
      nStartFileNum = s;
      nStopFileNum  = s + nNumIntervals - 1;

      stringstream ss;
      ss << "ERROR!  Intervals do not seem to be equally spaced!\n"
         << " nNumIntervals  fWidth  h  l " << nNumIntervals << " " << fWidth << " " << h << " " << l << std::endl; 
      RUNTIME_ASSERT( fabs( nNumIntervals * fWidth - (h - l ) ) < 0.001,
                      ss.str() );
      vRangeList = vRange;
      InitializeIndexList( vRange );
    }
    
    //---------------------------------------------------------
    // InitializeIndexList  
    //  Note:  IndexList is basically a lookup table of the
    //         file index given a angular range.  Everyone is i
    //         initialized as NoMatch initially until proven
    //         otherwise
    //---------------------------------------------------------
    void InitializeIndexList( const vector<SRange> &vRange )
    {
      vIndexList.resize( nNumIntervals );

      for( Size_Type i = 0; i < vIndexList.size(); i++ )
        vIndexList[i] = NoMatch;
      
      for( Size_Type i = 0; i < vRangeList.size(); i ++ )
      {
        Float fCenter = ( vRangeList[i].fHigh + vRangeList[i].fLow ) / Float(2) ;
        vIndexList[ AngleToIndex( fCenter ) ] = i;
      }
    }
    
    //---------------------------------------------------------
    // maps fAngle -> index (which fits into a range)
    // note that fAngle is in radian
    // 
    // 
    //---------------------------------------------------------
    Size_Type ToFileNumber ( Float fAngle ) const
    {
      Int n = AngleToIndex( fAngle );

      if ( n > nStopFileNum || n < 0  )
        return NoMatch;
      else
        return n + nStartFileNum;
    }

    //---------------------------------------------------------
    //
    //   NEED input sanity check
    //
    //---------------------------------------------------------
    inline Float IndexToIntervalCenter( Int nIndex ) const
    {
      Float fAngle = nIndex * fWidth + fLow + Float( 0.5 ) * fWidth;
      return fAngle;
    }

    //---------------------------------------------------------
    //
    //   NEED input sanity check
    //
    //---------------------------------------------------------
    inline SRange IndexToInterval( Int nIndex ) const
    {
      DEBUG_ASSERT( vRangeList[nIndex].fLow < vRangeList[nIndex].fHigh,
                    "[CSimulationRange::IndexToInterval]:  Interval swapped!\n");
      return vRangeList[ nIndex ];
    }
    
    //---------------------------------------------------------
    // overloaded operator () as function object for indexing
    // purpose
    //
    // maps fAngle -> index (which fits into a range)
    // note that fAngle is in radian
    //---------------------------------------------------------
    Size_Type operator() ( Float fAngle ) const
    {
      Int n = AngleToIndex( fAngle );

      if ( n >= nNumIntervals || n < 0 )
        return NoMatch;
      else
        return vIndexList[ n ];
    }


    //---------------------------------------------------------
    //  checks to see if interval is valid
    //---------------------------------------------------------
    bool ValidIndex( Size_Type n ) const
    {
      return ( n < nNumIntervals && n >= 0 );
    }
    
    //---------------------------------------------------------
    //  GetIntervalWidth
    //---------------------------------------------------------
    Float GetIntervalWidth() const
    {
      return fWidth;
    }
    
    //---------------------------------------------------------
    // overloaded operator () as function object for indexing
    // purpose
    //
    // maps [ fAngleMin, fAngleMax ] -> [ IterBegin, IterEnd );
    // note that fAngle is in radian
    //
    //  WARNING:  This function currently only works with continuous ranges
    //  Note:     The list of omega returned may still contain NoMatch.
    //            All we are gauranteeing is that it is "check-able."
    //---------------------------------------------------------
    std::pair< IndexIter, IndexIter >
    operator() ( Float fAngleMin, Float fAngleMax ) const
    {
      Int nMin = AngleToIndex( fAngleMin );
      Int nMax = AngleToIndex( fAngleMax );

      // ensure ordering
      if( nMin > nMax )
        std::swap( nMin, nMax );
      
      // clamp to limit
      nMin = std::max( 0, nMin );
      nMin = std::min( nMin, nNumIntervals - 1 );

      nMax = std::max( 0, nMax );
      nMax = std::min( nMax, nNumIntervals - 1 );           
      IndexIter pFirst = vIndexList.begin()
                       + static_cast<Size_Type>( nMin ) ;

      if ( pFirst == vIndexList.end() )         // clamp to ends
        return std::make_pair( pFirst, pFirst );

      IndexIter pEnd = vIndexList.begin()
                     + static_cast<Size_Type>( nMax + 1 ) ;
      
      
      return std::make_pair( pFirst, pEnd  );
    }

    //---------------------------------------------------------
    //  Accessor:  GetnumIntervals
    //---------------------------------------------------------
    Int GetNumIntervals() const
    {
      return vRangeList.size();
    }
    
    //---------------------------------------------
    //  Save
    //---------------------------------------------
    Bool Save( CSerializer & oSerialBuf ) const
    {
      Bool bSuccess = true;
      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fLow);
      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fHigh );
      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fWidth );
      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nNumIntervals );

      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nStartFileNum );
      bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nStopFileNum );

      bSuccess = bSuccess && oSerialBuf.InsertCompactVector( vRangeList );
      bSuccess = bSuccess && oSerialBuf.InsertCompactVector( vIndexList );
      
      return bSuccess;
    }

    //---------------------------------------------
    //  Restore
    //---------------------------------------------
    Bool Restore( CDeserializer & oSerialBuf )
    {
      Bool bSuccess = true;
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fLow);
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fHigh );
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fWidth );
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nNumIntervals );
      
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nStartFileNum );
      bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nStopFileNum );
      
      bSuccess = bSuccess && oSerialBuf.GetCompactVector( vRangeList );
      bSuccess = bSuccess && oSerialBuf.GetCompactVector( vIndexList );
      
      return bSuccess;
      
    }

    
    
    
  };

  typedef SimulationData< int > CSimulationData;
}

#endif
