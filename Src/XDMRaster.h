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
//  Filename:  XDMRaster.h
//  Author:    Frankie Li
//
//
//  Purpose:   Specialization of Raster.h for the purpose of
//             XDM simulation application.
//
////////////////////////////////////////////////////////////


#ifndef _XDMRASTER_H_
#define _XDMRASTER_H_

#include "Raster.h"

#include <limits>
#include <vector>
#include <string>
#include <fstream>
#include "Types.h"
#include "Peak.h"
#include "Pixel.h"
#include "BBox.h"
#include "Error.h"
#include "SimpleGraph.h"
#include "Serializer.h"
#include <boost/pending/disjoint_sets.hpp>
#include <boost/function.hpp>

#include <boost/numeric/ublas/vector_of_vector.hpp>

#include "BitMatrix.h"

namespace GeneralLib
{


  //------------------------------------------------------------------------------------
  //
  //  class  CRaster
  //
  //  An abstraction of a raster to hide the detailed implementation, i.e., how pixels
  //  are stored and shapes rasterized.
  //
  //------------------------------------------------------------------------------------
  class CXDMRaster : public CRaster<IntensityT, boost::numeric::ublas::compressed_matrix<IntensityT> >   // using default compressed matrix
   //  class CXDMRaster : public CRaster<IntensityT, BitMatrix >   // using default compressed matrix
  {
  
  private:

    vector<Point> oTriangle;
    vector<Point> RectList;
    typedef CRaster<IntensityT,boost::numeric::ublas::compressed_matrix<IntensityT>  > Base;
    //typedef CRaster<IntensityT, BitMatrix  > Base;
    typedef IntensityT  BaseIntensityType;  
    struct SCompactCoordinate
    {
    
      typedef short int CoordType;
      CoordType nX;
      CoordType nY;
      BaseIntensityType fVal;
      SCompactCoordinate( Int _x, Int _y, BaseIntensityType _f): nX( _x ),
                                                                 nY( _y ),
                                                                 fVal( _f )
      {
        DEBUG_ASSERT( _x >= std::numeric_limits<CoordType>::min()
                      && _x <= std::numeric_limits<CoordType>::max()
                      && _y >= std::numeric_limits<CoordType>::min()
                      && _y <= std::numeric_limits<CoordType>::max()  );
      }
    
    };
  
    //------------------------------------------
    //  Specialized Template to count the number
    //  of pixels to be rastered on the detector
    //------------------------------------------
    struct SPixelRecordProcessor
    {
      Int nPixelOnDetector;
      Int nPixelOverlap;
      inline void operator() ( MatrixElementT_ConstRef p, IntensityT f)
      {
        nPixelOnDetector ++;
        if ( p > 0 )
          nPixelOverlap ++;
      }
    };

    //----------------------------
    //  Convert Pixel numbers to ID
    //----------------------------
    Int PixelToID( Int x, Int y ) const
    {
      DEBUG_ASSERT( IsInBound( x, y ), "Error:  PixelToID (x, y) out of bound\n" );
      return x * nRasterHeight + y;
    }
  
    //----------------------------
    //  Convert ID back to Pixel numbers
    //----------------------------
    std::pair<Int, Int> IDToPixel( Int nID ) const
    {
      DEBUG_ASSERT( nID < nRasterWidth * nRasterHeight,
                    "Error:  IDToPixel (nID) out of bound\n" );
      Int x = floor( nID / nRasterHeight );
      Int y = nID - (x * nRasterHeight);
      return std::make_pair( x, y );
    }

    //----------------------------
    // Special options for client-sever apps
    //
    // bServerMode is used with oRawPixelList.  When bServerMode is true,
    // the raw data read from UFF format will not be saved to the
    // sparse matrix in the server  This leads to a factor of 2 saving
    // in RAM usage.  
    //
    //----------------------------
    bool bSerializeMode;
    vector<SCompactCoordinate>  oRawPixelList;
    
  public:

    CXDMRaster() :
      Base(), oTriangle(3), RectList(4), bSerializeMode( false ), oRawPixelList() {}
    
    CXDMRaster( Int nWidth, Int nHeight ) :
      Base( nWidth, nHeight ),
      oTriangle(3), RectList(4), 
      bSerializeMode(false), oRawPixelList()  {}
    
    //------------------------------------------
    // GetPolygonOverlap
    // Find number of overlapping pixels
    //------------------------------------------
    Int GetPolygonOverlap( const vector<Point> & oVertices );

    //------------------------------------------
    //  GetTriangleOverlap
    //  Find number of overlapping pixels
    //------------------------------------------
    Int GetTriangleOverlap( const Point & v0, const Point & v1, const Point &v2 );

    //------------------------------------------
    //  GetTriangleOverlapProperty
    //------------------------------------------
    std::pair<Int, Int> GetTriangleOverlapProperty( const Point & v0,
                                                    const Point & v1,
                                                    const Point &v2 );
    std::pair<Int, Int>
    GetSquareOverlapProperty( const Point & v0, const Point & v1,
                              const Point & v2, const Point & v3  );
    //------------------------------------------
    //  SerializedAdd
    //------------------------------------------
    void SerializedAdd( Int nJPixel, Int nKPixel, Float fValue );
    void EnableSerializeMode();
    //------------------------------------------
    // Extract the list of connected pixels
    //------------------------------------------
    vector<CDetectorPeak> GetPeakList( bool bPurgePixels = true ) const;
    
    //------------------------------------------
    //  Save and Restore feature
    //------------------------------------------
    Bool Save   ( CSerializer & oSerialBuf ) const;
    Bool Restore( CDeserializer & oSerialBuf );


    //------------------------------------------
    //  For debugging only
    //------------------------------------------
    template< typename CustomRasterType >
    Bool Overwrite( const CustomRasterType & Image )
    {
      pImage.clear();
      for( int n = 0; n < Image.GetWidth(); n ++)
        for( int m = 0; m < Image.GetHeight(); m ++)
        {
          if( Image(n, m) > 0 )
            pImage( n, m ) = Image(n, m);
        }
      return true;
    }
    
  };

  

} // end Namespace

#endif
