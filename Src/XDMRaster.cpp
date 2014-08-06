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
///////////////////////////////////////////////////////////
//
//  Filename:  XDMRaster.cpp
//  Author:    Frankie Li
//
//
//  Purpose:  Implementation of a specialization of CRaster for
//            XDM simulations.
//
////////////////////////////////////////////////////////////


#include "XDMRaster.h"

namespace GeneralLib
{

//-------------------------------------------------------------------------------
//  GetPolygonOverlap
//
//  Count the number of pixels the inputted polygons overlaping with current raster
//
//  Precondition:  The list of vertices, oPolygon, must be in winding order (right hand rule.)
//
//  Note that majority of this code is EXACTLY the same as RasterizePolygon, with the
//  exception of the pixel processor used.
//
//-------------------------------------------------------------------------------
Int CXDMRaster::GetPolygonOverlap( const vector<Point> & oPolygon )
{
  OverlapCounterT oProcessor;
  
  oProcessor.nPixelOverlap = 0;
  GeneralRasterizePolygon( oPolygon, 0, oProcessor );

  return oProcessor.nPixelOverlap;
}

//-------------------------------------------------------------------------------
//  GetTriangleOverlap
//
//  Find the overlap between the triangle from the input with the current raster.
//
//  Precondition:  The list of vertices, oVertices, must be in winding order (right hand rule.)
//-------------------------------------------------------------------------------
Int CXDMRaster::GetTriangleOverlap( const Point & v0, const Point & v1, const Point &v2 )
{
  //  oTriangle.resize(3);
  oTriangle[0] = v0 ;
  oTriangle[1] = v1 ;
  oTriangle[2] = v2 ;

  return GetPolygonOverlap( oTriangle );
}

//-------------------------------------------------------------------------------
//  GetTriangleOverlapProperty
//
//  
//  Precondition:  The list of vertices, oVertices, must be in winding order (right hand rule.)
//-------------------------------------------------------------------------------
std::pair<Int, Int>
CXDMRaster::GetTriangleOverlapProperty( const Point & v0, const Point & v1, const Point &v2 )
{
  SPixelRecordProcessor oProcessor;
  oProcessor.nPixelOnDetector = 0;
  oProcessor.nPixelOverlap = 0 ;

  //  vector<Point> oTriangle( 3 );
  // oTriangle.resize( 3 );
  oTriangle[0] = v0 ;
  oTriangle[1] = v1 ;
  oTriangle[2] = v2 ;
  
  GeneralRasterizePolygon( oTriangle, 0, oProcessor );
  std::pair<Int, Int> oRes = std::make_pair( oProcessor.nPixelOverlap,
                                             oProcessor.nPixelOnDetector );
  return oRes;
}

//-------------------------------------------------------------------------------
//  GetSquareOverlapProperty
//-------------------------------------------------------------------------------
std::pair<Int, Int>
CXDMRaster::GetSquareOverlapProperty( const Point & v0, const Point & v1,
                                      const Point & v2, const Point & v3  )
{
  SPixelRecordProcessor oProcessor;
  oProcessor.nPixelOnDetector = 0;
  oProcessor.nPixelOverlap = 0 ;


  RectList[0] = v0 ;
  RectList[1] = v1 ;
  RectList[2] = v2 ;
  RectList[3] = v3 ;
  
  GeneralRasterizePolygon( RectList, 0, oProcessor );
  std::pair<Int, Int> oRes = std::make_pair( oProcessor.nPixelOverlap,
                                             oProcessor.nPixelOnDetector );
  return oRes;
}

//----------------------------------------------------------------------------------
//
// Public: ExtractPeak
//
// Note:  This isn't the most efficient want to do things, but it really doesn't
//        matter all that much since this function will only be called exactly once.
// Note:  This function exists in the Raster level because the convention of X and Y
//        pixel is specified by CXDMRaster.
//
// HACK!!!  The pixel list is emptied at the end because it isn't being used.  This is
//          done to reduce memory consumption significantly.
//----------------------------------------------------------------------------------
vector<CDetectorPeak> CXDMRaster::GetPeakList( bool bPurgePixels ) const
{
  // using a the ID of each pixel as the element ID
  boost::disjoint_sets_with_storage<> oDSet( nRasterHeight * nRasterWidth );
  for ( SparseIterator1 iter1 = pImage.begin1();
        iter1 != pImage.end1(); ++ iter1 )
    for ( SparseIterator2 iter2  = iter1.begin();
          iter2 != iter1.end(); ++ iter2)
      if( *iter2 > 0 )
        oDSet.make_set( PixelToID( iter2.index1(), iter2.index2() ) );
  
  //  Note that iter1, iter2 goes through the dimensions
  for ( SparseIterator1 iter1 = pImage.begin1();
        iter1 != pImage.end1(); ++ iter1 )
  {
    for ( SparseIterator2 iter2  = iter1.begin();
          iter2 != iter1.end(); ++ iter2)
    {
      Int nX      = iter2.index1();
      Int nY      = iter2.index2();
      Int nCenter = PixelToID( nX, nY );
      for ( Int nDx = 0; nDx <= 1; nDx ++ )
      {
        for ( Int nDy = 0; nDy <= 1; nDy ++ )
        {
          if ( nDx != 0 || nDy != 0 )
          {
            Int nNewX = nX + nDx;
            Int nNewY = nY + nDy;
            if ( IsInBound( nNewX, nNewY  ) && pImage( nNewX, nNewY ) > 0 )
            {
              oDSet.link( nCenter, PixelToID( nNewX, nNewY ) );
            }
          }
        } 
      }
    }
  }

  vector< CDetectorPeak > oPeakList;
  std::map< Int, Int> oIDToIndexMap;
  Int nCompIndex = 0;
  for ( SparseIterator1 iter1 = pImage.begin1();
        iter1 != pImage.end1(); ++ iter1 )
    for ( SparseIterator2 iter2  = iter1.begin();
          iter2 != iter1.end(); ++ iter2)
    {
      if( *iter2 > 0 )
      {
        Int nRep = oDSet.find_set( PixelToID( iter2.index1(), iter2.index2() ) );
        std::map< Int, Int >::iterator pFound = oIDToIndexMap.find( nRep );
        
        Pixel p;
        p.x          = iter2.index1();
        p.y          = iter2.index2();
        p.fIntensity = *iter2;
        
        if( pFound  == oIDToIndexMap.end() )
        {
          oIDToIndexMap[ nRep ] = nCompIndex;
          oPeakList.push_back( CDetectorPeak() );
          oPeakList[ nCompIndex ].vPixelList.push_back( p );
          nCompIndex ++;
        }
        else
        {
          oPeakList[ pFound->second ].vPixelList.push_back( p );
        }
      }
    }
  
  for( Size_Type i = 0; i < oPeakList.size(); i ++ )
  {
    oPeakList[i].CalculateBoundingBox();
    if( bPurgePixels )
      oPeakList[i].vPixelList.clear();    // HACK -- removing the pixel list to reduce memory consumption!!!!
  }
  return oPeakList;
}

//------------------------------------------
//  SerializedAdd
//------------------------------------------
void CXDMRaster::SerializedAdd( Int nJPixel, Int nKPixel, Float fValue )
{
  RUNTIME_ASSERT( bSerializeMode,
                  "[CXDMRaster]: SerializeAdd should not be used when not in serialization mode\n");
  oRawPixelList.push_back( SCompactCoordinate( nJPixel, nKPixel, fValue ) );
}
//------------------------------------------
// EnableSerializeMode
//------------------------------------------
void CXDMRaster::EnableSerializeMode()
{
  bSerializeMode = true;
}
//------------------------------------------
//  Save
//------------------------------------------
Bool CXDMRaster::Save   ( CSerializer & oSerialBuf ) const
{
  Bool bSuccess;
  vector<SCompactCoordinate> vElementList;
  bSuccess = oSerialBuf.InsertCompactObj( nRasterWidth  );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nRasterHeight );

  if( bSerializeMode )
  {
    bSuccess = bSuccess && oSerialBuf.InsertCompactVector( oRawPixelList );
  }
  else
  {
    for ( SparseIterator1 iter1 = pImage.begin1();
          iter1 != pImage.end1(); ++ iter1 )
      for ( SparseIterator2 iter2  = iter1.begin();
            iter2 != iter1.end(); ++ iter2)
        vElementList.push_back( SCompactCoordinate( iter2.index1(),
                                                    iter2.index2(),
                                                    *iter2 ) );
    bSuccess = bSuccess && oSerialBuf.InsertCompactVector( vElementList );
  }
  return bSuccess;
}

//------------------------------------------
//  Restore
//------------------------------------------
Bool CXDMRaster::Restore( CDeserializer & oSerialBuf )
{
  Int nWidth, nHeight;
  Bool bSuccess;
  bSuccess = oSerialBuf.GetCompactObj( &nWidth  );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nHeight );
  Base::Resize( nWidth, nHeight );
  vector<SCompactCoordinate> vElementList;
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vElementList );

  RUNTIME_ASSERT( bSuccess, "Restore of raster unsuccessful!\n" );
  for( Size_Type i = 0; i < vElementList.size(); i ++ )
  {
    Int nX = vElementList[i].nX;
    Int nY = vElementList[i].nY;
    pImage( nX, nY ) = vElementList[i].fVal;
  }

  return bSuccess;
}
  

}// end namespace
