////////////////////////////////////////////////////////////
//
//  Filename:  Raster.tmpl.cpp
//  Author:    Frankie Li
//
//
//  Purpose:   Templated  implementation of rastering for convex
//             polygons.
//
//
//  NOTE:      This is file was copied from the IceNine branch manually.
// - reverted
////////////////////////////////////////////////////////////


namespace GeneralLib
{
  //-------------------------------------------------------------------------------
  //
  //  CRaster
  //
  //  Pixels go from 0 -> nWidth -1 and 0 -> nHeight -1, but the setup of oClipper is not correct
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  CRaster< DataT, MatrixT >::CRaster( Int nWidth, Int nHeight ): oEdgeTable(),
                                                                 nRasterHeight(nHeight),
                                                                 nRasterWidth(nWidth)
                                                     
  {
    InitializeEdgeBuffer( nRasterHeight );
    InitializeImageBuffer( nWidth, nHeight );
  }
  //-------------------------------------------------------------------------------
  //  CRaster ( copy constructor )
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  CRaster< DataT, MatrixT >::CRaster( const CRaster< DataT, MatrixT > &  oRHS ): oEdgeTable(),
                                                                                 nRasterHeight( 0 ),
                                                                                 nRasterWidth( 0 )
                                                              
  {
    Copy( oRHS );
  }
  
  //-------------------------------------------------------------------------------
  //
  //  ~CRaster
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  CRaster< DataT, MatrixT >::~CRaster()
  {
    DeleteEdgeBuffer( nRasterHeight );
    DeleteImageBuffer();
  }

  //-------------------------------------------------------------------------------
  //
  //  DeleteImageBuffer
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::DeleteImageBuffer()
  {
    pImage.clear();
  }

  //-------------------------------------------------------------------------------
  //
  //  InitializeImageBuffer
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::InitializeImageBuffer( Int nWidth, Int nHeight )
  {
    if ( nWidth > 0 && nHeight > 0 )
      pImage.resize( nWidth, nHeight, false );
  }

  //-------------------------------------------------------------------------------
  //
  //   InitializeEdgeBuffer
  //   Allocate the memory required for oEdgeTable based on the hight
  //   of the detector.  Return false if the buffer is already allocated
  //   or if there is not enough memory.
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::InitializeEdgeBuffer( Int nRasterHeight )
  {
    DEBUG_ASSERT( nRasterHeight >= 0, " InitializeEdgeBuffer:  negative raster height detected!!" );
    if( nRasterHeight <= 0 )
    {
      return;
    }
    oEdgeTable = vector< vector<Int> >(  nRasterHeight, vector<Int>(2) );
  }

  //-------------------------------------------------------------------------------
  //
  //  Safely remove pEdgeTable and delete its memory
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  Bool CRaster< DataT, MatrixT >::DeleteEdgeBuffer( Int nRasterHeight )
  {
    if( nRasterHeight == 0 )
      return false;

    if( oEdgeTable.size() == 0 )
      return false;
    
    oEdgeTable.clear();
    return true;
  }

  //-------------------------------------------------------------------------------
  //
  //  ResetEdgeBuffer from [nYMin, nYMax], inclusively
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::ResetEdgeBuffer( Int nYMin, Int nYMax )
  {
    for( Int i = nYMin; i <= nYMax; i ++ )
    {
      oEdgeTable[i][0] = oEdgeTable[i][1] = -1; 
    }
  }
  //-------------------------------------------------------------------------------
  //  Copy, make *this the same as oRhs
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster<DataT, MatrixT>::Copy( const CRaster< DataT, MatrixT > & oRHS )  
  {
    this->DeleteEdgeBuffer( nRasterHeight );
    this->DeleteImageBuffer();
    this->InitializeEdgeBuffer( oRHS.nRasterHeight );
    this->InitializeImageBuffer( oRHS.nRasterWidth, oRHS.nRasterHeight );
    nRasterHeight = oRHS.nRasterHeight;
    nRasterWidth  = oRHS.nRasterWidth;
    // Copy data via sparse
    for ( SparseIterator1 iter1 = oRHS.pImage.begin1();
          iter1 != oRHS.pImage.end1(); ++ iter1 )
    {
      for ( SparseIterator2 iter2  = iter1.begin();
            iter2 != iter1.end(); ++ iter2)
      {
        pImage( iter2.index1(), iter2.index2() ) = * iter2;
      }
    }
  }
  //-------------------------------------------------------------------------------
  //  operator =
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  CRaster<DataT, MatrixT> & CRaster<DataT, MatrixT>::operator= ( const CRaster<DataT, MatrixT>  &oRHS )
  {
    Copy( oRHS );
    return *this;
  }
  
  //-------------------------------------------------------------------------------
  //  CalculateScanline
  //
  //  Save endpoints of the scanline to pEdgeTable
  //
  //  NOTE:  A single point will be saved as both start and end point in pEdgeTable
  //         if it is the first line drawn and recorded.  In another words, if pEdgeTable
  //         has just been reset, then all position of a sloped line will be recorded 
  //         as both the start and end points.
  //
  //  Precondition:  pEdgeTable initialized to the size of the screen.  oVertex0, oVertex1
  //  will not have pixels outside of pEdgeTable.  Vertex0 and Vertex1 are also within the
  //  screen.  Therefore recording the position of the lines that's drawn between v0 and v1
  //  will not result in segfault.
  //
  //  This implmentation is from wikipedia:
  //  http://en.wikipedia.org/wiki/Bresenhams_line_algorithm
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::CalculateScanline( const Pixel &v0, const Pixel & v1 )
  {

    Int v0x = v0.x;
    Int v0y = v0.y;
    Int v1x = v1.x;
    Int v1y = v1.y;
  
    bool bSteep = abs( v1y - v0y) > abs( v1x - v0x );

    if ( bSteep )
    {
      std::swap( v0x, v0y );
      std::swap( v1x, v1y );
    }

    if ( v0x > v1x )
    {
      std::swap( v0x, v1x );
      std::swap( v0y, v1y );
    }

    Int nDeltaX = v1x - v0x;
    Int nDeltaY = abs( v1y - v0y );
    Int nError = nDeltaX;
    Int nYStep;
    Int y = v0y;
  
    if ( v0y < v1y )
      nYStep = 1;
    else
      nYStep = -1;

    for ( Int x = v0x; x <= v1x; x ++ )
    {
      if ( bSteep )   // plot( y, x), which means oEdgeTable[x][ 0 or 1 ] will be assigned a value
      {
        if ( oEdgeTable[x][0] < 0 )
        {
          oEdgeTable[x][0] = y;
        }
        else if ( oEdgeTable[x][0] > y )   // find left extreme
        {
          oEdgeTable[x][0] = y;
        }

        if ( oEdgeTable[x][1] < y )   // find right extreme (True for uninitizalized (-1) or any value)
        {
          oEdgeTable[x][1] = y;
        }
      }
      else            // plot( y, x), which means oEdgeTable[y][ 0 or 1 ] will be assigned a value
      {
      
        if ( oEdgeTable[y][0] < 0 )
        {
          oEdgeTable[y][0] = x;
        }
        else if ( oEdgeTable[y][0] > x )   // find left extreme
        {
          oEdgeTable[y][0] = x;
        }

        if ( oEdgeTable[y][1] < x )   // find right extreme  (True for uninitizalized (-1) or any value)
        {                             // Note that this is checked even if x has been found to be the left extreme.
          oEdgeTable[y][1] = x;       // The reason is that the first pass requires the line to be both left and right
        }                             // extreme.
      }

      nError -=  2 * nDeltaY;
      if (nError < 0)
      {
        y = y + nYStep;
        nError += 2 * nDeltaX;
      }
    
    }
  }

  //-------------------------------------------------------------------------------
  //  Private: GeneralRasterizePolygon
  //
  //
  //  TODO:  Change this to use iterators instead of vector<Point>
  //-------------------------------------------------------------------------------
  //template <  typename InputParmT, typename OutputParmT >
  template < typename DataT, typename MatrixT >
  template< typename HelperFnT >
  void CRaster< DataT, MatrixT >::GeneralRasterizePolygon( const vector<Point> & oPolygon,
                                                           DataT fFillValue,
                                                           HelperFnT &ProcessorFn )
  {
    vector< Pixel > oVertexList = ClipPolygon( oPolygon );

    if ( oVertexList.size() <= 0 )
      return ;
  
    // find range of vertices
  
    Int nMaxY = 0;
    Int nMinY = MAX_INT;
  
    for ( Size_Type i = 0; i < oVertexList.size(); i ++)
    {
      if ( oVertexList[i].y > nMaxY )
        nMaxY = oVertexList[i].y;

      if ( oVertexList[i].y < nMinY )
        nMinY = oVertexList[i].y;
    }

    ResetEdgeBuffer( nMinY, nMaxY );
  
    // calculate scanlines
    for ( Size_Type i = 0; i < oVertexList.size(); i ++)
    {
      if ( i == 0 )   // wrap around
      {
        CalculateScanline( oVertexList[ oVertexList.size() - 1 ],  oVertexList[ i ] );
      }
      else
      {
        CalculateScanline( oVertexList[ i - 1 ], oVertexList[ i ] );
      }
    }

    // fill polygon
    for ( int y = nMinY; y <= nMaxY; y ++)
    {
      if ( oEdgeTable[y][0] >= 0 || oEdgeTable[y][1] >= 0 )   // at least one pixel is lit
      {
        if ( oEdgeTable[y][0] < 0 )                          // single pixel
        {
          ProcessorFn( pImage( oEdgeTable[y][1], y ), fFillValue );
        }
        else if ( oEdgeTable[y][1] < 0 )                     // single pixel
        {
          ProcessorFn( pImage( oEdgeTable[y][0] , y ), fFillValue );
        }
        else
        {
          if ( oEdgeTable[y][0] > oEdgeTable[y][1] )           // order reversed
            std::swap( oEdgeTable[y][0], oEdgeTable[y][1] );
        
          for ( int x = oEdgeTable[y][0]; x <= oEdgeTable[y][1]; x ++ )
          {
            ProcessorFn( pImage( x, y ), fFillValue );
          }
        }
      }
    }  
  
  }

  //-------------------------------------------------------------------------------
  //  RasterizePolygon
  //
  //  Rasterize, or fill the polygon specified by vertices.
  //
  //  Precondition:  The list of vertices, oPolygon, must be in winding order (right hand rule.)
  //
  //  Note that majority of this code is EXACTLY the same as GetPolygonOverlap, with the
  //  exception of the type of pixel processor used.
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::RasterizePolygon( const vector<Point> & oVertices, DataT fFillValue )
  {
    AccumulatorT oAccumulator;
    GeneralRasterizePolygon( oVertices, fFillValue, oAccumulator );
  }


  //-------------------------------------------------------------------------------
  //  RasterizePolygon
  //
  //  Rasterize, or fill the polygon specified by vertices.
  //
  //
  //  Maybe a specialized version if the Polygon version is not fast enough
  //
  //  Precondition:  The list of vertices, oVertices, must be in winding order (right hand rule.)
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::RasterizeTriangle( const Point & v0,
                                            const Point & v1, const Point &v2,
                                            DataT fFillValue )
  {
    vector<Point> oTriangle;

    oTriangle.resize(3);
    oTriangle[0] = v0 ;
    oTriangle[1] = v1 ;
    oTriangle[2] = v2 ;

    RasterizePolygon( oTriangle, fFillValue);
  }

  //-------------------------------------------------------------------------------
  //  Resize
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::Resize( Int nWidth, Int nHeight )
  {
    DeleteEdgeBuffer( nRasterHeight );
    DeleteImageBuffer();

    nRasterHeight = nHeight;
    nRasterWidth = nWidth;
    InitializeEdgeBuffer( nHeight );
    InitializeImageBuffer( nWidth, nHeight );
  
  }

  //-------------------------------------------------------------------------------
  //  Resize
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::Fill( DataT value )
  {
    for ( SparseIterator1 iter1 = pImage.begin1();
          iter1 != pImage.end1(); ++ iter1 )
    {
      for ( SparseIterator2 iter2  = iter1.begin();
            iter2 != iter1.end(); ++ iter2)
      {
        pImage( iter2.index1(), iter2.index2() ) = value;
      }
    }
  }

  
  //-------------------------------------------------------------------------------
  //  SetPixel
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::SetPixel( Int nI, Int nJ, DataT fValue )
  {
    pImage( nI, nJ ) = fValue;
  }

  //-------------------------------------------------------------------------------
  //  AddToPixel
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::AddToPixel( Int nI, Int nJ, DataT fValue )
  {
    pImage( nI, nJ ) += fValue;
  }
  //-------------------------------------------------------------------------------
  //
  //  ClearImage
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::ClearImage()
  {
    pImage.clear();
  }

  //-------------------------------------------------------------------------------
  //
  //  PrintRaster
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  void CRaster< DataT, MatrixT >::PrintRaster( string sFilename ) const
  {
    ofstream outfile;
    outfile.open( sFilename.c_str() );
    if(!outfile){
      RUNTIME_ASSERT(0, "[PrintRaster]: Cannot Open File\n");
      return;
    }
  
    for ( Int i = 0; i < nRasterWidth; i ++ )
      for ( Int j = 0; j < nRasterHeight; j ++ )
      {
      
        if ( pImage( i, j ) > 0 )
          outfile << i << ", " 
                  << j << ", "
                  << pImage( i, j )
                  << endl;
      
      }

    outfile.close();
  
  }

  //-------------------------------------------------------------------------------
  //
  //  GetHeight
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  Int CRaster< DataT, MatrixT >::GetHeight() const
  {
    return nRasterHeight;
  }

  //-------------------------------------------------------------------------------
  //
  //  GetWidth
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  Int CRaster< DataT, MatrixT >::GetWidth() const
  {
    return nRasterWidth;
  }


  //-------------------------------------------------------------------------------
  //  GetBegin1()
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  typename CRaster< DataT, MatrixT >::SparseIterator1 CRaster< DataT, MatrixT >::GetBegin1( ) const
  {
    return pImage.begin1();
  }

  //-------------------------------------------------------------------------------
  //  GetEnd1()
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  typename CRaster< DataT, MatrixT >::SparseIterator1 CRaster< DataT, MatrixT >::GetEnd1( ) const
  {
    return pImage.end1();
  }

  //-------------------------------------------------------------------------------
  //
  //  operator ()
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  typename MatrixT::const_reference
  CRaster< DataT, MatrixT >::operator() ( Int nI, Int nJ ) const
  {
    return pImage( nI, nJ );
  }
  
  //-------------------------------------------------------------------------------
  //
  //  operator ()
  //
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  typename MatrixT::reference
  CRaster< DataT, MatrixT >::operator() ( Int nI, Int nJ )
  {
    return pImage( nI, nJ );
  }
  //-------------------------------------------------------------------------------
  //  InBound
  //-------------------------------------------------------------------------------
  template < typename DataT, typename MatrixT >
  Bool CRaster< DataT, MatrixT >::IsInBound( Int nI, Int nJ ) const
  {
    return ( ( nI >= 0 ) && ( nJ >= 0 ) && ( nI < nRasterWidth ) && ( nJ < nRasterHeight ) );
  }



}// end of GeneralLib namespace
