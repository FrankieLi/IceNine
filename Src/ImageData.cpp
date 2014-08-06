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
//
//   ImageData.cpp
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:  Implementation of ImageData.h
//
////////////////////////////////////////////////////////////

#include "ImageData.h"

//----------------------------------------------------------------------------------
//  CImageData
//----------------------------------------------------------------------------------
CImageData::CImageData() : oImage( 0, 0 )
{
}

//----------------------------------------------------------------------------------
//  CImageData( nNumJPixels, nNumKPixels )
//----------------------------------------------------------------------------------
CImageData::CImageData( Int nWidth, Int nHeight ) : oImage( nWidth, nHeight )
{
}

//----------------------------------------------------------------------------------
//  CImageData( const CImageData & ), copy constructor
//----------------------------------------------------------------------------------
CImageData::CImageData( const CImageData & oRHS ) : oImage(  oRHS.oImage )
{
  
}

//----------------------------------------------------------------------------------
//  operator =
//----------------------------------------------------------------------------------
CImageData & CImageData::operator=( const CImageData & oRHS )
{
  oImage = oRHS.oImage;
  return *this;
}
//----------------------------------------------------------------------------------
//  ~CImageData( )
//----------------------------------------------------------------------------------
CImageData::~CImageData()
{
}

//----------------------------------------------------------------------------------
//  Resize
//----------------------------------------------------------------------------------
void CImageData::Resize( Int nWidth, Int nHeight )
{
  oImage.Resize( nWidth, nHeight );
}

//----------------------------------------------------------------------------------
//   IsDark
//
//   Return true of the pixel at ( nJpixel, nKPixel ) is dark
//----------------------------------------------------------------------------------
Bool CImageData::IsDark( Int nJPixel, Int nKPixel ) const
{
  DEBUG_ASSERT( IsInBound( nJPixel, nKPixel ), " Failed at IsDark\n" );
  return ( oImage( nJPixel, nKPixel ) <= 0 );
}

//----------------------------------------------------------------------------------
//  GetImageBBox
//----------------------------------------------------------------------------------
BBox2D CImageData::GetImageBBox( ) const
{
  BBox2D oRes( Point(0, 0), Point( oImage.GetWidth() -1, oImage.GetHeight() -1 ) );
  return oRes;
}

//----------------------------------------------------------------------------------
//   IsBright
//
//   Return true of the pixel at ( nJpixel, nKPixel ) is bright
//----------------------------------------------------------------------------------
Bool CImageData::IsBright( Int nJPixel, Int nKPixel ) const
{
  DEBUG_ASSERT( IsInBound( nJPixel, nKPixel), " Failed at IsBright\n" );
  return ( oImage( nJPixel, nKPixel ) > 0 ) ;
}

//----------------------------------------------------------------------------------
//   IsInBound
//
//   Return true of the pixel at ( nJpixel, nKPixel ) is in bound
//----------------------------------------------------------------------------------
Bool CImageData::IsInBound( Int nJPixel, Int nKPixel ) const
{
  return oImage.IsInBound( nJPixel, nKPixel );
}

//----------------------------------------------------------------------------------
//  SetPixel:   Set the value of the specified pixel to fValue
//----------------------------------------------------------------------------------
void CImageData::SetPixel( Int nJPixel, Int nKPixel, Float fValue )
{
  DEBUG_ASSERT( IsInBound( nJPixel, nKPixel ), " Failed at SetPixel\n" );
  oImage.SetPixel( nJPixel, nKPixel, fValue );
}

//----------------------------------------------------------------------------------
//  AddToImage:  Add pixel to image
//----------------------------------------------------------------------------------
void CImageData::AddToPixel( Int nJPixel, Int nKPixel, Float fValue )
{
  DEBUG_ASSERT( IsInBound( nJPixel, nKPixel ), " Failed at AddToPixel\n" );
  oImage.AddToPixel( nJPixel, nKPixel, fValue );
}

//----------------------------------------------------------------------------------
//   Add Triangle
//----------------------------------------------------------------------------------
void CImageData::AddTriangle( const Point & v0, const Point & v1,
                              const Point & v2, Float fIntensity )
{
  oImage.RasterizeTriangle( v0, v1, v2, fIntensity );   
}

//----------------------------------------------------------------------------------
//   Add Polygon
//----------------------------------------------------------------------------------
void CImageData::AddPolygon( const vector<Point> & oVertices, Float fFillValue )
{
  oImage.RasterizePolygon( oVertices, fFillValue );
}

//----------------------------------------------------------------------------------
//   GetTriangleOverlap
//----------------------------------------------------------------------------------
Int CImageData::GetTriangleOverlap( const Point &v0, const Point & v1, const Point & v2 )
{
  return oImage.GetTriangleOverlap( v0, v1, v2 );
}

//----------------------------------------------------------------------------------
//  GetTriangleOverlapProperty
//----------------------------------------------------------------------------------
std::pair<Int, Int>
CImageData::GetTriangleOverlapProperty( const Point &v0, const Point & v1, const Point & v2 )
{
  return oImage.GetTriangleOverlapProperty( v0, v1, v2 );
}

//----------------------------------------------------------------------------------
//  GetTriangleOverlapProperty
//----------------------------------------------------------------------------------
std::pair<Int, Int>
CImageData::GetSquareOverlapProperty( const Point &v0, const Point & v1,
                                      const Point & v2, const Point & v3 )
{
  return oImage.GetSquareOverlapProperty( v0, v1, v2, v3 );
}

//----------------------------------------------------------------------------------
//  GetNumPixelsLit
//
//  Purpose:  count the number of pixels that the triangle specified by
//            the three verticies lit on the detector.
//----------------------------------------------------------------------------------
Int CImageData::GetNumPixelsLit( const Point &v0, const Point & v1, const Point & v2 )
{
  CXDMRaster::HitCounterT oPixelCounter;
  oPixelCounter.nPixelOnDetector = 0;
  vector<Point> oPolygon(3);
  oPolygon[0] = v0;
  oPolygon[1] = v1;
  oPolygon[2] = v2;

  Float fPlaceHolder = 0;
  oImage.GeneralRasterizePolygon( oPolygon, fPlaceHolder, oPixelCounter );
  
  return oPixelCounter.nPixelOnDetector;
}

//----------------------------------------------------------------------------------
//   ClearImage
//----------------------------------------------------------------------------------
void CImageData::ClearImage()
{
  oImage.ClearImage();
}

//--------------------------
// Save and Restore
//--------------------------
Bool CImageData::Save   ( CSerializer & oSerialBuf ) const
{
  return oImage.Save( oSerialBuf );
}

//--------------------------
//  Restore
//--------------------------
Bool CImageData::Restore( CDeserializer & oSerialBuf )
{
  return oImage.Restore( oSerialBuf );
}

//----------------------------------------------------------------------------------
//
//
//  F I L E   I / O
//
//
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
//   ReadCXDMSimulationPeakFile
//----------------------------------------------------------------------------------
bool CImageData::ReadCXDMSimulationDataFile( const string &filename )
{
  char *pBuffer = NULL;
  Size_Type nBufSize = 0;
  vector< vector<string> > vsTokens;
	
  pBuffer = InitFileIO::ReadFileToBuf(nBufSize, filename);

  if(pBuffer == NULL)
  {
    // test to see if file exist, but is empty
    ifstream oInfile;
    oInfile.open( filename.c_str() );
    Bool bFail = oInfile.fail();
    oInfile.close();
    if ( bFail )
    {
      cerr << "[CImageData::ReadCXDMSimulationDataFile] Error:  File: " << filename << "  cannot be opened" << std::endl;
      return false;
    }
    else
    {
      return true;
    }
  }
  
  GeneralLib::Tokenize( vsTokens, string( pBuffer, nBufSize ), ",# \t\n");
    
  for( Size_Type i = 0; i < vsTokens.size(); i ++)
  {
    Int nX, nY;
    Float fIntensity;

    if ( vsTokens[i].size() == 3 )
    {
      nX = atoi( vsTokens[i][0].c_str() );
      nY = atoi( vsTokens[i][1].c_str() );
      fIntensity = atof( vsTokens[i][2].c_str() );
			
      RUNTIME_ASSERT( nX >= 0 && nX < oImage.GetWidth() && nY >= 0 && nY < oImage.GetHeight(),
                      " [ReadCXDMSimulationDataFile] ERROR:  Index out of bounds!\n");
      SetPixel( nX, nY, fIntensity);
    }
    else if( vsTokens[i].size() != 0 )   // 0 is empty line, ignore
    {
      cerr << "[CImageData::ReadCXDMSimulationDataFile] Error: Incorrect number of columns, line: "
           << i  << " in file: " << filename << endl;
      return false;
    }
  }

  delete[] pBuffer;
  
  return true;
}

//----------------------------------------------------------------------------------
//   WritePixelsToASCII
//----------------------------------------------------------------------------------
bool CImageData::WritePixelsToASCII( const string & filename )
{
  oImage.PrintRaster( filename );
  return true;
}

//----------------------------------------------------------------------------------
//
//
//  Class   CSearchableImageData
//
//
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
// BuildSearhTree
//----------------------------------------------------------------------------------
void CSearchableImageData::BuildSearchTree( )
{
  //  std::cout << "Building Search Tree " << std::endl;
  for( Size_Type i = 0; i < vDetectorPeakList.size(); i ++ )
  {
    BBox2D oBBox  = vDetectorPeakList[i].GetBoundingBox( );
    RangeSearch::CPeakPtr pPeak( &vDetectorPeakList[i], Null_Deleter() );
    oPeakTree.Add( pPeak, oBBox );
  }
}

//-------------------------------------------------------------------------
//
//  Override the read
//
//-------------------------------------------------------------------------
bool CSearchableImageData::ReadCXDMSimulationDataFile( const string &sFilename )
{
  bool bSuccess  = CImageData::ReadCXDMSimulationDataFile( sFilename );
  if ( !bSuccess )
    return false;
  InitializePeakList();
  return true;
}

//----------------------------------------------------------------------------------
//  ReadCXDMSimulationUFFFile - binary format
//----------------------------------------------------------------------------------
bool CSearchableImageData::ReadCXDMSimulationUFFFile( const string &filename,
                                                      bool bServerMode )
{
  CInputPeakFile oInputFile;
  bool bSuccess = oInputFile.Read( filename );

  if( !bSuccess )
    return false;

  if( bServerMode )
    oImage.EnableSerializeMode();
  
  if( oInputFile.NumPixels() <= 0 )   // no pixel -- it happens
    return true;
  
  // read to raster first
  Pixel oTmpPixel;
  Int   nPeakID;
  std::map<Int, Int> nIDMap;    // remap the ID into a continuous strip
  Int   nCurID = 0;
  
  vDetectorPeakList.clear();
  vDetectorPeakList.resize( oInputFile.NumPeaks() );   // nCurID = number of peaks
  
  vector< Pixel > oPixelList;                          // buffer pixels
  oPixelList.reserve( oInputFile.NumPixels() );
  while( oInputFile.GetNextPixel( oTmpPixel, nPeakID ) )
  {
    std::map<Int, Int>::iterator pLoc = nIDMap.find( nPeakID );
    if( pLoc == nIDMap.end() )
    {
      nIDMap[ nPeakID ] = nCurID;
      nCurID ++;
    }
    RUNTIME_ASSERT( oTmpPixel.x >= 0 && oTmpPixel.x < oImage.GetWidth()
                    && oTmpPixel.y >= 0 && oTmpPixel.y < oImage.GetHeight(),
                    " [ReadCXDMSimulationDataFile] ERROR:  Index out of bounds!\n");

    oPixelList.push_back( oTmpPixel );
    vDetectorPeakList[ nIDMap[ nPeakID ] ].AddPointToBBox( oTmpPixel );
    //vDetectorPeakList[ nIDMap[ nPeakID ] ].vPixelList.push_back( oTmpPixel );   // used for calculating bbox, may be changed later
  }
  
  std::sort( oPixelList.begin(), oPixelList.end(),
             SPixelCoordCmp( oImage.GetWidth(), oImage.GetHeight() ) );         // Sort for fast insertion into compressed matrix
  
  for( Size_Type i = 0; i < oPixelList.size(); i ++ )
  {
    if( bServerMode )
      oImage.SerializedAdd( oPixelList[i].x, oPixelList[i].y, oPixelList[i].fIntensity );
    else
      SetPixel( oPixelList[i].x, oPixelList[i].y, oPixelList[i].fIntensity );
  }
  for( Size_Type i = 0; i < vDetectorPeakList.size(); i ++ )  // WARNING:  pixels removed after bbox calculation 
  {
    vDetectorPeakList[ i ].CalculateBoundingBox();
    vDetectorPeakList[ i ].vPixelList.clear();                // since they are not being used now (will be needed later)
  }
  if( !bServerMode )
    BuildSearchTree();
  return true;
}

//----------------------------------------------------------------------------------
//
//   GetNumOverlapPeaks
//
//
//----------------------------------------------------------------------------------
Int  CSearchableImageData::GetNumOverlapPeaks( const Point &v0, const Point & v1, const Point & v2 ) const
{
  
  BBox2D oTestRegion;
  
  oTestRegion = Union( oTestRegion, v0 );
  oTestRegion = Union( oTestRegion, v1 );
  oTestRegion = Union( oTestRegion, v2 );

  return oPeakTree.NumPeakOverlaps( oTestRegion );
}

//----------------------------------------------------------------------------------
//  HasOverlapPeak
//
//  Return true if there exist one peak that overlaps the bounding
//  box oBox
//----------------------------------------------------------------------------------
Bool CSearchableImageData::HasOverlapPeak( const BBox2D & oTestBox ) const
{
  return oPeakTree.Overlaps( oTestBox );
}

//----------------------------------------------------------------------------------
//   InitializePeakList
//----------------------------------------------------------------------------------
void CSearchableImageData::InitializePeakList()
{
  vDetectorPeakList.clear();
  vDetectorPeakList = oImage.GetPeakList();
  BuildSearchTree();
}

//----------------------------------------------------------------------------------
//
//   GetNumOverlapPeaks
//
//
//----------------------------------------------------------------------------------
Int  CSearchableImageData::GetNumOverlapPeaks( const vector<Point> & vPoints ) const
{
  BBox2D oTestRegion;
  for( Size_Type i = 0; i < vPoints.size(); i ++ )
    oTestRegion = Union( oTestRegion, vPoints[i] );

  return  oPeakTree.NumPeakOverlaps( oTestRegion );
}

//----------------------------------------------------------------------------------
//   CheckOverlap  --
//
//   Create bounding box out of the points given in vPoints.  Calculate overlap
//   between bounding boxes.
//
//----------------------------------------------------------------------------------
Bool CSearchableImageData::HasOverlap( const vector<Point> & vPoints ) const
{
  BBox2D oTestRegion;
  for( Size_Type i = 0; i < vPoints.size(); i ++ )
  {
    oTestRegion = Union( oTestRegion, vPoints[i] );   // This could be optimized to use O( log(vPoints.size() ) )comparison
  }

  return oPeakTree.Overlaps( oTestRegion );
}

//----------------------------------------------------------------------------------
//   CheckOverlap  --
//
//   TODO: evolve this function into returning area of overlap?
//
//----------------------------------------------------------------------------------
Bool CSearchableImageData::HasOverlap( const Point &v0, const Point & v1, const Point & v2 ) const
{
  Point pMax;
  Point pMin;
  
  pMax.x = std::max( v0.x, v1.x);
  pMax.x = std::max( v2.x, pMax.x);
  
  pMax.y = std::max( v0.y, v1.y);
  pMax.y = std::max( v2.y, pMax.y);
    
  pMin.x = std::min( v0.x, v1.x);
  pMin.x = std::min( v2.x, pMin.x);
  
  pMin.y = std::min( v0.y, v1.y);
  pMin.y = std::min( v2.y, pMin.y);
  
  BBox2D oVoxelBBox( pMin, pMax );
  

  return oPeakTree.Overlaps( oVoxelBBox );
}

//----------------------------------------------------------------------------------
//   GetNumPeaksOverlap  --
//
//
//----------------------------------------------------------------------------------
Int CSearchableImageData::GetNumPeaksOverlap( const Point &v0, const Point & v1, const Point & v2 ) const
{
  Int nPeakOverlap = 0;

  for ( Size_Type i = 0; i < vDetectorPeakList.size(); i ++ )
  {
    BBox2D oBBox  =  vDetectorPeakList[i].GetBoundingBox( );
    Bool bOverlap = Geometry2D::TriangleAABBOverlapEstimate( oBBox, v0, v1, v2 );
    
    if ( bOverlap )
      nPeakOverlap ++;
  }
  
  return nPeakOverlap;
}

//--------------------------
// Save
//--------------------------
Bool CSearchableImageData::Save   ( CSerializer & oSerialBuf ) const
{
  Bool bSuccess;
  bSuccess = CImageData::Save( oSerialBuf );
  Size_Type nVectorSize = vDetectorPeakList.size();
  bSuccess = oSerialBuf.InsertCompactObj( nVectorSize );
  for( Size_Type i = 0; i < nVectorSize; i ++ )
  {
    bSuccess = bSuccess && vDetectorPeakList[i].Save( oSerialBuf );
  }
  return bSuccess;
}

//--------------------------
//  Restore
//--------------------------
Bool CSearchableImageData::Restore( CDeserializer & oSerialBuf,
                                    bool bBuildSearchTree )
{
  Bool bSuccess;
  bSuccess = CImageData::Restore( oSerialBuf );

  Size_Type nVectorSize;
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nVectorSize );
  vDetectorPeakList.clear();
  vDetectorPeakList.resize( nVectorSize );

  for ( Size_Type i = 0; i < nVectorSize; i ++ )
  {
    bSuccess = bSuccess && vDetectorPeakList[i].Restore( oSerialBuf );
  }

  BuildSearchTree();      // construct secondary structure for quick searching
  
  DEBUG_ASSERT( bSuccess, " Restore failed at CSearchableImageData::Restore \n" );
  return bSuccess;
}

//-------------------------------------------------------------------------------------
//
//  2D geometry namespace for simple things like triangle bounding box intersection,
//  separation axis test  (TODO:  move to another file in the future?) 
//
//-------------------------------------------------------------------------------------
namespace Geometry2D
{
  //-------------------------------------------------------------------------------------
  //
  //  CheckSeparatingAxis
  //
  //  Verticse of a rectangle is given as four vertices. 
  //
  //  Given separating axis specified by vertices v1 and v2 containing only the X and Y
  //  coordinates,  (Note that this is 2D)
  //
  //  Return true if a separating axis exists
  //
  //-------------------------------------------------------------------------------------
  bool CheckSeparatingAxis( const BBox2D & oBox, const SVoxel &oVoxel,
                            const SVector3 & v1, const SVector3 & v2 )
  {

    SVector3 oDir = v1 - v2;
    //   oDir.Normalize(); // don't have to normalize when every vertex is being dotted with the same vector
    
    // rotate 90 degress y <- x, x <- -y
    Float tmp = oDir.m_fY;
    oDir.m_fY = oDir.m_fX;
    oDir.m_fX = -tmp;

    Float fMin = MAX_FLOAT;
    Float fMax = MIN_FLOAT;
    
    Float fAxisProjection;

    //---------------------
    // All 4 vertices of the bounding box
    //
    //  1 ----- 2
    //  |       |
    //  |       |
    //  4 ----- 3   where pMin = vertex4,  pMax = vertex2
    //
    //  vertex1 = ( pMin.x, pMax.y )
    //  vertex2 = ( pMax.x, pMax.y )
    //  vertex3 = ( pMax.x, pMin.y )
    //  vertex4 = ( pMin.x, pMin.y )
    //---------------------

    // vertex 1
    fAxisProjection =  oBox.pMin.x * oDir.m_fX + oBox.pMax.y * oDir.m_fY ;
    fMin = std::min( fAxisProjection,  fMin );
    fMax = std::max( fAxisProjection,  fMax );

    // vertex 2
    fAxisProjection =  oBox.pMax.x * oDir.m_fX + oBox.pMax.y * oDir.m_fY ;
    fMin = std::min( fAxisProjection,  fMin );
    fMax = std::max( fAxisProjection,  fMax );

    // vertex 3
    fAxisProjection =  oBox.pMax.x * oDir.m_fX + oBox.pMin.y * oDir.m_fY ;
    fMin = std::min( fAxisProjection,  fMin );
    fMax = std::max( fAxisProjection,  fMax );

    // vertex 4
    fAxisProjection =  oBox.pMin.x * oDir.m_fX + oBox.pMin.y * oDir.m_fY ;
    fMin = std::min( fAxisProjection,  fMin );
    fMax = std::max( fAxisProjection,  fMax );


    // iterate through vertices of voxel
    for ( int i = 0; i < 3; i ++ )
    {
      Float fVertex = Dot( oVoxel.pVertex[i], oDir );  
      if( fVertex <= fMax  && fVertex >= fMin )   // if intersecting
      {
        return false;
      }
    }
    return true;
  }
  
  //-------------------------------------------------------------------------------------
  //
  //  VoxelAABBOverlap()
  //
  //  Return:    True     if the AABB (axis aligned bounding box) overlaps with the voxel
  //             False    otherwise  (This condition can be met if one separating axis is found)
  //
  //  NOTE:  This function only works in 2D.  Specifically, it only works with the x and y
  //         coordinates.
  //
  //
  //-------------------------------------------------------------------------------------
  bool VoxelAABBOverlap( const BBox2D &oBox, const SVoxel &oVoxel )
  {
    Float fVoxelMinX = MAX_FLOAT;
    Float fVoxelMinY = MAX_FLOAT;
    Float fVoxelMaxX = MIN_FLOAT;
    Float fVoxelMaxY = MIN_FLOAT;
    
    // find separating axis for the bounding box manuall.y
    for ( Int i = 0; i < 3; i ++ )
    {
      fVoxelMinX = std::min( oVoxel.pVertex[i].m_fX, fVoxelMinX );
      fVoxelMinY = std::min( oVoxel.pVertex[i].m_fY, fVoxelMinY );
      fVoxelMaxX = std::max( oVoxel.pVertex[i].m_fX, fVoxelMaxX );
      fVoxelMaxY = std::max( oVoxel.pVertex[i].m_fY, fVoxelMaxY );
    }

    bool bXAxisNoOverlap = ( fVoxelMinX > oBox.pMax.x ) || ( oBox.pMin.x > fVoxelMaxX );
    bool bYAxisNoOverlap = ( fVoxelMinY > oBox.pMax.y ) || ( oBox.pMin.y > fVoxelMaxY );

    if ( bXAxisNoOverlap || bYAxisNoOverlap )
      return false;

    // else, both axis from the bounding box overlap.  check to see if axis from
    // voxel is overlaping as well.

    bool bV01AxisExists = CheckSeparatingAxis( oBox, oVoxel, oVoxel.pVertex[0], oVoxel.pVertex[1] );
    bool bV12AxisExists = CheckSeparatingAxis( oBox, oVoxel, oVoxel.pVertex[1], oVoxel.pVertex[2] );
    bool bV20AxisExists = CheckSeparatingAxis( oBox, oVoxel, oVoxel.pVertex[2], oVoxel.pVertex[0] );
    
    if( bV01AxisExists || bV12AxisExists || bV20AxisExists  )
    {
      return false;
    }
    return true;      // no separating axis found
  }
  
  //-------------------------------------------------------------------------------------
  //
  //  TriangleAABBOverlap()
  //
  //  Return:    True     if the AABB (axis aligned bounding box) overlaps with the voxel
  //             False    otherwise  (This condition can be met if one separating axis is found)
  //
  //  NOTE:  This function only works in 2D.  Specifically, it only works with the x and y
  //         coordinates.
  //
  //
  //-------------------------------------------------------------------------------------
  bool TriangleAABBOverlap( const BBox2D &oBox, const SVector3 &p0, const SVector3 &p1, const SVector3 &p2 )
  {
    SVoxel oVoxel;
    oVoxel.pVertex[0] = p0;
    oVoxel.pVertex[1] = p1;
    oVoxel.pVertex[2] = p2;
    
    return  VoxelAABBOverlap( oBox, oVoxel );
    
  }

  //-------------------------------------------------------------------------------------
  //
  //  A point wrapper for TriangleAABBOverlap
  //
  //-------------------------------------------------------------------------------------
  bool TriangleAABBOverlap( const BBox2D &oBox, const Point & p1, const Point & p2, const Point & p3)
  {
    return TriangleAABBOverlap( oBox,
                                SVector3( p1.x, p1.y, 0),
                                SVector3( p2.x, p2.y, 0),
                                SVector3( p3.x, p3.y, 0)  );
  }

  //-------------------------------------------------------------------------------------
  //
  //  Calculate a bounding box for the triangle specified by the points, then
  //  check for overlap
  //
  //-------------------------------------------------------------------------------------
  bool TriangleAABBOverlapEstimate( const BBox2D &oBox, const Point & p1, const Point & p2, const Point & p3)
  {
    Point pMax;
    Point pMin;

    pMax.x = std::max( p1.x, p2.x);
    pMax.x = std::max( p3.x, pMax.x);

    pMax.y = std::max( p1.y, p2.y);
    pMax.y = std::max( p3.y, pMax.y);


    pMin.x = std::min( p1.x, p2.x);
    pMin.x = std::min( p3.x, pMin.x);

    pMin.y = std::min( p1.y, p2.y);
    pMin.y = std::min( p3.y, pMin.y);
        
    BBox2D oVoxelBBox( pMin, pMax );
    return oVoxelBBox.Overlaps( oBox );
  }
}
