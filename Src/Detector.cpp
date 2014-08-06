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
//   Detector.cpp
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:
//    Implementation of detector geometry and physical properties.
//
//
////////////////////////////////////////////////////////////

#include "Detector.h"


///////////////////////////////////////////////////////////////////////////
//
//
//
//                      Class  CDetector
//
//
//
///////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
//
//  CDetector::Constructor:  Note that parameters are required
//
//-------------------------------------------------------------------------
CDetector::CDetector(const BBox2D &range,
                     Float fBeamCenterX, Float fBeamCenterY,
                     Float fPixelXLength, Float fPixelYLength,
                     UInt nRows, UInt nCols,
                     const SVector3 & oImageBasisJ,
                     const SVector3 & oImageBasisK,
                     const SVector3 & oImageCoordinateOrigin)
{
  
  SetImageParameter( range, fBeamCenterX, fBeamCenterY, fPixelXLength,
                     fPixelYLength, nRows, nCols, oImageBasisJ,
                     oImageBasisK, oImageCoordinateOrigin );
  
}

//-------------------------------------------------------------------------
//  SetImageParameter
//-------------------------------------------------------------------------
void CDetector::SetImageParameter(const BBox2D &range,
                                  Float fBeamCenterX, Float fBeamCenterY,
                                  Float fPixelXLength, Float fPixelYLength,
                                  UInt nCols, UInt nRows,
                                  const SVector3 & oImageBasisJ,
                                  const SVector3 & oImageBasisK,
                                  const SVector3 & oImageCoordinateOrigin )
{
  nNumRows = nRows;
  nNumCols = nCols;
  fPixelHeight = fPixelYLength;
  fPixelWidth = fPixelXLength;
  fPixelHalfHeight = fPixelHeight / (Float)2.0;
  fPixelHalfWidth  = fPixelWidth / (Float)2.0;;

  // setting beam origin
  // this is for the YZ configuration

  this->fBeamCenterJ = fBeamCenterX;  // fBeamCenterX and Y are measured in pixels
  this->fBeamCenterK = fBeamCenterY;
  
  this->oImageBasisJ = oImageBasisJ;
  this->oImageBasisK = oImageBasisK;
  this->oImageCoordinateOrigin = oImageCoordinateOrigin;
  
  oDetFrameImageBasisJ = oImageBasisJ;
  oDetFrameImageBasisK = oImageBasisK;
  oDetFrameImageCoordinateOrigin = oImageCoordinateOrigin;
	
  // Setting Detector location and orientation
  oLabToDetectorMatrix.SetIdentity();
  oPrevOrientationMatrix.SetIdentity();

  oDetectorRotationCenter.Set( 0, 0, 0 );    // beam origin in lab frame
  // v3Orientation.Set(0, 0, 0);

  // TODO:  obsolete this with the oImageBasisJ and oImageBasisK  (Need a point from the detector, i.e., beam origin)
  //
  // changed from x-y to y-z plane
  // detector plane is always the y-z plane (by convention)
  {  // set detector plane
    SVector3 v1, v2, v3;
    v1.Set(1, 0, 0);
    v2.Set(0, 1, 0);
    v3.Set(0, 0, 0);
    SetDetectorPlane(v1, v2, v3);
  }
  CalculateImagePlane();
}
	
//-------------------------------------------------------------------------
//
//  SetLocation -- only sets the location, doesn't do anything else
//
//-------------------------------------------------------------------------
void CDetector::SetLocation(const SVector3 &loc)
{
  oDetectorRotationCenter = loc;
  CalculateImagePlane();
}

//-------------------------------------------------------------------------
//
//  SetOrientation --  Only sets the orientation, doesn't do anything else.
//
//-------------------------------------------------------------------------
void CDetector::SetOrientation( const SMatrix3x3 &oOrientation )
{
  oLabToDetectorMatrix = oOrientation;

  // Rotate the basis vectors
  oImageBasisJ = oLabToDetectorMatrix * oDetFrameImageBasisJ;
  oImageBasisK = oLabToDetectorMatrix * oDetFrameImageBasisK;

  oImageCoordinateOrigin = oLabToDetectorMatrix * oDetFrameImageCoordinateOrigin;


  CalculateImagePlane();
}

//-------------------------------------------------------------------------
//
//  SetOrientation --  Only sets the orientation, doesn't do anything else.
//
//-------------------------------------------------------------------------
void CDetector::SetOrientation(Float phi, Float theta, Float psi)
{
  oLabToDetectorMatrix.BuildActiveEulerMatrix( phi, theta, psi );
  SetOrientation( oLabToDetectorMatrix );
}

//-------------------------------------------------------------------------
//
//  Translate
//
//-------------------------------------------------------------------------
void CDetector::Translate(const SVector3 &loc)
{
  oDetectorRotationCenter += loc;
  CalculateImagePlane();
}

//-------------------------------------------------------------------------
//
//  Rotate
//
//-------------------------------------------------------------------------
void CDetector::Rotate(Float phi, Float theta, Float psi)
{
  SMatrix3x3 tmpMatrix;
  
  tmpMatrix.BuildActiveEulerMatrix( phi, theta, psi );
  oLabToDetectorMatrix = tmpMatrix * oLabToDetectorMatrix;

  // Rotate the basis vectors
  oImageBasisJ = oLabToDetectorMatrix * oImageBasisJ;
  oImageBasisK = oLabToDetectorMatrix * oImageBasisK;
  
  oImageCoordinateOrigin = oLabToDetectorMatrix * oImageCoordinateOrigin;  
  CalculateImagePlane();
}

//-------------------------------------------------------------------------
//
//  CalculateImagePlane
//
//-------------------------------------------------------------------------
void CDetector::CalculateImagePlane()
{
  // Rotate
  //
  SVector3 p1 = oLabToDetectorMatrix * oDetectorPt1;
  SVector3 p2 = oLabToDetectorMatrix * oDetectorPt2;
  SVector3 p3 = oLabToDetectorMatrix * oDetectorPt3;

  // Move the points to the proper location  
  //
  p1 = p1 + oDetectorRotationCenter;
  p2 = p2 + oDetectorRotationCenter;
  p3 = p3 + oDetectorRotationCenter;
  
  SVector3 oEdge1 = p3 - p1;
  SVector3 oEdge2 = p2 - p1;
  SVector3 oNormal = Cross( oEdge2, oEdge1 );
  
  oNormal.Normalize();
  
  SVector4 oPlaneNorm;

  // Calculate triangle plane
  oPlaneNorm.m_fX = oNormal.m_fX;
  oPlaneNorm.m_fY = oNormal.m_fY;
  oPlaneNorm.m_fZ = oNormal.m_fZ;
  oPlaneNorm.m_fW = -Dot( oNormal, p1 );  // This determines the location


  oLabFrameImagePlane.m_fA = oPlaneNorm.m_fX;
  oLabFrameImagePlane.m_fB = oPlaneNorm.m_fY;
  oLabFrameImagePlane.m_fC = oPlaneNorm.m_fZ;
  oLabFrameImagePlane.m_fD = oPlaneNorm.m_fW;
}

//-------------------------------------------------------------------------
//
//  Intersect:  produces oRes, which is
//              the vector at the point of
//              intersection in the detector
//              coordinates
//
//-------------------------------------------------------------------------
bool CDetector::Intersects(const CRay &oRay,
                                Float &fT) const
{
  // find the intersection
  return Collision::Intersects(oLabFrameImagePlane, oRay, fT);
}

//-------------------------------------------------------------------------
//
//  LabToPixelCoordinate
//
//-------------------------------------------------------------------------
void CDetector::LabToDetectorCoordinate( Float &fJ, Float &fK, 
                                         const SVector3 &oPos ) const
{
  // THE CORRECT VERSION
  SVector3 oRayDetIntersect = oPos - oDetectorRotationCenter; 
  SVector3 oPixelLoc = oRayDetIntersect - oImageCoordinateOrigin; 
  fJ = Dot( oPixelLoc, oImageBasisJ );
  fK = Dot( oPixelLoc, oImageBasisK );  
}

//-------------------------------------------------------------------------
//
//  DetectorToLabCoordinate
//
//-------------------------------------------------------------------------
SVector3 CDetector::DetectorToLabCoordinate( Float fJ, Float fK ) const
{
  //
  // [ NOT CHECKED ]
  //
  SVector3 oRes = fJ * oImageBasisJ + fK * oImageBasisK;
  oRes = oRes + oDetectorRotationCenter + oImageCoordinateOrigin;
  return oRes;
}

//-------------------------------------------------------------------------
//
//  LabToPixel
//
//------------------------------------------------------------------------
void CDetector::LabToPixel( Float &row, Float &col, const SVector3 &oPos ) const
{
  Float fK, fJ;
  LabToDetectorCoordinate( fJ, fK, oPos );

  //------------------------------------------------------------
  // Note:  Bob's coordinate system
  // Float fJ = oPos.m_fY + oImageCoordinateOrigin.m_fY;
  // Float fK = -( oPos.m_fZ - oImageCoordinateOrigin.m_fZ );
  //
  //  This system is now generalized
  //--------------------------------------------------------------

  row = ToRowPixel( fK );
  col = ToColPixel( fJ );
}
//-------------------------------------------------------------------------
//
//  PixelToLabCoordinate
//  Col = X pixel, Row = Y pixel
//-------------------------------------------------------------------------
SVector3 CDetector::PixelToLabCoordinate( const Float & fXPixel, const Float & fYPixel ) const
{
  SVector3 oRes;

  Float fJ = ColPixelToImageJ( fXPixel );
  Float fK = RowPixelToImageK( fYPixel );

  oRes = DetectorToLabCoordinate( fJ, fK );
  return oRes;
}

//-------------------------------------------------------------------------
//
//  AddDirectBeam - add direct beam into the current detector image
//  
//-------------------------------------------------------------------------
void CDetector::AddDirectBeam( CImageData &oOutputImage, Float fBeamHeight, Float fBeamWidth )
{

  SVector3 oPixelLoc = oDetectorRotationCenter - oImageCoordinateOrigin;
  Float oBeamCenterJ = Dot( oPixelLoc, oImageBasisJ );  // TODO:  Warning, may need to be generalize
  Float oBeamCenterK = Dot( oPixelLoc, oImageBasisK );

  vector<Point> p;
  p.resize(4);
  
  p[0].x = ToRowPixel( oBeamCenterJ + fBeamWidth / 2.0 );
  p[0].y = ToColPixel( oBeamCenterK + fBeamHeight / 2.0 );

  p[1].x = ToRowPixel( oBeamCenterJ - fBeamWidth / 2.0 );
  p[1].y = ToColPixel( oBeamCenterK + fBeamHeight / 2.0 );

  p[2].x = ToRowPixel( oBeamCenterJ - fBeamWidth / 2.0 );
  p[2].y = ToColPixel( oBeamCenterK - fBeamHeight / 2.0 );
  
  p[3].x = ToRowPixel( oBeamCenterJ + fBeamWidth / 2.0 );
  p[3].y = ToColPixel( oBeamCenterK - fBeamHeight / 2.0 );

  oOutputImage.AddPolygon( p, 10000 ); // fill with arbitary value.

}


///////////////////////////////////////////////////////////////////////////
//
//
//
//                     A C C E S S O R S
//
//
//
///////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------
//  GetLocation  (Location is the same as detector center)
//------------------------------------------------------------
SVector3 CDetector::GetLocation( ) const
{
  SVector3 oRet = oDetectorRotationCenter;
  return oRet;
}

//------------------------------------------------------------
//  GetOrientationMatrix
//------------------------------------------------------------
SMatrix3x3 CDetector::GetOrientationMatrix( ) const
{
  SMatrix3x3 oRes =  oLabToDetectorMatrix;
  return oRes;
}
//------------------------------------------------------------
//  GetRotationCenter
//------------------------------------------------------------
SVector3 CDetector::GetRotationCenter( ) const
{
  SVector3 oRet = oDetectorRotationCenter;
  return oRet;
}

//------------------------------------------------------------
//  GetRotationCenter
//------------------------------------------------------------
Float CDetector::GetBeamCenterJ () const
{
  return fBeamCenterJ;
}

//------------------------------------------------------------
//  GetRotationCenter
//------------------------------------------------------------
Float CDetector::GetBeamCenterK () const
{
  return fBeamCenterK;
}

//------------------------------------------------------------
//  GetPixelWidth
//------------------------------------------------------------
Float CDetector::GetPixelWidth  () const
{
  return fPixelWidth;
}

//------------------------------------------------------------
//  GetPixelHeight
//------------------------------------------------------------
Float  CDetector::GetPixelHeight () const
{
  return fPixelHeight;
}

//------------------------------------------------------------
//  GetDetectorCoordinateOrigin -- return coordinate origin of the detector
//  in the lab frame
//------------------------------------------------------------
SVector3 CDetector::GetDetectorCoordinateOrigin( ) const
{
  SVector3 oRet = oImageCoordinateOrigin + oDetectorRotationCenter;
  return oRet;
}

//------------------------------------------------------------
//  GetDetFrameCoordinateOrigin()
//------------------------------------------------------------
SVector3 CDetector::GetDetFrameCoordinateOrigin() const
{
  return oImageCoordinateOrigin;
}

//------------------------------------------------------------
//  GetDetFrameCoordBasis
//------------------------------------------------------------
void CDetector::GetDetFrameCoordBasis( SVector3 & oJBasis, SVector3 & oKBasis ) const
{
  oJBasis = oDetFrameImageBasisJ;
  oKBasis = oDetFrameImageBasisK;
} 

//------------------------------------------------------------
//  GetCoordinateBasis
//------------------------------------------------------------
void CDetector::GetCoordinateBasis( SVector3 & oJBasis, SVector3 & oKBasis ) const
{
  oJBasis = oImageBasisJ;
  oKBasis = oImageBasisK;
}

//------------------------------------------------------------
//  GetDetectorWidth
//------------------------------------------------------------
Float CDetector::GetDetectorWidth( ) const
{
  return (Float) nNumCols * fPixelWidth; 
}

//------------------------------------------------------------
//  GetDetectorHeight
//------------------------------------------------------------
Float CDetector::GetDetectorHeight( ) const
{
  return (Float) nNumRows * fPixelHeight;
}

//------------------------------------------------------------
//  GetNumCols
//------------------------------------------------------------
Int CDetector::GetNumCols( ) const
{
  return nNumCols;
}

//------------------------------------------------------------
//  GetNumRows
//------------------------------------------------------------
Int CDetector::GetNumRows( ) const
{
  return nNumRows;
}

//------------------------------------------------------------
//  GetLabDetectorPlane
//------------------------------------------------------------
CPlane CDetector::GetLabDetectorPlane( )
{
  return oLabFrameImagePlane;
}

//------------------------------------------------------------
//  RowPixelToImageY -- inverse of ToRowPixel
//------------------------------------------------------------
Float CDetector::RowPixelToImageK( Float p ) const
{
  return p * fPixelHeight;
}

//------------------------------------------------------------
//  ColPixelToImageX  -- inverse of ToColPixel
//------------------------------------------------------------
Float CDetector::ColPixelToImageJ( Float p ) const
{
  return p * fPixelWidth;
}

//------------------------------------------------------------
//  ToRowPixel
//
//  Given a coordinate d from 0, convert from d to pixels
//
//  Post condition:  Negative pixels does not exist
//
//------------------------------------------------------------
Int CDetector::ToRowPixel(Float d)  const
{
  if ( d < 0 )
    return -1;    // there should be negative pixel
  else
   return (Int)( (d +   fPixelHalfHeight )/ fPixelHeight );  // use positive truncation instead

}

//------------------------------------------------------------
//  ToColPixel
//
//  Given a coordinate d from 0, convert from d to pixels
//------------------------------------------------------------
Int CDetector::ToColPixel(Float d) const
{
  if ( d < 0 )
    return -1;    // there should be negative pixel
  else
    return (Int)( (d +   fPixelHalfWidth )/ fPixelWidth );  // use positive truncation instead
}


///////////////////////////////////////////////////////////////////////////
//
//
//
//                     U T I L I T I E S
//
//
//
///////////////////////////////////////////////////////////////////////////


//------------------------------------------------------------
//  GetPixelExtent
//
//  Return the width of the voxel
//  in number of pixels
//------------------------------------------------------------
Int CDetector::GetPixelExtent( const SVoxel & oVoxel ) const
{

  SVector3 oMax(0, 0, 0);
  for( Int i = 0; i < 3; i ++ )
  {
    SVector3 oDiff = oVoxel.pVertex[i] - oVoxel.pVertex[ (i + 1) % 3 ];
    
    oMax.m_fX = std::max( oMax.m_fX, fabs( oDiff.m_fX ) );
    oMax.m_fY = std::max( oMax.m_fY, fabs( oDiff.m_fY ) );
    oMax.m_fZ = std::max( oMax.m_fZ, fabs( oDiff.m_fZ ) );
  }

  Float fMax = oMax.m_fX;
  fMax = std::max( fMax, oMax.m_fY );
  fMax = std::max( fMax, oMax.m_fZ );

  Float fDim = std::min( fPixelWidth, fPixelHeight );
  return Int( fMax / fDim );
}



///////////////////////////////////////////////////////////////////////////
//                  CXDMDetectorFactory
///////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------
//  GetCoordOrigin
//------------------------------------------------------------
SVector3 CXDMDetectorFactory::GetCoordOrigin( Float fBeamCenterJ, Float fBeamCenterK,
                                              Float fPixelWidth,  Float fPixelHeight,
                                              const SVector3 & oJUnitVector,
                                              const SVector3 & oKUnitVector )
{
  SVector3 oCoordOrigin;
  oCoordOrigin  = - (Float) fBeamCenterJ * fPixelWidth  * oJUnitVector;
  oCoordOrigin += - (Float) fBeamCenterK * fPixelHeight * oKUnitVector;
  return oCoordOrigin;
}

//------------------------------------------------------------
// MakeDetector
//
//------------------------------------------------------------
CDetector CXDMDetectorFactory::MakeDetector( Int nNumJPixels, Int nNumKPixels,
                                             const SVector3 & oDetLabFrameLoc,
                                             Float fBeamCenterJ, Float fBeamCenterK,
                                             const SVector3 &oJUnitVector,
                                             const SVector3 &oKUnitVector,
                                             Float fPixelWidth, Float fPixelHeight,
                                             const SMatrix3x3 & oDetLabFrameOrientation )
{
  BBox2D range( Point( 0, 0 ),
                Point( nNumJPixels, nNumKPixels ) );

  SVector3 oCoordOrigin = GetCoordOrigin( fBeamCenterJ, fBeamCenterK,
                                          fPixelWidth, fPixelHeight,
                                          oJUnitVector, oKUnitVector );
  
  // note that the col <-> x, row <-> y; hence the flip in the last two parameters
  CDetector oDetector( range, fBeamCenterJ, 
                       fBeamCenterK, fPixelWidth,
                       fPixelHeight, nNumJPixels, nNumKPixels,
                       oJUnitVector, oKUnitVector, oCoordOrigin );
  
  oDetector.Translate( oDetLabFrameLoc );
  oDetector.SetOrientation( oDetLabFrameOrientation );
  
  return oDetector;
}

//------------------------------------------------------------
//  ModifyImageParameters
//------------------------------------------------------------
void CXDMDetectorFactory::ModifyImageParameters( CDetector & oDetector,
                                                 const SVector3 & oDetLabFrameLoc,
                                                 Float fBeamCenterJ, Float fBeamCenterK,
                                                 Float fPixelWidth, Float fPixelHeight,
                                                 const SMatrix3x3 & oDetLabFrameOrientation )
{

  Int nNumJPixels = oDetector.GetNumCols();
  Int nNumKPixels  = oDetector.GetNumRows();
  
  BBox2D range( Point(0, 0),
                Point( nNumJPixels, nNumKPixels) );
  SVector3 oJUnitVector, oKUnitVector;
  oDetector.GetDetFrameCoordBasis( oJUnitVector, oKUnitVector );
  SVector3 oNewCoordOrigin = GetCoordOrigin( fBeamCenterJ, fBeamCenterK,
                                             fPixelWidth, fPixelHeight,
                                             oJUnitVector, oKUnitVector );
  oDetector.SetImageParameter( range, fBeamCenterJ, 
                               fBeamCenterK, fPixelWidth,
                               fPixelHeight, nNumJPixels, nNumKPixels,
                               oJUnitVector, oKUnitVector, oNewCoordOrigin );

  oDetector.SetLocation( oDetLabFrameLoc );
  oDetector.SetOrientation( oDetLabFrameOrientation );
}

//------------------------------------------------------------
//  GetImageParameters
//------------------------------------------------------------
CXDMDetectorFactory::SDetParameters
CXDMDetectorFactory::GetImageParameters( const CDetector & oDetector )
{
  SDetParameters oRes;

  oRes.oOrientation = oDetector.GetOrientationMatrix();
  oRes.oPosition    = oDetector.GetRotationCenter   ();
  oRes.fBeamCenterJ = oDetector.GetBeamCenterJ      ();
  oRes.fBeamCenterK = oDetector.GetBeamCenterK      ();
  oRes.fPixelWidth  = oDetector.GetPixelWidth       ();
  oRes.fPixelHeight = oDetector.GetPixelHeight      ();
  oRes.nNumJPixels  = oDetector.GetNumCols          ();
  oRes.nNumKPixels  = oDetector.GetNumRows          ();

  return oRes;
}

//------------------------------------------------------------
//  GetImageInfo
//------------------------------------------------------------
InitFileIO::CDetectorInfo
CXDMDetectorFactory::GetImageInfo( const CDetector & oDetector )
{

  InitFileIO::CDetectorInfo oRes;
  
  oDetector.GetDetFrameCoordBasis( oRes.vJUnitVector, oRes.vKUnitVector );
  
  oRes.fBeamCenterJ          = oDetector.GetBeamCenterJ      ();
  oRes.fBeamCenterK          = oDetector.GetBeamCenterK      ();
  oRes.oLabFrameLocation     = oDetector.GetRotationCenter   ();
  oRes.oLabFrameOrientMatrix = oDetector.GetOrientationMatrix();
  
  oRes.vLabFrameOrientation  = oDetector.GetOrientationMatrix().GetEulerAngles();     

  oRes.nNumJPixels           = oDetector.GetNumCols          ();
  oRes.nNumKPixels           = oDetector.GetNumRows          (); 
  oRes.fPixelWidth           = oDetector.GetPixelWidth       ();
  oRes.fPixelHeight          = oDetector.GetPixelHeight      ();

  return oRes;
}
