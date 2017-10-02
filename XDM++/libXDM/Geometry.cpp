/////////////////////////////////////////////////////////////////
//
//  File:    Geometry.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//
//  Purpose: Implementations of geometric primitives.
//
/////////////////////////////////////////////////////////////////


#include "Geometry.h"

//-----------------------------------------------------------------
//
// Useful mesh geometric operations
//
//-----------------------------------------------------------------
namespace MeshGeometry
{

  //-----------------------------------------------------------------
  //  MidPoint
  //-----------------------------------------------------------------
  SVector3 MidPoint(const SVector3 & p1, const SVector3 & p2)
  {
    SVector3 oMidPoint;

    oMidPoint = p1 + p2;
    oMidPoint /= (Float)2.0;

    return oMidPoint;
  }
  
  //---------------------------------------------------------------
  //  IsLeft - Basic line side test.  Given a line that goes from p1 to p2,
  //            is p3 on the left or the right side of this line
  //
  //  This can be computed by looking at the cross product of ( p2 - p1 ) x (p3 - p1) 
  //
  //  Returns true of point3 is to the left of OR ON the line p1 -> p2.
  //  
  //
  //-----------------------------------------------------------------
  bool IsLeft(const SVector3 &p1, const SVector3 &p2, const SVector3 &p3, Float fError)
  {

    // in the 2D case of MIC file, we only use the first 2 coordinates

    Float fUx = p2.m_fX - p1.m_fX;
    Float fUy = p2.m_fY - p1.m_fY;
    
    Float fVx = p3.m_fX - p1.m_fX;
    Float fVy = p3.m_fY - p1.m_fY;
    return ( ( fUx * fVy - fUy * fVx )  >  - fError );
    
  }

  //-----------------------------------
  //  Intersections
  //  Purpose: Line segment to line segment intersection
  //
  //  Input:   (u0, u1), (v0, v1) form two line segments
  //           fError is the relative error (fractional)
  //
  //  Return:  t1, t2, such that the two lines intersect
  //           intersect with the segment within epsilon specified.
  //  Valididty:  t1, t2 must be within the range [0, 1] if an         
  //              intersection has occored.
  //
  //-----------------------------------
  std::pair<Float, Float> Intersects( const SVector3 & u0, const SVector3 & u1,
                                      const SVector3 & v0, const SVector3 & v1,
                                      Float fRelError )
  {
    // choose u0, v0 as p0, p0'
    SVector3 oDiff = u0 - v0;
    SVector3 R0    = u1 - u0;
    SVector3 R1    = v1 - v0;
    Float Norm0 = R0.GetLength();
    Float Norm1 = R1.GetLength();
    Float fAbsError = fRelError * ( Norm0 + Norm1 );
    R0 /= Norm0;
    R1 /= Norm1;

    Float fDenom0 = R0.m_fX * R1.m_fY - R1.m_fX * R0.m_fY;
    Float fDenom1 = R1.m_fX * R0.m_fY - R0.m_fX * R1.m_fY ;
    if(    fabs( fDenom0 ) > fAbsError
        && fabs( fDenom1 ) >  fAbsError )
    {
      Float t0 = fabs( oDiff.m_fX * R0.m_fY - oDiff.m_fY * R0.m_fX );
      Float t1 = fabs( oDiff.m_fX * R1.m_fY - oDiff.m_fY * R1.m_fX );
      t0 /= fDenom0;
      t1 /= fDenom1;
      return std::pair<Float, Float>( t0, t1 );
    }

    return std::pair<Float, Float>( Float(-1), Float(-1) );
  }
 

  
} // end MeshGeometry
