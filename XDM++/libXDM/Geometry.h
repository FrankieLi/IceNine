/////////////////////////////////////////////////////////////////
//
//  File:    Geometry.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//
//  Purpose: A collection of geometric primitives.
//
/////////////////////////////////////////////////////////////////

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "3dMath.h"
#include <math.h>
#include <vector>
using namespace GeneralLib;

//-----------------------------------
//
// Useful mesh geometric operations
//
//-----------------------------------
namespace MeshGeometry
{
  //-----------------------------------
  //  MidPoint
  //-----------------------------------
  SVector3 MidPoint( const SVector3 & p1, const SVector3 & p2 );
  
  //-----------------------------------
  //  IsLeft
  //-----------------------------------
  bool IsLeft( const SVector3 &p1, const SVector3 &p2, const SVector3 &p3, Float fError );
  
  //-----------------------------------
  //  GetPlanarOverlappingEdge
  //
  //  Purpose:  Given two overlapping planar shapes, return the edge that's shared 
  //            by the two objects.  Note that the smallest overlapping edge will be
  //            returned.  In another words, line returned will be formally the
  //            intersection of v1, v1
  //  Return:   True if such a line exists.  False is returned otherwise.
  //-----------------------------------
  template< typename ShapeT, typename EdgeT > 
  Bool GetOverlappingEdge( EdgeT & oRes, ShapeT oVBound, ShapeT oNgb, Float fError )
  {
    if( oVBound.nGeneration < oNgb.nGeneration )
      std::swap( oVBound, oNgb );
    
    // find connecting vertices;
    std::vector< typename ShapeT::Trait::VertexT > oConnectingVertices;
    for( Int i = 0; i < ShapeT::NumVertices; i ++ )
      if( ShapeT::Connected( oVBound.pVertex[i], oNgb, fError ) )
        oConnectingVertices.push_back( oVBound.pVertex[i] );

    if( oConnectingVertices.size() == 2 )
    {
      oRes = EdgeT( oConnectingVertices[0], oConnectingVertices[1] );
      return true;
    }
    else
    {
      return false;
    }
  }


  //-----------------------------------
  //  Intersections
  //  Purpose: Line segment to line segment intersection
  //
  //  fError is the RELATIVE error (fractional)
  //
  //  Input:   (u0, u1), (v0, v1) form two line segments
  //  Return:  t1, t2, such that the two lines intersect
  //           intersect with the segment within epsilon specified.
  //  Valididty:  t1, t2 must be within the range [0, 1] if an         
  //              intersection has occored.
  //
  //  NOT TESTED
  //
  //-----------------------------------
  std::pair<Float, Float> Intersects( const SVector3 & u0,
                                      const SVector3 & u1,
                                      const SVector3 & v0,
                                      const SVector3 & v1,
                                      Float fRelError = 0.001 );
 
  
}

#endif
