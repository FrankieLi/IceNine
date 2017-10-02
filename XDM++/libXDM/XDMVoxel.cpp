////////////////////////////////////////////////////////////////
//
//  File:    XDMVoxel.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Specialization of voxel class
//
//
/////////////////////////////////////////////////////////////////

#include "XDMVoxel.h"

namespace XDMUtility
{
 
  //---------------------------------------------------------------------------------------------------
  //  Public: Set
  //---------------------------------------------------------------------------------------------------
  void CXDMVoxel::Set( const SVector3 &v1, const SVector3 &v2,
                       const SVector3 & v3,
                       bool bIsUp, Float fSW, Int nGen )
  {
    pVertex[0] = v1;
    pVertex[1] = v2;
    pVertex[2] = v3;
    bVoxelPointsUp = bIsUp;
    fSideLength = fSW;
    nGeneration = nGen;
  }
  
  //---------------------------------------------------------------------------------------------------
  //  Make the vertices follow a counter clockwise fashion.
  //  This voxel must be on the x-y plane
  //---------------------------------------------------------------------------------------------------
  void CXDMVoxel::MakeXYCounterClockwise()
  {
    SVector3 v01 = pVertex[0] - pVertex[1];
    SVector3 v21 = pVertex[2] - pVertex[1];
    SVector3 crossRes;
    crossRes = Cross( v21, v01 );
    
    // result is in z direction only, if not CCW, swap
    if ( crossRes.m_fZ < 0 )
    {
      SVector3 tmp = pVertex[2];
      pVertex[2] = pVertex[0];
      pVertex[0] = tmp;
    }
  }

  //---------------------------------------------------------------------------------------------------
  //  GetBoundingBox  (note: only on xy plane)
  //---------------------------------------------------------------------------------------------------
  BBox2D CXDMVoxel::GetBoundingBox ( ) const
  {
    BBox2D oRes;
    Float fMinX = std::min( pVertex[0].m_fX, pVertex[1].m_fX );
    fMinX       = std::min( fMinX, pVertex[2].m_fX );
    Float fMinY = std::min( pVertex[0].m_fY, pVertex[1].m_fY );
    fMinY       = std::min( fMinY, pVertex[2].m_fY );

    Float fMaxX = std::max( pVertex[0].m_fX, pVertex[1].m_fX );
    fMaxX       = std::max( fMaxX, pVertex[2].m_fX );
    Float fMaxY = std::max( pVertex[0].m_fY, pVertex[1].m_fY );
    fMaxY       = std::max( fMaxY, pVertex[2].m_fY );

    oRes.pMin = Point( fMinX, fMinY );
    oRes.pMax = Point( fMaxX, fMaxY );

    return oRes;
  }

  //---------------------------------------------------------------------------------------------------
  //  GetGeneration
  //---------------------------------------------------------------------------------------------------
  Int CXDMVoxel::GetGeneration() const
  {
    return nGeneration;
  }
  
  //---------------------------------------------------------------------------------------------------
  //  Vertex
  //---------------------------------------------------------------------------------------------------
  SVector3 & CXDMVoxel::Vertex( Int n )
  {
    DEBUG_ASSERT( n <= Int( Polygon::NUM_VERTICES ), "GetVertex: index exceeds Polygon::NUM_VERTICES\n");
    return pVertex[n];
  }

  //---------------------------------------------------------------------------------------------------
  //  GetVertex
  //---------------------------------------------------------------------------------------------------
  const SVector3 & CXDMVoxel::Vertex( Int n ) const
  {
    DEBUG_ASSERT( n <= Int( Polygon::NUM_VERTICES ), "GetVertex: index exceeds Polygon::NUM_VERTICES\n");
    return pVertex[n];
  }

  //---------------------------------------------------------------------------------------------------
  //
  //  Private: CheckAxis
  //
  //  Return true if the axis perpendicular to u = v2 - v1 is a separating axis  (i.e., no overlap)
  //  Therefore, we project to this perpendicular axis
  //
  //  Currently only works with 2D objects on the x-y plane
  //
  //
  //---------------------------------------------------------------------------------------------------
  bool CXDMVoxel::CheckSeparatingAxis( const CXDMVoxel & oVoxel, const CXDMVoxel & oTestVoxel,
                                       const SVector3 &oV1, const SVector3 &oV2) const
  {
    SVector3 oDir = oV2 - oV1;
    oDir.Normalize();

    // rotate 90 degress y <- x, x <- -y
    Float tmp = oDir.m_fY;
    oDir.m_fY = oDir.m_fX;
    oDir.m_fX = - tmp;

    Float fMin = MAX_FLOAT;
    Float fMax = MIN_FLOAT;

    for ( int i = 0; i < 3; i ++ )
    {
      Float fAxisProjection = Dot( oTestVoxel.pVertex[i], oDir );

      if ( fAxisProjection < fMin )
        fMin = fAxisProjection;
      if ( fAxisProjection > fMax )
        fMax = fAxisProjection;
    }
  
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

  //---------------------------------------------------------------------------------------------------
  //  MinDistance calculation
  //  --  clearly this could be optimized.  Let's get it right first.
  //
  //  Precondition:  All of the Z must be exactly the same!
  //
  //---------------------------------------------------------------------------------------------------
  bool CXDMVoxel::Connected( const SVector3 & oV1, const SVector3 & oV2,
                             const SVector3 & p, Float fError )
  {
    SVector3 oDist_Top = p - oV2;
    SVector3 oDist = p - oV1;
    oDist_Top.m_fZ = 0;   // zeroing out the out of plane components
    oDist.m_fZ     = 0;
    // check for vertex overlap
    if ( (    fabs( oDist_Top.m_fX ) < fError
           && fabs( oDist_Top.m_fY ) < fError
         )
      ||
         (    fabs( oDist.m_fX ) < fError
           && fabs( oDist.m_fY ) < fError
         )
       )
    {
      return true;
    }
    
    SVector3 oLine = oV2 - oV1;
    oLine.m_fZ = 0;
    Float    fLineLength = oLine.GetLength();
    DEBUG_ASSERT( fLineLength > std::numeric_limits<Float>::epsilon(),
                  "CXDMVoxel:: MinDistance: numerical limit reached for this floaing point numbers\n" ); 

    SVector3 oLineDir = oLine / fLineLength;
    Float    fProjected = Dot( oDist, oLineDir );
    if( fProjected < 0 || fProjected > fLineLength )
      return false;
    
    SVector3 oPerpDir( -oLineDir.m_fY, oLineDir.m_fX, oLineDir.m_fZ );
    Float fPerpDist = Dot( oDist, oPerpDir );
    return ( fabs(fPerpDist) < fError );

//----------------------------------------
//  -- Also correct, but I'm afraid of dividing and sqrt-ing very small numbers
//     SVector3 oDistPerp = oDist - fProjected * oLineDir;
//     return ( oDistPerp.GetLength() < fError );
//----------------------------------------
  }
  
  //---------------------------------------------------------------------------------------------------
  //
  //  Public: Contains
  //
  //  Use barycentric coordinate to define "in triangle"
  //
  //---------------------------------------------------------------------------------------------------
  bool CXDMVoxel::Contains( const SVector3 &p ) const
  {
  
    SVector3 oB = pVertex[1] - pVertex[0];
    SVector3 oC = pVertex[2] - pVertex[0];
    SVector3 oP = p - pVertex[0];

    // should check close-ness to zero
    Float fU = ( oP.m_fY * oC.m_fX - oP.m_fX * oC.m_fY ) / ( oB.m_fY * oC.m_fX - oB.m_fX * oC.m_fY );
    Float fV = ( oP.m_fY * oB.m_fX - oP.m_fX * oB.m_fY ) / ( oC.m_fY * oB.m_fX - oC.m_fX * oB.m_fY );
    if ( fU >= Float( 0.0 ) && fV >= Float( 0.0 ) && ( fU + fV ) <= Float( 1.0 ) )
      return true;
    else
      return false;
  }

  //---------------------------------------------------------------------------------------------------
  //  Public:   PerfectMatch
  //  Purpose:  Handles the degenerate case where v and this voxel
  //            have the same vertices.
  //---------------------------------------------------------------------------------------------------
  bool CXDMVoxel::PerfectMatch( const CXDMVoxel & v ) const
  {

    if ( v.nGeneration != this->nGeneration )
      return false;
    if ( v.bVoxelPointsUp != this->bVoxelPointsUp )
      return false;
  
    Float fError = fSideLength / Float( 10.0 ); 
  
    Bool pVertexMatch[3];

    for( Int i = 0; i < 3; i ++ )
    {
      pVertexMatch[i] = false;
      for( Int j = 0; j < 3; j ++ )
      {
        SVector3 oDisp = v.pVertex[i] - this->pVertex[j];
        if ( oDisp.GetLength() < fError )
          pVertexMatch[i] = true;
      }
    }
  
    return ( pVertexMatch[0] && pVertexMatch[1] && pVertexMatch[2] );
  }
  
}

