/////////////////////////////////////////////////////////////////
//
//  File:    Voxel.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Contains basic voxel and Voxel iterators
//
//
/////////////////////////////////////////////////////////////////

#ifndef _VOXEL_H_
#define _VOXEL_H_

#include "3dMath.h"
#include "Types.h"
#include <vector>
#include "BBox.h"
#include "Geometry.h"
#include <utility>
#include "XDMVoxel.h"

using std::vector;
using namespace XDMUtility;

///////////////////////////////////////////////////////////////////
//
//  Class CEverythingVoxel
//
//  General Coordinate Voxel
//
//
///////////////////////////////////////////////////////////////////

class CEverythingVoxel : public CXDMVoxel
{
  
private:  
  //----------------------------------------------
  //  CopyData
  //----------------------------------------------
  void CopyData( CEverythingVoxel & oLHS, const CEverythingVoxel & oRHS ) const
  {
    oLHS.oDeformation    = oRHS.oDeformation;
    oLHS.oOrientMatrix   = oRHS.oOrientMatrix;
    oLHS.fConfidence     = oRHS.fConfidence;
    oLHS.fCost           = oRHS.fCost;
    oLHS.fPixelOverlapRatio = oRHS.fPixelOverlapRatio;
    oLHS.nPhase          = oRHS.nPhase;
    oLHS.nID             = oRHS.nID;
    oLHS.fFittingTime    = oRHS.fFittingTime;
  }
  
public:
  
  SMatrix3x3 oDeformation;   // This isn't exactly strain.  More like a metric modifier
  SMatrix3x3 oOrientMatrix;   
  Float      fConfidence;
  Float      fCost;
  Float      fPixelOverlapRatio;
  Float      fFittingTime;
  Int        nPhase;
  Int        nID;

  CEverythingVoxel( ){ oDeformation.SetIdentity(); }
  CEverythingVoxel( const CXDMVoxel & oRHS ):CXDMVoxel( oRHS ){}
          
  CEverythingVoxel( const CEverythingVoxel & oRHS ):
    CXDMVoxel( oRHS ),
    oDeformation ( oRHS.oDeformation ),
    oOrientMatrix( oRHS.oOrientMatrix ),
    fConfidence  ( oRHS.fConfidence   ),
    fCost( oRHS.fCost ),
    fPixelOverlapRatio( oRHS.fPixelOverlapRatio),
    fFittingTime( oRHS.fFittingTime ),
    nPhase       ( oRHS.nPhase ),
    nID          ( oRHS.nID )
  {}

    
  // subdivde
  void Subdivide( vector<CEverythingVoxel> & vChildren ) const
  {
    CXDMVoxel::Subdivide( vChildren );
    for( Size_Type i = 0; i < vChildren.size(); i ++ )
      CopyData( vChildren[i], *this );
  }
  void Subdivide( vector<CEverythingVoxel> & vChildren, Int nMaxGen ) const
  {
    CXDMVoxel::Subdivide( vChildren, nMaxGen );
    for( Size_Type i = 0; i < vChildren.size(); i ++ )
      CopyData( vChildren[i], *this );
  }    
};



struct SquareVoxel
{
  SVector3 pVertex[4];
  SVector3 Center;
  
  Float    fSideLength;
  //  Bool     bVoxelPointsUp;  // up or down
  short    nPhase;
  Int      nID;
  
  //  SMatrix3x3 oDeformation;   // This isn't exactly strain.  More like a metric modifier
  SMatrix3x3 oOrientMatrix;   
  Float fConfidence;
  Float fCost;
  Float fPixelOverlapRatio;
  Float fFittingTime;

  //  default ctor
  SquareVoxel()
  {}

  // type cast
  //
  //   THIS IS A HACK - this is done so that I don't have to make a new
  //   file type just now.  Will have to be changed soon.
  SquareVoxel( const CEverythingVoxel & v )
  {
    fSideLength        = v.fSideLength;
    nPhase             = v.nPhase;
    nID                = v.nID;
    //    oDeformation       = v.oDeformation;
    oOrientMatrix      = v.oOrientMatrix; 
    fConfidence        = v.fConfidence;
    fCost              = v.fCost;
    fPixelOverlapRatio = v.fPixelOverlapRatio;
    fFittingTime       = v.fFittingTime;
    SetCenter( CEverythingVoxel::LeftVertex( v ), fSideLength );
  }

  void SetCenter( const SVector3 & v, float fSideLength_ )
  {
    fSideLength = fSideLength_;
    SetCenter(v);
  }
  void SetCenter( const SVector3 & v)
  {
    pVertex[0] = pVertex[1] = pVertex[2] = pVertex[3] = Center = v;
    pVertex[0].m_fX -= fSideLength / static_cast< Float > ( 2 );
    pVertex[1].m_fX -= fSideLength / static_cast< Float > ( 2 );
    pVertex[2].m_fX += fSideLength / static_cast< Float > ( 2 );
    pVertex[3].m_fX += fSideLength / static_cast< Float > ( 2 );

    pVertex[0].m_fY -= fSideLength / static_cast< Float > ( 2 );
    pVertex[1].m_fY += fSideLength / static_cast< Float > ( 2 );
    pVertex[2].m_fY -= fSideLength / static_cast< Float > ( 2 );
    pVertex[3].m_fY += fSideLength / static_cast< Float > ( 2 );
  }
  
  SVector3 GetCenter() const
  {
    return Center;
  }
  
  //-------------------------------
  //  Vertex()
  //-------------------------------
//   SVector3 & Vertex( Int n )
//   {
//     return pVertex[n];
//   }
    
  const SVector3 & Vertex( Int n ) const
  {
    return pVertex[n];
  }
};


typedef CEverythingVoxel SVoxel;



#endif
