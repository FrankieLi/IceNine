/////////////////////////////////////////////////////////////////
//
//  File:    XDMVoxel.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Specialization of voxel class
//
//
/////////////////////////////////////////////////////////////////

#ifndef _XDM_VOXEL_H_
#define _XDM_VOXEL_H_

#include "Polygon.h"
#include <vector>
#include "BBox.h"
#include "Geometry.h"
#include <limits>

namespace XDMUtility
{
  using std::vector;
  using namespace PBRMath;
  using namespace GeneralLib;
  

  //-------------------------------
  //  MIC file / XDM specific
  //-------------------------------
  enum ObjLabel {
    UP = 1,
    DOWN = 2,
    UNFITTED = 0
  };

  struct SFloatError
  {
    Float fError;
    SFloatError( Float fError_ ):fError( fError_ ) {}
    Float operator()() const { return fError;  }
  };
  
  //-------------------------------
  //  Trait of CXDMVoxel
  //-------------------------------
  struct CXDMVoxelTrait
  {
    typedef BBox2D   BBoxT;
    typedef SVector3 VertexT;
    typedef Int      VertexIndexT;
  };
  
  //-------------------------------
  //  Generalized
  //-------------------------------
  class CXDMVoxel : public Polygon::CTriangle< SVector3 >,
                    public Polygon::CShape< CXDMVoxelTrait >
  {
  private:

    //------------------------------------
    //  CheckSeparatingAxis
    //  Return true if the axis perpendicular to u = v2 - v1 is a separating axis  (i.e., no overlap)
    //  Therefore, we project to this perpendicular axis
    //
    //  Currently only works with 2D objects on the x-y plane
    //
    //------------------------------------
    bool CheckSeparatingAxis( const CXDMVoxel & oVoxel, const CXDMVoxel &oTestVoxel, 
                              const SVector3 &oV1, const SVector3 &oV2 ) const;

    //------------------------------------
    //  Connected
    //  Purpose:   Return true if the line described by v1, v2 is "connected"
    //             to the point p within error fError.  This is done by calculating
    //             the minium distance between the line specified by
    //             the first two vertices and the point on the third one.
    //             Minimum distance is defined to be the perpendicular distance
    //             of the point to the line.
    //
    //------------------------------------
    static bool Connected( const SVector3 & oV1, const SVector3 & oV2,
                           const SVector3 & p, Float fError );

  public:
    
    //------------------------------------
    typedef CXDMVoxelTrait Trait;
    //------------------------------------
    
    Int nGeneration;
    Bool bVoxelPointsUp;  // up or down
    Float fSideLength;

    //------------------------------------
    //  Using standard default constructors
    //------------------------------------
    CXDMVoxel(): nGeneration(-1), bVoxelPointsUp( false ), fSideLength(0){}
    CXDMVoxel( const CXDMVoxel & v_ ):
      Polygon::CTriangle<SVector3>( v_ ),
      nGeneration   ( v_.nGeneration ),
      bVoxelPointsUp( v_.bVoxelPointsUp ),
      fSideLength   ( v_.fSideLength ){}
    
    //-------------------------------
    //  Set
    //-------------------------------
    void Set( const SVector3 &v1, const SVector3 &v2,
              const SVector3 & v3,
              bool bIsUp, Float fSW, Int nGen );

    //-------------------------------
    //  Make the vertices follow a counter clockwise fashion.
    //  This voxel must be on the x-y plane
    //-------------------------------
    void MakeXYCounterClockwise();
     
    //----------------------------
    //  Geometry
    //----------------------------
  
    //---------------------------------
    //  Contains(const SVector3 &p)
    //  Action:  Return true if p is
    //           within this voxel
    //---------------------------------
    bool Contains( const SVector3 &p ) const;

    //---------------------------------
    //  PerfectMatch
    //  Purpose:  Handles the degenerate case where v and this voxel
    //            have the same vertices.
    //
    //  Action:   Return true if v is overlapping exactly to this voxel.
    //---------------------------------
    bool PerfectMatch( const CXDMVoxel & v ) const;
    
    //----------------------------
    //  Accessors
    //----------------------------
    
    //-------------------------------
    //  GetBoundingBox
    //
    //  Action:  Return in plane bounding box.
    //           i.e., bounding box in x-y plane
    //           only.
    //-------------------------------
    BBox2D GetBoundingBox ( ) const;

    //-------------------------------
    //  GetGeneration
    //-------------------------------
    Int GetGeneration() const;

    //-------------------------------
    //  Vertex()
    //-------------------------------
    SVector3 & Vertex( Int n );
    const SVector3 & Vertex( Int n ) const;
    
    //-------------------------------
    //  GetCenter - return the center of the voxel
    //-------------------------------
    SVector3 GetCenter() const
    {
      SVector3 oRet;
      oRet = pVertex[0] + pVertex[1];
      oRet += pVertex[2];
      oRet /= 3.0;
      return oRet;
    }

    //-------------------------------
    //  ResizeAndCenter
    //
    //  Purpose:  Resize the voxel and recenter to the previous
    //            center location.
    //
    //  Main usage:  Applying ResizeAndCenter with a new sidelength
    //               of ( 1 - delta) * originalSideLength will lead
    //               a better resolution between "overlap" of two
    //               voxels vs. voxels that share edges.
    //-------------------------------
    void ResizeAndCenter( Float fNewSideLength )
    {
      SVector3 oCenter = GetCenter();
      Float    fScale  = fNewSideLength / fSideLength;

      for( Int i = 0; i < 3; i ++ )
      {
        SVector3 oDiff = pVertex[i] - oCenter;
        oDiff         *= fScale;
        pVertex[i]     = oCenter + oDiff;
      }
      fSideLength = fNewSideLength;
    }
    
    //---------------------------------------------------------------
    //  Templated functions required from Shape Concept
    //---------------------------------------------------------------
    template< typename ObjT >
    static SVector3 LeftVertex( const ObjT & oObj )
    {
      Float fMinX = std::numeric_limits<Float>::max();
      Int nLeftMost = 0;
      for ( int i = 0; i < 3; i ++ )
      {
        if( oObj.Vertex(i).m_fX < fMinX)
        {
          fMinX = oObj.Vertex(i).m_fX;
          nLeftMost = i;
        }
      }
      return oObj.Vertex( nLeftMost );
    }

    //---------------------------------------------------------------
    //  Return direction of object (if there is one)
    //---------------------------------------------------------------
    template< typename ObjT >
    static ObjLabel Label( const ObjT & oObj )
    {
      ObjLabel nDir = UNFITTED;
      if ( oObj.bVoxelPointsUp )
      {
        nDir = UP;
      }
      else
      {
        nDir = DOWN;
      }
      return nDir;
    }
    
    //-------------------------------
    //  Connected
    //   -- Test to see if the point p is connected to the object's
    //      *boundary*
    //-------------------------------
    template< typename ObjT >
    static bool Connected( const SVector3 & p, const ObjT & oObj, Float fEpsilon )
    {
      for( Int j = 0; j < 3; j ++ )
      {
        bool bVertexConnected = Connected( oObj.pVertex[j],
                                           oObj.pVertex[ (j + 1) % 3 ],
                                           p, fEpsilon );
        if( bVertexConnected )
          return true;
      }
      return false;
    }
    
    //-------------------------------
    //  Connected
    //   -- test for voxel connectedness
    //-------------------------------
    template < typename ObjT, typename FNumericLimit >
    bool Connected( const ObjT & oRHS, FNumericLimit FError,
                    Int nVerticesToConnect = 2 )  const
    {
      Float fEpsilon = FError();
      Int nVertexConnected = 0;
      for( Int i = 0; i < 3; i ++ )
      {
        if( Connected( pVertex[i], oRHS, fEpsilon ) )
          nVertexConnected ++;
        if( nVertexConnected >= nVerticesToConnect )
          return true;
      }

      nVertexConnected = 0;
      for( Int i = 0; i < 3; i ++ )
      {
        if( Connected( oRHS.pVertex[i], *this, fEpsilon ) )
          nVertexConnected ++;
        if( nVertexConnected >= nVerticesToConnect )
          return true;
      }
      return false;
    }
    
    //-------------------------------
    //  ShareEdge
    //   -- test for voxel connectedness
    //-------------------------------
    template < typename ObjT, typename FNumericLimit >
    bool ShareEdge( const ObjT & oRHS, FNumericLimit FError )  const
    {
      Float fEpsilon = FError();
      Int nVertexConnected = 0;
      for( Int i = 0; i < 3; i ++ )
        for( Int j = 0; j < 3; j ++ )
        {
          Bool bSuccess = Connected( pVertex[i],
                                     pVertex[ (i + 1) % 3 ],
                                     oRHS.pVertex[j], fEpsilon );
          bSuccess = bSuccess && Connected( pVertex[i],
                                            pVertex[ (i + 1) % 3 ],
                                            oRHS.pVertex[ (j + 1) % 3 ], fEpsilon );
          if( bSuccess )
            return true;
        }
      return false;
    }

    //-------------------------------
    //  Overlap
    //-------------------------------
    template < typename ObjT >
    bool Overlaps( const ObjT & oRHS ) const
    {
      bool bThisV01 = CheckSeparatingAxis( oRHS, *this, pVertex[0], pVertex[1] );
      bool bThisV12 = CheckSeparatingAxis( oRHS, *this, pVertex[1], pVertex[2] );
      bool bThisV20 = CheckSeparatingAxis( oRHS, *this, pVertex[2], pVertex[0] );
      
      bool bQueryV01 = CheckSeparatingAxis( *this, oRHS, oRHS.pVertex[0], oRHS.pVertex[1] ) ;
      bool bQueryV12 = CheckSeparatingAxis( *this, oRHS, oRHS.pVertex[1], oRHS.pVertex[2] ) ;
      bool bQueryV20 = CheckSeparatingAxis( *this, oRHS, oRHS.pVertex[2], oRHS.pVertex[0] ) ;
      
      return ( ( ! bThisV01  && ! bThisV12 && ! bThisV20  )
               || ( ! bQueryV01 && ! bQueryV12 && ! bQueryV20 )  );
    }

    
    //------------------------------------
    //  Subdivide
    //
    //  Split current voxel into 4 new voxels
    //
    //------------------------------------
    template< typename VoxelT >
    void Subdivide( vector<VoxelT> & vChildren ) const
    {
      vChildren = vector<VoxelT>( 4, static_cast<const VoxelT > (*this)   );   // will this ever be problematic?
      // midpoint between v0 and v1
      SVector3 oV01 = MeshGeometry::MidPoint( pVertex[0], pVertex[1] );
      SVector3 oV12 = MeshGeometry::MidPoint( pVertex[1], pVertex[2] );
      SVector3 oV20 = MeshGeometry::MidPoint( pVertex[2], pVertex[0] );
      Float fNewSideLength = fSideLength / (Float) 2.0;
      
      // check which orientation THIS triangle is
      if( bVoxelPointsUp )
      {
        vChildren[0].Set( oV01, oV12, oV20, false, fNewSideLength, nGeneration + 1 );
        vChildren[1].Set( oV20, oV12, pVertex[2], true, fNewSideLength, nGeneration + 1  );
        vChildren[2].Set( pVertex[0], oV01, oV20, true, fNewSideLength, nGeneration + 1  );
        vChildren[3].Set( oV01, pVertex[1], oV12, true, fNewSideLength, nGeneration + 1  );
        
      }else{
        
        vChildren[0].Set( oV20, oV01, oV12, true, fNewSideLength, nGeneration + 1);
        vChildren[1].Set( oV20, pVertex[0], oV01, false, fNewSideLength, nGeneration + 1 );
        vChildren[2].Set( oV01, pVertex[1], oV12, false, fNewSideLength, nGeneration + 1 );
        vChildren[3].Set( oV20, oV12, pVertex[2], false, fNewSideLength, nGeneration + 1 );
      }

    }
    
    //------------------------------------
    //  Subdivide
    //  Divide current voxel to generation specified
    //------------------------------------
    template< typename VoxelT >
    void Subdivide( vector<VoxelT> & vChildren, Int nMaxGen ) const
    {
      if( nMaxGen <= nGeneration )
        return;
      else
      {
        vector<VoxelT> vNewChildren;
        Subdivide( vNewChildren );
        if( nMaxGen - 1 == nGeneration )
        {
          vChildren.insert( vChildren.end(), vNewChildren.begin(), vNewChildren.end() );
        }
        else
        {
          for( Int i = 0; i < 4; i ++ )
            vNewChildren[i].Subdivide( vChildren, nMaxGen );
        }
      }
    }
  };
  
};

#endif
