/////////////////////////////////////////////////////////////////
//
//  File:    MicMesh.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Implementation of MicMesh.h
//
//
/////////////////////////////////////////////////////////////////


#include "MicMesh.h"


namespace XDMUtility
{
  
  //----------------------------------------------------------------------
  //  ResetErrorLimit
  //----------------------------------------------------------------------
  void CMicMesh::UpdateErrorLimit( const BBox oBoundingBox )
  {
    Float fDX = oBoundingBox.pMax.x - oBoundingBox.pMin.x;
    Float fDY = oBoundingBox.pMax.y - oBoundingBox.pMin.y;
    Float fNewError = std::max( fDX, fDY );
    fNewError /= 10;
    fNewError = std::max( std::numeric_limits<Float>::epsilon(), fNewError );
    FErrorEstimator.fError = std::max( FErrorEstimator.fError, fNewError ); 
  }
  
  //----------------------------------------------------------------------
  //  Public: Insert
  //----------------------------------------------------------------------
  CMicMesh::ShapePtr CMicMesh::Insert( const ShapeT & oShape )
  {
    BBox oBoundingBox = oShape.GetBoundingBox();
    UpdateErrorLimit( oBoundingBox );
    ShapePtr pNewVertex = ShapePtr( new ShapeT( oShape ) );
    oShapeLocator.Add( pNewVertex, oBoundingBox );
    oShapePtrList.push_back( pNewVertex );
    return pNewVertex;
  }
  
  //----------------------------------------------------------------------
  //  Public: Insert
  //----------------------------------------------------------------------
  CMicMesh::ShapePtr CMicMesh::Insert( const ShapePtr & pShape )
  {
    BBox oBoundingBox = pShape->GetBoundingBox();
    UpdateErrorLimit( oBoundingBox );
    oShapeLocator.Add( pShape, oBoundingBox );
    oShapePtrList.push_back( pShape );
    return pShape;
  }

  //----------------------------------------------------------------------
  //  Delete
  //----------------------------------------------------------------------
  void CMicMesh::Delete( const ShapeT & oPoint )
  {
    SDeleteCmp FToDelete( oPoint );
    SVector3 oCenter = oPoint.GetCenter();
    oShapeLocator.Delete( Point( oCenter.m_fX, oCenter.m_fY ),
                          FToDelete );
    ShapePtrIter pLastToKeep =
      std::remove_if( oShapePtrList.begin(), oShapePtrList.end(), FToDelete );
    oShapePtrList.erase( pLastToKeep, oShapePtrList.end() );
  }

  //----------------------------------------------------------------------
  //  Refine
  //----------------------------------------------------------------------
  void CMicMesh::Refine( const ShapePtr & pVoxel )
  {
    SDeleteCmp FToDelete( pVoxel );
    
    ShapePtrIter pLastToKeep =
      std::remove_if( oShapePtrList.begin(), oShapePtrList.end(), FToDelete );
    oShapePtrList.erase( pLastToKeep, oShapePtrList.end() );
    RUNTIME_ASSERT( oShapePtrList.end() - pLastToKeep == 1,
                    "More than one voxel perfectly matching...  \n" );
    SVector3 oCenter = pVoxel->GetCenter();
    oShapeLocator.Delete( Point( oCenter.m_fX, oCenter.m_fY ),
                          FToDelete );

    vector<ShapeT> oNewVoxels;
    pVoxel->Subdivide( oNewVoxels );
    for( Size_Type i = 0; i < oNewVoxels.size(); i ++ )
      Insert( oNewVoxels[i] );
  }
  
  //-------------------------------------------------------------------------------
  //  FindOverlap
  //-------------------------------------------------------------------------------
  vector<CMicMesh::ShapePtr>
  CMicMesh::FindOverlap( const ShapeT &oCenter, Float fRelError ) const
  {
    vector<ShapePtr> oNgbs;
    BBox oCenterBB = oCenter.GetBoundingBox();
    Float fError = FErrorEstimator();      // extend Bounding box by epsilon of this mesh
    oCenterBB.pMin.x -= fError;            //  (unrelated to fRelError)
    oCenterBB.pMin.y -= fError;

    oCenterBB.pMax.x += fError;
    oCenterBB.pMax.y += fError;
    Find( oNgbs, oCenterBB );

    if( oNgbs.size() == 0 )
      return oNgbs;
    
    // need to filter out triangles that are not overlaping oCenter
    // don't want anything that's overlapping and connected and
    // the vertices are not in the center
    //
    //  ( as definition of connected is just side/vertex connected )
    //
    Float fMinSideLength = fError;
    for( Size_Type i = 0; i < oNgbs.size(); i ++ )
      fMinSideLength = std::min( oNgbs[i]->fSideLength , fMinSideLength );        
    
    Float fAbsError = fRelError * fMinSideLength;
    ShapeT oTestCenter = oCenter;
    oTestCenter.ResizeAndCenter( oCenter.fSideLength - fAbsError );

    std::sort( oNgbs.begin(), oNgbs.end(),
               XDMUtility::FShapePtrLessCmp<ShapePtr>() );
    
    vector<ShapePtr>::iterator
      pLast = std::unique( oNgbs.begin(),
                           oNgbs.end(),
                           XDMUtility::FShapePtrEqCmp<ShapePtr>() );
    oNgbs.erase( pLast, oNgbs.end() );
    pLast = std::remove_if( oNgbs.begin(), oNgbs.end(),
                            FDisconnectedCmp<ShapePtr> ( oTestCenter ) );
    oNgbs.erase( pLast, oNgbs.end() );
    return oNgbs;
  }
  
  //-------------------------------------------------------------------------------
  //  FindOverlap
  //-------------------------------------------------------------------------------
  vector<CMicMesh::ShapePtr>
  CMicMesh::FindOverlap( const BBox & oCenter ) const
  {
    vector<ShapePtr> oNgbs;
    Find( oNgbs, oCenter );
    return oNgbs;
  }
  
  //-------------------------------------------------------------------------------
  //  GetNeighbors
  //-------------------------------------------------------------------------------
  vector<CMicMesh::ShapePtr>
  CMicMesh::GetNeighbors( const ShapeT &oCenter ) const
  {
    vector<ShapePtr> oNgbs;
    BBox oCenterBB = oCenter.GetBoundingBox();
    Float fError = FErrorEstimator();          // extend Bounding box by epsilon
    oCenterBB.pMin.x -= fError;
    oCenterBB.pMin.y -= fError;
    oCenterBB.pMax.x += fError;
    oCenterBB.pMax.y += fError;
    
    vector< ShapePtr > oCandidates = FindOverlap( oCenterBB );
    
    for( Size_Type i = 0; i < oCandidates.size(); i ++ )   // have to exclude self
    {
      SVoxel oTmp = *oCandidates[i];
      oTmp.Vertex(0).m_fZ = oCenter.Vertex(0).m_fZ;
      oTmp.Vertex(1).m_fZ = oCenter.Vertex(1).m_fZ;
      oTmp.Vertex(2).m_fZ = oCenter.Vertex(2).m_fZ;
      
      if( oCenter.Connected( oTmp, FErrorEstimator ) )  // look for edge connection
        if( ! oTmp.PerfectMatch( oCenter ) )
          oNgbs.push_back( oCandidates[i] );
    }
    return oNgbs;
  }
  
  //-------------------------------------------------------------------------------
  //  Find
  //  Precondition:  vertices of Voxel oCenter must already exist.
  //  Purpose: GetNeighbors by voxel
  //-------------------------------------------------------------------------------
  void CMicMesh::Find( vector<CMicMesh::ShapePtr> & oRes, const BBox & oSearchBox ) const
  {
    FNgbLookup ProcessFn;
    oShapeLocator.RangeSearch( oSearchBox, ProcessFn );
    oRes.reserve( ProcessFn.oNgbList.size() );
    for( Size_Type i = 0; i < ProcessFn.oNgbList.size(); i ++ )
    {
      FDuplicateCheck<ShapePtr> oEqCmp( ProcessFn.oNgbList[i] );
      vector<ShapePtr>::iterator pDuplicate = std::find_if( oRes.begin(),
                                                            oRes.end(),
                                                            oEqCmp );
      if( pDuplicate == oRes.end() )
        oRes.push_back( ProcessFn.oNgbList[i] );
    }
  }
  
  //-------------------------------------------------------------------------------
  //  Exists
  //
  //  Purpose:  Check to see if oCenter is in the tree already.  The existance
  //            requires that all of the vertices of the ShapeT to match.
  //-------------------------------------------------------------------------------
  bool CMicMesh::Exists( const ShapeT & oCenter ) const 
  {
    // not the smartest way - but will have a better method with the improved quadtree
    BBox oCenterBB = oCenter.GetBoundingBox();
    Float fError = FErrorEstimator();          // extend Bounding box by epsilon
    oCenterBB.pMin.x -= fError;
    oCenterBB.pMin.y -= fError;

    oCenterBB.pMax.x += fError;
    oCenterBB.pMax.y += fError;
    vector< ShapePtr > oCandidates;
    Find( oCandidates, oCenterBB );

    for( Size_Type i = 0; i < oCandidates.size(); i ++ )   // have to exclude self
    {
      SVoxel oTmp = *oCandidates[i];
      oTmp.Vertex(0).m_fZ = oCenter.Vertex(0).m_fZ;
      oTmp.Vertex(1).m_fZ = oCenter.Vertex(1).m_fZ;
      oTmp.Vertex(2).m_fZ = oCenter.Vertex(2).m_fZ;
      if( oTmp.PerfectMatch( oCenter ) )
        return true;
    }
    return false;
  }
  
}; // end namespace



