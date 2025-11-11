/////////////////////////////////////////////////////////////////
//
//  File:     MicMesh.tmpl.cpp
//  Authors:  Frankie Li
//  e-mail:   sfli@cmu.edu
//
//
//  Implementations of templated MicMesh functions and classes
//
//
/////////////////////////////////////////////////////////////////


namespace XDMUtility
{
  //----------------------------------------------------------------------
  //  ResetErrorLimit
  //----------------------------------------------------------------------
  template< class ShapeContainer >
  void CGeneralShapeLocator<ShapeContainer>::
  UpdateErrorLimit( const BBox2D oBoundingBox )
  {
    Float fDX = oBoundingBox.pMax.x - oBoundingBox.pMin.x;
    Float fDY = oBoundingBox.pMax.y - oBoundingBox.pMin.y;
    Float fNewError = std::max( fDX, fDY );
    fNewError /= 10;
    fNewError = std::max( std::numeric_limits<Float>::epsilon(), fNewError );
    FErrorEstimator.fError = std::max( FErrorEstimator.fError, fNewError ); 
  }

  //-------------------------------------------------------------------------------
  //    Insert
  //-------------------------------------------------------------------------------
  template< class ShapeContainer >
  typename CGeneralShapeLocator<ShapeContainer>::ShapePtr
  CGeneralShapeLocator<ShapeContainer>::Insert( const ShapeT & oShape )
  {
    UpdateErrorLimit( oShape.GetBoundingBox() );
    Float fZ = oShape.Vertex(0).m_fZ;
    typename ShapeTree3D::iterator pIter = oObjLocator.find( fZ );
    if( pIter == oObjLocator.end() )
    {
      ShapeTreePtr pNewTree = ShapeTreePtr( new ShapeTreeT( oPlanarBBox, nMaxDepth, FErrorEstimator() ) );
      bool bSuccess;
      std::tie( pIter, bSuccess) = oObjLocator.insert( std::pair<Float, ShapeTreePtr>( fZ, pNewTree ) );
    }
    return pIter->second->Insert( oShape );
  }

  //-------------------------------------------------------------------------------
  //    Insert
  //-------------------------------------------------------------------------------
  template< class ShapeContainer >
  typename CGeneralShapeLocator<ShapeContainer>::ShapePtr
  CGeneralShapeLocator<ShapeContainer>::Insert( const ShapePtr & pShape )
  {
    UpdateErrorLimit( pShape->GetBoundingBox() );
    Float fZ = pShape->Vertex(0).m_fZ;
    typename ShapeTree3D::iterator pIter = oObjLocator.find( fZ );
    if( pIter == oObjLocator.end() )
    {
      ShapeTreePtr pNewTree = ShapeTreePtr( new ShapeTreeT( oPlanarBBox, nMaxDepth, FErrorEstimator() ) );
      bool bSuccess;
      std::tie( pIter, bSuccess) = oObjLocator.insert( std::pair<Float, ShapeTreePtr>( fZ, pNewTree ) );
    }
    pIter->second->Insert( pShape );
    return pShape;
  }

  //-------------------------------------------------------------------------------
  //    Find
  //-------------------------------------------------------------------------------
  template< class ShapeContainer >
  vector< typename CGeneralShapeLocator<ShapeContainer>::ShapePtr >
  CGeneralShapeLocator<ShapeContainer>::Find( const BBox3D & oBox ) const
  {
    Float fMinZ = oBox.m_oBoxMin.m_fZ;
    Float fMaxZ = oBox.m_oBoxMax.m_fZ;
    
    ConstLayerIter pLower = oObjLocator.lower_bound( fMinZ );
    ConstLayerIter pUpper = oObjLocator.upper_bound( fMaxZ );   // inclusive bound [ fMinZ, fMax ]
    
    BBox2D oLayerBox( Point( oBox.m_oBoxMin.m_fX, oBox.m_oBoxMin.m_fY ),
                      Point( oBox.m_oBoxMax.m_fX, oBox.m_oBoxMax.m_fY ) );
    
    vector< ShapePtr > oShapeList;
    for( ; pLower != pUpper; ++ pLower )
    {
      vector<ShapePtr> oLocalList = pLower->second->FindOverlap( oLayerBox );
      oShapeList.insert( oShapeList.end(), oLocalList.begin(), oLocalList.end() );
    }
    return oShapeList;
  }
  
  //-------------------------------------------------------------------------------
  //    Find
  //-------------------------------------------------------------------------------
  template< class ShapeContainer >
  vector< typename CGeneralShapeLocator<ShapeContainer>::ShapePtr >
  CGeneralShapeLocator<ShapeContainer>::Find( const BBox2D & oBox,
                                              Float fLayerLoc ) const
  {
    ConstLayerIter pIter = oObjLocator.find( fLayerLoc );
    vector<ShapePtr> oShapeList = pIter->second->FindOverlap( oBox );
    return oShapeList;
  }

  //-------------------------------------------------------------------------------
  //    GetNeighbors
  //    Purpose:  Get neighbors in 3D, i.e., from in layer and neighboring
  //              layers.
  //-------------------------------------------------------------------------------
  template< class ShapeContainer >
  vector< typename CGeneralShapeLocator<ShapeContainer>::ShapePtr >
  CGeneralShapeLocator<ShapeContainer>::GetNeighbors( const ShapePtr & pShape ) const
  {
    Float fCurZ =  pShape->Vertex(0).m_fZ;
    ConstLayerIter pCenter = oObjLocator.find( fCurZ );

    vector< ShapePtr > oShapeList;
     
    if( pCenter == oObjLocator.end() )
      return oShapeList;
   
    oShapeList = pCenter->second->GetNeighbors( *pShape );

    ConstLayerIter pUpper = pCenter;
    ++ pUpper;
    if( pUpper != oObjLocator.end() )
    {
      vector<ShapePtr> oLocalList = pUpper->second->FindOverlap( *pShape );
      oShapeList.insert( oShapeList.begin(), oLocalList.begin(), oLocalList.end() );
    }
        
    ConstLayerIter pLower = oObjLocator.begin();
    if( pLower == pCenter )
      return oShapeList;
    ConstLayerIter pNext = pLower;
    ++ pNext;
    while( pNext != pCenter )
    {
      ++ pNext;
      ++ pLower;
    }
    
    if( pLower != oObjLocator.end() )
    {
      vector<ShapePtr> oLocalList = pLower->second->FindOverlap( *pShape );
      oShapeList.insert( oShapeList.begin(), oLocalList.begin(), oLocalList.end() );
    }
    return oShapeList;
  }


  
} // end namespace
