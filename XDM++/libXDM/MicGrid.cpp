/////////////////////////////////////////////////////////////////
//
//  File:    MicGrid.cp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//
/////////////////////////////////////////////////////////////////
#include "MicGrid.h"


namespace MicAnalysis
{
  
  //-------------------------------------------------------
  //  operator() ( nI, nJ )
  //-------------------------------------------------------
  CMicGrid::ShapePtr CMicGrid::operator() ( Int nI, Int nJ ) const
  {
    DEBUG_ASSERT( InRange( nI, nJ ), "Requested Triangle Out of Bound\n");
    return oGrid[ nI ][ nJ ];
  }
  
  //-------------------------------------------------------
  //  operator() ( SVector3 )
  //-------------------------------------------------------
  CMicGrid::ShapePtr CMicGrid::operator() ( const SVector3 & oPos ) const
  {
    Int nI, nJ;
    PositionToIndices( nI, nJ, oPos );
    return operator()( nI, nJ );
  }
  //-------------------------------------------------------
  //-------------------------------------------------------
  CMicGrid::ShapePtr CMicGrid::operator( ) ( const CMicGrid::ShapeT & v  ) const
  {
    Int nI, nJ;
    GetIndices( nI, nJ, v );
    return operator()( nI, nJ );
  }
  //-------------------------------------------------------
  //-------------------------------------------------------
  CMicGrid::ShapePtr CMicGrid::operator( ) ( const CMicGrid::ShapePtr & v  ) const
  {
    return operator()( *v );
  }

  //-------------------------------------------------------
  //  PositionToIndices
  //-------------------------------------------------------
  void CMicGrid::PositionToIndices( Int & nI, Int & nJ, const SVector3 & oLoc ) const
  {
    SVector3 oPos = oLoc - oOrigin;
    Float fI = Dot( oPos, oBasisI ) / ( fMinSideLength / Float(2) );   // gauranteed to be positve numbers 
    Float fJ = Dot( oPos, oBasisJ ) / ( fMinSideLength * sqrt(3)/Float(2) );
    
    // fI and fJ should be close to integer.  This is a convention
    nI = floor( fI + 0.1 );
    nJ = floor( fJ + 0.1 );

  }
  
  //-------------------------------------------------------
  //  GetIndicies
  //-------------------------------------------------------
  void CMicGrid::GetIndices( Int & nI, Int & nJ, const Shape & oShape ) const
  {
    SVector3 oShiftedCenter;  // this is important for the grid property
    
    Float fMinY = oShape.Vertex(0).m_fY;
    Float fMaxY = fMinY;
    for ( int i = 1; i < 3; i ++ )
    {
      fMinY = std::min( fMinY, oShape.Vertex(i).m_fY );
      fMaxY = std::max( fMaxY, oShape.Vertex(i).m_fY );
    }
    oShiftedCenter = oShape.GetCenter();
    oShiftedCenter.m_fX = oShape.GetCenter().m_fX;
    oShiftedCenter.m_fY = (fMinY + fMaxY) / static_cast<Float>( 2 );
    //    std::cout << " ---> " << (fMinY + fMaxY)/float(2) << " " <<  oShiftedCenter.m_fY << std::endl;
    PositionToIndices( nI, nJ, oShiftedCenter );
  }

  //-------------------------------------------------------
  //  Insert
  //-------------------------------------------------------
  void CMicGrid::Insert( const Mic & oMic )
  {
    InsertReplace( oMic.GetInitialSideLength(), oMic.VoxelListBegin(), oMic.VoxelListEnd() );
  }

  //-------------------------------------------------------
  //  InsertWithChildren
  //-------------------------------------------------------
  void CMicGrid::InsertChildren( const ShapePtr & pShape )
  {
    vector< Shape > vChildren;
    pShape->Subdivide( vChildren , nMaxGeneration );
    for( Size_Type i = 0; i < vChildren.size(); i ++ )
      Insert( vChildren[i] );
    
  }

  //-------------------------------------------------------
  //  Insert
  //-------------------------------------------------------
  CMicGrid::ShapePtr CMicGrid::Insert( const ShapePtr & pShape )
  {
    if( pShape->nGeneration < nMaxGeneration )
    {
      InsertChildren( pShape );
      return pShape;   // not really sensible at this point
    }
    
    Int nI, nJ;
    GetIndices( nI, nJ, *pShape );

    //-------------------------------------------------
    {

      
//       SVector3 oPos = oLoc - oOrigin;
//       Float fI = Dot( oPos, oBasisI ) / ( fMinSideLength / Float(2) );   // gauranteed to be positve numbers 
//       Float fJ = Dot( oPos, oBasisJ ) / ( fMinSideLength * sqrt(3)/Float(2) );
      
//       // fI and fJ should be close to integer.
//       nI = floor( fI + 0.1 );
//       nJ = floor( fJ + 0.1 );
    }
    
     std::stringstream ss;
     ss << nI << " " << nJ << " " << pShape->GetCenter() << std::endl
        << " Origin " << oOrigin << std::endl;    
    RUNTIME_ASSERT( InRange( nI, nJ ), "PositionToIndices: Triangle does not belong in the Mic file\n" + ss.str() );
    //    RUNTIME_ASSERT( InRange( nI, nJ ), "PositionToIndices: Triangle does not belong in the Mic file\n"  );
    if( oGrid[ nI ][ nJ ] != ShapePtr() )
    {
      std::cout << "Warning: Original Mic file has Duplicates, trying to insert: " << std::endl
                << pShape->Vertex(0) << std::endl
                << pShape->Vertex(1) << std::endl
                << pShape->Vertex(2) << std::endl;
      std::cout << "Original: " << std::endl
                << oGrid[ nI ][ nJ ]->Vertex(0) << std::endl
                << oGrid[ nI ][ nJ ]->Vertex(1) << std::endl
                << oGrid[ nI ][ nJ ]->Vertex(2) << std::endl;
    }
    oGrid[ nI ][ nJ ] = pShape;
    return pShape;
  }

  //-------------------------------------------------------
  //  Insert
  //-------------------------------------------------------
  CMicGrid::ShapePtr CMicGrid::Insert( const Shape & oShape )
  {
    ShapePtr pShape = ShapePtr( new Shape( oShape ) );
    return Insert( pShape );
  }
  
  //-------------------------------------------------------
  //  GetNeighbors (formally, only first neighbors)
  //-------------------------------------------------------
  vector<CMicGrid::ShapePtr> CMicGrid::GetNeighbors( const Shape & oShape ) const
  {
    Int nI, nJ;
    vector<ShapePtr> oNgbs;
    GetIndices( nI, nJ, oShape );
    if( ! InRange(nI, nJ) ||  ! IsValid( oGrid[ nI ][ nJ ] ) )
      return oNgbs;
    
    DEBUG_ASSERT( InRange( nI, nJ ), "PositionToIndices: Triangle does not belong in the Mic file\n" );
    if( InRange( nI - 1, nJ ) && IsValid( oGrid[ nI - 1 ][ nJ ] ) )
      oNgbs.push_back( oGrid[ nI - 1 ][ nJ ] );
    if( InRange( nI + 1, nJ ) && IsValid( oGrid[ nI + 1 ][ nJ ] ) )
      oNgbs.push_back( oGrid[ nI + 1 ][ nJ ] );
    
    if( oShape.bVoxelPointsUp )
    {
      // neighbors are ( x - 1, y ), ( x + 1, y ), ( x, y - 1 )
      if( InRange( nI , nJ - 1 ) && IsValid( oGrid[ nI ][ nJ - 1 ] ) )
        oNgbs.push_back( oGrid[ nI ][ nJ - 1 ] );
    }
    else
    {
      // neighbors are ( x - 1, y ), ( x + 1, y ), ( x, y + 1)
      if( InRange( nI, nJ + 1 )  && IsValid( oGrid[ nI ][ nJ + 1 ] ) )
        oNgbs.push_back( oGrid[ nI ][ nJ + 1 ] );
    }
    
    return oNgbs;
  }

  //-------------------------------------------------------
  //  GetNeighbors (formally, only first neighbors)
  //-------------------------------------------------------
  vector<CMicGrid::ShapePtr> CMicGrid::GetNeighbors( const ShapePtr & pShape ) const  
  {
    return GetNeighbors( *pShape );
  }

  //-----------------------------------
  //  FindOverlap(const &SVoxel)
  //
  //  Purpose:  Find all voxels within the voxel tree that overlaps the voxel, oCenter 
  //            All (Data) voxels in oRes must overlap oCenter.
  //
  //-----------------------------------
  CMicGrid::ShapePtrList CMicGrid::FindOverlap( const BBox & oCenter ) const
  {
    ShapePtrList oRes;
    Find( oRes, oCenter );
    return oRes;
  }
  
  //-----------------------------------
  //  FindOverlap
  //  Purpose:  Find all overlapping voxel with oCenter, given a relative
  //            error limit.  (Fractional)   Note that this limits the resolution
  //            to fRelError * fSideLength
  //-----------------------------------
  CMicGrid::ShapePtrList
  CMicGrid::FindOverlap( const ShapeT & oCenter, Float fRelError ) const
  {
    using namespace boost::lambda;
    vector<ShapePtr> oNgbs;
    BBox oCenterBB = oCenter.GetBoundingBox();
    Find( oNgbs, oCenterBB );
    
    if( oNgbs.size() == 0 )
      return oNgbs;
    
    Float fAbsError = fRelError * oCenter.fSideLength;
    ShapeT oTestCenter = oCenter;
    oTestCenter.ResizeAndCenter( oCenter.fSideLength - fAbsError );
    
    vector<ShapePtr>::iterator pLast
      = std::remove_if( oNgbs.begin(), oNgbs.end(),
                        XDMUtility::FDisconnectedCmp<ShapePtr> ( oTestCenter ) );
    oNgbs.erase( pLast, oNgbs.end() );
    return oNgbs;
  }

  //-----------------------------------
  // Find
  //-----------------------------------
  void CMicGrid::Find( vector<ShapePtr> & oRes, const BBox & oSearchBox ) const
  {
    SVector3 oMin( oSearchBox.pMin.x, oSearchBox.pMin.y, 0 );
    SVector3 oMax( oSearchBox.pMax.x, oSearchBox.pMax.y, 0 );
    
    Int nMinI, nMaxI, nMinJ, nMaxJ;
    
    PositionToIndices( nMinI, nMinJ, oMin );
    PositionToIndices( nMaxI, nMaxJ, oMax );

    for( Int i = nMinI; i <= nMaxI; i ++ )
      for( Int j = nMinJ; j <= nMaxJ; j ++ )
        if( InRange( i, j ) && IsValid( oGrid[i][j] ) )
          oRes.push_back( oGrid[i][j] );
  }
  
  //-----------------------------------------------
  //
  //  GetMic
  //
  //-----------------------------------------------
  CMicGrid::Mic CMicGrid::GetMic() const
  {
    Mic oRes;
    oRes.SetInitialSideLength( fInitialSideLength );
    
    for( Int i = 0; i < Size1(); i ++ )
      for( Int j = 0; j < Size2(); j ++ )
        if( IsValid( oGrid[i][j] ) )
          oRes.AddVoxel( *oGrid[i][j] );
    
    return oRes;
  }

  //---------------------------------------------------------------------------------------
  //  CalculateNyeTensorField
  //---------------------------------------------------------------------------------------
  TensorField3 CalculateNyeTensorField( const SMicGridVolume & oMicVolume,
                                        const SVector3 & oStepSize )
  {
    TensorField3 oRes;
    
    SVector3 oOrigin(1, 1, 1);
    oOrigin *= std::numeric_limits<Float>::max();
    Int nXSize = 0;
    Int nYSize = 0;
    
    for( Size_Type i = 0; i < oMicVolume.oGridList.size(); i ++ )
    {
      SVector3 oLayerOrigin = oMicVolume.oGridList[i].Origin();
      oOrigin.m_fX = std::min( oLayerOrigin.m_fX, oOrigin.m_fX );
      oOrigin.m_fY = std::min( oLayerOrigin.m_fY, oOrigin.m_fY );
      oOrigin.m_fZ = std::min( oLayerOrigin.m_fZ, oOrigin.m_fZ );
      nXSize       = std::max( nXSize, oMicVolume.oGridList[i].Size1() );
      nYSize       = std::max( nYSize, oMicVolume.oGridList[i].Size2() );
    }

    oRes.resize( boost::extents[ oMicVolume.oGridList.size() ][ nXSize ][ nYSize ] );

    // z varying slowest
    for( Int k = 0; k < Int( oMicVolume.oGridList.size() ); k ++ )
    {
      for( Int i = 2; i < nXSize - 2; i ++ )
      {
        for( Int j = 2; j < nYSize - 2; j ++ )
        {
          SVector3 oPos = oOrigin;
          oPos.m_fX += oStepSize.m_fX * Float( i );
          oPos.m_fY += oStepSize.m_fY * Float( j );
          oPos.m_fZ += oStepSize.m_fZ * Float( k );
          Int oIndices[3];
          oIndices[0] = i;
          oIndices[1] = j;
          oIndices[2] = k;
          
          oRes[k][i][j] = Details::CalculateNyeTensor( oMicVolume, oIndices, oStepSize );
        }
      }
    }
    return oRes;
  }

  //------------------------------------------------------------
  //  Finite difference
  //
  //  This is method one - literally calculating the differential property
  //
  //    Method 2 is to calculate misorientation over distance
  //
  //------------------------------------------------------------
  Float SMicGridVolume::Differentiate( Int oIndices[3],
                                       const SVector3 & oDiffStepSize,
                                       Int nSpatialInd, Int l, Int m ) const
  {      
    int oStep[3]; oStep[0] = 0; oStep[1] = 0; oStep[2] = 0;
    oStep[ nSpatialInd ] = 1;

    int oUpIndices[3];
    int oDownIndices[3];
    for( int i = 0; i < 3; i ++ )
    {
      oUpIndices[i]   = oIndices[i] + oStep[i];
      oDownIndices[i] = oIndices[i];
    }

    if( ! InRange( oUpIndices ) || ! InRange( oDownIndices ) )
      return 0;
    
    ShapePtr pDown = oGridList[ oDownIndices[2] ]( oDownIndices[0], oDownIndices[1] );
    ShapePtr pUp   = oGridList[ oUpIndices[2] ]  ( oUpIndices[0], oUpIndices[1] );

    // need to consider boost::optional
    if( pDown != ShapePtr() && pUp != ShapePtr() )
    {
      SVector3 oDiff = pUp->GetCenter() - pDown->GetCenter();
      Float fDelta = oDiff[ nSpatialInd ];
      //std::cout << oDiff << std::endl;
    
      if( fDelta < 0.001 )
        return 0;
          
      const SMatrix3x3 & g0 = pDown->oOrientMatrix; 
      const SMatrix3x3 & g1 = pUp->oOrientMatrix;
      return ( g1.m[l][m] - g0.m[l][m] ) / fDelta;
    }
    else
    {
      return 0;
    }
  }

  
 
  
  namespace Details
  {
    //---------------------------------------------------------------------------------------
    //
    //  CalculateNyeTensor
    //
    //  StepSize is currently not used due to lack of resampling
    //
    //---------------------------------------------------------------------------------------
    SMatrix3x3 CalculateNyeTensor( const SMicGridVolume & oMicVolume,
                                   Int oIndices[3], const SVector3 & oDiffStepSize )
    {
      SMatrix3x3 oAlpha;

      for( Int i = 0; i < 3; i ++ )
      {
        for( Int j = 0; j < 3; j ++ )
        {
          Int l = (i + 1) % 3;
          Int k = (i + 2) % 3;
          oAlpha.m[i][j] = oMicVolume.Differentiate( oIndices, oDiffStepSize, k, l, j ) 
                         - oMicVolume.Differentiate( oIndices, oDiffStepSize, l, k, j );
        }
      }
      return oAlpha;
    }



  }
  
}
