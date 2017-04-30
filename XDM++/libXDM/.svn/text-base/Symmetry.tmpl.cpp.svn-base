////////////////////////////////////////////////////////////////
//
//  File:    Symmetry.tmpl.cpp
//  Authors:  Frankie Li, Jon Lind
//  e-mail:  sfli@cmu.edu
//  
//
// Implementation of template functions used in Symmetry
//
/////////////////////////////////////////////////////////////////


namespace LatticeSymmetry
{

  //--------------------------------------------------------------------------------------------------------
  //
  //  GetMisorientation
  ///   Calculate the misorientation between input oEulerOrient1 and oEulerOrient2
  //
  //--------------------------------------------------------------------------------------------------------
  template< class CSymmetryType >
  Float GetMisorientation( const CSymmetryType & oSymOps, 
                           const SMatrix3x3 &m1, 
                           const SMatrix3x3 &m2 )
  {
    SQuaternion q1;
    SQuaternion q2;
    q1.Set( m1 );
    q2.Set( m2 );
    Float fMinAngle = GetMisorientation( oSymOps, q1, q2 );
    return fMinAngle;
  }
  
  //--------------------------------------------------------------------------------------------------------
  //
  //  GetMisorientation
  ///   Calculate the misorientation between input oEulerOrient1 and oEulerOrient2
  //
  //--------------------------------------------------------------------------------------------------------
  template< class CSymmetryType >
  Float GetMisorientation( const CSymmetryType & oSymOps, 
                           const SVector3 &oEulerOrient1, 
                           const SVector3 &oEulerOrient2 )
  {

    SMatrix3x3 m1, m2;
    m1.BuildActiveEulerMatrix( oEulerOrient1.m_fX, oEulerOrient1.m_fY, oEulerOrient1.m_fZ );
    m2.BuildActiveEulerMatrix( oEulerOrient2.m_fX, oEulerOrient2.m_fY, oEulerOrient2.m_fZ );
    
    Float fMinAngle = GetMisorientation( oSymOps, m1, m2 );
    return fMinAngle ;
  }

  //--------------------------------------------------------------------------------
  //
  //  Misorientation  (Quanterion)
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Float GetMisorientation( const CSymmetryType & oSymOps, 
                           const SQuaternion &q1, 
                           const SQuaternion &q2 )
  {
    Float fAngle;
    SQuaternion qProd = q1.Inverse() * q2;
    qProd = ReduceToFundamentalZone( oSymOps, qProd ) ;
    
    fAngle = Float( 2 ) * acos( std::min( Float( 1 ), qProd.m_fW ) );
    return fAngle;
  }
  
  //--------------------------------------------------------------------------------
  //
  //  Symmetry reduction
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  SQuaternion ReduceToFundamentalZone( const CSymmetryType & oSymOps, const SQuaternion & q )
  {

    SQuaternion oRes = q;
    Float fMaxAbsCosComponent = fabs(q.m_fW); 

    const vector<SQuaternion> &vOperatorList = oSymOps.GetQuatOperatorList();
    
    for ( Size_Type i = 0; i < vOperatorList.size(); i ++ )
    {
      SQuaternion oTmp = q * vOperatorList[i];

      if ( fabs( oTmp.m_fW ) > fMaxAbsCosComponent )
      {
        oRes = oTmp;
        fMaxAbsCosComponent = fabs( oRes.m_fW );
      }
    }
    return oRes;
  }

  //--------------------------------------------------------------------------------
  //  IsInFundamentalZone
  //  Return true if the quaternion specified is in the fundamental zone.
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Bool IsInFundamentalZone( const CSymmetryType & oSym, const SQuaternion & q )
  {
    SQuaternion oQFz = ReduceToFundamentalZone( oSym, q );
    SQuaternion oQDiff = q - oQFz;

    // check for difference against something reasonable
    if ( oQDiff.EuclideanNorm() < 0.01 )
    {
      return true;
    }
    return false;
  }

  //--------------------------------------------------------------------------------
  //
  //  Return true of the two vectors are equivilent under symmetry operations
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Bool Equivilent( const CSymmetryType & oSymOps, const SVector3 &v1, const SVector3 &v2, Float fError )
  {
    const vector<SMatrix3x3> &vOperatorList = oSymOps.GetOperatorList();
    
    for( Size_Type i = 0; i < vOperatorList.size(); i ++ )
    {
      SVector3 oTmp = vOperatorList[i] * v1;
      SVector3 oDiff = oTmp - v2;
      Float fDiff = oDiff.GetLength();
      
      if( fDiff < fError )
      {
        return true;
      }
    }
    return false;
  }

  //--------------------------------------------------------------------------------
  //
  //  Reducing the list of vectors into the unique set
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  vector<SVector3> GetUniqueVectors( const CSymmetryType & oSymOps, const vector<SVector3> &oVectorList)
  {
    vector<SVector3> oUniqueList;
    
    for( Size_Type i = 0; i < oVectorList.size(); i ++ )
    {
      SVector3 oCurrentVec = oVectorList[i];
      Bool bUnique = true;
      for( Size_Type j = 0; j < oUniqueList.size(); j ++ )
      {      
        if ( Equivilent( oSymOps, oCurrentVec, oUniqueList[ j ] ) )
        {
          bUnique = false;
        }
      }
      if ( bUnique )
        oUniqueList.push_back( oCurrentVec );
    }
    
    return oUniqueList;
  }
  
  //--------------------------------------------------------------------------------
  //
  //  Reducing the list of recipocal vectors into the unique set
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  void GetUniqueRecpVectors(  vector<CRecpVector> &oUniqueList,
                              const CSymmetryType & oSymOps, const vector< CRecpVector > &oVectorList )
  {
    
    for( Size_Type i = 0; i < oVectorList.size(); i ++ )
    {
      CRecpVector oCurrentVec = oVectorList[i];
      Bool bUnique = true;
      for( Size_Type j = 0; j < oUniqueList.size(); j ++ )
      {      
        if ( Equivilent( oSymOps, oCurrentVec.v, oUniqueList[ j ].v ) )
        {
          bUnique = false;
        }
      }
      if ( bUnique )
        oUniqueList.push_back( oCurrentVec );
    }
  }

  //--------------------------------------------------------------------------------
  //
  //  Average  (repeated code so that it is less memory intensive)
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  SQuaternion Average( const CSymmetryType & oSym, const vector< SVector3 > &vAngleList )
  {
    return Average( oSym, vAngleList.begin(), vAngleList.end(),
                    Utilities::AngleToQuaternion() );
  }
  
  //--------------------------------------------------------------------------------
  //
  //  Average  (Quaternion)
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  SQuaternion Average( const CSymmetryType & oSym, const vector< SQuaternion > &vQuatList )
  {
    return  Average( oSym, vQuatList.begin(), vQuatList.end(),
                     Utilities::QuaternionToQuaternion() );
  }

  //--------------------------------------------------------------------------------
  //  Generalized Average
  //  Author:  Jon Lind
  //  Added quaternion cloud sorting that was missing in the previous version somehow...
  //  
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType, typename Iterator, typename PropertyMap >
  SQuaternion Average( const CSymmetryType & oSym, Iterator pBegin, Iterator pEnd,
                       PropertyMap oMap )
  {
    vector<SQuaternion> oSymOpList = oSym.GetQuatOperatorList();
    vector<SQuaternion> oCloudAverage( oSymOpList.size()  );
    for( Size_Type n = 0; n < oSymOpList.size(); n ++ )
    {
      oCloudAverage[n] = oMap.GetQuaternion( *pBegin ) * oSymOpList[n];
      oCloudAverage[n].ToConvention();      
    }

    vector<Int> oCloudPointsList( oSymOpList.size(), 1 );
    for( Iterator pCur = (pBegin + 1); pCur != pEnd; ++ pCur )
    {
      for( Size_Type n = 0; n < oSymOpList.size(); n ++ )
      {
        Int nCloudIndex = 0;
        SQuaternion q = oMap.GetQuaternion( *pCur ) * oSymOpList[ n ];;
        Float fMinDist = 200;
        for( Size_Type m = 0; m < oSymOpList.size(); m ++ )  // find cloud
        {
          SQuaternion qProd = q.Inverse() * oCloudAverage[m];
          Float fDist =  Float( 2 ) * acos( std::min( Float( 1 ), qProd.m_fW ) );
          if( fDist < fMinDist )
          {
            fMinDist = fDist;
            nCloudIndex = m;
          }
        }

        Int nPoints = oCloudPointsList[ nCloudIndex ];
        SQuaternion & qAve = oCloudAverage[ nCloudIndex ];
        q.ToConvention();
        
        qAve.m_fX = ( ( qAve.m_fX * nPoints  ) + q.m_fX );
        qAve.m_fY = ( ( qAve.m_fY * nPoints  ) + q.m_fY );
        qAve.m_fZ = ( ( qAve.m_fZ * nPoints  ) + q.m_fZ );
        qAve.m_fW = ( ( qAve.m_fW * nPoints  ) + q.m_fW );
        qAve = qAve / qAve.EuclideanNorm();
        qAve.ToConvention();
       oCloudPointsList[ nCloudIndex ] ++;
      }
      
    }
    
    SQuaternion oOrientationAverage;
    oOrientationAverage.Set( 0, 0, 0, 0 );
    for( Size_Type n = 0; n < oSymOpList.size(); n ++ )
    {
      oCloudAverage[n] = ReduceToFundamentalZone( oSym, oCloudAverage[n] );
      oCloudAverage[n].ToConvention();
      oOrientationAverage += oCloudAverage[n];
    }

    oOrientationAverage = oOrientationAverage / oOrientationAverage.EuclideanNorm();
    
    return oOrientationAverage;
  }

  
} // end namespace
