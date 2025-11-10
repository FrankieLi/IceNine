////////////////////////////////////////////////////////////////
//
//  File:    Sampling.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Implementation of Sampling methods
//
/////////////////////////////////////////////////////////////////

#include "Sampling.h"

namespace GeneralLib
{
  //--------------------------------------------------------------------------------------------------------
  //
  //  namespace UniformGrid
  //
  //--------------------------------------------------------------------------------------------------------
  namespace UniformGrid
  {

  //------------------------------------------------------------------------
  //
  //  MakeSukarevGridPoints
  //
  //  Parameter:  nLevel specifies the resolution, where resolution is 
  //  specified by a L_infinite metric. (i.e., bounding box.)  nLevel of 0
  //  implies that there is exactly 1 point in the box.  Given a bounding
  //  box of length L each side, then the L_infinite 'radius' of nLevel = 0
  //  is 2.5.  Note that the L_2 distance is sqrt(3) * r_infinite   
  //
  //  Note that the bounding box is always of length L = 1 in our example
  //
  //  TODO:  Make this into something that takes an approximate resolution
  //
  //  Grid samples ( [0, 1], [0, 1], [0, 1] ) x fSideLength
  //
  //------------------------------------------------------------------------
  vector<SVector3> MakeSukarevGridPoints( Int nLevel, Float fSideWidth )
  {
    vector<SVector3> oGrid;
    Int nDiv = pow( 2, nLevel );
    Float fScale = fSideWidth / Float( nDiv );

    SVector3 oXGen( fScale, 0, 0 );
    SVector3 oYGen( 0, fScale, 0 );
    SVector3 oZGen( 0, 0, fScale );

    SVector3 oOffset( 0.5 * fScale , 0.5 * fScale, 0.5 * fScale );
    Int nPointsPerAxis = round( pow( 2, nLevel ) );

    for ( Int nX = 0; nX < nPointsPerAxis; nX ++  )
    {
      for ( Int nY = 0; nY < nPointsPerAxis; nY ++  )
      {
        for ( Int nZ = 0; nZ < nPointsPerAxis; nZ ++  )
        {
          SVector3 oNewPoint = nX * oXGen + nY * oYGen + nZ * oZGen;
          oNewPoint += oOffset;
          oGrid.push_back( oNewPoint );
        }
      }
    }
    
    return oGrid;
  }


    
    //------------------------------------------------------------------------
    //
    //  GenerateGrid
    //
    //  A helper function for GetSukharevGridPoint
    //!-----------------------Sukharev grid sequence-----------------------------//
    //! Written by Steven Lindemann in the context of----------------------------//
    //!"Incremental low-discrepancy lattice methods for motion planning"---------//
    //! by S. R. Lindemann and S. M. LaValle-------------------------------------//
    //! In Proc. IEEE Internation Conference on Robotics and Automation, 2003----//
    //------------------------------------------------------------------------
    SVector3 GetLayeredSukharevGridPoint( Int nState, Int nDim )
    {
      //calculate offset based on index
      nState++;
      Int level = 0;
      Float offset  = 0.5;
      Float sampIndex = pow( 2.0, nDim * level);
      
      while ( nState > sampIndex )
      {
        nState -= ( Int ) sampIndex;
        offset = offset / 2.0;
        level++;
        sampIndex = pow( 2.0, nDim * level);
      }
      nState--;   //now i is the index into the grid, offset is the offset
            
      SVector3 blah;
      for (int j = 0; j < nDim; j++)
      {
        blah[j] = offset;
      }
      
      SVector3 sample = blah + GetSukarevGridPoint( nState, nDim );
      return sample;
    }
    
    
    //------------------------------------------------------------------------
    //
    //  GenerateGrid
    //
    //  A helper function for GetLayeredSukharevGridPoint
    //
    //!-----------------------Sukharev grid sequence-----------------------------//
    //! Written by Steven Lindemann in the context of----------------------------//
    //!"Incremental low-discrepancy lattice methods for motion planning"---------//
    //! by S. R. Lindemann and S. M. LaValle-------------------------------------//
    //! In Proc. IEEE Internation Conference on Robotics and Automation, 2003----//
    //
    //------------------------------------------------------------------------
    SVector3 GetSukarevGridPoint( Int nState, Int nDim )
    {
      
      int numberCells = nVertices;
      SVector3 sample(0, 0, 0);
      Float currentFactor = 0.5;
      Int currentIndex = nState % numberCells;
      Int blah = nState / numberCells;
      
      while (blah > 0)
      {
        for (Int j = 0; j < nDim; j++)
        {
          sample[j] += currentFactor * VertexOrderGrayCode[currentIndex][j];
        }
        currentFactor *= 0.5;
        currentIndex = blah % numberCells;
        blah = blah / numberCells;
      }
      
      for (Int j = 0; j < nDim; j++)
      {
        sample[j] += currentFactor * VertexOrderGrayCode[currentIndex][j];
      }
      
      return sample;
    }

    //----------------------------------------------------------------------------------------------
    //
    //
    //                              C Q u a t e r n i o n G r i d
    //
    //
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    // CQuaternionGrid
    //----------------------------------------------------------------------------------------------
    CQuaternionGrid::CQuaternionGrid()
    {
      InitializeFaces();
    }
    
    //----------------------------------------------------------------------------------------------
    //  Private:  InitializeFaces
    //----------------------------------------------------------------------------------------------
    void CQuaternionGrid::InitializeFaces()
    {
      
      // HyperFace 0
      pHyperFaceVertices[0][0].Set(  0.5, -0.5, -0.5, -0.5 ); 
      pHyperFaceVertices[0][1].Set(  0.5,  0.5, -0.5, -0.5 ); 
      pHyperFaceVertices[0][2].Set(  0.5, -0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[0][3].Set(  0.5,  0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[0][4].Set(  0.5, -0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[0][5].Set(  0.5,  0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[0][6].Set(  0.5, -0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[0][7].Set(  0.5,  0.5,  0.5,  0.5 );
    
      // HyperFace 1
      pHyperFaceVertices[1][0].Set( -0.5,  0.5, -0.5, -0.5 ); 
      pHyperFaceVertices[1][1].Set(  0.5,  0.5, -0.5, -0.5 ); 
      pHyperFaceVertices[1][2].Set( -0.5,  0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[1][3].Set(  0.5,  0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[1][4].Set( -0.5,  0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[1][5].Set(  0.5,  0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[1][6].Set( -0.5,  0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[1][7].Set(  0.5,  0.5,  0.5,  0.5 );
    
      // HyperFace 2
      pHyperFaceVertices[2][0].Set( -0.5, -0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[2][1].Set(  0.5, -0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[2][2].Set( -0.5, -0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[2][3].Set(  0.5, -0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[2][4].Set( -0.5,  0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[2][5].Set(  0.5,  0.5,  0.5, -0.5 ); 
      pHyperFaceVertices[2][6].Set( -0.5,  0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[2][7].Set(  0.5,  0.5,  0.5,  0.5 );
    
      // HyperFace 3
      pHyperFaceVertices[3][0].Set( -0.5, -0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[3][1].Set(  0.5, -0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[3][2].Set( -0.5, -0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[3][3].Set(  0.5, -0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[3][4].Set( -0.5,  0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[3][5].Set(  0.5,  0.5, -0.5,  0.5 ); 
      pHyperFaceVertices[3][6].Set( -0.5,  0.5,  0.5,  0.5 ); 
      pHyperFaceVertices[3][7].Set(  0.5,  0.5,  0.5,  0.5 ); 
    
    }

    
    
    //----------------------------------------------------------------------------------------------
    //  Public:  Interpolate
    //
    //  Using spherical linear interpolation to interpolate across the cube (hyper-face).
    //----------------------------------------------------------------------------------------------
    SQuaternion CQuaternionGrid::BarycentricToQuaternion( Int nFaceIndex, Float fAlpha, Float fBeta, Float fGamma) const
    {
      SQuaternion oRes;
    
      SQuaternion x1 = Interpolate( pHyperFaceVertices[nFaceIndex][0], pHyperFaceVertices[nFaceIndex][1], fAlpha );
      SQuaternion x2 = Interpolate( pHyperFaceVertices[nFaceIndex][2], pHyperFaceVertices[nFaceIndex][3], fAlpha );
      SQuaternion x3 = Interpolate( pHyperFaceVertices[nFaceIndex][4], pHyperFaceVertices[nFaceIndex][5], fAlpha );
      SQuaternion x4 = Interpolate( pHyperFaceVertices[nFaceIndex][6], pHyperFaceVertices[nFaceIndex][7], fAlpha );
    
      SQuaternion y1 = Interpolate( x1, x3, fBeta );
      SQuaternion y2 = Interpolate( x2, x4, fBeta );
    
      oRes = Interpolate(y1, y2, fGamma );

      oRes.ToConvention();
    
      return oRes;
    }

    //---------------------------------------------------------------------------------------
    //  Private: GetViolatingPointPairs
    //
    //  Purpose:  Given oQuat, a point in oOrientationGrid, find the nearest
    //            neighbor (with least misorientation) and return it.  If no such
    //            element exists (i.e., oSamplePoint does not violate the fMinMisorientation
    //            criterion), false is returned in the first element of the returned
    //            pair.
    //
    //  Parameters:  oSamplePoint - the sample point of interest that lies in the oOrientationGrid
    //               fMinMisorientation - minimum misorientation required in the grid criterion
    //
    //---------------------------------------------------------------------------------------
    std::pair<Bool, SQuaternion> CQuaternionGrid::GetViolatingPoint( const vector<SQuaternion>::iterator & pSamplePoint,
                                                                     const vector<SQuaternion> & oOrientationGrid,
                                                                     const LatticeSymmetry::CSymmetry & oSym,
                                                                     Float fMinMisorientCriterion )
    {
      using namespace LatticeSymmetry;
      SQuaternion oViolatingPoint;
      Bool bSamplePointViolation = true;
      Float fMinMisorientFound = Float(2.0) * PI;
      
      for( vector<SQuaternion>::const_iterator pGridPtr = oOrientationGrid.begin();
          pGridPtr != oOrientationGrid.end(); ++pGridPtr )
      {
        if( pGridPtr != pSamplePoint )
        {
          Float fCurrentMisorient = GetMisorientation( oSym, *pSamplePoint, *pGridPtr );

          if( fCurrentMisorient < fMinMisorientFound )
          {
            fMinMisorientFound = fCurrentMisorient;
            oViolatingPoint = *pGridPtr;
          }
        }
      }
      
      if( fMinMisorientFound > fMinMisorientCriterion )
        bSamplePointViolation = true;
      else
        bSamplePointViolation = false;
        
      return std::make_pair( bSamplePointViolation, oViolatingPoint );
    }

    //---------------------------------------------------------------------------------------
    //
    //  Private:   GetQuatPertubationRadius
    //
    //  Purpose:   Given oQuat, calculate the radius of perturbation required to generate
    //             random quaternion points needed to cover the area of fMisorientation.
    //
    //  Return:    The radius on S^3 that'd satistfy the condition above.
    //
    //  Parameters:  fMisorientation
    //               oQuat -- the orientation point
    //
    //---------------------------------------------------------------------------------------
    Float CQuaternionGrid::GetQuatPerturbationRadius( const SQuaternion & oCenter, Float fMisorientation )
    {
      Float fDelta;
      fDelta  =  Float( 2 ) * cos( fMisorientation / Float( 2)  ) + Float( 2 ) ;
      return fDelta;
    }

    //---------------------------------------------------------------------------------------
    //
    // PerturbQuaternionPoint
    //
    //---------------------------------------------------------------------------------------
    void CQuaternionGrid::PerturbQuatSampPoints( vector<SQuaternion> & oOrientationGrid, Float fDelta )
    {
      DEBUG_ASSERT(0, " NOT YET IMPLEMENTED ");
    }
    
    //---------------------------------------------------------------------------------------
    //  Private:  PatchOrientationGrid
    // 
    //  Purpose:  When we want just a grid in the fundamental zone under a certain symmetry,
    //            it is possible to break the min misorientation criterion.  This is because
    //            The fundamental zone boundary disallows some of the points from the sample.
    //            A simple solution is to preturb these points with large misorientation or
    //            to add new points between these points and their nearest neighbor(s).
    //
    //  Parameters:  oOrientationGrid  -- a set of sample points that we would like to patch
    //               to meet the minimum misorientation requirement.
    //               fMinMisorientation -- specification of minimum misorientation between points
    //
    //---------------------------------------------------------------------------------------
    void CQuaternionGrid::PatchOrientationGrid( vector<SQuaternion> & oOrientationGrid,
                                                const LatticeSymmetry::CSymmetry & oSym,
                                                Float fMinMisorientation )
    {
     
      for( Size_Type i = 0; i < oOrientationGrid.size(); i ++ )
      {
        Bool bSamplePointViolation;
        SQuaternion oViolatingNgb;
        vector<SQuaternion>::iterator pCurPointPtr = oOrientationGrid.begin()  + i;
        boost::tie( bSamplePointViolation, oViolatingNgb ) = GetViolatingPoint( pCurPointPtr, oOrientationGrid,
                                                                                oSym, fMinMisorientation );

        //----------------------------------------------------
        //  COMMENT:   Add a new point "half way" between the
        //             pair of points that are violating the
        //             minimum misorientation criterion.  Note
        //             that since most of these violating points
        //             are near the boundary of the fundamental
        //             zone, one must allow the interpolation to
        //             go across zone boundary. (otherwise, we'll
        //             end up interpolating between two points in
        //             the "wrong direction" in the 3-sphere.
        //-----------------------------------------------------
        if( bSamplePointViolation )
        {
          SQuaternion oBestPatch;
          Float fMinMisorientSum = MAX_FLOAT;
          const vector<SQuaternion> & oSymOpList=  oSym.GetQuatOperatorList();
          for ( Size_Type nOpIndex = 0; nOpIndex < oSymOpList.size(); nOpIndex ++ )
          {
            SQuaternion oSymEquivSampPoint = *pCurPointPtr * oSymOpList[ nOpIndex ];
            SQuaternion oMidPoint = Interpolate( oViolatingNgb, oSymEquivSampPoint, Float( 0.5 ) );

            Float fMisorientSum = GetMisorientation( oSym, oMidPoint, oSymEquivSampPoint ) +
                                  GetMisorientation( oSym, oMidPoint, oViolatingNgb );

            if( fMisorientSum < fMinMisorientSum )
            {
              fMinMisorientSum = fMisorientSum;
              oBestPatch = oMidPoint;
            }
          }
          
          oBestPatch =  LatticeSymmetry::ReduceToFundamentalZone( oSym, oBestPatch );
          oOrientationGrid.push_back( oBestPatch );
        }
      }
    }
    
    //---------------------------------------------------------------------------------------
    //
    //  GetGrid
    //
    //  Return a grid with nNumGridPoints uniformly distributed across S^3 with unit 1. (i.e., uniform
    //  rotation that's evenly spaced)
    //
    //  This function is from sample.C of 
    //  Anna Yershova and Steven M. LaValle, "Deterministic Sampling Methods for Spheres and SO(3)"
    //  
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> CQuaternionGrid::GetGrid( Int nNumGridPoints )
    {
      vector<SQuaternion> oOrientationGrid;
      for ( Int i = 0; i < nNumGridPoints; i ++ )
      {
        Int nFaceIndex = i % nNumHyperFaces;
        int nIndex = (Int) floor( (Float) i / (Float) nNumHyperFaces );
        SVector3 vInterp = UniformGrid::GetLayeredSukharevGridPoint( nIndex );
        SQuaternion oGridPoint = BarycentricToQuaternion( nFaceIndex,
                                                          vInterp.m_fX,
                                                          vInterp.m_fY,
                                                          vInterp.m_fZ );
        oOrientationGrid.push_back( oGridPoint );
      }
    
      return oOrientationGrid;
    }
  
    //---------------------------------------------------------------------------------------
    //
    //  GetGrid
    //
    //  Return a grid with uniformly distribued across S^3 such that no two points are more than
    //  nMaxDistance apart
    //
    //  Note:  fMaxDistance is the dispersion, or the radius of the largest empty ball in the sample
    //         grid.  Note also that fMaxDistance is the Eucledian distance in S^3
    //
    //         This is due to Proposition 4.5 of  Anna Yershova and Steven M. LaValle,
    //         "Deterministic Sampling Methods for Spheres and SO(3)"
    //
    //         d_rho( T ) <=  2 Pi / (n ( 2^d - 1 )/( 2 ( d + 1 )) + 1 ) ^(1/d )
    //
    //          n  ~ 1/2 * [(2 * pi / d_rho)^3 - 1 ]
    //
    //  TODO:  Make sure that d_rho is measured in angle (i.e., misorientation) rather than
    //         Euclidean norm (i.e., |q1 - q2|) or geodesic distance.
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> CQuaternionGrid::GetGrid( Float nMaxDistance )
    {
      vector<SQuaternion> oOrientationGrid;

      // Calculate number of points required based on criterion specified
      Int nNumPoints;
      nNumPoints = ceil( Float( 1 ) / Float( 2 ) * ( pow( Float(2) * PI / nMaxDistance, Float(3.) ) - Float( 1 ) ) );

      //      std::cout << nNumPoints << std::endl;
      oOrientationGrid = GetGrid( nNumPoints );
      return oOrientationGrid;
    }

    //---------------------------------------------------------------------------------------
    //  GetNearIdentityPoint
    //  Purpose:  Return a point of the quanterion space that's interpolated on the
    //            identity face around the origin.
    //---------------------------------------------------------------------------------------
    SQuaternion CQuaternionGrid::GetNearIdentityPoint( Float fX, Float fY, Float fZ ) const
    {
      const Int nIdentityFace = 0;
      SVector3 oOrigin(0.5, 0.5, 0.5);
      SVector3 oNewPos( fX, fY, fZ );
      oNewPos += oOrigin;
      SQuaternion oGridPoint = BarycentricToQuaternion( nIdentityFace,
                                                        oNewPos.m_fX,
                                                        oNewPos.m_fY,
                                                        oNewPos.m_fZ );
      return oGridPoint;
    }
    
    //---------------------------------------------------------------------------------------
    //
    //  GetRandomLocalGrid
    //
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> CQuaternionGrid::GetRandomLocalGrid( Float fMaxDistance, Int nGridPoints ) const
    {
      boost::mt19937 oRngEngine;
      RandomRealT oRandomReal( oRngEngine , boost::uniform_real<>( -fMaxDistance/2.0, fMaxDistance/2.0 ) );
      vector<SQuaternion> oLocalGrid;      
      for ( Int i = 0; i < nGridPoints; i ++ )
      {
        SQuaternion oGridPoint = GetNearIdentityPoint( oRandomReal(), oRandomReal(), oRandomReal() );
        oLocalGrid.push_back( oGridPoint );
      }
      return oLocalGrid;
    }
    
    //---------------------------------------------------------------------------------------
    //  GetRandomLocalPoint
    //---------------------------------------------------------------------------------------
    SQuaternion CQuaternionGrid::GetRandomLocalPoint( RandomRealT & oRandomReal ) const
    {
      return GetNearIdentityPoint( oRandomReal(), oRandomReal(), oRandomReal() );
    }
    
    //---------------------------------------------------------------------------------------
    //
    //  GetStructuredLocalGrid
    //
    //  Purpose:    A grid is to be generated with "side" of fSideLength in quaternion
    //              representation of SO(3).  Therefore, given a small fSideLength, we're
    //              generating a locally uniform, Eulerian grid.  Note that the distortion
    //              from this grid increases as a function of fSideLength.
    //
    //  Note:       To generate this local grid, we must convert fSideLength to a "arc length"
    //              on S^3.  This is done by mapping [0, pi) -> [0, 1) 
    //
    //  Parameter:  fSideLength is a distance defined by misorientation,
    //              or d = d(q1, q2) = acos(dot(q1, q2))
    //
    //              nLevel is the resolution, where grid points are expected to have neighbors
    //              approximately d/2^nLevel away.
    //
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> CQuaternionGrid::GetStructuredLocalGrid( Float fSideLength,
                                                                 Int nLevel )
    {
      // Project back from angular space to parametric space
      fSideLength = Float(2) * sin( PI/ Float(8)) * tan( fSideLength / Float(2) );        
      vector<SQuaternion> oLocalGrid;
      vector<SVector3> vOffsetList = MakeSukarevGridPoints( nLevel, fSideLength );
      SVector3 oSukarevGridCenter( fSideLength / Float(2), fSideLength / Float(2), fSideLength / Float(2) );
      SVector3 oOrigin( 0.5, 0.5, 0.5 );
      oOrigin -= oSukarevGridCenter;
      const Int nIdentityFace = 0;
      for ( Size_Type i = 0; i < vOffsetList.size(); i ++ )
      {
        SVector3 oNewPos = vOffsetList[i] + oOrigin;
        SQuaternion oGridPoint = BarycentricToQuaternion( nIdentityFace,
                                                          oNewPos.m_fX,
                                                          oNewPos.m_fY,
                                                          oNewPos.m_fZ );
        oLocalGrid.push_back( oGridPoint );
      }
      return oLocalGrid;
    }
    
    //---------------------------------------------------------------------------------------
    //
    //  GetMinimalFZGrid
    //
    //  Return minimal uniform grid distributed across S^3 that lies inside the fundamental zone
    //  of oSym.
    //
    //  Note that minimal uniform grid defined here to be a uniform grid point generated with
    //  upperbound of fMaxDistance such that no two points are closer than fMinDistance.
    //  Clearly this operation is extremely slow.
    //
    //  WARNING:  Currently, we're simply finding points that are in the fundamental zone.  Note
    //            that this will lead to a possible increase of discrepency by a factor of two
    //            near the boundaries.  A simple solution is to reduce everything back to the fundamental
    //            zone, then remove them one at a time according to fMinDistance.  Alternative is
    //            to calculate the apply stratified sampling.  Finally, it's also possible to
    //            reduce the original grid generated by first figuring out which of the hyperfaces
    //            does not contain any point in the fundamental zone.
    //
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> CQuaternionGrid::GetMinimalFZGrid( Float fMinDistance,
                                                           Float fMaxDistance,
                                                           const LatticeSymmetry::CSymmetry & oSym )
    {
      using namespace LatticeSymmetry;
      vector<SQuaternion> oOrientationGrid;
      vector<SQuaternion> oFZOrientationGrid;
      oOrientationGrid = GetGrid( fMaxDistance );

      for( Size_Type i = 0; i < oOrientationGrid.size(); i ++ )
      {
        if( IsInFundamentalZone( oSym, oOrientationGrid[i] ) )
        {
          oFZOrientationGrid.push_back( oOrientationGrid[i] );
        }
      }

      //      std::cout << "Before patching " << oFZOrientationGrid.size() << std::endl;
      PatchOrientationGrid( oFZOrientationGrid, oSym, fMaxDistance );
      return oFZOrientationGrid;
    }
  
    
  }
}
