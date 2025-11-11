////////////////////////////////////////////////////////////////
//
//  File:    Sampling.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  A class with implementations of sampling methods.  (i.e., uniform
//  grid on SO(3), stratified sampling, and so forth.
//
/////////////////////////////////////////////////////////////////


#ifndef SAMPLING_H_
#define SAMPLING_H_

#include "3dMath.h"
#include "Quaternion.h"
#include "Symmetry.h"
#include <vector>
#include <random>
#include <tuple>

using std::vector;

namespace GeneralLib
{

  //---------------------------------------------------------------------
  //  CUniformRandomReal
  //
  //  A simple random wrapper
  //---------------------------------------------------------------------
  class CUniformRandomReal
  {
  private:
    std::mt19937  oRngEngine;
    std::uniform_real_distribution<Float> oDistribution;
  public:

    //---------------------------------------------------------------------
    //  Constructor specifying a range
    //---------------------------------------------------------------------
    CUniformRandomReal( Float fMin, Float fMax ): oRngEngine(),
      oDistribution( fMin, fMax ) {}

    //---------------------------------------------------------------------
    //  Default constructor (range is [0, 1) )
    //---------------------------------------------------------------------
    CUniformRandomReal(): oRngEngine(),
      oDistribution( 0, 1 ) {}

    //---------------------------------------------------------------------
    //  GetRandomVariable
    //---------------------------------------------------------------------
    Float GetRandomVariable( Float fMin, Float fMax )
    {
      std::uniform_real_distribution<Float> dist( fMin, fMax );
      return dist( oRngEngine );
    }

    //---------------------------------------------------------------------
    //  GetRandomVariable
    //   This returns the natural range of the constructed random number
    //   generator.  Note that the default range is [0, 1)
    //---------------------------------------------------------------------
    Float GetRandomVariable( )
    {
      return oDistribution( oRngEngine );
    }

    //---------------------------------------------------------------------
    //  operator() - to be used as a functor
    //---------------------------------------------------------------------
    Float operator()()
    {
      return oDistribution( oRngEngine );
    }

    //---------------------------------------------------------------------
    //  operator() - to be used as a functor
    //---------------------------------------------------------------------
    Float operator()( Float fMin, Float fMax )
    {
      return GetRandomVariable( fMin, fMax );
    }
  };

  //--------------------------------------------------------------------------------------------------------
  //  SingletonUniformRandomReal
  //
  //  A singleton wrapper around the random number generator.  Sometimes
  //  a presistent random number generator is needed for the entire life time
  //  of a program.  Therefore a singleton, or a global object with only one
  //  instance may be more appropriate.
  //--------------------------------------------------------------------------------------------------------
  class CSingletonUniformRandomReal : public CUniformRandomReal
  {
  public:
    static CSingletonUniformRandomReal & Get()
    {
      static CSingletonUniformRandomReal SingletonRandomRealGenerator;
      return SingletonRandomRealGenerator;
    }
    
  private:
    CSingletonUniformRandomReal() : CUniformRandomReal() {}
    CSingletonUniformRandomReal(const CSingletonUniformRandomReal &c);  // prevent copy construction
    const CSingletonUniformRandomReal & operator= (const CSingletonUniformRandomReal &c);
  };
 
//--------------------------------------------------------------------------------------------------------
//
//  namespace UniformGrid
//
//--------------------------------------------------------------------------------------------------------
namespace UniformGrid
{
  static const Int nVertices = 8;
  static const Int nCubeDimension = 3;
  static const Int VertexOrderGrayCode[ nVertices ][ nCubeDimension ] =
    {
      { 0, 0, 0 },
      { 1, 1, 1 },
      { 0, 0, 1 },
      { 1, 1, 0 },
      { 0, 1, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 1, 0, 0 }      
    };

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
  //------------------------------------------------------------------------
  vector<SVector3> MakeSukarevGridPoints( Int nLevel, Float fSideWidth = 1.0 );
 
  //------------------------------------------------------------------------
  //
  //  GenerateGrid
  //
  //  A helper function for GetSukharevGridPoint
  //
  //------------------------------------------------------------------------
  SVector3 GetSukarevGridPoint( Int nState, Int nDim  = nCubeDimension);
    
    
  //------------------------------------------------------------------------
  //
  //  Code from Steven Lindemann, UIUC for Layered Sukharev Grid Sequence 
  //
  //  Generates a Layered Sukharev Grid sequence for d-dimensional cube.
  //  This function maps: int -> ( real, real, real... )  Grid generated 
  //  by this function is incremental and uniform.
  //
  //------------------------------------------------------------------------
  SVector3 GetLayeredSukharevGridPoint( Int nState, Int nDim  = nCubeDimension );


  //----------------------------------------------------------
  //
  //  CQuaternionGrid
  //
  //  Purpose:  Produce an adaptive grid in SO(3) using quaternions
  //----------------------------------------------------------
  class CQuaternionGrid
  {
    
  private:

    static const int nNumHyperFaces = 4;
    static const int nNumVertices = 8;

    //-----------------------------
    // Copy constructor - Not allowed
    //-----------------------------
    CQuaternionGrid( const CQuaternionGrid * oGrid );

    //-----------------------------
    // Initialize vertices of each face
    //-----------------------------
    void InitializeFaces();

    SQuaternion pHyperFaceVertices[ nNumHyperFaces ][ nNumVertices ];


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
    //  Warning:   This is definitely a slow operation!
    //
    //---------------------------------------------------------------------------------------
    void PatchOrientationGrid( vector<SQuaternion> & oOrientationGrid,
                               const LatticeSymmetry::CSymmetry & oSym,
                               Float fMinMisorientation );

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
    //               oSym               - Symmetry group that this orientation grid belongs to.
    //---------------------------------------------------------------------------------------
    std::pair<Bool, SQuaternion> GetViolatingPoint( const vector<SQuaternion>::iterator & pSamplePoint,
                                                    const vector<SQuaternion> & oOrientationGrid,
                                                    const LatticeSymmetry::CSymmetry & oSym,
                                                    Float fMinMisorientation );
    
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
    Float GetQuatPerturbationRadius( const SQuaternion & oCenter, Float fMisorientation );


    //---------------------------------------------------------------------------------------
    //
    // PerturbQuaternionPoint
    //
    //---------------------------------------------------------------------------------------
    void PerturbQuatSampPoints( vector<SQuaternion> & oOrientationGrid, Float fDelta );
    
  public:

    //-----------------------------
    // Default constructor
    //-----------------------------
    CQuaternionGrid();

    //---------------------------------------------------------------------------------------
    // BarycentricToQuaternion
    //
    // Purpose:  Using the defined function to get a point on the S^3 sphere parameterized by
    //           alpha, beta, and gamma.  A face index is required to identify which of the
    //           4 faces in the upper hemi-hypersphere is to be considered.  In another words,
    //           all points on the S^3 is parameterized by the tuple, ( n, Alpha, Beta, Gamma)
    //
    //
    // Parameters:  nFaceIndex -  [0, 3]
    //              fAlpha, fBeta, fGamma - [0, 1]
    //
    // Note:     Return value will be restricted to the fZ >= 0 hemi-hypersphere.
    //
    //---------------------------------------------------------------------------------------
    SQuaternion BarycentricToQuaternion( Int nFaceIndex, Float fAlpha, Float fBeta, Float fGamma ) const;
    
    //---------------------------------------------------------------------------------------
    //
    //  GetGrid
    //
    //  Return a grid with nNumPoints uniformly distributed across S^3 with unit 1. (i.e., uniform
    //  rotation that's evenly spaced)
    //
    //
    //---------------------------------------------------------------------------------------
    vector< SQuaternion > GetGrid( Int nNumPoints );
    
    //---------------------------------------------------------------------------------------
    //
    //  GetGrid
    //
    //  Return a grid with uniformly distribued across S^3 such that no two points are more than
    //  fMaxDistance apart.
    //
    //  Note:  fMaxDistance is the dispersion, or the radius of the largest empty ball in the sample
    //         grid.  Note also that fMaxDistance is the Eucledian distance in S^3
    //
    //         This is due to Proposition 4.5 of  Anna Yershova and Steven M. LaValle,
    //         "Deterministic Sampling Methods for Spheres and SO(3)"
    //
    //         d_rho( T ) <=  2 Pi / (n ( 2^d - 1 )/( 2 ( d + 1 )) + 1 ) ^(1/d )
    //
    //         Where T is the sequence (or grid) of points generated, d_rho is the dispersion.
    //         d is the dimension of the sphere.  In our case S^3, d = 3.
    //
    //         This works out to give us 
    //
    //          n  ~ 1/2 * [(2 * pi / d_rho)^3 - 1 ]
    //
    //  TODO:  Make sure that d_rho is measured in angle (i.e., misorientation) rather than
    //         Euclidean norm (i.e., |q1 - q2|) or geodesic distance.
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> GetGrid( Float fMaxDistance );

    SQuaternion GetNearIdentityPoint( Float fX, Float fY, Float fZ ) const;
    //---------------------------------------------------------------------------------------
    //
    //  GetRandomLocalGrid
    //
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> GetRandomLocalGrid( Float fMaxDistance, Int nGridPoints ) const;

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
    vector<SQuaternion> GetStructuredLocalGrid( Float fMaxDistance, Int nLevel = 2 );
    
    //---------------------------------------------------------------------------------------
    //
    //  GetMinimalFZGrid
    //
    //  Return minimal uniform grid distributed across S^3 that lies inside the fundamental zone
    //  of oSym.
    //
    //  Note:  Minimal uniform grid defined here to be a uniform grid point generated with
    //         upperbound of fMaxDistance.  Note that fMinDistance is really a guide rather
    //         than a gaurantee.  It provides the following constraint:
    //    
    //     1.  \forall p_i, p_j \in S, the sampling grid of SO(3), dist( p_i, p_j ) >= d_l, where
    //         d_l <= fMinDistance, d_l is specified by Proposition 4.5 of [Yershova and LaValle]
    //    
    //     2.  \forall p_i, p_j \in S, the sampling grid of SO(3), dist( p_i, p_j ) <= fMaxDistance.        
    //         This is strictly enforced.
    //
    //
    //         The definition of fMaxDistance is the same as above.  fMinDistance is measured
    //         using misorientation.  
    //
    //  WARNING:  Clearly this operation is extremely slow.
    // 
    //  WARNING:  Dispersion gaurantee may be violated around the edges of the fundamental zone!!!
    //            (It is possible to fix it, but it will be much more costly.)
    //
    //  TODO:    Resolve boundary gaurantee problem.  [ Done April 27, 2009.  Not tested ]
    //
    //---------------------------------------------------------------------------------------
    vector<SQuaternion> GetMinimalFZGrid( Float fMinDistance,
                                          Float fMaxDistance,
                                          const LatticeSymmetry::CSymmetry & oSym );
    
  };

} // UniformGrid
} // Namespace GeneralLib


#endif
