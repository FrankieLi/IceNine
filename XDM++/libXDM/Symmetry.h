/////////////////////////////////////////////////////////////////
//
//  File:    Symmetry.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Collection of all possible symmetry operators in singletons
//
//
/////////////////////////////////////////////////////////////////
#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "3dMath.h"
#include "Debug.h"
#include "Error.h"
#include <math.h>
#include <vector>
#include "Quaternion.h"
#include "CrystalPrimitives.h"


using std::vector;
using namespace GeneralLib;
namespace LatticeSymmetry
{

  enum ESymmetryT{
    eCubic, eHexagonal, eTetragonal, eNone, eNumSymmetry
  };

  //--------------------------------------------------------------------------------
  //  ParseSymmetry
  //--------------------------------------------------------------------------------
  ESymmetryT ParseSymmetry( const string & oTag );
  
  //--------------------------------------------------------------------------------
  //
  //  General Symmetry Class
  //
  //  NOTE:  By convention, the identity element is the last element in the operator list!
  //
  //
  //--------------------------------------------------------------------------------
  class CSymmetry
  {
  public:
    const vector<SMatrix3x3>  & GetOperatorList() const { return oOperatorList;};
    const vector<SQuaternion> & GetQuatOperatorList() const { return oQuatOperatorList;};
  protected:
    vector<SMatrix3x3> oOperatorList;
    vector<SQuaternion> oQuatOperatorList;
    void BuildQuaternionOperatorList( const vector<SMatrix3x3> &oMatrixOperatorList );
  };

  //--------------------------------------------------------------------------------
  //  TODO: Make this into one of those recurring templates in the near future to
  //        accomodate for multipel type  of crystal symmetries.  All symmetry classes
  //        should be singleton
  //  
  //--------------------------------------------------------------------------------
  class CCubicSymmetry : public CSymmetry
  {
    
  public:
    static CCubicSymmetry & Get();

  private:
    CCubicSymmetry();
    CCubicSymmetry(const CCubicSymmetry &c);  // prevent copy construction
    const CCubicSymmetry & operator= (const CCubicSymmetry &c);
    
  };

  //--------------------------------------------------------------------------------
  //  TODO:  Make this into the base class of Symmetry - just like CubicSymmetry
  //
  //--------------------------------------------------------------------------------
  class CHexagonalSymmetry : public CSymmetry
  {
  public:
    static CHexagonalSymmetry & Get();

  private:
    CHexagonalSymmetry();
    CHexagonalSymmetry(const CHexagonalSymmetry &c);  // prevent copy construction
    const CHexagonalSymmetry & operator= (const CHexagonalSymmetry &c);
    
  };

  //--------------------------------------------------------------------------------
  //  Tetragonal Symmetry
  //--------------------------------------------------------------------------------
  class CTetragonalSymmetry : public CSymmetry
  {
    
  public:
    static CTetragonalSymmetry & Get();

  private:
    CTetragonalSymmetry();
    CTetragonalSymmetry(const CTetragonalSymmetry &c);  // prevent copy construction
    const CTetragonalSymmetry & operator= (const CTetragonalSymmetry &c);
  };
  
  //--------------------------------------------------------------------------------
  //  Misorientation - Computes based on two given orientations matrices
  // 
  //  
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Float GetMisorientation( const CSymmetryType & oSymOps, 
                           const SMatrix3x3 &oM1,
                           const SMatrix3x3 &oM2 );


  
  //--------------------------------------------------------------------------------
  //  Misorientation - Computes based on two given orientations e1, e2, in SVector3
  // 
  //  NOTE  - Computations are done in radians
  //
  //  
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Float GetMisorientation( const CSymmetryType & oSymOps, 
                           const SVector3 &oEulerOrient1, 
                           const SVector3 &oEulerOrient2 );


  //--------------------------------------------------------------------------------
  //
  //  Misorientation  (Quanterion)
  //
  //  Return angle that takes  q1 -> q2, i.e, the angle that the operator inverse( q1 ) * q2
  //  is away from the "reference" frame.
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Float GetMisorientation( const CSymmetryType & oSymOps, 
                           const SQuaternion &q1, 
                           const SQuaternion &q2 );

  //--------------------------------------------------------------------------------
  //
  //  Symmetry reduction
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  SQuaternion ReduceToFundamentalZone( const CSymmetryType & oSym, const SQuaternion & q );

  //--------------------------------------------------------------------------------
  //  IsInFundamentalZone
  //  Return true if the quaternion specified is in the fundamental zone.
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Bool IsInFundamentalZone( const CSymmetryType & oSym, const SQuaternion & q );

  
  //--------------------------------------------------------------------------------
  //
  //  Return true of the two vectors are equivilent under symmetry operations
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  Bool Equivilent( const CSymmetryType & oSym, const SVector3 &v1, const SVector3 &v2, Float fError = EPSILON );


  //--------------------------------------------------------------------------------
  //
  //  Average -- Orientation average with respect to the rotational symmetry of the crystal.
  //
  //  NOTE:  This is an approixmation to the true average.  (See Rollett's paper on orientation
  //         averaging:  "Determination of a Mean Orientation in EBSD Measurements", Cho, Rollett, Oh )
  //
  //
  //  Precondition:  The set of orientation input must already be reduced to a fundamental zone.
  //  
  //
  //  Return the average of the orientations, given either a list of Euler angles in radians,
  //  quaternions, or Euler matrices
  //
  //  The return value is in quanterion.  (it just makes more sense that way)
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  SQuaternion Average( const CSymmetryType & oSym, const vector< SVector3 > &vAngleList );

  //--------------------------------------------------------------------------------
  //  Average
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  SQuaternion Average( const CSymmetryType & oSym, const vector< SQuaternion > &vQuatList );

  //--------------------------------------------------------------------------------
  //  Average
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType, typename Iterator, typename PropertyMap >
  SQuaternion Average( const CSymmetryType & oSym, Iterator pBegin, Iterator pEnd,
                       PropertyMap oMap );
  
  //--------------------------------------------------------------------------------
  //
  //  Reducing the list of vectors into the unique set
  //
  //  TODO:  template this, and use a default functor for vector<SVector3>
  //
  //--------------------------------------------------------------------------------
  template< typename CSymmetryType >
  vector<SVector3> GetUniqueVectors( const CSymmetryType & oSym,
                                     const vector<SVector3> &oVectorList);


  template< typename CSymmetryType >
  void GetUniqueRecpVectors( vector<CRecpVector> & oUniqueList, const CSymmetryType & oSym,
                             const vector< CRecpVector > &oRecpVectorList);



  // Utilities function
  namespace Utilities
  {
    CSymmetry* GetSymmetryObj( ESymmetryT eSym );

    struct AngleToQuaternion
    {
      SQuaternion GetQuaternion( const SVector3 & v )
      {
        SQuaternion q;
        q.Set( v.m_fX, v.m_fY, v.m_fZ );
        return q;
      }
    };

    struct QuaternionToQuaternion
    {
      SQuaternion GetQuaternion( const SQuaternion & q )
      {
        return q;
      }
    };


  }
  
}

#include "Symmetry.tmpl.cpp"
#endif
