/////////////////////////////////////////////////////////////////
//
//  File:    Symmetry.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
/////////////////////////////////////////////////////////////////
#include "Symmetry.h"
#include <iostream>
using namespace std;

namespace LatticeSymmetry
{
  const char* pSymmetryKeyword[] =  {
    "Cubic", "Hexagonal", "Tetragonal", "None"
  }; 

  //-------------------------------------------------------------------------------------------
  //  ParseSymmetry
  //-------------------------------------------------------------------------------------------
  ESymmetryT ParseSymmetry( const string & oTag )
  {
    for( UInt j = 0; j < eNumSymmetry; j++)
      if ( oTag == pSymmetryKeyword[ j ] )
        return ESymmetryT(j);

    return eNumSymmetry;
  }
  
  //-------------------------------------------------------------------------------------------
  //
  //
  //
  //-------------------------------------------------------------------------------------------
  void CSymmetry::BuildQuaternionOperatorList( const vector<SMatrix3x3> &oMatrixOperatorList )
  {
    for ( Size_Type i = 0; i < oMatrixOperatorList.size(); i ++ )
    {
      SQuaternion q;
      q.Set( oMatrixOperatorList[i] );

      // round to zero (there should be no values that would be approximately < 1e-8
      if( fabs( q.m_fX ) < EPSILON )
        q.m_fX = 0;
      if( fabs( q.m_fY ) < EPSILON )
        q.m_fY = 0;
      if( fabs( q.m_fZ ) < EPSILON )
        q.m_fZ = 0;
      if( fabs( q.m_fW ) < EPSILON )
        q.m_fW = 0;
      
      oQuatOperatorList.push_back( q );
    }
  }
  


  CCubicSymmetry & CCubicSymmetry::Get()
  {
    static CCubicSymmetry oCubicSymmetry;
    return oCubicSymmetry;
  }

 
  //-------------------------------------------------------------------------------------------
  //
  //
  //  Initialization of the Cubic symmetry list
  //
  //
  //-------------------------------------------------------------------------------------------
  CCubicSymmetry::CCubicSymmetry()
  {
    SMatrix3x3 oSymMatrix;
    oSymMatrix.SetZero();
    
    oSymMatrix.m[0][1] = 1;
    oSymMatrix.m[1][0] = -1;
    oSymMatrix.m[2][2] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // L_{001}^2
    oSymMatrix.m[0][0] = -1;
    oSymMatrix.m[1][1] = -1;
    oSymMatrix.m[2][2] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();


    // L_{00-1}^4
    oSymMatrix.m[0][1] = -1;
    oSymMatrix.m[1][0] = 1;
    oSymMatrix.m[2][2] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();


    // L_{010}^4
    oSymMatrix.m[0][2] = -1;
    oSymMatrix.m[1][1] = 1;
    oSymMatrix.m[2][0] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // L_{010}^2
    oSymMatrix.m[0][0] = -1;
    oSymMatrix.m[1][1] = 1;
    oSymMatrix.m[2][2] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{0-10}^4
    oSymMatrix.m[0][2] = 1;
    oSymMatrix.m[1][1] = 1;
    oSymMatrix.m[2][0] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{100}^4
    oSymMatrix.m[0][0] = 1;
    oSymMatrix.m[1][2] = 1;
    oSymMatrix.m[2][1] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{100}^2
    oSymMatrix.m[0][0] = 1;
    oSymMatrix.m[1][1] = -1;
    oSymMatrix.m[2][2] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{-100}^4
    oSymMatrix.m[0][0] = 1;
    oSymMatrix.m[1][2] = -1;
    oSymMatrix.m[2][1] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();



    //
    //  120 degrees about {111}s's 8 operations
    // % L_{111}^3
    oSymMatrix.m[0][1] = 1;
    oSymMatrix.m[1][2] = 1;
    oSymMatrix.m[2][0] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();
    // % L_{-1-1-1}^3
    oSymMatrix.m[0][2] = 1;
    oSymMatrix.m[1][0] = 1;
    oSymMatrix.m[2][1] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{11-1}^3
    oSymMatrix.m[0][2] = -1;
    oSymMatrix.m[1][0] = 1;
    oSymMatrix.m[2][1] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();
    // % 	L_{-1-11}^3;
    oSymMatrix.m[0][1] = 1;
    oSymMatrix.m[1][2] = -1;
    oSymMatrix.m[2][0] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();
    // % L_{1-11}^3
    oSymMatrix.m[0][2] = 1;
    oSymMatrix.m[1][0] = -1;
    oSymMatrix.m[2][1] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{-11-1}^3
    oSymMatrix.m[0][1] = -1;
    oSymMatrix.m[1][2] = -1;
    oSymMatrix.m[2][0] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{-111}^3
    oSymMatrix.m[0][2] = -1;
    oSymMatrix.m[1][0] = -1;
    oSymMatrix.m[2][1] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{1-1-1}^3
    oSymMatrix.m[0][1] = -1;
    oSymMatrix.m[1][2] = 1;
    oSymMatrix.m[2][0] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    //////////////////////////////////////////////////////////////////
    //       % 180 degrees about {110}'s: 6 operations

    // % L_{011}^2
    oSymMatrix.m[0][0] = -1;
    oSymMatrix.m[1][2] = 1;
    oSymMatrix.m[2][1] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{101}^2
    oSymMatrix.m[0][2] = 1;
    oSymMatrix.m[1][1] = -1;
    oSymMatrix.m[2][0] = 1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{110}^2
    oSymMatrix.m[0][1] = 1;
    oSymMatrix.m[1][0] = 1;
    oSymMatrix.m[2][2] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{01-1}^2
    oSymMatrix.m[0][0] = -1;
    oSymMatrix.m[1][2] = -1;
    oSymMatrix.m[2][1] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{10-1}^2
    oSymMatrix.m[0][2] = -1;
    oSymMatrix.m[1][1] = -1;
    oSymMatrix.m[2][0] = -1;

    oOperatorList.push_back(oSymMatrix);
    oSymMatrix.SetZero();

    // % L_{-110}^2
    oSymMatrix.m[0][1] = -1;
    oSymMatrix.m[1][0] = -1;
    oSymMatrix.m[2][2] = -1;

    oOperatorList.push_back(oSymMatrix);


    ////////////////////////////////
    //      % I: identity

    oSymMatrix.SetIdentity();
    oOperatorList.push_back(oSymMatrix);

    DEBUG_ASSERT(oOperatorList.size() == 24, "ERROR:  Missing Cubic Symmetry Operator\n");

    BuildQuaternionOperatorList( oOperatorList );
  }

  //-------------------------------------------------------------------------------------------
  //
  //
  //  Public:    const CHexagonalSymmetry & CHexagonalSymmetry::Get()
  //
  //
  //-------------------------------------------------------------------------------------------
  CHexagonalSymmetry & CHexagonalSymmetry::Get()
  {

    static CHexagonalSymmetry oHexagonalSymmetry;;
    return oHexagonalSymmetry;
  }


  //-------------------------------------------------------------------------------------------
  //
  //
  //  Public:  CHexagonalSymmetry()
  //
  //
  //-------------------------------------------------------------------------------------------
  CHexagonalSymmetry::CHexagonalSymmetry()
  {

    SMatrix3x3 oSymMatrix;
    
    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = 1;
    oSymMatrix.m[0][1] = - Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] =   Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);
    
    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = 1;
    oSymMatrix.m[0][1] = - Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] =   Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = - 1;
    oSymMatrix.m[1][1] = - 1;
    oSymMatrix.m[2][2] = 1;
    oOperatorList.push_back(oSymMatrix);

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = 1;
    oSymMatrix.m[0][1] =   Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] = - Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = 1;
    oSymMatrix.m[0][1] =   Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] = - Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = 1;
    oSymMatrix.m[1][1] = - 1;
    oSymMatrix.m[2][2] = - 1;
    oOperatorList.push_back(oSymMatrix);

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = - 1;
    oSymMatrix.m[0][1] =  Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] =  Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);    

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = - 1;
    oSymMatrix.m[0][1] =  Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] =  Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);    

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = - 1;
    oSymMatrix.m[1][1] =   1;
    oSymMatrix.m[2][2] = - 1;
    oOperatorList.push_back(oSymMatrix);

    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = - 1;
    oSymMatrix.m[0][1] = - Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] = - Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);
    
    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = Float( 1.0 / 2.0 );
    oSymMatrix.m[1][1] = - Float( 1.0 / 2.0 );
    oSymMatrix.m[2][2] = - 1;
    oSymMatrix.m[0][1] = - Float( sqrt(3.0) / 2.0 );
    oSymMatrix.m[1][0] = - Float( sqrt(3.0) / 2.0 );
    oOperatorList.push_back(oSymMatrix);
    


    oSymMatrix.SetIdentity();
    oOperatorList.push_back(oSymMatrix);

    
    DEBUG_ASSERT(oOperatorList.size() == 12, "ERROR:  Missing Hexagonal Symmetry Operator\n");

    BuildQuaternionOperatorList( oOperatorList );
  }


  //-------------------------------------------------------------------------------------------
  //
  //
  //  Public:    const CHexagonalSymmetry & CHexagonalSymmetry::Get()
  //
  //
  //-------------------------------------------------------------------------------------------
  CTetragonalSymmetry & CTetragonalSymmetry::Get()
  {
    static CTetragonalSymmetry oTetragonalSymmetry;
    return oTetragonalSymmetry;
  }

  //-------------------------------------------------------------------------------------------
  //
  //  Public:  CTetragonalSymmetry
  //
  //-------------------------------------------------------------------------------------------
  CTetragonalSymmetry::CTetragonalSymmetry()
  {
    SMatrix3x3 oSymMatrix;
    // 1
    //(001,0);
    oSymMatrix.SetIdentity();
    oOperatorList.push_back(oSymMatrix);

    // 2
    //(001,180);
    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = -1;
    oSymMatrix.m[1][1] = -1;
    oSymMatrix.m[2][2] =  1;
    oOperatorList.push_back(oSymMatrix);

    // 3
    //(001,90);
    oSymMatrix.SetZero();
    oSymMatrix.m[0][1] =  1;
    oSymMatrix.m[1][0] = -1;
    oSymMatrix.m[2][2] =  1;
    oOperatorList.push_back(oSymMatrix);
 
    // 4
    //(001,270);
    oSymMatrix.SetZero();
    oSymMatrix.m[0][1] = -1;
    oSymMatrix.m[1][0] =  1;
    oSymMatrix.m[2][2] =  1;
    oOperatorList.push_back(oSymMatrix);

    // 5
    //(100,180);
    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] =  1;
    oSymMatrix.m[1][1] = -1;
    oSymMatrix.m[2][2] = -1;
    oOperatorList.push_back(oSymMatrix);

    // 6
    //(010,180)
    oSymMatrix.SetZero();
    oSymMatrix.m[0][0] = -1;
    oSymMatrix.m[1][1] =  1;
    oSymMatrix.m[2][2] = -1;
    oOperatorList.push_back(oSymMatrix);

    // 7
    //(1-10,180)
    oSymMatrix.SetZero();
    oSymMatrix.m[0][1] = -1;
    oSymMatrix.m[1][0] = -1;
    oSymMatrix.m[2][2] = -1;
    oOperatorList.push_back(oSymMatrix);

    // 8
    //(110,180)
    oSymMatrix.SetZero();
    oSymMatrix.m[0][1] =  1;
    oSymMatrix.m[1][0] =  1;
    oSymMatrix.m[2][2] = -1;
    oOperatorList.push_back(oSymMatrix);

    
    DEBUG_ASSERT(oOperatorList.size() == 8, "ERROR:  Missing Cubic Symmetry Operator\n");

    BuildQuaternionOperatorList( oOperatorList );
  }

  // Utilities function
  namespace Utilities
  {

    //---------------------------------------------------
    //  GetSymmmetryObj
    //
    //  The only reason why returning address here works
    //  is because these objects are singletons.
    //
    //  
    //---------------------------------------------------
    CSymmetry* GetSymmetryObj( ESymmetryT eSym )
    {
      switch( eSym )
      {
        case LatticeSymmetry::eCubic:
          return &( LatticeSymmetry::CCubicSymmetry::Get() );
        case LatticeSymmetry::eHexagonal:
          return & ( LatticeSymmetry::CHexagonalSymmetry::Get() );
        case LatticeSymmetry::eTetragonal:
          return & ( LatticeSymmetry::CTetragonalSymmetry::Get() );
        default:
          cerr << "WARNING!!!!  Not using any symmetry!!" << std::endl;
          return 0;
      }
    }
  }
  
  
}// end namespace
