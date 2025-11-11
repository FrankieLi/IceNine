////////////////////////////////////////////////////////////////
//
//  File:    Quaternion.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Quaternion class for basic SO(3), S^4 math
//  
//
/////////////////////////////////////////////////////////////////


#ifndef QUATERNION_H_
#define QUATERNION_H_

#include "3dMath.h"
#include <limits>
#include <math.h>
#include <random>

typedef Float FLOAT;

namespace GeneralLib
{
  
  // Forward declearation
  class SQuaternion;

  
  //----------------------------------------------------------
  //
  //  General lib functions
  //
  //----------------------------------------------------------
  
  //-----------------------------
  // Dot product
  //-----------------------------
  Float Dot( const SQuaternion & q1, const SQuaternion & q2 );
  
  //-----------------------------
  // class CRandomQuaternian
  //
  // Initializes a random number generator to enable random quanterion
  // generation.
  //
  //-----------------------------
  class CRandomRotationGenerator
  {
  private:
    std::mt19937 oRngEngine;
    std::uniform_real_distribution<Float> oDistribution;

  public:

    //----------------------------------------------------------
    //  Default constructor:  initializes the random number generator
    //                        automatically
    //----------------------------------------------------------
    CRandomRotationGenerator(): oRngEngine(),
                                oDistribution( 0, 1 )
    {
    };
    


    //----------------------------------------------------------
    //  Return a random quanternion using the random number generator
    //  specified.
    //----------------------------------------------------------
    SQuaternion GetRandomQuaternion( );

    //----------------------------------------------------------
    //  Return a random quanternion using the random number generator
    //  specified, with the range scaled by the parameter.
    //
    //  This is an ad hoc way to produce local random orientations
    //
    //----------------------------------------------------------
    SQuaternion GetRandomQuaternion( Float fScale );

    
  };


  
  //----------------------------------------------------------
  // Interpolate (SLERP) -- spherical linear interpolation
  //----------------------------------------------------------
  SQuaternion Interpolate(const SQuaternion &oParam1, const SQuaternion &oParam2, FLOAT fInterp );
  

  //----------------------------------------------------------
  //
  //  SQuaternion
  //
  //----------------------------------------------------------
  class SQuaternion
  {
  public:
    Float m_fW;
    Float m_fX;
    Float m_fY;
    Float m_fZ;

    SQuaternion();
    SQuaternion( const SVector3 &oAxis, FLOAT fAngle );
    ~SQuaternion();
    
    void CreateFromAxisAngle( const SVector3 &oAxis, FLOAT fAngle );
    SMatrix4x4 BuildTransformationMatrix( const SVector3 &oScale,
                                          const SVector3 &oTranslate );

    //-----------------------------
    //  Set ( from Euler angle of an active matrix )
    //  
    //  Note:  this is in z-x-z convention
    //-----------------------------
    void Set( Float fPhi, Float fTheta, Float fPsi );

    //-----------------------------
    //  Set ( from rotation matrix)
    //-----------------------------
    void Set( const SMatrix4x4 & mActiveMatrix );
    void Set( const SMatrix3x3 & mActiveMatrix );
    void Set( FLOAT fW , FLOAT fX, FLOAT fY, FLOAT fZ);
    

    
    //-----------------------------
    // EuclideanNorm
    //-----------------------------
    Float EuclideanNorm() const
    {
      return sqrt( Dot( *this, *this ) );
    };

    //-----------------------------
    // ToConvention  -- make sure that
    // the quaternion is in the positive hemi-sphere
    //-----------------------------
    void ToConvention()
    {
      if ( m_fW < 0 )
        *this = - *this;
    };
    
    //-----------------------------
    //  Convert quaternion to active rotation
    //  matrix
    //-----------------------------
    SMatrix4x4 GetRotationMatrix4x4() const;
    SMatrix3x3 GetRotationMatrix3x3() const;

    //-----------------------------
    //  Get the axis angle representation of
    //  quaternion
    //-----------------------------
    std::pair<SVector3, Float> GetAxisAngle() const;

    
    // Return the inverse
    SQuaternion Inverse() const;
    
    //
    //  operator*  -- Hamiltonian Product
    //  (For composition of rotations)
    //
    SQuaternion operator *( const SQuaternion &oQuat ) const;
    SQuaternion operator -( ) const;
    
    //
    //  Operations that thinks of quaternions as R^4
    //  Use with caution  (i.e., it doesn't make sense to add
    //  rotations, but it makes perfect sense to add 2 quaternions
    //  together.
    SQuaternion operator - ( const SQuaternion & oRHS ) const;
    SQuaternion operator + ( const SQuaternion & oRHS ) const;
    SQuaternion & operator +=( const SQuaternion & oRHS );

    SQuaternion operator/ ( Float fScale ) const;
  };

  //----------------------------------------------------------
  //  Stream operator
  //----------------------------------------------------------
  std::ostream & operator<< ( std::ostream & os, const SQuaternion & oRHS ); 
 

}



#endif
