////////////////////////////////////////////////////////////////
//
//  File:    Quaternion.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Implementation of quaternion class for basic SO(3), S^4 math
//  
//
/////////////////////////////////////////////////////////////////

#include "Quaternion.h"

namespace GeneralLib
{
  
  //----------------------------------------------------------------------------------------------
  // Public : SQuaternion
  //----------------------------------------------------------------------------------------------
  SQuaternion::SQuaternion() : m_fW(1.0f), m_fX(0.0f), m_fY(0.0f), m_fZ(0.0f)
  {

  }

  //----------------------------------------------------------------------------------------------
  // Public : SQuaternion
  //----------------------------------------------------------------------------------------------
  SQuaternion::SQuaternion( const SVector3 &oAxis, FLOAT fAngle )
  {
    CreateFromAxisAngle( oAxis, fAngle );
  }

  //----------------------------------------------------------------------------------------------
  // Public : ~SQuaternion
  //----------------------------------------------------------------------------------------------
  SQuaternion::~SQuaternion()
  {

  }

  //----------------------------------------------------------------------------------------------
  //  Public:  Set ( from Euler angle of an active matrix )
  //----------------------------------------------------------------------------------------------
  void SQuaternion::Set( Float fPhi, Float fTheta, Float fPsi )
  {
    SMatrix3x3 oMat;
    oMat.BuildActiveEulerMatrix( fPhi, fTheta, fPsi );

    Set( oMat );
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Set
  //
  // This method is adapation of Donald Boyce's OdfPf.  It is more numerically stable than the other
  // methods
  //----------------------------------------------------------------------------------------------  
  void SQuaternion::Set( const SMatrix3x3 & oMat )
  {
    
    Float fTrace = oMat.Trace() - Float( 1.0 );
    Float fCosAngle = 0.5 * fTrace;
    //  const Float fEpsilon = Float(10) * std::numeric_limits<Float>::epsilon();
    const Float fEpsilon = 1e-4;
  
    fCosAngle = std::min( fCosAngle, Float(  1 ) );    // clamp cos(angle) to [-1, 1]
    fCosAngle = std::max( fCosAngle, Float( -1 ) );
    Float fAngle = acos( fCosAngle );

    SVector3 oRotAxis;
    oRotAxis.Set( oMat.m[2][1] - oMat.m[1][2],
                  oMat.m[0][2] - oMat.m[2][0],
                  oMat.m[1][0] - oMat.m[0][1] );
        
    if ( fAngle < fEpsilon )
    {
      fAngle = 0;
      oRotAxis.Set( 1, 1, 1 );
    }

    if( fAngle > (PI - fEpsilon) )   // near PI, find axis of rotation by
    {                                // looking for the largest component
      SMatrix3x3 oTmp;
      oTmp.SetIdentity();
      oTmp += oMat;
      
      Int nMaxIndex = 0;
      Float fMaxDot = -1;   // the dot product computed will always be positive here
      for( Int i = 0; i < 3; i ++ )
      {
        Float fDot =   oTmp.m[i][0] * oTmp.m[i][0]    // magnitude is always positive
                     + oTmp.m[i][1] * oTmp.m[i][1]
                     + oTmp.m[i][2] * oTmp.m[i][2];

        if ( fDot > fMaxDot )
        {
          fMaxDot = fDot;
          nMaxIndex = i;
        }
      }
      oRotAxis.Set( oTmp.m[ nMaxIndex ][0],
                    oTmp.m[ nMaxIndex ][1],
                    oTmp.m[ nMaxIndex ][2] );
    }
    oRotAxis.Normalize();
    CreateFromAxisAngle( oRotAxis, fAngle );
  }
  //----------------------------------------------------------------------------------------------
  // Public : Set
  //----------------------------------------------------------------------------------------------
  void SQuaternion::Set( FLOAT fW, FLOAT fX, FLOAT fY, FLOAT fZ)
  {
    m_fW = fW;
    m_fX = fX;
    m_fY = fY;
    m_fZ = fZ;
  }

  //----------------------------------------------------------------------------------------------
  // Public : CreateFromAxisAngle
  //
  //  Note Angle is in Radian
  //----------------------------------------------------------------------------------------------
  void SQuaternion::CreateFromAxisAngle( const SVector3 &oAxis, FLOAT fAngle )
  {
    FLOAT fSinAngle = sin( fAngle/2.0f ) ;
    m_fW = FLOAT( cos( fAngle/2.0f ) );

    // Calculate the x, y and z of the quaternion
    m_fX = FLOAT( oAxis.m_fX * fSinAngle );
    m_fY = FLOAT( oAxis.m_fY * fSinAngle );
    m_fZ = FLOAT( oAxis.m_fZ * fSinAngle );
  }

  //----------------------------------------------------------------------------------------------
  //  Inverse Operator
  //----------------------------------------------------------------------------------------------
  SQuaternion SQuaternion::Inverse() const
  {
    SQuaternion oRes;
    oRes.m_fW = - m_fW;
    oRes.m_fX =   m_fX;
    oRes.m_fY =   m_fY;
    oRes.m_fZ =   m_fZ;
    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SQuaternion SQuaternion::operator*( const SQuaternion  &oQuat ) const
  {
    SQuaternion oRes;
    oRes.m_fW = m_fW * oQuat.m_fW - m_fX * oQuat.m_fX - m_fY * oQuat.m_fY - m_fZ * oQuat.m_fZ;
    oRes.m_fX = m_fW * oQuat.m_fX + m_fX * oQuat.m_fW + m_fY * oQuat.m_fZ - m_fZ * oQuat.m_fY;
    oRes.m_fY = m_fW * oQuat.m_fY + m_fY * oQuat.m_fW + m_fZ * oQuat.m_fX - m_fX * oQuat.m_fZ;
    oRes.m_fZ = m_fW * oQuat.m_fZ + m_fZ * oQuat.m_fW + m_fX * oQuat.m_fY - m_fY * oQuat.m_fX;

    oRes.ToConvention( );
    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator-
  //----------------------------------------------------------------------------------------------
  SQuaternion SQuaternion::operator - ( const SQuaternion & oRHS ) const
  {
    SQuaternion oRes;
    oRes.m_fW =  m_fW - oRHS.m_fW;
    oRes.m_fX =  m_fX - oRHS.m_fX;
    oRes.m_fY =  m_fY - oRHS.m_fY;
    oRes.m_fZ =  m_fZ - oRHS.m_fZ;
    return oRes;
  }
  //----------------------------------------------------------------------------------------------
  // Public : operator-
  //----------------------------------------------------------------------------------------------
  SQuaternion SQuaternion::operator-() const
  {
    SQuaternion oRes;
    oRes.Set( -m_fW, -m_fX, -m_fY, -m_fZ );
    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator+
  //----------------------------------------------------------------------------------------------
  SQuaternion SQuaternion::operator+( const SQuaternion & oRHS  ) const
  {
    SQuaternion oRes;
    oRes.m_fW =  m_fW + oRHS.m_fW;
    oRes.m_fX =  m_fX + oRHS.m_fX;
    oRes.m_fY =  m_fY + oRHS.m_fY;
    oRes.m_fZ =  m_fZ + oRHS.m_fZ;

    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator+
  //----------------------------------------------------------------------------------------------
  SQuaternion & SQuaternion::operator+=( const SQuaternion & oRHS  ) 
  {
    m_fW += oRHS.m_fW;
    m_fX += oRHS.m_fX;
    m_fY += oRHS.m_fY;
    m_fZ += oRHS.m_fZ;
    return *this;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator/
  //----------------------------------------------------------------------------------------------
  SQuaternion SQuaternion::operator/ ( Float fScale ) const
  {
    SQuaternion oRes;

    oRes.m_fW = m_fW / fScale;
    oRes.m_fX = m_fX / fScale;
    oRes.m_fY = m_fY / fScale;
    oRes.m_fZ = m_fZ / fScale;
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  //  Stream operator
  //----------------------------------------------------------------------------------------------
  std::ostream & operator<< ( std::ostream & os, const SQuaternion & oRHS )
  {
    os << oRHS.m_fW << " " << oRHS.m_fX << " " << oRHS.m_fY << " " << oRHS.m_fZ << " ";
    return os;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : CreateMatrix
  //----------------------------------------------------------------------------------------------
  SMatrix4x4 SQuaternion::GetRotationMatrix4x4() const
  {
    SMatrix4x4 oRes;
    //oRes.m[0][0] = (m_fW * m_fW) + (m_fX * m_fX) - (m_fY * m_fY) - (m_fZ * m_fZ); <-- Alternative form
    oRes.m[0][0] = 1.0f - 2.0f * ( m_fY * m_fY + m_fZ * m_fZ );  
    oRes.m[0][1] = 2.0f * ( m_fX * m_fY - m_fW * m_fZ );  
    oRes.m[0][2] = 2.0f * ( m_fX * m_fZ + m_fW * m_fY );  
    oRes.m[0][3] = 0.0f;  
    oRes.m[1][0] = 2.0f * ( m_fX * m_fY + m_fW * m_fZ );  
    //oRes.m[1][1] = (m_fW * m_fW) - (m_fX * m_fX) + (m_fY * m_fY) - (m_fZ * m_fZ); <-- Alternative form
    oRes.m[1][1] = 1.0f - 2.0f * ( m_fX * m_fX + m_fZ * m_fZ );  
    oRes.m[1][2] = 2.0f * ( m_fY * m_fZ - m_fW * m_fX );  
    oRes.m[1][3] = 0.0f;  
    oRes.m[2][0] = 2.0f * ( m_fX * m_fZ - m_fW * m_fY );  
    oRes.m[2][1] = 2.0f * ( m_fY * m_fZ + m_fW * m_fX );  
    //oRes.m[2][2] = (m_fW * m_fW) - (m_fX * m_fX) - (m_fY * m_fY) + (m_fZ * m_fZ); <-- Alternative form
    oRes.m[2][2] = 1.0f - 2.0f * ( m_fX * m_fX + m_fY * m_fY );  
    oRes.m[2][3] = 0.0f;  
    oRes.m[3][0] = 0;  
    oRes.m[3][1] = 0;  
    oRes.m[3][2] = 0;  
    oRes.m[3][3] = 1.0f;
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : CreateMatrix
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SQuaternion::GetRotationMatrix3x3() const
  {
    SMatrix3x3 oRes;
    //oRes.m[0][0] = (m_fW * m_fW) + (m_fX * m_fX) - (m_fY * m_fY) - (m_fZ * m_fZ); <-- Alternative form
    oRes.m[0][0] = 1.0f - 2.0f * ( m_fY * m_fY + m_fZ * m_fZ );  
    oRes.m[0][1] = 2.0f * ( m_fX * m_fY - m_fW * m_fZ );  
    oRes.m[0][2] = 2.0f * ( m_fX * m_fZ + m_fW * m_fY );  

    oRes.m[1][0] = 2.0f * ( m_fX * m_fY + m_fW * m_fZ );  
    //oRes.m[1][1] = (m_fW * m_fW) - (m_fX * m_fX) + (m_fY * m_fY) - (m_fZ * m_fZ); <-- Alternative form
    oRes.m[1][1] = 1.0f - 2.0f * ( m_fX * m_fX + m_fZ * m_fZ );  
    oRes.m[1][2] = 2.0f * ( m_fY * m_fZ - m_fW * m_fX );  

    oRes.m[2][0] = 2.0f * ( m_fX * m_fZ - m_fW * m_fY );  
    oRes.m[2][1] = 2.0f * ( m_fY * m_fZ + m_fW * m_fX );  
    //oRes.m[2][2] = (m_fW * m_fW) - (m_fX * m_fX) - (m_fY * m_fY) + (m_fZ * m_fZ); <-- Alternative form
    oRes.m[2][2] = 1.0f - 2.0f * ( m_fX * m_fX + m_fY * m_fY );  

    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public:  GetAxisAngle
  //----------------------------------------------------------------------------------------------
  std::pair<SVector3, Float> SQuaternion::GetAxisAngle() const
  {
    Float fAngle = 2 * acos( m_fW );
    Float fDenom = sqrt( Float(1) - m_fW * m_fW );
    const Float epsilon = 1e-4;
    std::pair<SVector3, Float> oRes;
    if( fabs( fAngle ) < epsilon )
      oRes.first.Set( 0, 0, 0 );
    else
      oRes.first.Set( m_fX / fDenom, m_fY / fDenom, m_fZ / fDenom );

    oRes.second = fAngle;

    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : BuildTransformationMatrix
  //----------------------------------------------------------------------------------------------
  SMatrix4x4 SQuaternion::BuildTransformationMatrix( const SVector3 &oScale,
                                                     const SVector3 &oTranslate )
  {
    SMatrix4x4 oRotateMat = GetRotationMatrix4x4();
    SMatrix4x4 oRes;
    
    // Mt * Mr * Ms.  Have to scale before rotate before translate. 
    oRes.m[0][0] = oRotateMat.m[0][0] * oScale.m_fX;
    oRes.m[0][1] = oRotateMat.m[0][1] * oScale.m_fY;
    oRes.m[0][2] = oRotateMat.m[0][2] * oScale.m_fZ;
    oRes.m[0][3] = oTranslate.m_fX;
    oRes.m[1][0] = oRotateMat.m[1][0] * oScale.m_fX;
    oRes.m[1][1] = oRotateMat.m[1][1] * oScale.m_fY;
    oRes.m[1][2] = oRotateMat.m[1][2] * oScale.m_fZ;
    oRes.m[1][3] = oTranslate.m_fY;
    oRes.m[2][0] = oRotateMat.m[2][0] * oScale.m_fX;
    oRes.m[2][1] = oRotateMat.m[2][1] * oScale.m_fY;
    oRes.m[2][2] = oRotateMat.m[2][2] * oScale.m_fZ;
    oRes.m[2][3] = oTranslate.m_fZ;
    oRes.m[3][0] = oRes.m[3][1] = oRes.m[3][2] = 0;
    oRes.m[3][3] = 1;

    return oRes;
  }


  //----------------------------------------------------------------------------------------------
  // Dot product
  //----------------------------------------------------------------------------------------------
  Float Dot( const SQuaternion & q1, const SQuaternion & q2 )
  {
    return q1.m_fW * q2.m_fW + q1.m_fX * q2.m_fX + q1.m_fY * q2.m_fY + q1.m_fZ * q2.m_fZ ; 
  }
  
  //----------------------------------------------------------------------------------------------
  // CRandomQuaternion::RandomQuaternion
  //----------------------------------------------------------------------------------------------
  SQuaternion CRandomRotationGenerator::GetRandomQuaternion( )
  {
  
    Float s = oRandomReal();
    Float theta1 = 2 * PI * oRandomReal();
    Float theta2 = 2 * PI * oRandomReal();
    Float sigma1 = sqrt( 1 - s );
    Float sigma2 = sqrt( s );
  
    Float w = cos( theta2 ) * sigma2 ;
    Float x = sin( theta1 ) * sigma1 ;
    Float y = cos( theta1 ) * sigma1 ;
    Float z = sin( theta2 ) * sigma2 ;
  
    SQuaternion qRet;
    qRet.m_fW = w;
    qRet.m_fX = x;
    qRet.m_fY = y;
    qRet.m_fZ = z;
  
    return qRet;
  
  }

  //----------------------------------------------------------
  //  Return a random quanternion using the random number generator
  //  specified, with the range scaled by the parameter.
  //
  //  This is an ad hoc way to produce local random orientations
  //
  //----------------------------------------------------------
  SQuaternion CRandomRotationGenerator::GetRandomQuaternion( Float fScale )
  {
    Float s = oRandomReal() * fScale;
    Float theta1 = 2 * PI * oRandomReal() * fScale;
    Float theta2 = 2 * PI * oRandomReal() * fScale;
    Float sigma1 = sqrt( 1 - s );
    Float sigma2 = sqrt( s );
  
    Float w = cos( theta2 ) * sigma2 ;
    Float x = sin( theta1 ) * sigma1 ;
    Float y = cos( theta1 ) * sigma1 ;
    Float z = sin( theta2 ) * sigma2 ;
  
    SQuaternion qRet;
    qRet.m_fW = w;
    qRet.m_fX = x;
    qRet.m_fY = y;
    qRet.m_fZ = z;
  
    return qRet;
  }
   
  //----------------------------------------------------------------------------------------------
  // Interpolate -- SLERP (Spherical Linear Interpolation)
  //
  // [  Not completely tested?!  (Previous comment indicates that there are bugs) ]
  //----------------------------------------------------------------------------------------------
  SQuaternion Interpolate( const SQuaternion &oParam1, const SQuaternion &oParam2, FLOAT fInterp )
  {
    FLOAT fTheta, fSinTheta, fCosTheta;
    SQuaternion oRes;

    // cosine theta = dot product of A and B with Euclidean metric
    fCosTheta = Dot( oParam1, oParam2 ); 

    // if B is on opposite hemisphere from A, use -B instead
    bool bFlip = true;
    if ( fCosTheta < FLOAT(0.0) )
      fCosTheta = -fCosTheta;
    else
      bFlip = false;

    //----------------------------------------------------------------------
    //
    // COMMENT:
    // limit theta -> 0 of  sin [ ( 1- t) theta ]/ sin( theta ) p1 + sin (t theta) / sin (theta) p2
    // =  ( 1 - t ) p1 + t p2.  (do Talyor expansion and divide the series or look at Wikipedia)
    //
    // (note though that lim theta -> 0 of sin ( theta ) = 0, which means a floating point error
    // will occur unless a linear interpolation case is checked.  Since sin(x) = x for small x
    // (1 part per 100 will do), the check for this zero condition can be leaning toward
    // not dividing by zero.  Therefore we will check the fCosTheta against 0.01 OR
    // (1 - fCosTheta) against 0.01
    //
    //----------------------------------------------------------------------
    FLOAT fAlpha = fInterp, fBeta;
    if ( FLOAT(1.0) - fCosTheta < 0.01 )
    {
      fBeta = FLOAT(1.0) - fInterp;  // linear interpolation
    }	// normal case
    else
    {
      fTheta = acos( fCosTheta );
      fSinTheta = sin( fTheta );
      fBeta = sin( fTheta - fInterp * fTheta ) / fSinTheta;
      fAlpha = sin( fInterp*fTheta ) / fSinTheta;
    }

    if ( bFlip )
      fAlpha = -fAlpha;

    // Interpolate
    oRes.m_fX = fBeta * oParam1.m_fX + fAlpha * oParam2.m_fX;
    oRes.m_fY = fBeta * oParam1.m_fY + fAlpha * oParam2.m_fY;
    oRes.m_fZ = fBeta * oParam1.m_fZ + fAlpha * oParam2.m_fZ;
    oRes.m_fW = fBeta * oParam1.m_fW + fAlpha * oParam2.m_fW;

    return oRes;
  }
  
}// end GeneralLib
