/********************************************************************
	created:	2004/09/25
	authors:		Shan-Min Chao & Frankie Li
	
	purpose:	General purpose math class for standard 3d math.
*********************************************************************/

#include "Debug.h"
#include <math.h>
#include <stdio.h>
#include <string>
#include "3dMath.h"

namespace GeneralLib
{

  

  //----------------------------------------------------------------------------------------------
  // Public : SMatrix3x3
  //----------------------------------------------------------------------------------------------

  SMatrix3x3::SMatrix3x3()
  {}

  //----------------------------------------------------------------------------------------------
  // Public : SMatrix3x3
  //----------------------------------------------------------------------------------------------
  SMatrix3x3::SMatrix3x3(FLOAT pMatrix[3][3])
  {
    m[0][0] = pMatrix[0][0];
    m[1][0] = pMatrix[1][0];
    m[2][0] = pMatrix[2][0];

    m[0][1] = pMatrix[0][1];
    m[1][1] = pMatrix[1][1];
    m[2][1] = pMatrix[2][1];


    m[0][2] = pMatrix[0][2];
    m[1][2] = pMatrix[1][2];
    m[2][2] = pMatrix[2][2];

  }


  //----------------------------------------------------------------------------------------------
  // Public : SetIdentity
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::SetIdentity()
  {
    memset( m, 0, sizeof(FLOAT) * 9 );
    m[0][0] = m[1][1] = m[2][2] = (FLOAT)1.0;
  }

  //----------------------------------------------------------------------------------------------
  // Public : SetZero
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::SetZero()
  {
    memset( m, 0, sizeof(FLOAT) * 9 );
  }

 //----------------------------------------------------------------------------------------------
  // Public : BuildRotationAboutAxis.
  //
  // Description:  An active transformation matrix that will perform the rotation in a positive sence.
  // (i.e., follows the right hand rule.)
  //
  // Precondition:  oAxis must be a UNIT VECTOR
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildRotationAboutAxis( const SVector3 &oAxis, FLOAT fAngle )
  {
    // 1. Perform transformations which align rotation axis with one of coordinate axis (x, y, z).
    //	  This is essentially a transformation matrix with the given vector as one of the rows, and
    //    orthogonal vectors as the other rows.
    // 2. Perform rotation about the axis
    // 3. Do inverse of (1)

    // Alternatively, use the method described on page 79-80 of Eric Lengyel's "Mathematics for 3d game programming
    // and computer graphics". The core idea involves finding basis vectors and writing the rotation as a
    // linear combination of these vectors.  The equation used here is a transpose of the one derived in Lengyel's
    // book, since multiplication is done here with the vector on the LHS and the matrix on the RHS.

    // Need proper epsilon
    DEBUG_ASSERT( fabs( oAxis.GetLength() - Float(1.0) ) < 1E-6,
                  "[SMatrix3x3::BuildRotationAboutAxis]: ERROR, oAxis needs to be a UNIT VECTOR\n" );
    
    FLOAT fCosTheta = cosf( fAngle );
    FLOAT fSinTheta = sinf( fAngle );
    m[0][0] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fX * oAxis.m_fX);
    m[1][0] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fY) + (fSinTheta * oAxis.m_fZ);
    m[2][0] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fZ) - (fSinTheta * oAxis.m_fY);

    m[0][1] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fY) - (fSinTheta * oAxis.m_fZ);
    m[1][1] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fY * oAxis.m_fY);
    m[2][1] = (1 - fCosTheta)*(oAxis.m_fY * oAxis.m_fZ) + (fSinTheta * oAxis.m_fX);

    m[0][2] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fZ) + (fSinTheta * oAxis.m_fY);
    m[1][2] = (1 - fCosTheta)*(oAxis.m_fY * oAxis.m_fZ) - (fSinTheta * oAxis.m_fX);
    m[2][2] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fZ * oAxis.m_fZ);
  }
  //----------------------------------------------------------------------------------------------
  // Public : BuildPassiveEulerMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildPassiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {
    FLOAT fCosX = cosf( fPhi ), fCosY = cosf( fTheta ), fCosZ = cosf( fPsi );
    FLOAT fSinX = sinf( fPhi ), fSinY = sinf( fTheta ), fSinZ = sinf( fPsi );

    DEBUG_ASSERT(fTheta <= PI, "BuildEulerRotationMatrix: Theta out of range");
    
    m[0][0] = fCosZ * fCosX - fCosY * fSinX * fSinZ;
    m[1][0] = -fSinZ * fCosX - fCosY * fSinX * fCosZ;
    m[2][0] = fSinY * fSinX;

    m[0][1] = fCosZ * fSinX + fCosY * fCosX * fSinZ;
    m[1][1] = -fSinZ * fSinX + fCosY * fCosX * fCosZ;
    m[2][1] = -fSinY * fCosX;

    m[0][2] = fSinZ * fSinY;
    m[1][2] = fCosZ * fSinY;
    m[2][2] = fCosY;


  }

  //----------------------------------------------------------------------------------------------
  // BuildActiveSmallRotation
  //  -- this is made for infinitesimal approximation
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildActiveSmallRotation(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {
    m[0][0] =   Float( 1 );
    m[1][0] =   fPsi      ;
    m[2][0] = - fTheta    ;

    m[0][1] = - fPsi      ;
    m[1][1] =   Float( 1 );
    m[2][1] =   fPhi      ;

    m[0][2] =   fTheta    ;
    m[1][2] = - fPhi      ;
    m[2][2] =   Float( 1 );
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : BuildActiveEulerMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildActiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {

    FLOAT fCosPhi = cosf( fPhi ), fCosTheta = cosf( fTheta ), fCosPsi = cosf( fPsi );
    FLOAT fSinPhi = sinf( fPhi ), fSinTheta = sinf( fTheta ), fSinPsi = sinf( fPsi );

    DEBUG_ASSERT(fTheta <= PI, "BuildEulerRotationMatrix: Theta out of range");

    m[0][0] = fCosPhi * fCosPsi - fSinPhi * fCosTheta * fSinPsi;
    m[1][0] = fSinPhi * fCosPsi + fCosPhi * fCosTheta * fSinPsi;
    m[2][0] = fSinTheta * fSinPsi;

    m[0][1] = -fCosPhi * fSinPsi - fSinPhi * fCosTheta * fCosPsi;
    m[1][1] = -fSinPhi * fSinPsi + fCosPhi * fCosTheta * fCosPsi;
    m[2][1] = fSinTheta * fCosPsi;

    m[0][2] = fSinPhi * fSinTheta;
    m[1][2] = -fCosPhi * fSinTheta;
    m[2][2] = fCosTheta;

  }
  
  //----------------------------------------------------------------------------------------------
  // Public : GetEulerAngles
  // return a vector of Euler angles in ZYZ convention in radians
  //
  //  This is taken from Bob's code
  //----------------------------------------------------------------------------------------------
  SVector3 SMatrix3x3::GetEulerAngles() const
  {
    SVector3 oEulerAngles;

    Float fCosThresh = 0.999999; 
    

    if( m[2][2] > fCosThresh )	                    //	 is chi approx 0.?
    {
      oEulerAngles.m_fX = Float( 0 );				// set omega and chi to zero
      oEulerAngles.m_fY = Float( 0 );
      oEulerAngles.m_fZ = atan2( m[1][0], m[0][0] );
    }
    else if ( m[2][2] < - fCosThresh )              //  is chi approx pi?
    {
      oEulerAngles.m_fX = Float( 0 );
      oEulerAngles.m_fY = PI;
      oEulerAngles.m_fZ = atan2( m[0][1], m[0][0] );
    }
    else                                            //  chi is not zero or pi
    { 
      oEulerAngles.m_fX = atan2( m[0][2], -m[1][2] );
      oEulerAngles.m_fY = atan2( sqrt( m[2][0] * m[2][0] +m[2][1] * m[2][1] ), m[2][2] );
      oEulerAngles.m_fZ = atan2( m[2][0], m[2][1] );
    }

    // bring back to proper region
    for (Int i = 0; i < 3; i ++ ) //			  % atan2 returns in [-pi:pi], we want [0:2pi]
      if( oEulerAngles[i] < Float (0 ) )
        oEulerAngles[i]+= 2 * PI;
    
 
    
    return oEulerAngles;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : BuildProjectionMatrix  -- oProjDir is a unit vector 
  //----------------------------------------------------------------------------------------------  
  void SMatrix3x3::BuildProjectionMatrix(const SVector3 &oProjDir )
  {
    m[0][0] = oProjDir.m_fX * oProjDir.m_fX;
    m[1][1] = oProjDir.m_fY * oProjDir.m_fY;
    m[2][2] = oProjDir.m_fZ * oProjDir.m_fZ;

    m[0][1] = m[1][0] = oProjDir.m_fX * oProjDir.m_fY;
    m[0][2] = m[2][0] = oProjDir.m_fX * oProjDir.m_fZ;
    m[1][2] = m[2][1] = oProjDir.m_fY * oProjDir.m_fZ;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Trace
  //----------------------------------------------------------------------------------------------  
  FLOAT SMatrix3x3::Trace() const
  {
    FLOAT fRet;
    fRet = m[0][0] + m[1][1] + m[2][2];
    return fRet;
  }

  //----------------------------------------------------------------------------------------------
  // Public : Transpose
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::Transpose()
  {
    // diagonals do not change, and we only have to swap
    // the upper triangle
    for ( int nRow = 0; nRow < 2; nRow++ )
    {
      for ( int nCol = nRow +1; nCol < 3; nCol++ )
      {
        // Swap
        FLOAT fTemp = m[nRow][nCol];
        m[nRow][nCol] = m[nCol][nRow];
        m[nCol][nRow] = fTemp;
      }
    }
  }

  //----------------------------------------------------------------------------------------------
  //  operator* float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator*( Float f ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] * f;
    oRes.m[0][1] =  m[0][1] * f;
    oRes.m[0][2] =  m[0][2] * f;
    
    oRes.m[1][0] =  m[1][0] * f;
    oRes.m[1][1] =  m[1][1] * f;
    oRes.m[1][2] =  m[1][2] * f;
  
    oRes.m[2][0] =  m[2][0] * f;
    oRes.m[2][1] =  m[2][1] * f;
    oRes.m[2][2] =  m[2][2] * f;
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  //  operator *= float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 & SMatrix3x3::operator*= ( Float f )
  {
    m[0][0] *= f;
    m[0][1] *= f;
    m[0][2] *= f;
   
    m[1][0] *= f;
    m[1][1] *= f;
    m[1][2] *= f;
  
    m[2][0] *= f;
    m[2][1] *= f;
    m[2][2] *= f;
    
    return *this; 
  }

    //----------------------------------------------------------------------------------------------
  //  operator/ float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator/( Float f ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] / f;
    oRes.m[0][1] =  m[0][1] / f;
    oRes.m[0][2] =  m[0][2] / f;
    
    oRes.m[1][0] =  m[1][0] / f;
    oRes.m[1][1] =  m[1][1] / f;
    oRes.m[1][2] =  m[1][2] / f;
 
    oRes.m[2][0] =  m[2][0] / f;
    oRes.m[2][1] =  m[2][1] / f;
    oRes.m[2][2] =  m[2][2] / f;
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  //  operator /= float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 & SMatrix3x3::operator/= ( Float f )
  {
    m[0][0] /= f;
    m[0][1] /= f;
    m[0][2] /= f;
   
    m[1][0] /= f;
    m[1][1] /= f;
    m[1][2] /= f;
  
    m[2][0] /= f;
    m[2][1] /= f;
    m[2][2] /= f;
    
    return *this; 
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator+
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator+( const SMatrix3x3& oRHS ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] + oRHS.m[0][0];
    oRes.m[0][1] =  m[0][1] + oRHS.m[0][1];
    oRes.m[0][2] =  m[0][2] + oRHS.m[0][2];
   
    oRes.m[1][0] =  m[1][0] + oRHS.m[1][0];
    oRes.m[1][1] =  m[1][1] + oRHS.m[1][1];
    oRes.m[1][2] =  m[1][2] + oRHS.m[1][2];
  
    oRes.m[2][0] =  m[2][0] + oRHS.m[2][0];
    oRes.m[2][1] =  m[2][1] + oRHS.m[2][1];
    oRes.m[2][2] =  m[2][2] + oRHS.m[2][2];
    
    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator+=
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 & SMatrix3x3::operator+=( const SMatrix3x3& oRHS ) 
  {
    m[0][0] +=  oRHS.m[0][0];
    m[0][1] +=  oRHS.m[0][1];
    m[0][2] +=  oRHS.m[0][2];
   
    m[1][0] +=  oRHS.m[1][0];
    m[1][1] +=  oRHS.m[1][1];
    m[1][2] +=  oRHS.m[1][2];
  
    m[2][0] +=  oRHS.m[2][0];
    m[2][1] +=  oRHS.m[2][1];
    m[2][2] +=  oRHS.m[2][2];
    
    return *this;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator-
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator-( const SMatrix3x3& oRHS ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] - oRHS.m[0][0];
    oRes.m[0][1] =  m[0][1] - oRHS.m[0][1];
    oRes.m[0][2] =  m[0][2] - oRHS.m[0][2];
   
    oRes.m[1][0] =  m[1][0] - oRHS.m[1][0];
    oRes.m[1][1] =  m[1][1] - oRHS.m[1][1];
    oRes.m[1][2] =  m[1][2] - oRHS.m[1][2];
  
    oRes.m[2][0] =  m[2][0] - oRHS.m[2][0];
    oRes.m[2][1] =  m[2][1] - oRHS.m[2][1];
    oRes.m[2][2] =  m[2][2] - oRHS.m[2][2];
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator*( const SMatrix3x3& oRHS ) const
  {
    SMatrix3x3 oRes;

    // Matrix multiplication
    oRes.m[0][0] = m[0][0] * oRHS.m[0][0] + m[0][1] * oRHS.m[1][0] + m[0][2] * oRHS.m[2][0];
    oRes.m[0][1] = m[0][0] * oRHS.m[0][1] + m[0][1] * oRHS.m[1][1] + m[0][2] * oRHS.m[2][1];
    oRes.m[0][2] = m[0][0] * oRHS.m[0][2] + m[0][1] * oRHS.m[1][2] + m[0][2] * oRHS.m[2][2];
   

    oRes.m[1][0] = m[1][0] * oRHS.m[0][0] + m[1][1] * oRHS.m[1][0] + m[1][2] * oRHS.m[2][0];
    oRes.m[1][1] = m[1][0] * oRHS.m[0][1] + m[1][1] * oRHS.m[1][1] + m[1][2] * oRHS.m[2][1];
    oRes.m[1][2] = m[1][0] * oRHS.m[0][2] + m[1][1] * oRHS.m[1][2] + m[1][2] * oRHS.m[2][2];
  

    oRes.m[2][0] = m[2][0] * oRHS.m[0][0] + m[2][1] * oRHS.m[1][0] + m[2][2] * oRHS.m[2][0];
    oRes.m[2][1] = m[2][0] * oRHS.m[0][1] + m[2][1] * oRHS.m[1][1] + m[2][2] * oRHS.m[2][1];
    oRes.m[2][2] = m[2][0] * oRHS.m[0][2] + m[2][1] * oRHS.m[1][2] + m[2][2] * oRHS.m[2][2];
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SVector3 SMatrix3x3::operator*( const SVector3 &oRHS ) const
  {
    SVector3 oRes;

    oRes.m_fX = m[0][0] * oRHS.m_fX + m[0][1] * oRHS.m_fY + m[0][2] * oRHS.m_fZ;
    oRes.m_fY = m[1][0] * oRHS.m_fX + m[1][1] * oRHS.m_fY + m[1][2] * oRHS.m_fZ;
    oRes.m_fZ = m[2][0] * oRHS.m_fX + m[2][1] * oRHS.m_fY + m[2][2] * oRHS.m_fZ;
    
    return oRes;
    
  }

  //----------------------------------------------------------------------------------------------
  // operator<<
  //----------------------------------------------------------------------------------------------
  std::ostream & operator<< ( std::ostream & os, const SMatrix3x3 &m )
  {
    for ( int i = 0; i < 3; i ++ )
    {
      for ( int j = 0; j < 3; j ++ )
      {
        if( fabs( m.m[i][j] ) > 0.00001 )
          os << m.m[i][j] << " ";
        else
          os << 0 << " ";
      }
      os << std::endl;
    }
    return os;
  }
  
  //----------------------------------------------------------------------------------------------
  // operator<<
  //----------------------------------------------------------------------------------------------
  std::ostream & operator<< ( std::ostream & os, const SVector3 &v )
  {
    os << v.m_fX << " " << v.m_fY << " " << v.m_fZ;
    return os;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : SMatrix4x4
  //----------------------------------------------------------------------------------------------
  SMatrix4x4::SMatrix4x4()
  {
  }

  //----------------------------------------------------------------------------------------------
  // Public : SMatrix4x4
  //----------------------------------------------------------------------------------------------
  SMatrix4x4::SMatrix4x4( FLOAT f11, FLOAT f12, FLOAT f13, FLOAT f14,
                          FLOAT f21, FLOAT f22, FLOAT f23, FLOAT f24,
                          FLOAT f31, FLOAT f32, FLOAT f33, FLOAT f34,
                          FLOAT f41, FLOAT f42, FLOAT f43, FLOAT f44 )
  {
  }

  //----------------------------------------------------------------------------------------------
  // Public : SetIdentity
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::SetIdentity()
  {
    memset( m, 0, sizeof(FLOAT) * 16 );
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
  }

  //----------------------------------------------------------------------------------------------
  // Public : Transpose
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::Transpose()
  {
    // diagonals do not change, and we only have to swap
    // the upper triangle
    for ( int nRow = 0; nRow < 3; nRow++ )
    {
      for ( int nCol = nRow +1; nCol < 4; nCol++ )
      {
        // Swap
        FLOAT fTemp = m[nRow][nCol];
        m[nRow][nCol] = m[nCol][nRow];
        m[nCol][nRow] = fTemp;
      }
    }
  }

  //----------------------------------------------------------------------------------------------
  // Public : Inverse
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::Inverse()
  {
    // Gauss-Jordan elimination routine from "Physically Based Rendering"
    int indxc[4], indxr[4];
    int ipiv[4] = { 0, 0, 0, 0 };
    //float minv[4][4];
    //memcpy(minv, m, 4*4*sizeof(float));
    for (int i = 0; i < 4; i++)
    {
      int irow = -1, icol = -1;
      float big = 0.;
      // Choose pivot
      for (int j = 0; j < 4; j++)
      {
        if (ipiv[j] != 1)
        {
          for (int k = 0; k < 4; k++)
          {
            if (ipiv[k] == 0)
            {
              if (fabsf(m[j][k]) >= big)
              {
                big = float(fabsf(m[j][k]));
                irow = j;
                icol = k;
              }
            }
            else if (ipiv[k] > 1)
              DEBUG_ASSERT(0, "Singular matrix in MatrixInvert");
          }
        }
      }
      ++ipiv[icol];
      // Swap rows _irow_ and _icol_ for pivot
      if (irow != icol)
      {
        for (int k = 0; k < 4; ++k)
        {
          // swap
          FLOAT tmp = m[irow][k];
          m[irow][k] = m[icol][k];
          m[icol][k] = tmp;

        }
      }
      indxr[i] = irow;
      indxc[i] = icol;
      if ( fabs( m[icol][icol] ) < EPSILON )
        DEBUG_ASSERT(0, "Singular matrix in MatrixInvert");
      // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
      float pivinv = 1.f / m[icol][icol];
      m[icol][icol] = 1.f;
      for (int j = 0; j < 4; j++)
        m[icol][j] *= pivinv;
      // Subtract this row from others to zero out their columns
      for (int j = 0; j < 4; j++)
      {
        if (j != icol)
        {
          float save = m[j][icol];
          m[j][icol] = 0;
          for (int k = 0; k < 4; k++)
            m[j][k] -= m[icol][k]*save;
        }
      }
    }
    // Swap columns to reflect permutation
    for (int j = 3; j >= 0; j--)
    {
      if (indxr[j] != indxc[j])
      {
        for (int k = 0; k < 4; k++)
        {
          // swap
          FLOAT tmp = m[k][indxr[j]];
          m[k][indxr[j]] = m[k][indxc[j]];
          m[k][indxc[j]] = tmp;
        }
      }
    }
    //return new Matrix4x4(minv);
  }

  //----------------------------------------------------------------------------------------------
  // Public : BuildRotationAboutAxis.
  //
  // Description:  An active transformation matrix that will perform the rotation in a positive sence.
  // (i.e., follows the right hand rule.)
  //
  // Precondition:  oAxis must be a UNIT VECTOR
  //
  //
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::BuildRotationAboutAxis( const SVector3 &oAxis, FLOAT fAngle )
  {
    // 1. Perform transformations which align rotation axis with one of coordinate axis (x, y, z).
    //	  This is essentially a transformation matrix with the given vector as one of the rows, and
    //    orthogonal vectors as the other rows.
    // 2. Perform rotation about the axis
    // 3. Do inverse of (1)

    // Alternatively, use the method described on page 79-80 of Eric Lengyel's "Mathematics for 3d game programming
    // and computer graphics". The core idea involves finding basis vectors and writing the rotation as a
    // linear combination of these vectors.  The equation used here is a transpose of the one derived in Lengyel's
    // book, since multiplication is done here with the vector on the LHS and the matrix on the RHS.

    // Need proper epsilon
    DEBUG_ASSERT( fabs( oAxis.GetLength() - Float(1.0) ) < 1e-6,
                  "[SMatrix4x4::BuildRotationAboutAxis]: ERROR, oAxis needs to be a UNIT VECTOR\n" );
    
    
    FLOAT fCosTheta = cosf( fAngle );
    FLOAT fSinTheta = sinf( fAngle );

    m[0][0] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fX * oAxis.m_fX);
    m[1][0] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fY) + (fSinTheta * oAxis.m_fZ);
    m[2][0] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fZ) - (fSinTheta * oAxis.m_fY);

    m[0][1] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fY) - (fSinTheta * oAxis.m_fZ);
    m[1][1] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fY * oAxis.m_fY);
    m[2][1] = (1 - fCosTheta)*(oAxis.m_fY * oAxis.m_fZ) + (fSinTheta * oAxis.m_fX);

    m[0][2] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fZ) + (fSinTheta * oAxis.m_fY);
    m[1][2] = (1 - fCosTheta)*(oAxis.m_fY * oAxis.m_fZ) - (fSinTheta * oAxis.m_fX);
    m[2][2] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fZ * oAxis.m_fZ);

    m[0][3] = m[1][3] = m[2][3] = m[3][0] = m[3][1] = m[3][2] = 0;
    m[3][3] = 1;
  }

  //----------------------------------------------------------------------------------------------
  // Public : BuildRotationRollPitchYaw
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::BuildRotationRollPitchYaw( FLOAT fRoll, FLOAT fPitch, FLOAT fYaw )
  {
    // Same as Rotation(Z) * Rotation(Y) * Rotation(X)
    // Roll = Z, Pitch = X, Yaw = Y
    FLOAT fCosX = cosf( fPitch ), fCosY = cosf( fYaw ), fCosZ = cosf( fRoll );
    FLOAT fSinX = sinf( fPitch ), fSinY = sinf( fYaw ), fSinZ = sinf( fRoll );
    m[0][0] = fCosZ * fCosY + fSinZ * fSinX * fSinY;
    m[1][0] = fSinZ * fCosX;
    m[2][0] = fCosZ * -fSinY + fSinZ * fSinX * fCosY;
    m[0][1] = -fSinZ * fCosY + fCosZ * fSinX * fSinY;
    m[1][1] = fCosZ * fCosX;
    m[2][1] = fSinZ * fSinY + fCosZ * fSinX * fCosY;
    m[0][2] = fCosX * fSinY;
    m[1][2] = -fSinX;
    m[2][2] = fCosX * fCosY;
    m[3][0] = m[1][3] = m[2][3] = m[0][3] = m[3][1] = m[3][2] = 0;
    m[3][3] = 1;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : SetPassiveEulerMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::SetPassiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi )
  {
    FLOAT fCosX = cosf( fPhi ), fCosY = cosf( fTheta ), fCosZ = cosf( fPsi );
    FLOAT fSinX = sinf( fPhi ), fSinY = sinf( fTheta ), fSinZ = sinf( fPsi );


    DEBUG_ASSERT(fTheta <= PI, "SetPassiveEulerMatrix: Theta out of range");

    m[0][0] = fCosZ * fCosX - fCosY * fSinX * fSinZ;
    m[1][0] = -fSinZ * fCosX - fCosY * fSinX * fCosZ;
    m[2][0] = fSinY * fSinX;

    m[0][1] = fCosZ * fSinX + fCosY * fCosX * fSinZ;
    m[1][1] = -fSinZ * fSinX + fCosY * fCosX * fCosZ;
    m[2][1] = -fSinY * fCosX;

    m[0][2] = fSinZ * fSinY;
    m[1][2] = fCosZ * fSinY;
    m[2][2] = fCosY;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : BuildEulerRotationMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::BuildPassiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {
    SetPassiveEulerMatrix( fPhi, fTheta, fPsi );
    m[3][0] = m[1][3] = m[2][3] = m[0][3] = m[3][1] = m[3][2] = 0;
    m[3][3] = 1;
  }

  //----------------------------------------------------------------------------------------------
  // Public : SetActiveEulerMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::SetActiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi )
  {
    FLOAT fCosPhi = cosf( fPhi ), fCosTheta = cosf( fTheta ), fCosPsi = cosf( fPsi );
    FLOAT fSinPhi = sinf( fPhi ), fSinTheta = sinf( fTheta ), fSinPsi = sinf( fPsi );

    DEBUG_ASSERT(fTheta <= PI, "BuildEulerRotationMatrix: Theta out of range");

    m[0][0] = fCosPhi * fCosPsi - fSinPhi * fCosTheta * fSinPsi;
    m[1][0] = fSinPhi * fCosPsi + fCosPhi * fCosTheta * fSinPsi;
    m[2][0] = fSinTheta * fSinPsi;

    m[0][1] = -fCosPhi * fSinPsi - fSinPhi * fCosTheta * fCosPsi;
    m[1][1] = -fSinPhi * fSinPsi + fCosPhi * fCosTheta * fCosPsi;
    m[2][1] = fSinTheta * fCosPsi;

    m[0][2] = fSinPhi * fSinTheta;
    m[1][2] = -fCosPhi * fSinTheta;
    m[2][2] = fCosTheta;
  }

  
  //----------------------------------------------------------------------------------------------
  // Public : BuildActiveEulerMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::BuildActiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {
    SetActiveEulerMatrix( fPhi, fTheta, fPsi );
    m[3][0] = m[1][3] = m[2][3] = m[0][3] = m[3][1] = m[3][2] = 0;
    m[3][3] = 1;

  }

  //----------------------------------------------------------------------------------------------
  // Public : BuildTranslation
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::BuildTranslation( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    memset( m, 0, sizeof(FLOAT) * 16 );
    m[0][0] = m[1][1] = m[2][2] = 1;
    SetTranslation( fX, fY, fZ );
  }
  //----------------------------------------------------------------------------------------------
  // Public : SetTranslation
  //
  // !!This only sets the translation, but nothing else!!
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::SetTranslation( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    m[0][3] = fX;
    m[1][3] = fY;
    m[2][3] = fZ;
    m[3][3] = 1;
  }

  //----------------------------------------------------------------------------------------------
  //  Set rotation component of this matrix
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::SetRotation( const SMatrix3x3 & oRHS)
  {
    m[0][0] = oRHS.m[0][0];
    m[0][1] = oRHS.m[0][1];
    m[0][2] = oRHS.m[0][2];

    m[1][0] = oRHS.m[1][0];
    m[1][1] = oRHS.m[1][1];
    m[1][2] = oRHS.m[1][2];

    m[2][0] = oRHS.m[2][0];
    m[2][1] = oRHS.m[2][1];
    m[2][2] = oRHS.m[2][2];
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : BuildScale
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::BuildScale( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    memset( m, 0, sizeof(FLOAT) * 16 );
    m[0][0] = fX;
    m[1][1] = fY;
    m[2][2] = fZ;
    m[3][3] = 1;
  }


  //----------------------------------------------------------------------------------------------
  // Public : Translate
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::Translate( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    m[0][3] += fX;
    m[1][3] += fY;
    m[2][3] += fZ;
  }

  //----------------------------------------------------------------------------------------------
  // Public : rotate
  //----------------------------------------------------------------------------------------------
  void SMatrix4x4::Rotate( const SMatrix3x3 & oRHS )
  {
    m[0][0] = m[0][0]*oRHS.m[0][0] + m[0][1]*oRHS.m[1][0] + m[0][2]*oRHS.m[2][0];
    m[0][1] = m[0][0]*oRHS.m[0][1] + m[0][1]*oRHS.m[1][1] + m[0][2]*oRHS.m[2][1];
    m[0][2] = m[0][0]*oRHS.m[0][2] + m[0][1]*oRHS.m[1][2] + m[0][2]*oRHS.m[2][2]; 
        
    m[1][0] = m[1][0]*oRHS.m[0][0] + m[1][1]*oRHS.m[1][0] + m[1][2]*oRHS.m[2][0];
    m[1][1] = m[1][0]*oRHS.m[0][1] + m[1][1]*oRHS.m[1][1] + m[1][2]*oRHS.m[2][1]; 
    m[1][2] = m[1][0]*oRHS.m[0][2] + m[1][1]*oRHS.m[1][2] + m[1][2]*oRHS.m[2][2];
    
    m[2][0] = m[2][0]*oRHS.m[0][0] + m[2][1]*oRHS.m[1][0] + m[2][2]*oRHS.m[2][0];
    m[2][1] = m[2][0]*oRHS.m[0][1] + m[2][1]*oRHS.m[1][1] + m[2][2]*oRHS.m[2][1]; 
    m[2][2] = m[2][0]*oRHS.m[0][2] + m[2][1]*oRHS.m[1][2] + m[2][2]*oRHS.m[2][2]; 
  }

  
  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SMatrix4x4 SMatrix4x4::operator*( const SMatrix4x4& oRHS ) const
  {
    SMatrix4x4 oRes;

    // Matrix multiplication
    oRes.m[0][0] = m[0][0]*oRHS.m[0][0] + m[0][1]*oRHS.m[1][0] + m[0][2]*oRHS.m[2][0] + m[0][3]*oRHS.m[3][0];
    oRes.m[0][1] = m[0][0]*oRHS.m[0][1] + m[0][1]*oRHS.m[1][1] + m[0][2]*oRHS.m[2][1] + m[0][3]*oRHS.m[3][1];
    oRes.m[0][2] = m[0][0]*oRHS.m[0][2] + m[0][1]*oRHS.m[1][2] + m[0][2]*oRHS.m[2][2] + m[0][3]*oRHS.m[3][2];
    oRes.m[0][3] = m[0][0]*oRHS.m[0][3] + m[0][1]*oRHS.m[1][3] + m[0][2]*oRHS.m[2][3] + m[0][3]*oRHS.m[3][3];

    oRes.m[1][0] = m[1][0]*oRHS.m[0][0] + m[1][1]*oRHS.m[1][0] + m[1][2]*oRHS.m[2][0] + m[1][3]*oRHS.m[3][0];
    oRes.m[1][1] = m[1][0]*oRHS.m[0][1] + m[1][1]*oRHS.m[1][1] + m[1][2]*oRHS.m[2][1] + m[1][3]*oRHS.m[3][1];
    oRes.m[1][2] = m[1][0]*oRHS.m[0][2] + m[1][1]*oRHS.m[1][2] + m[1][2]*oRHS.m[2][2] + m[1][3]*oRHS.m[3][2];
    oRes.m[1][3] = m[1][0]*oRHS.m[0][3] + m[1][1]*oRHS.m[1][3] + m[1][2]*oRHS.m[2][3] + m[1][3]*oRHS.m[3][3];

    oRes.m[2][0] = m[2][0]*oRHS.m[0][0] + m[2][1]*oRHS.m[1][0] + m[2][2]*oRHS.m[2][0] + m[2][3]*oRHS.m[3][0];
    oRes.m[2][1] = m[2][0]*oRHS.m[0][1] + m[2][1]*oRHS.m[1][1] + m[2][2]*oRHS.m[2][1] + m[2][3]*oRHS.m[3][1];
    oRes.m[2][2] = m[2][0]*oRHS.m[0][2] + m[2][1]*oRHS.m[1][2] + m[2][2]*oRHS.m[2][2] + m[2][3]*oRHS.m[3][2];
    oRes.m[2][3] = m[2][0]*oRHS.m[0][3] + m[2][1]*oRHS.m[1][3] + m[2][2]*oRHS.m[2][3] + m[2][3]*oRHS.m[3][3];

    oRes.m[3][0] = m[3][0]*oRHS.m[0][0] + m[3][1]*oRHS.m[1][0] + m[3][2]*oRHS.m[2][0] + m[3][3]*oRHS.m[3][0];
    oRes.m[3][1] = m[3][0]*oRHS.m[0][1] + m[3][1]*oRHS.m[1][1] + m[3][2]*oRHS.m[2][1] + m[3][3]*oRHS.m[3][1];
    oRes.m[3][2] = m[3][0]*oRHS.m[0][2] + m[3][1]*oRHS.m[1][2] + m[3][2]*oRHS.m[2][2] + m[3][3]*oRHS.m[3][2];
    oRes.m[3][3] = m[3][0]*oRHS.m[0][3] + m[3][1]*oRHS.m[1][3] + m[3][2]*oRHS.m[2][3] + m[3][3]*oRHS.m[3][3];

    return oRes;
  }
 

  //----------------------------------------------------------------------------------------------
  // Public : SVector2
  //----------------------------------------------------------------------------------------------
  SVector2::SVector2()
  {}

  //----------------------------------------------------------------------------------------------
  // Public : SVector2
  //----------------------------------------------------------------------------------------------
  SVector2::SVector2( FLOAT fX, FLOAT fY )
  {
    this->m_fX = fX;
    this->m_fY = fY;
  }

  //----------------------------------------------------------------------------------------------
  // Public : SVector3
  //----------------------------------------------------------------------------------------------
  SVector3::SVector3()
  {}

  //----------------------------------------------------------------------------------------------
  // Public : SVector3
  //----------------------------------------------------------------------------------------------
  SVector3::SVector3( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    Set( fX, fY, fZ );
  }

  //----------------------------------------------------------------------------------------------
  // Public : Set
  //----------------------------------------------------------------------------------------------
  void SVector3::Set( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    this->m_fX = fX;
    this->m_fY = fY;
    this->m_fZ = fZ;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Get
  //----------------------------------------------------------------------------------------------
  FLOAT SVector3::Get( UInt nAxis ) const
  {
    if ( nAxis == 0 )
      return this->m_fX;
    else if ( nAxis == 1 )
      return this->m_fY;
    else
      return this->m_fZ;
  }

  //----------------------------------------------------------------------------------------------
  // Public : Set
  //----------------------------------------------------------------------------------------------
  void SVector3::Set( UInt nAxis, FLOAT fVal )
  {
    if ( nAxis == 0 )
      this->m_fX = fVal;
    else if ( nAxis == 1 )
      this->m_fY = fVal;
    else
      this->m_fZ = fVal;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : GetLength
  //----------------------------------------------------------------------------------------------
  FLOAT SVector3::GetLength() const
  {
    return sqrtf(m_fX * m_fX + m_fY * m_fY + m_fZ * m_fZ);
  }

  //----------------------------------------------------------------------------------------------
  // Public : TransformNormal
  //----------------------------------------------------------------------------------------------
  void SVector3::TransformNormal( const SMatrix4x4 &oMatrix )
  {
    FLOAT fX = m_fX;
    FLOAT fY = m_fY;
    FLOAT fZ = m_fZ;
    m_fX = oMatrix.m[0][0] * fX + oMatrix.m[0][1] * fY + oMatrix.m[0][2] * fZ;
    m_fY = oMatrix.m[1][0] * fX + oMatrix.m[1][1] * fY + oMatrix.m[1][2] * fZ;
    m_fZ = oMatrix.m[2][0] * fX + oMatrix.m[2][1] * fY + oMatrix.m[2][2] * fZ;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Transform
  //----------------------------------------------------------------------------------------------
  void SVector3::Transform( const SMatrix4x4 &oMatrix )
  {
    FLOAT fX = m_fX;
    FLOAT fY = m_fY;
    FLOAT fZ = m_fZ;
    m_fX = fX * oMatrix.m[0][0] + fY * oMatrix.m[0][1] + fZ * oMatrix.m[0][2] + oMatrix.m[0][3];
    m_fY = fX * oMatrix.m[1][0] + fY * oMatrix.m[1][1] + fZ * oMatrix.m[1][2] + oMatrix.m[1][3];
    m_fZ = fX * oMatrix.m[2][0] + fY * oMatrix.m[2][1] + fZ * oMatrix.m[2][2] + oMatrix.m[2][3];
  }

  //----------------------------------------------------------------------------------------------
  // Public : Transform
  //
  //  Transform via oMatrix * thisVector
  //----------------------------------------------------------------------------------------------
  void SVector3::Transform( const SMatrix3x3 &oMatrix )
  {
    FLOAT fX = m_fX;
    FLOAT fY = m_fY;
    FLOAT fZ = m_fZ;

    m_fX = oMatrix.m[0][0] * fX + oMatrix.m[0][1] * fY + oMatrix.m[0][2] * fZ;
    m_fY = oMatrix.m[1][0] * fX + oMatrix.m[1][1] * fY + oMatrix.m[1][2] * fZ ; 
    m_fZ = oMatrix.m[2][0] * fX + oMatrix.m[2][1] * fY + oMatrix.m[2][2] * fZ ; 
  }

  
  //----------------------------------------------------------------------------------------------
  // Public : Normalize
  //----------------------------------------------------------------------------------------------
  void SVector3::Normalize()
  {
    FLOAT fLength = GetLength();
    m_fX /= fLength;
    m_fY /= fLength;
    m_fZ /= fLength;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator=
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator=( const SVector3 &oRHS )
  {
    m_fX = oRHS.m_fX;
    m_fY = oRHS.m_fY;
    m_fZ = oRHS.m_fZ;
    return *this;
  }

  //----------------------------------------------------------------------------------------------
  //  Public : operator+
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator+( const SVector3& oRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX + oRHS.m_fX;
    oRetObj.m_fY = m_fY + oRHS.m_fY;
    oRetObj.m_fZ = m_fZ + oRHS.m_fZ;
    return oRetObj;
  }

  //----------------------------------------------------------------------------------------------
  //  Public : operator+= 
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator+=( const SVector3 &oRHS )
  {
    m_fX += oRHS.m_fX;
    m_fY += oRHS.m_fY;
    m_fZ += oRHS.m_fZ;
    return *this;
  }
  
  //----------------------------------------------------------------------------------------------
  //  operator -
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator-( ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = - m_fX;
    oRetObj.m_fY = - m_fY;
    oRetObj.m_fZ = - m_fZ;
    return oRetObj;
  }
  //----------------------------------------------------------------------------------------------
  // Public : operator-
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator-( const SVector3& oRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX - oRHS.m_fX;
    oRetObj.m_fY = m_fY - oRHS.m_fY;
    oRetObj.m_fZ = m_fZ - oRHS.m_fZ;
    return oRetObj;
  }

  //----------------------------------------------------------------------------------------------
  //  Public : operator-= 
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator-=( const SVector3 &oRHS )
  {
    m_fX -= oRHS.m_fX;
    m_fY -= oRHS.m_fY;
    m_fZ -= oRHS.m_fZ;
    return *this;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator*( FLOAT fRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX * fRHS;
    oRetObj.m_fY = m_fY * fRHS;
    oRetObj.m_fZ = m_fZ * fRHS;
    return oRetObj;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator*  (left multiply)
  //----------------------------------------------------------------------------------------------
  SVector3 operator*( FLOAT fLHS, const SVector3 &oRHS )
  {
    SVector3 oRetObj;
    oRetObj.m_fX = fLHS * oRHS.m_fX;
    oRetObj.m_fY = fLHS * oRHS.m_fY;
    oRetObj.m_fZ = fLHS * oRHS.m_fZ;
    return oRetObj;
  }


  //----------------------------------------------------------------------------------------------
  // Public : operator*=
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator*=( FLOAT fRHS )
  {
    m_fX *= fRHS;
    m_fY*= fRHS;
    m_fZ *= fRHS;
    return *this;
  }


  //----------------------------------------------------------------------------------------------
  // Public : operator/
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator/( FLOAT fRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX / fRHS;
    oRetObj.m_fY = m_fY / fRHS;
    oRetObj.m_fZ = m_fZ / fRHS;
    return oRetObj;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator/=
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator/=( FLOAT fRHS )
  {
    m_fX /= fRHS;
    m_fY /= fRHS;
    m_fZ /= fRHS;
    return *this;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator[]
  //----------------------------------------------------------------------------------------------
  FLOAT SVector3::operator[]( Int i) const
  {

    DEBUG_ASSERT( i >= 0 && i < 3, "ERROR:  index out of range for SVector3\n");
    switch(i)
    {
      case 0:
        return m_fX;
      case 1:
        return m_fY;
      case 2:
        return m_fZ;
    };
    // error branch  (throw exception in the future)
    return std::numeric_limits<Float>::signaling_NaN();
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator[]
  //----------------------------------------------------------------------------------------------
  FLOAT& SVector3::operator[](Int i)
  {
    DEBUG_ASSERT( i >= 0 && i < 3, "ERROR:  index out of range for SVector3\n");
    switch(i)
    {
      case 0:
        return m_fX;
      case 1:
        return m_fY;
      case 2:
        return m_fZ;
    };
    // error branch  (throw exception in the future)
    return m_fX;
  }

  //----------------------------------------------------------------------------------------------
  // Cross
  //----------------------------------------------------------------------------------------------
  SVector3 Cross( const SVector3 &oLHS, const SVector3 &oRHS )
  {
    SVector3 oRes;

    oRes.m_fX = oLHS.m_fY * oRHS.m_fZ - oLHS.m_fZ * oRHS.m_fY;
    oRes.m_fY = oLHS.m_fZ * oRHS.m_fX - oLHS.m_fX * oRHS.m_fZ;
    oRes.m_fZ = oLHS.m_fX * oRHS.m_fY - oLHS.m_fY * oRHS.m_fX;

    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Dot
  //----------------------------------------------------------------------------------------------
  FLOAT Dot( const SVector3 &oLHS, const SVector3 &oRHS )
  {
    return oLHS.m_fX * oRHS.m_fX + oLHS.m_fY * oRHS.m_fY + oLHS.m_fZ * oRHS.m_fZ;
  }

  //----------------------------------------------------------------------------------------------
  //  OuterProduct
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 OuterProduct( const SVector3 & oLHS, const SVector3 & oRHS )
  {
    SMatrix3x3 oRes;
    oRes.m[0][0] = oLHS.m_fX * oRHS.m_fX;
    oRes.m[0][1] = oLHS.m_fX * oRHS.m_fY;
    oRes.m[0][2] = oLHS.m_fX * oRHS.m_fZ;

    oRes.m[1][0] = oLHS.m_fY * oRHS.m_fX;
    oRes.m[1][1] = oLHS.m_fY * oRHS.m_fY;
    oRes.m[1][2] = oLHS.m_fY * oRHS.m_fZ;
    
    oRes.m[2][0] = oLHS.m_fZ * oRHS.m_fX;
    oRes.m[2][1] = oLHS.m_fZ * oRHS.m_fY;
    oRes.m[2][2] = oLHS.m_fZ * oRHS.m_fZ;

    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // DegreeToRadian
  //----------------------------------------------------------------------------------------------
  SVector3 DegreeToRadian( const SVector3 & oRHS )
  {
    SVector3 oRes;
    oRes.m_fX =  DEGREE_TO_RADIAN( oRHS.m_fX );
    oRes.m_fY =  DEGREE_TO_RADIAN( oRHS.m_fY );
    oRes.m_fZ =  DEGREE_TO_RADIAN( oRHS.m_fZ );
    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // RadianToDegree
  //----------------------------------------------------------------------------------------------
  SVector3 RadianToDegree( const SVector3 & oRHS )
  {
    SVector3 oRes;
    oRes.m_fX =  RADIAN_TO_DEGREE( oRHS.m_fX );
    oRes.m_fY =  RADIAN_TO_DEGREE( oRHS.m_fY );
    oRes.m_fZ =  RADIAN_TO_DEGREE( oRHS.m_fZ );
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : SVector4
  //----------------------------------------------------------------------------------------------
  SVector4::SVector4()
  {}

  //----------------------------------------------------------------------------------------------
  // Public : SVector4
  //----------------------------------------------------------------------------------------------
  SVector4::SVector4( FLOAT fX, FLOAT fY, FLOAT fZ, FLOAT fW )
  {
    Set( fX, fY, fZ, fW );
  }

  //----------------------------------------------------------------------------------------------
  // Public : Set
  //----------------------------------------------------------------------------------------------
  void SVector4::Set( FLOAT fX, FLOAT fY, FLOAT fZ, FLOAT fW )
  {
    this->m_fX = fX;
    this->m_fY = fY;
    this->m_fZ = fZ;
    this->m_fW = fW;
  }

  //----------------------------------------------------------------------------------------------
  // Public : Transform
  //----------------------------------------------------------------------------------------------
  void SVector4::Transform( const SMatrix4x4 &oMatrix )
  {
    FLOAT fX = m_fX;
    FLOAT fY = m_fY;
    FLOAT fZ = m_fZ;
    FLOAT fW = m_fW;
    m_fX = fX * oMatrix.m[0][0] + fY * oMatrix.m[0][1] + fZ * oMatrix.m[0][2] + fW * oMatrix.m[0][3];
    m_fY = fX * oMatrix.m[1][0] + fY * oMatrix.m[1][1] + fZ * oMatrix.m[1][2] + fW * oMatrix.m[1][3];
    m_fZ = fX * oMatrix.m[2][0] + fY * oMatrix.m[2][1] + fZ * oMatrix.m[2][2] + fW * oMatrix.m[2][3];
    m_fW = fX * oMatrix.m[3][0] + fY * oMatrix.m[3][1] + fZ * oMatrix.m[3][2] + fW * oMatrix.m[3][3];
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Overlaps
  //----------------------------------------------------------------------------------------------
  bool SBoundingBox::Overlaps( const SBoundingBox &oBox ) const
  {
    // Calculate the center point of each bounding box
    SVector3 oCenterPointA, oCenterPointB;
    oCenterPointA.m_fX = (m_oBoxMax.m_fX + m_oBoxMin.m_fX) / 2.0f;
    oCenterPointA.m_fY = (m_oBoxMax.m_fY + m_oBoxMin.m_fY) / 2.0f;
    oCenterPointA.m_fZ = (m_oBoxMax.m_fZ + m_oBoxMin.m_fZ) / 2.0f;
    oCenterPointB.m_fX = (oBox.m_oBoxMax.m_fX + oBox.m_oBoxMin.m_fX) / 2.0f;
    oCenterPointB.m_fY = (oBox.m_oBoxMax.m_fY + oBox.m_oBoxMin.m_fY) / 2.0f;
    oCenterPointB.m_fZ = (oBox.m_oBoxMax.m_fZ + oBox.m_oBoxMin.m_fZ) / 2.0f;

    // Calculate distance between center point and the side of the box
    SVector3 oExtentA, oExtentB;
    oExtentA.m_fX = fabs( oCenterPointA.m_fX - m_oBoxMin.m_fX );
    oExtentA.m_fY = fabs( oCenterPointA.m_fY - m_oBoxMin.m_fY );
    oExtentA.m_fZ = fabs( oCenterPointA.m_fZ - m_oBoxMin.m_fZ );
    oExtentB.m_fX = fabs( oCenterPointB.m_fX - oBox.m_oBoxMin.m_fX );
    oExtentB.m_fY = fabs( oCenterPointB.m_fY - oBox.m_oBoxMin.m_fY );
    oExtentB.m_fZ = fabs( oCenterPointB.m_fZ - oBox.m_oBoxMin.m_fZ );

    //vector from A to B
    SVector3 oSeparatingAxis;
    oSeparatingAxis.m_fX = oCenterPointB.m_fX - oCenterPointA.m_fX;
    oSeparatingAxis.m_fY = oCenterPointB.m_fY - oCenterPointA.m_fY;
    oSeparatingAxis.m_fZ = oCenterPointB.m_fZ - oCenterPointA.m_fZ;

    // If the two boxes are disjoint (do not overlap), then at least one
    // of these will form a separating axis.
    if ( fabs(oSeparatingAxis.m_fX) <= (oExtentA.m_fX + oExtentB.m_fX) )
    {
      if ( fabs(oSeparatingAxis.m_fY) <= (oExtentA.m_fY + oExtentB.m_fY) )
      {
        if ( fabs(oSeparatingAxis.m_fZ) <= (oExtentA.m_fZ + oExtentB.m_fZ) )
        {
          return true;
        }
      }
    }

    return false;
  }


  //----------------------------------------------------------------------------------------------
  // Public : CalculateNormal
  //----------------------------------------------------------------------------------------------
  void STriangle::CalculateNormal( SVector3 &oNormal ) const
  {
    SVector3 oEdge1 = m_oPt2 - m_oPt1;
    SVector3 oEdge2 = m_oPt3 - m_oPt1;
    oNormal = Cross( oEdge1, oEdge2 );
  }

  //----------------------------------------------------------------------------------------------
  // Public : CalculateNormalizedNormal
  //----------------------------------------------------------------------------------------------
  void STriangle::CalculateNormalizedNormal( SVector3 &oNormal ) const
  {
    SVector3 oEdge1 = m_oPt2 - m_oPt1;
    SVector3 oEdge2 = m_oPt3 - m_oPt1;
    oNormal = Cross( oEdge1, oEdge2 );
    oNormal.Normalize();
  }

  //----------------------------------------------------------------------------------------------
  // Public : CalculateBoundingBox
  //----------------------------------------------------------------------------------------------
  void STriangle::CalculateBoundingBox( SBoundingBox &oBox ) const
  {
    SVector3 &oBoxMin = oBox.m_oBoxMin;
    oBoxMin.m_fX = std::min( m_oPt1.m_fX, m_oPt2.m_fX );
    oBoxMin.m_fX = std::min( oBoxMin.m_fX, m_oPt3.m_fX );
    oBoxMin.m_fY = std::min( m_oPt1.m_fY, m_oPt2.m_fY );
    oBoxMin.m_fY = std::min( oBoxMin.m_fY, m_oPt3.m_fY );
    oBoxMin.m_fZ = std::min( m_oPt1.m_fZ, m_oPt2.m_fZ );
    oBoxMin.m_fZ = std::min( oBoxMin.m_fZ, m_oPt3.m_fZ );
    SVector3 &oBoxMax = oBox.m_oBoxMax;
    oBoxMax.m_fX = std::max( m_oPt1.m_fX, m_oPt2.m_fX );
    oBoxMax.m_fX = std::max( oBoxMax.m_fX, m_oPt3.m_fX );
    oBoxMax.m_fY = std::max( m_oPt1.m_fY, m_oPt2.m_fY );
    oBoxMax.m_fY = std::max( oBoxMax.m_fY, m_oPt3.m_fY );
    oBoxMax.m_fZ = std::max( m_oPt1.m_fZ, m_oPt2.m_fZ );
    oBoxMax.m_fZ = std::max( oBoxMax.m_fZ, m_oPt3.m_fZ );
  }

  //----------------------------------------------------------------------------------------------
  // Public : Normalize
  //----------------------------------------------------------------------------------------------
  void CPlane::Normalize()
  {
    FLOAT fLength = sqrtf(m_fA * m_fA + m_fB * m_fB + m_fC * m_fC);
    m_fA /= fLength;
    m_fB /= fLength;
    m_fC /= fLength;
  }
  //----------------------------------------------------------------------------------------------
  // Public : CalculateBoundingBox
  //----------------------------------------------------------------------------------------------
  void CSphere::CalculateBoundingBox( SBoundingBox &oBox ) const
  {
    oBox.m_oBoxMin.Set( m_oCenter.m_fX-m_fRadius, m_oCenter.m_fY-m_fRadius, m_oCenter.m_fZ-m_fRadius );
    oBox.m_oBoxMax.Set( m_oCenter.m_fX+m_fRadius, m_oCenter.m_fY+m_fRadius, m_oCenter.m_fZ+m_fRadius );
  }

  
  //----------------------------------------------------------------------------------------------
  // QuadraticFormula
  //----------------------------------------------------------------------------------------------
  bool QuadraticFormula( FLOAT fA, FLOAT fB, FLOAT fC, FLOAT &fRes1, FLOAT &fRes2 )
  {
    FLOAT fRootSquared = fB*fB - 4*fA*fC;
    if ( fRootSquared >= 0 && fabs(fA) > EPSILON )
    {
      FLOAT fRoot = sqrtf( fRootSquared );
      FLOAT fD = 1 / (2 * fA);
      fRes1 = ( -fB + fRoot ) * fD;
      fRes2 = ( -fB - fRoot ) * fD;

      // real roots found
      return true;
    }

    // no real roots found (complex roots)
    return false;
  }

  //----------------------------------------------------------------------------------------------
  // Barycentric
  //----------------------------------------------------------------------------------------------
  void Barycentric( const STriangle &oTriangle, const SVector3 &oPoint, FLOAT &fU, FLOAT &fV, FLOAT &fW )
  {
    const SVector3 &oPt1 = oTriangle.m_oPt1;
    const SVector3 &oPt2 = oTriangle.m_oPt2;
    const SVector3 &oPt3 = oTriangle.m_oPt3;

    {
      // Code from "Real time collision detection", page 47
      // Explanation "Mathematics for 3d game programming & computer graphics", page 143
      SVector3 oB = oPt2 - oPt1;
      SVector3 oC = oPt3 - oPt1;
      SVector3 oP = oPoint - oPt1;
      FLOAT fD00 = Dot(oB,oB);
      FLOAT fD01 = Dot(oB,oC);
      FLOAT fD11 = Dot(oC,oC);
      FLOAT fD20 = Dot(oP,oB);
      FLOAT fD21 = Dot(oP,oC);
      FLOAT fDenom = fD00 * fD11 - fD01 * fD01;
      fV = (fD11*fD20-fD01*fD21) / fDenom;
      fW = (fD00*fD21-fD01*fD20) / fDenom;
      fU = 1.0f - fV - fW;
    }

    //{
    //	// WARNING: The below function suffers from divide-by-zero if y-coordinate is 0 for all points
    //	// Use the barycentric equation to determine whether point lies in the triangle or not
    //	// Recall barycentric equation is: u(B-A) + v(C-A) = (P-A), and it is inside if u>=0 & v>=0 & u+v<=1
    //	//	There's 3 equations (for each x,y,z), and 2 unknowns.  Substitute: b=(B-A), c=(C-A), p=(P-A)
    //	//	Solving algebraically, we get u=(py*cx-px*cy)/(by*cx-bx*cy), v=(py*bx-px*by)/(cy*bx-cx*by)
    //	SVector3 oB, oC, oP;
    //	Vector3Subtract( oB, oPt2, oPt1 );
    //	Vector3Subtract( oC, oPt3, oPt1 );
    //	Vector3Subtract( oP, oPoint, oPt1 );
    //	fU = (oP.m_fY * oC.m_fX - oP.m_fX * oC.m_fY) / (oB.m_fY*oC.m_fX-oB.m_fX*oC.m_fY);
    //	fV = (oP.m_fY * oB.m_fX - oP.m_fX * oB.m_fY) / (oC.m_fY*oB.m_fX-oC.m_fX*oB.m_fY);
    //	fW = 1.0f - fU - fV;
    //}
  }

  //----------------------------------------------------------------------------------------------
  // Interpolate
  //----------------------------------------------------------------------------------------------
  void Interpolate( SVector3 &oRes, const SVector3 &oParam1, const SVector3 &oParam2, const SVector3 &oParam3, FLOAT fU, FLOAT fV )
  {
    // Barycentric equation: 
    //		P = uA + vB + wC  ==> P = C + u(A-C) + v(B-C)
    //		Substitute a=u(A-C), b=v(B-C)
    SVector3 oA = (oParam1 - oParam3) * fU;
    SVector3 oB = (oParam2 - oParam3) * fV;
    oRes = oParam3 + oA + oB;
  }

  //----------------------------------------------------------------------------------------------
  // Interpolate
  //----------------------------------------------------------------------------------------------
  SVector3 Interpolate( const SVector3 &oParam1, const SVector3 &oParam2, FLOAT fInterp )
  {
    SVector3 oRes;
    oRes = (oParam2 - oParam1) * fInterp;
    oRes += oParam1;
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------
  // Collision namespace
  //----------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------

  namespace Collision
  {

    //----------------------------------------------------------------------------------------------
    // Intersects
    //----------------------------------------------------------------------------------------------
    bool Intersects( const CSphere &oSphereA, const SVector3 &oVelocityA,
                     const CSphere &oSphereB,
                     FLOAT &fCollisionTime1, FLOAT &fCollisionTime2 )
    {
      // Adapted from http://www.gamasutra.com/features/19991018/Gomez_2.htm

      SVector3 oPrevPosA = oSphereA.GetCenter();
      SVector3 oPrevPosB = oSphereB.GetCenter();

      //vector from A0 to A1
      const SVector3 va = oVelocityA;

      //vector from A0 to B0
      SVector3 AB = oPrevPosB - oPrevPosA;

      //relative velocity (in normalized time)
      SVector3 vab = oVelocityA * (FLOAT) ( -1 ) ;

      const FLOAT rab = oSphereA.GetRadius() + oSphereB.GetRadius();

      //u*u coefficient
      const FLOAT a = Dot( vab, vab );

      //u coefficient
      const FLOAT b = 2 * Dot( vab, AB );

      //constant term
      float fABDotAB = Dot( AB, AB );
      float fRabRab = ( rab * rab );
      const FLOAT c = fABDotAB - fRabRab;

      //check if they're currently overlapping
      if( Dot(AB, AB) <= fRabRab )
      {
        fCollisionTime1 = 0;
        fCollisionTime2 = 0;
        return true;
      }

      //check if they hit each other
      // during the frame
      if( QuadraticFormula( a, b, c, fCollisionTime1, fCollisionTime2 ) )
      {
        if( fCollisionTime1 > fCollisionTime2 )
        {
          // Swap
          FLOAT fTemp = fCollisionTime1;
          fCollisionTime1 = fCollisionTime2;
          fCollisionTime2 = fTemp;
        }

        return true;
      }

      return false;
    }

    //----------------------------------------------------------------------------------------------
    // Intersects
    //----------------------------------------------------------------------------------------------
    bool Intersects( const CSphere &oSphereA, const SVector3 &oVelocityA,
                     const CSphere &oSphereB, const SVector3 &oVelocityB,
                     FLOAT &fCollisionTime1, FLOAT &fCollisionTime2 )
    {
      SVector3 oVel = oVelocityA - oVelocityB;
      return Intersects( oSphereA, oVel, oSphereB, fCollisionTime1, fCollisionTime2 );
    }

    //----------------------------------------------------------------------------------------------
    // Public : Intersects
    //----------------------------------------------------------------------------------------------
    bool Intersects( const CSphere &oSphere, const CRay &oRay, FLOAT &fT )
    {
      FLOAT fRes1, fRes2;
      bool bFound = Intersects( oSphere, oRay, fRes1, fRes2 );
      
      // Return the closest intersection point
      fT = fRes1;
      if ( fRes1 < 0.0f || ((fRes2 < fRes1) && (fRes2 > 0.0f)) )
        fT = fRes2;
      
      return bFound;
    }

    //----------------------------------------------------------------------------------------------
    // Public : Intersects
    //----------------------------------------------------------------------------------------------
    bool Intersects( const CSphere &oSphere, const CRay &oRay, FLOAT &fT1, FLOAT &fT2 )
    {
      // Equation of sphere: (x - c)^2 - r^2 = 0
      // Equation of ray: p0 + td (p0 is starting point, d is direction)
      // Plug equation of ray into x of sphere equation.  Using quadratic equation, solve and you get:
      // a = d^2, b = 2*d*(p0-c), c=(p0-c)(p0-c)-r^2

      const SVector3 &oCenter = oSphere.GetCenter();
      const SVector3 &oStart = oRay.GetStart();
      const SVector3 &oDir = oRay.GetDirection();

      // Calculate A
      FLOAT fA = Dot( oDir, oDir );

      // Calculate B
      SVector3 oStartMinusCenter = oStart - oCenter;
      FLOAT fB = Dot( oDir, oStartMinusCenter );
      fB *= 2.0f; 

      // Calculate C
      FLOAT fC = Dot( oStartMinusCenter, oStartMinusCenter );
      fC -= oSphere.GetRadius() * oSphere.GetRadius();

      bool bFound = QuadraticFormula( fA, fB, fC, fT1, fT2 );

      // Intersection behind the ray
      if ( fT1 < 0.0f && fT2 < 0.0f )
        bFound = false;

      return bFound;
    }


    //----------------------------------------------------------------------------------------------
    // Public : Intersects
    //----------------------------------------------------------------------------------------------
    bool Intersects( const SVector3 &oPt1, const SVector3 &oPt2, const SVector3 &oPt3, const CRay &oRay, FLOAT &fT )
    {
      // First, calculate normal of triangle. Assume triangle given in clockwise winding order, and calculations
      // done in left-handed (directx) coordinate system (uses left-hand rule).
      SVector3 oEdge1 = oPt3 - oPt1;
      SVector3 oEdge2 = oPt2 - oPt1;
      SVector3 oNormal = Cross( oEdge2, oEdge1 );
      oNormal.Normalize();

      // Calculate triangle plane
      CPlane oPlane;
      oPlane.m_fA = oNormal.m_fX;
      oPlane.m_fB = oNormal.m_fY;
      oPlane.m_fC = oNormal.m_fZ;
      oPlane.m_fD = -Dot( oNormal, oPt1 );

      // Ray-plane intersection
      if ( !Intersects( oPlane, oRay, fT ) )
        return false;

      // Calculate ray-plane intersection point
      SVector3 oIntersectionPoint;
      SVector3 oTemp =  oRay.GetDirection() * fT;
      oIntersectionPoint = oTemp + oRay.GetStart();

      // Use the barycentric equation to determine whether point lies in the triangle or not
      // Recall barycentric equation is: u(B-A) + v(C-A) = (P-A), and it is inside if u>=0 & v>=0 & u+v<=1
      //	There's 3 equations (for each x,y,z), and 2 unknowns.  Substitute: b=(B-A), c=(C-A), p=(P-A)
      //	Solving algebraically, we get u=(py*cx-px*cy)/(by*cx-bx*cy), v=(py*bx-px*by)/(cy*bx-cx*by)
      SVector3 oB = oPt2 - oPt1;
      SVector3 oC = oPt3 - oPt1;
      SVector3 oP = oIntersectionPoint - oPt1;
      FLOAT fU = (oP.m_fY * oC.m_fX - oP.m_fX * oC.m_fY) / (oB.m_fY*oC.m_fX-oB.m_fX*oC.m_fY);
      FLOAT fV = (oP.m_fY * oB.m_fX - oP.m_fX * oB.m_fY) / (oC.m_fY*oB.m_fX-oC.m_fX*oB.m_fY);

      if ( fU >= 0.0f && fV >= 0.0f && (fU+fV) <= 1.0f )
        return true;

      return false;
    }

    //----------------------------------------------------------------------------------------------
    // Public : Intersects
    //----------------------------------------------------------------------------------------------
    bool Intersects( const CPlane &oPlane, const CRay &oRay, FLOAT &fT )
    {
      // To solve for intersection, substitute parametric ray equation into the plane equation
      // and solve for t
      //	ie: Plane: NX+d=0, Ray: P+tD=y
      //			N(P+tD)+d=0, t=(-d-NP)/(ND)

      SVector3 oPlaneNormal( oPlane.m_fA, oPlane.m_fB, oPlane.m_fC );
      FLOAT fND = Dot( oPlaneNormal, oRay.GetDirection() );

      // TODO: Use proper epsilon
      // This happens if the plane and ray are perpendicular to each other
      if ( fabs(fND) < 0.001f )
        return false;

      FLOAT fNP = Dot( oPlaneNormal, oRay.GetStart() );
      fT = (-oPlane.m_fD - fNP) / (fND);

      if ( fT < 0.0f )
        return false;

      return true;
    }


    //----------------------------------------------------------------------------------------------
    // Public : Intersects
    //----------------------------------------------------------------------------------------------
    bool Intersects( const SBoundingBox &oAABB, const CRay &oRay, FLOAT &fNearT, FLOAT &fFarT )
    {
      fNearT = MIN_FLOAT;
      fFarT = MAX_FLOAT;
      FLOAT fTempNearT, fTempFarT;
      for ( UInt nAxis = 0; nAxis < 3; nAxis++ )
      {
        // Projected ray length (Same as dot(normal-axis, ray-direction))
        // Ray length is 1.  Thus fNearT/fFarT == 1 when the intersection
        // is right at the point where the ray length is 1.
        FLOAT fProjRayLen = oRay.GetDirection().Get( nAxis );

        if ( fabs( fProjRayLen ) > EPSILON )
        {
          // Projected ray length inverse
          FLOAT fInvProjRayLen = 1.0f / fProjRayLen;

          // Distance to plane (this is implictly dot-product with the normal, which simply
          // is 1 for that axis)
          FLOAT fDistToNearPlane = oAABB.m_oBoxMin.Get( nAxis ) - oRay.GetStart().Get( nAxis );
          FLOAT fDistToFarPlane =  oAABB.m_oBoxMax.Get( nAxis ) - oRay.GetStart().Get( nAxis );

          // t-value is calculated by distance to plane divided by total ray length
          // The ratio between the projected distances is the same as the ratio between
          // the actual distances
          fTempNearT = fDistToNearPlane * fInvProjRayLen;
          fTempFarT =  fDistToFarPlane * fInvProjRayLen;

          if ( fTempNearT > fTempFarT ) 
            std::swap( fTempNearT, fTempFarT );

          fNearT = ( fTempNearT > fNearT ) ? fTempNearT : fNearT;
          fFarT  = ( fTempFarT  < fFarT  ) ? fTempFarT  : fFarT;

          if ( fFarT < fNearT ) 
            return false;
        }				
      }

      return true;
    }
    
  }
}
