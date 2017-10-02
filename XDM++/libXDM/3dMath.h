/********************************************************************
	created:	2004/09/25
	authors:    Shan-Min Chao & Frankie Li
	
	purpose:	General purpose math class for standard 3d math.
*********************************************************************/


#ifndef _3DMATH_H_ 
#define _3DMATH_H_

#include "Types.h"
#include "Debug.h"
#include <limits>

typedef Float FLOAT; 


#ifndef PI
#define PI ((FLOAT)  3.14159265358979323846)
#endif

#ifndef SQRT3
#define SQRT3 ((FLOAT)1.73205080756888f)
#endif

#define DEGREE_TO_RADIAN( degree ) ((degree) * (PI / (FLOAT)180.0f))
#define RADIAN_TO_DEGREE( radian ) ((radian) * ( (FLOAT) 180.0f / PI))



namespace GeneralLib
{
  class SVector3;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////

  class SMatrix3x3
  {
  public:
    SMatrix3x3();
    SMatrix3x3(FLOAT pMatrix[3][3]);
    void SetIdentity();
    void SetZero();

    //  BuildActiveSmalLRotation
    void BuildActiveSmallRotation(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);
    
    //   BuildActiveEulerMatrix - oEulerAngles are in radians 
    void BuildActiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);

    void BuildPassiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);

    // Builds matrix that rotates about the specified axis
    // Precondition:  oAxis must be a UNIT VECTOR
    void BuildRotationAboutAxis( const SVector3 &oAxis, FLOAT fAngle );

    
    //  Given a unit vector oProjDirection, build a projection matrix.
    void BuildProjectionMatrix(const SVector3 &oProjDirection );

      
    // return a vector of Euler angles in ZYZ convention in radians
    SVector3 GetEulerAngles() const;
    
    FLOAT Trace() const;
    
    void Transpose();

    SMatrix3x3 operator*( Float f ) const;
    SMatrix3x3 & operator*=( Float f );
    SMatrix3x3 operator/( Float f ) const;
    SMatrix3x3 & operator/=( Float f );


    
    SMatrix3x3 operator*( const SMatrix3x3& oRHS ) const;
    SMatrix3x3 operator+( const SMatrix3x3& oRHS ) const;
    SMatrix3x3 & operator+=( const SMatrix3x3& oRHS );
    SMatrix3x3 operator-( const SMatrix3x3& oRHS ) const;
    SVector3 operator*( const SVector3 &oRHS ) const;

    FLOAT m[3][3];
  };

  
  std::ostream & operator<< ( std::ostream & os, const SMatrix3x3 &m ); 
  
  ////////////////////////////////////////////////////////////////////////////////////////////////

  class SMatrix4x4
  {
  public:
    SMatrix4x4();

    SMatrix4x4( FLOAT f11, FLOAT f12, FLOAT f13, FLOAT f14,
                FLOAT f21, FLOAT f22, FLOAT f23, FLOAT f24,
                FLOAT f31, FLOAT f32, FLOAT f33, FLOAT f34,
                FLOAT f41, FLOAT f42, FLOAT f43, FLOAT f44 );

    // Set matrix to identity
    void SetIdentity();

    // Transpose the matrix
    void Transpose();

    // Inverse the matrix (SLOW)
    void Inverse();

    // Builds matrix that rotates about the specified axis
    // Precondition:  oAxis must be a UNIT VECTOR
    void BuildRotationAboutAxis( const SVector3 &oAxis, FLOAT fAngle );

    // Builds matrix that rotates about the specific roll pitch yaw angles (in radians)
    void BuildRotationRollPitchYaw( FLOAT fYaw, FLOAT fPitch, FLOAT fRoll );

    // Builds matrix that rotates about the specific Euler Angles, following
    // classical physics convention, aka the z-x-z convention
    void BuildPassiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);

    // Set active euler matrix (only rotational part)
    void SetPassiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi );

    // set passive euler matrix (only rotational part)
    void SetActiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi );

    // Set the translational part of the matrix
    void SetTranslation( FLOAT fX, FLOAT fY, FLOAT fZ );

    //  SetRotation
    void SetRotation( const SMatrix3x3 & oRotComponent);
    
    // Builds active transformat matrix using Euler rotation (inverse of
    // the coordinate transformation version)
    void BuildActiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);

    // Builds a matrix using the specified offsets.
    void BuildTranslation( FLOAT fX, FLOAT fY, FLOAT fZ );

    // Builds a scaling matrix
    void BuildScale( FLOAT fX, FLOAT fY, FLOAT fZ );

    // Add translation to the matrix
    void Translate( FLOAT fX, FLOAT fY, FLOAT fZ );

    //
    //  Perform right multiplication with oRot to rotate the rotational
    //  component of this matrix.
    //  i.e., this.m * oRot
    void Rotate(const SMatrix3x3 & oRot);
    
    // Add scale to the matrix
    void Scale( FLOAT fX, FLOAT fY, FLOAT fZ );

    SMatrix4x4 operator*( const SMatrix4x4& oRHS ) const;

    FLOAT m[4][4];
  };


  ////////////////////////////////////////////////////////////////////////////////////////////////

  class SVector2
  {
  public:
    SVector2();
    SVector2( FLOAT fX, FLOAT fY );
	void Set( FLOAT fX, FLOAT fY );

    FLOAT m_fX;
    FLOAT m_fY;
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////

  class SVector3
  {
  public:
    SVector3();
    SVector3( FLOAT fX, FLOAT fY, FLOAT fZ );

    void Set( FLOAT fX, FLOAT fY, FLOAT fZ );

    FLOAT Get( UInt nAxis ) const;
    void Set( UInt nAxis, FLOAT fVal );
    
    FLOAT GetLength() const;

    // Transforms the vector by the given matrix
    // The vector is treated as a normal in the form (x,y,z,0).  
    // Thus, the translation part of the matrix is ignored.
    void TransformNormal( const SMatrix4x4 &oMatrix );
	
	// Transforms the vector by the given matrix
	// The vector is treated as if in the form (x,y,z,1).  
	void Transform( const SMatrix4x4 &oMatrix );
    
	// Transforms the vector by the given matrix (normal 3x3 vector-matrix multiply)
    void Transform( const SMatrix3x3 &oMatrix );
      
    // Normalize the vector
    void Normalize();

    // Assignment & comparisons
    SVector3 &operator=( const SVector3 &oRHS );
        
    // Add & subtract
    SVector3 operator+( const SVector3 &oRHS ) const;
    SVector3 &operator+=( const SVector3 &oRHS );
    SVector3 operator-( const SVector3& oRHS ) const;
    SVector3 operator-( ) const;
    
    SVector3 &operator-=( const SVector3 &oRHS );
    
    // Multiply & divide
    SVector3 operator*( FLOAT fRHS ) const;
    SVector3 &operator*=( FLOAT fRHS );
    SVector3 operator/( FLOAT fRHS ) const;
    SVector3 &operator/=( FLOAT fRHS );		

    // Accessor

    FLOAT operator[]( Int i) const;
    Float& operator[](Int i);

    
    FLOAT m_fX;
    FLOAT m_fY;
    FLOAT m_fZ;
   
    
  };
  
  SVector3 Cross( const SVector3 &oLHS, const SVector3 &oRHS );
  FLOAT Dot( const SVector3 &oLHS, const SVector3 &oRHS );
  SMatrix3x3 OuterProduct( const SVector3 & oLHS, const SVector3 & oRHS );
  SVector3 operator*( FLOAT fLHS, const SVector3 &oRHS );
  std::ostream & operator<< ( std::ostream & os, const SVector3 &v ); 

  SVector3 DegreeToRadian( const SVector3 & oRHS );
  SVector3 RadianToDegree( const SVector3 & oRHS );


  
  ////////////////////////////////////////////////////////////////////////////////////////////////

  class SVector4
  {
  public:
    SVector4();
    SVector4( FLOAT fX, FLOAT fY, FLOAT fZ, FLOAT fW );

    void Set( FLOAT fX, FLOAT fY, FLOAT fZ, FLOAT fW );

    // Transforms the vector by the given matrix
    void Transform( const SMatrix4x4 &oMatrix );

    FLOAT m_fX;
    FLOAT m_fY;
    FLOAT m_fZ;
    FLOAT m_fW;
  };


  ////////////////////////////////////////////////////////////////////////////////////////////////

  // An Axis-aligned bounding box
  struct SBoundingBox
  {										
    SVector3 m_oBoxMin;					
    SVector3 m_oBoxMax;

    // Determines whether the given box overlaps with this one
    // Not a trivially-cheap calculation
    bool Overlaps( const SBoundingBox &oBox ) const;

    //
    //  Union this box with the point p
    void Union( const SVector3 & p )
    {
      for( Int i = 0; i < 3; i ++ )
      {
        m_oBoxMin[i] = std::min( m_oBoxMin[i], p[i] );
        m_oBoxMax[i] = std::max( m_oBoxMax[i], p[i] );
      }
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////

  struct STriangle
  {
    SVector3 m_oPt1;
    SVector3 m_oPt2;
    SVector3 m_oPt3;

    void CalculateNormal( SVector3 &oNormal ) const;
    void CalculateNormalizedNormal( SVector3 &oNormal ) const;
    void CalculateBoundingBox( SBoundingBox &oBox ) const;
  };



  ////////////////////////////////////////////////////////////////////////////////////////////////

  class CPlane
  {
  public:

    // Normalizes a plane so that |a,b,c|==1
    void Normalize();

    FLOAT m_fA;
    FLOAT m_fB;
    FLOAT m_fC;
    FLOAT m_fD;
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////


  struct CSphere
  {
  public:
		
    CSphere() {}
    CSphere( const SVector3 &oCenter, FLOAT fRadius ) { m_oCenter = oCenter; m_fRadius = fRadius; }

    // Accessors
    const SVector3 &GetCenter() const { return m_oCenter; }
    FLOAT GetRadius() const { return m_fRadius; }
    void SetCenter( const SVector3 &oCenter ) { m_oCenter = oCenter; }
    void SetRadius( FLOAT fRadius ) { m_fRadius = fRadius; }
	
	// Calculate an AABB for the sphere
	void CalculateBoundingBox( SBoundingBox &oBox ) const;

  private:

    SVector3 m_oCenter;
    FLOAT m_fRadius;
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////

  class CRay
  {
  public:

    CRay() {}
    CRay( const SVector3 &oStart, const SVector3 &oDir ) : m_oStart(oStart), m_oDir(oDir) {}

    // Accessors
    const SVector3 &GetStart() const { return m_oStart; }
    const SVector3 &GetDirection() const { return m_oDir; }
    void SetStart( const SVector3 &oStart ) { m_oStart = oStart; }
    void SetDirection( const SVector3 &oDir ) { m_oDir = oDir; }
		
		
    //
    //  oStart is a point, and w cannot be zero for a point in Affine Geometry.
    //
    void SetStart( const SVector4 &oStart ) { 
      DEBUG_ASSERT(oStart.m_fW > 0, "SetStart:Weights of a point must be nonzero.");
      m_oStart.m_fX = oStart.m_fX / oStart.m_fW; 
      m_oStart.m_fY = oStart.m_fY / oStart.m_fW; 
      m_oStart.m_fZ = oStart.m_fZ / oStart.m_fW; 

    }

    // Directions do not care about weights, as it is a vector, and
    // weight of a vector is 0
    void SetDirection( const SVector4 &oDir ) { 
      m_oDir.m_fX = oDir.m_fX; 
      m_oDir.m_fY = oDir.m_fY; 
      m_oDir.m_fZ = oDir.m_fZ; 
    }
		
		
  
  
      SVector3 Evaluate( FLOAT fT) const{
      SVector3 oRes;
      oRes = m_oDir * fT + m_oStart;
      return oRes;
    }


  private:

    SVector3 m_oStart;
    SVector3 m_oDir;
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Quadratic formula solver.  Return true if real solution found.
  bool QuadraticFormula( FLOAT fA, FLOAT fB, FLOAT fC, FLOAT &fRes1, FLOAT &fRes2 );

  // Compute barycentric coordinates (u,v,w) for point p with respect to the given triangle
  void Barycentric( const STriangle &oTriangle, const SVector3 &oPoint, FLOAT &fU, FLOAT &fV, FLOAT &fW );

  // Interpolate the given parameters given the barycentric coordinates.  The w-parmater is not necessary since u+v+w=1
  void Interpolate( SVector3 &oRes, const SVector3 &oParam1, const SVector3 &oParam2, const SVector3 &oParam3, FLOAT fU, FLOAT fV );

  // Linearly interpolate 2 vectors
  SVector3 Interpolate( const SVector3 &oParam1, const SVector3 &oParam2, FLOAT fInterp );

  // Spherical linear interpolation between two quaternions
  //  void Interpolate( SQuaternion &oRes, const SQuaternion &oParam1, const SQuaternion &oParam2, FLOAT fInterp );



  namespace Collision
  {
    //-------------------------------------------------------------------------------------------
    // Sphere-to-sphere dynamic intersection
    //-------------------------------------------------------------------------------------------
    // One moving sphere, one static sphere
    bool Intersects( const CSphere &oSphereA, const SVector3 &oVelocityA,		// Sphere A properties
                     const CSphere &oSphereB,									// Sphere B properties
                     FLOAT &fCollisionTime1, FLOAT &fCollisionTime2					// Normalized times for 1st & 2nd collision
                     );

    // Two moving spheres
    bool Intersects( const CSphere &oSphereA, const SVector3 &oVelocityA,		// Sphere A properties
                     const CSphere &oSphereB, const SVector3 &oVelocityB,		// Sphere B properties
                     FLOAT &fCollisionTime1, FLOAT &fCollisionTime2					// Normalized times for 1st & 2nd collision
                     );

    //-------------------------------------------------------------------------------------------
    // Ray  intersection
    //-------------------------------------------------------------------------------------------
    // Ray-sphere intersection
    bool Intersects( const CSphere &oSphere, const CRay &oRay, FLOAT &fT );

	// Ray-sphere intersection.  Returns both intersection t (in no particular order).
	bool Intersects( const CSphere &oSphere, const CRay &oRay, FLOAT &fT1, FLOAT &fT2 );


    // Ray-triangle intersection. Assume clockwise winding order for triangles.
    bool Intersects( const SVector3 &oPt1, const SVector3 &oPt2, const SVector3 &oPt3, const CRay &oRay, FLOAT &fT );

	// Ray-triangle intersection. Assume clockwise winding order for triangles.
	bool Intersects( const STriangle &oTriangle, const CRay &oRay, FLOAT &fT );

    // Ray-plane intersection
    bool Intersects( const CPlane &oPlane, const CRay &oRay, FLOAT &fT );

  }
  
  //-------------------------------------------------------------------------------------------
  //
  //  Class SRange
  //
  //  TODO:  Make this a template?
  //
  //-------------------------------------------------------------------------------------------
  class SIntRange
  {
  public:
    Int nHigh;
    Int nLow;
    
    bool Contains(Int n) const
    {
      return ( (n <= nHigh) && (n >= nLow) );
    };  
  };

  //-------------------------------------------------------------------------------------------
  // SRange
  //-------------------------------------------------------------------------------------------
  class SRange
  {
  public:
    FLOAT fLow;
    FLOAT fHigh;
    SRange(){};
    SRange( FLOAT fL, FLOAT fH):fLow(fL), fHigh(fH) {};
    
    bool Contains(FLOAT f) const
    {
      return ( (f <= fHigh) && (f >= fLow) );
    };
  };
  
}

#endif
