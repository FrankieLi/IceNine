//--------------------------------------------
//
//  ShapeGenerator for code testing
//
//--------------------------------------------


#ifndef SHAPE_GENERATOR_H
#define SHAPE_GENERATOR_H

#include "3dMath.h"

namespace GeometricValidation
{
  
  //-----------------------------------------------------------------------------
  //  GenerateSphere
  //-----------------------------------------------------------------------------
  template< class ImageData >
  void GenerateSphere( ImageData & oData,
                       int nX, int nY, int nZ,
                       const GeneralLib::SVector3 & PixelScale,
                       const GeneralLib::SVector3 & Center,
                       float Radius, int ObjectID )
  {
    typedef GeneralLib::SVector3 SVector3;
    oData.SetPixelSize( PixelScale.m_fX, PixelScale.m_fY, PixelScale.m_fZ );
    oData.Resize( nX, nY, nZ );
    for( int i = 0; i < nX; i ++ )
      for( int j = 0; j < nY; j ++ )
        for( int k = 0; k < nZ; k ++ )
          oData.Set( i, j, k, 1 );
  
    for( int k = 1; k < nZ-1; k ++ )
    {
      for( int j = 1; j < nY-1; j ++ )
      {
        for( int i = 1; i < nX-1; i ++ )     
        {
          SVector3 v = SVector3( i, j, k ) - Center;

          if( v.GetLength() <= Radius )
          {
            if( v.m_fX < 0 )
              oData.Set( i, j , k, ObjectID );
            else
              oData.Set( i, j , k, ObjectID -1 );
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------------
  //  GenerateCylinder
  //
  //  Always perpendicular to x-y plane
  //
  //-----------------------------------------------------------------------------
  template< class ImageData >
  void GenerateCylinder( ImageData & oData,
                         int nX, int nY, int nZ,
                         const GeneralLib::SVector3 & PixelScale,
                         const GeneralLib::SVector3 & Center,
                         float Radius, float Length, int ObjectID )
  {
    typedef GeneralLib::SVector3 SVector3;
    oData.SetPixelSize( PixelScale.m_fX, PixelScale.m_fY, PixelScale.m_fZ );
    oData.Resize( nX, nY, nZ );
    for( int i = 0; i < nX; i ++ )
      for( int j = 0; j < nY; j ++ )
        for( int k = 0; k < nZ; k ++ )
          oData.Set( i, j, k, -1 );
  
    for( int k = 1; k < nZ-1; k ++ )
    {
      for( int j = 1; j < nY-1; j ++ )
      {
        for( int i = 1; i < nX-1; i ++ )     
        {
          SVector3 v = SVector3( i, j, 0 ) - Center;
          if( std::abs( v.m_fZ ) < Length )
          {
            v.m_fZ = 0;
            if( v.GetLength() <= Radius )
              oData.Set( i, j , k, ObjectID );
          }
        }
      }
    }
  }
  
  //-----------------------------------------------------------------------------
  //  GenerateCube
  //
  //  Always perpendicular to x-y plane
  //
  //-----------------------------------------------------------------------------
  template< class ImageData >
  void GenerateCube( ImageData & oData,
                     int nX, int nY, int nZ,
                     const GeneralLib::SVector3 & PixelScale,
                     const GeneralLib::SVector3 & Center,
                     float fSideLength,
                     const GeneralLib::SMatrix3x3 & oRot,
                     int ObjectID )
  {
    typedef GeneralLib::SVector3 SVector3;
    oData.SetPixelSize( PixelScale.m_fX, PixelScale.m_fY, PixelScale.m_fZ );
    oData.Resize( nX, nY, nZ );
    for( int i = 0; i < nX; i ++ )
      for( int j = 0; j < nY; j ++ )
        for( int k = 0; k < nZ; k ++ )
          oData.Set( i, j, k, -1 );
  
    for( int k = 1; k < nZ-1; k ++ )
    {
      for( int j = 1; j < nY-1; j ++ )
      {
        for( int i = 1; i < nX-1; i ++ )     
        {
          SVector3 v(i, j, k);
          v -= Center;
          v = oRot * v;
          float dX = floor( fSideLength / 2. ) + 0.5;
          if( v.m_fX    <= dX && v.m_fX > -dX
              && v.m_fY <= dX && v.m_fY > -dX
              && v.m_fZ <= dX && v.m_fZ > -dX )
          {
            oData.Set( i, j , k, ObjectID );
          }
        }
      }
    }
  }


  //-----------------------------------------------------------------------------
  //  GenerateConstrainedCube
  //
  //  Generate a constrained cube center on in the region.  Constrain simply
  //  means that it has neighboring "grains" that are not empty space, and therefore
  //  we can test the effects of the triple line + quad point identifier.
  //
  //  Assumption:  nX == nY == nZ;
  //
  //-----------------------------------------------------------------------------
  template< class ImageData >
  void GenerateConstrainedCube( ImageData & oData,
                                int nX, int nY, int nZ,
                                const GeneralLib::SVector3 & PixelScale,
                                const GeneralLib::SVector3 & Center,
                                float fSideLength,
                                const GeneralLib::SMatrix3x3 & oRot,
                                int ObjectID )
  {
    typedef GeneralLib::SVector3 SVector3;
    oData.SetPixelSize( PixelScale.m_fX, PixelScale.m_fY, PixelScale.m_fZ );
    oData.Resize( nX, nY, nZ );
    for( int i = 0; i < nX; i ++ )
      for( int j = 0; j < nY; j ++ )
        for( int k = 0; k < nZ; k ++ )
          oData.Set( i, j, k, -1 );


    // 27 grains total, only the one at the center will recieve ObjectID.
    vector<float> fMin(3);
    vector<float> fMax(3);

    float SideGrainWidth = floor( ( nX - fSideLength) / 2. );
    fMin[0] = 0;
    fMin[1] = SideGrainWidth;
    fMin[2] = SideGrainWidth + fSideLength;

    fMax[0] = fMin[1];
    fMax[1] = fMin[2];
    fMax[2] = nX;
    
    for( int k = 1; k < nZ-1; k ++ )
    {
      for( int j = 1; j < nY-1; j ++ )
      {
        for( int i = 1; i < nX-1; i ++ )     
        {
          SVector3 v(i, j, k);
          v -= Center;
          v = oRot * v + Center;
          
          bool bOutside = false;
          vector<int> CubeCode(3);
          for( int n = 0; n < 3; n++)
          {
            bool bCodeFound = false;
            for( int nCode = 0; nCode < 3 && ! bCodeFound; nCode ++)
            {
              // std::cout << fMin[nCode] << " " << v[n] << " " << fMax[nCode] << std::endl;
              if( ( v[n] > fMin[nCode] ) && ( v[n] <= fMax[nCode] ) )
              {
                bCodeFound = true;
                CubeCode[n] = nCode;
              }
            }
            if( ! bCodeFound )
              bOutside = true;
          }
          if( !bOutside )  // note: [111] = center =  1 * 3^2 + 1*3^1 + 1*3^0 =  13
          {
            int ID = CubeCode[2] * 9 + CubeCode[1] * 3 + CubeCode[0];
            if( ID != 13 )
              oData.Set( i, j , k, ID + ObjectID + 1 );
            else
              oData.Set( i, j , k, ObjectID );
          }

        }
      }
    }
  }
  
}



#endif
