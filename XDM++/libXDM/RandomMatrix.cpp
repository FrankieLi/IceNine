//----------------------------------------------------------
//  RandomMatrix.cpp
//
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Implementation of Random Matrix generator
//----------------------------------------------------------

#include "RandomMatrix.h"
namespace GeneralLib
{
  namespace UniformGrid
  {
    //-------------------------------------
    //  DiagonalRandomMatrixGenerator::operator()
    //  Purpose:  Generator of random diagonal matrices
    //-------------------------------------
    SMatrix3x3 DiagonalRandomMatrixGenerator::operator() ( ) const
    {
      SMatrix3x3 oRes;
      oRes.SetZero();
      oRes.m[0][0] =
        RandomRealGen.GetRandomVariable( BndBox.m_oBoxMin.m_fX,
                                         BndBox.m_oBoxMax.m_fX );
      oRes.m[1][1] =
        RandomRealGen.GetRandomVariable( BndBox.m_oBoxMin.m_fY,
                                         BndBox.m_oBoxMax.m_fY );
      oRes.m[2][2] =
        RandomRealGen.GetRandomVariable( BndBox.m_oBoxMin.m_fZ,
                                         BndBox.m_oBoxMax.m_fZ );
      return oRes;
    }
    
    //-------------------------------------
    //  DiagonalRandomMatrixGenerator::SetVolume()
    //  Purpose:  Setvolume to the new bounding box
    //-------------------------------------
    void DiagonalRandomMatrixGenerator::SetVolume( const SBoundingBox & NewBBox )
    {
      BndBox = NewBBox;
    }
    
    //-------------------------------------
    //  SymmetricRandomMatrixGenerator::operator()
    //  Purpose:  Generator of random diagonal matrices
    //-------------------------------------
    SMatrix3x3 SymmetricRandomMatrixGenerator::operator() ( ) const
    {
      SMatrix3x3 oRes;

      oRes.m[0][1] =  RandomRealGen.GetRandomVariable( _min.m[0][1], _max.m[0][1] );
      oRes.m[0][2] =  RandomRealGen.GetRandomVariable( _min.m[0][2], _max.m[0][2] );
      oRes.m[1][2] =  RandomRealGen.GetRandomVariable( _min.m[1][2], _max.m[1][2] );

      oRes.m[0][0] =  RandomRealGen.GetRandomVariable( _min.m[0][0], _max.m[0][0] );
      oRes.m[1][1] =  RandomRealGen.GetRandomVariable( _min.m[1][1], _max.m[1][1] );
      oRes.m[2][2] =  RandomRealGen.GetRandomVariable( _min.m[2][2], _max.m[2][2] );

      oRes.m[1][0] =  oRes.m[0][1];
      oRes.m[2][0] =  oRes.m[0][2];
      oRes.m[2][1] =  oRes.m[1][2];
      
      return oRes;
    }
  }
}
