//----------------------------------------------------------
//
//  RandomMatrix.h
//
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Implementation of Random Matrix generator
//
//
//----------------------------------------------------------

#ifndef RANDOM_MATRIX_H
#define RANDOM_MATRIX_H

#include "3dMath.h"
#include "Sampling.h"
namespace GeneralLib
{
  namespace UniformGrid
  {
    class CRandomMatrixGenerator
    {
    protected:
      mutable CUniformRandomReal RandomRealGen;
    public:
      virtual SMatrix3x3 operator()() const = 0;      
    };

    //-------------------------------------
    //  DiagonalRandomMatrixGenerator
    //  - generate random matrices that
    //    are diagonal
    //-------------------------------------
    class DiagonalRandomMatrixGenerator
      : public CRandomMatrixGenerator
    {
    private:
      SBoundingBox BndBox;

      DiagonalRandomMatrixGenerator();
    public:

      DiagonalRandomMatrixGenerator( const SBoundingBox & oRangeBBox ):
        BndBox( oRangeBBox ) { }

      void SetVolume( const SBoundingBox & NewBBox );
      
      SMatrix3x3 operator()( ) const;
    };

    //-------------------------------------
    //  SymmetricRandomMatrixGenerator
    //  - generate random matrices that
    //    are symmetric
    //-------------------------------------
    class SymmetricRandomMatrixGenerator
      : public CRandomMatrixGenerator
    {
    private:
      SymmetricRandomMatrixGenerator();
      SMatrix3x3 _max;
      SMatrix3x3 _min;
    public:
      SymmetricRandomMatrixGenerator( const SMatrix3x3 & MaxValues,
                                      const SMatrix3x3 & MinValues ):
        _max( MaxValues ), _min( MinValues ) { }
      
      SMatrix3x3 operator()( ) const;
    };
    
  }
}
#endif
