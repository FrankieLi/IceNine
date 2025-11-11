//============================================================================== 
// Copyright (c) 2014, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by S. F. Li (li31@llnl.gov)
// LLNL-CODE-657639
// All rights reserved.
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the Lawrence Livermore National Lab nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL LAB BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//============================================================================== 

//------------------------------------------------------------------------------------
//  Author:  S. F. Li (Frankie)
//  e-mail:  li31@llnl.gov; sfli@cmu.edu 
//------------------------------------------------------------------------------------
////////////////////////////////////////////////////////////
//
//  OrientationSearch.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//  Purpose:  Experimental orientation search methods.
//            This is the place for the supporting routines
//
//  TODO:     Abstract orientation space search and place into this namespace
//
////////////////////////////////////////////////////////////


#ifndef _ORIENTATION_SEARCH_H_
#define _ORIENTATION_SEARCH_H_

#include "3dMath.h"
#include "Types.h"
#include "Debug.h"
#include "math.h"
#include <vector>
#include "Sampling.h"
#include <fstream>

//  BOOST Mersenne Twister and distribution generators



#include <boost/numeric/ublas/vector.hpp>

#include "Quaternion.h"

#include <tuple>

using namespace GeneralLib;

//--------------------------------------------------------------------------------------------------------
//
//  namespace  OrientationSearch
//
//--------------------------------------------------------------------------------------------------------
namespace OrientationSearch
{

  //--------------------------------------------------------------------------------------------------------
  //
  //  GetAlignRotation
  //  -- Given *unit vector* oRefDir, oVecDir, return a matrix that rotates oVec -> oRef
  //
  //--------------------------------------------------------------------------------------------------------
  SMatrix3x3 GetAlignmentRotation( const SVector3 &oVecDir, const SVector3 &oRefDir );
  
}

//--------------------------------------------------------------------------------------------------------
//
//  namespace  OrientationOptimization
//
//  Purpose:   A collection of optimization routines to find the orientation of the sample that'd
//             give a maximally overlapping simulated peaks.  Currently, the planned methods of
//             optimizations are:
//             
//             Zero temperature Monte Carlo
//             Simulated Annealing
//             Conjugate Gradient
//      
//             Note that parallelized versions may also be implemented in the future.
//
//--------------------------------------------------------------------------------------------------------
namespace OrientationOptimization
{

  //---------------------------------------------------------------------
  //
  //  A test class
  //
  //---------------------------------------------------------------------
  class COverlapOptimizerTrait
  {
  public:
    typedef Float ValueType;
    typedef SMatrix3x3 CoordinateType; 
  };
  
  //---------------------------------------------------------------------
  //
  //  ObjFnTraitT
  //
  //---------------------------------------------------------------------
  
  //---------------------------------------------------------------------
  //
  // CObjectiveFunction
  //
  // Purpose:   General definition of objective functions to be a function
  //            object that returns a value given a position.  Note that
  //            a class must be constructed and derived from CObjectiveFunction.
  //            The operator, (), must be implemented.  ObjFnTraitT has the
  //            requirement that it must contain a ValueType and CoordinateType.
  //
  //            Less than operator must be well defined for the value type.
  //
  //---------------------------------------------------------------------
  template< class ObjFnTraitT >
  class CObjectiveFunction
  {
  public:
    typedef typename ObjFnTraitT::ValueType ValueType;
    typedef typename ObjFnTraitT::CoordinateType CoordinateType;
    typedef typename std::less< ValueType > LessThanOp;
    virtual ValueType operator() ( const CoordinateType & oPos ) = 0;
  };

  
  //---------------------------------------------------------------------
  //  A test class for the monte carlo searh functions
  //
  //  Purpose:  This is to test and make sure that COrientationMC works
  //
  //  TODO:   Move this to UniformReconstruction.h
  //
  //
  //---------------------------------------------------------------------
  class COverlapFunction : public CObjectiveFunction< COverlapOptimizerTrait >
  {
  public:
    Float operator() ( const SMatrix3x3 & oPos )
    {
      RUNTIME_ASSERT( 0, "COverlapFunction is being called directly without derivation\n");
      return -1;
    };
  };
  
  //--------------------------------------------------------------------------------------------------------
  //
  //  class  COrientationMC
  //
  //  Purpose:  A class that holds all of the Monte Carlo optimization related routines
  //            along with their helper functions.
  //
  //--------------------------------------------------------------------------------------------------------
  class COrientationMC
  {
  public:

    //    std::ofstream oLogFile;

  private:

    CSingletonUniformRandomReal & oRandomReal;  // singleton object
    UniformGrid::CQuaternionGrid oUniformGridGen;

    //---------------------------------------------------------------------
    //  ScaleRandomVariable
    //---------------------------------------------------------------------
    Float GetRandomVariable( Float fMin, Float fMax );
    
  public:
    
    //---------------------------------------------------------------------
    //  Default Constructor
    //---------------------------------------------------------------------
    COrientationMC(): oRandomReal( CSingletonUniformRandomReal::Get() ),
                      oUniformGridGen()  {}

    //---------------------------------------------------------------------
    //  Public:   ZeroTemperatureOptimization
    //
    //  Purpose:  Using zero temperature Monte Carlo method to optimize the
    //            orientation, which is subjected tot he objective function,
    //            oObjectiveFunction.  (Rationale - this allows for a good amount
    //            of generalization with very minimal overhead.  In fact, an inline
    //            implementation of the oObjectiveFunction will basically be exactly
    //            the same as pasting the function into ZeroTemperatureOptimization.)
    //
    //
    //  Comments: The random number generator used is the Mersenne Twister.  This function
    //            is statistically sound when used in series.  However, a paralleled version
    //            *should require* randomized seeding to prevent correlation between successive
    //            elements.  Note that to a large extent, this is a limitation of the Monte
    //            Carlo method implemented.  An alternative would be to implementa a quicksort-like
    //            algorithm that randomly chooses a pivot.  Consequently, the solution will be
    //            deterministic even though the algorithm itself is random.  (To be implemented soon)
    //
    //
    //  Algorithm:  A state is randomly chosen within the oInitialOrientation + oVariationalRange.
    //              The state will be accepted if the objective function evaluated at the new orientation
    //              is lower than the previous state. 
    //
    //
    //  Parameters:  oInitialOrientation -  Initial orientation in euler angles
    //
    //               oVariationRange     -  The maximum and minium displacement of the orientation in
    //                                      Euler angle to be searched.  The angles are specified in radians.
    //                                      Clearly, the angles must be small.  Otherwise the non-Euclidiean
    //                                      behavior of the Euler space will become prominent.
    //
    //               oObjectiveFunction  -  A function that returns a cost for a given orientation using
    //                                      the overloaded operator().  This is the function that is being
    //                                      minimized by this optimization.
    //
    //
    //               nMaxMCStep          -  Maximum number of Monte Carlo Steps before convergence
    //
    //               fAngularRadius      -  Area that is covered in one evaluation of the cost function (voxel size) 
    //
    //              
    //  TODO:  -- Add convergence criterion as a function object?
    //
    //  Return:    An object, SVector3, containing the optimized orientation in Euler angles.  Return value
    //             optimization (RVO) is expected from the compiler.  (g++ 3.2 or above and icpc)
    //---------------------------------------------------------------------
    std::pair<SMatrix3x3, Float>
    ZeroTemperatureOptimization( const SMatrix3x3 & oInitialOrientation,
                                 Float fAngularRadius, 
                                 //const CObjectiveFunction< ObjFnTraitT > & oObjectiveFunction ); // To be used later
                                 COverlapFunction & oObjectiveFunction,
                                 Int nMaxMCStep );
    
    std::tuple<SMatrix3x3, Float, Float>
    ZeroTemperatureOptimizationWithVariance( const SMatrix3x3 & oInitialOrientation,
                                             Float fAngularRadius, 
                                             COverlapFunction & oObjectiveFunction,
                                             Int nMaxMCStep );

    //---------------------------------------------------------------------
    //
    //  RandomRestartZeropTemp
    //
    //  Purpose:  A somewhat "discretized version" of simulated annealing.
    //            Given an initial orientation, used as the starting point,
    //            a zero temperature Monte Carlo optimization is performed
    //            to optimize the cost function locally.  After a certain
    //            number of steps, n_ergodic, which is calculated based on
    //            the search space volume specified by fAngularBoxSideLength
    //            and fAngularStepSize, the optimization will restart with a
    //            new random point within the bounding box delinated by
    //            fAngularBoxSideLength will be selected.  A second optimization
    //            will be ran.  This continues until nMaxRestarts or convergence.
    //
    //  Parameters:  oInitialOrientation -- initial orientation
    //
    //               fAngularStepSize -- step size used in the Monte Carlo optimization
    //                                   measured in radians.
    //               fAngularBoxSideLength -- The side length of the bounding volume
    //                                        of the local search space, measured in radians.
    //               oObjectiveFunction -- Cost function to be minimized.
    //               nMaxMCStep         --  Maximum number of Monte Carlo steps allowed.
    //               nMaxFailedRestarts -- Maximum number of consecutively failed restarts
    //                                     allowed before quitting
    //               fMaxConvergenceCost -- Maximum cost allowed to be considered as converged.
    //                                      Note that this is an *early quit* criterion
    //
    //---------------------------------------------------------------------
    std::pair<SMatrix3x3, Float>
    RandomRestartZeroTemp( const SMatrix3x3 & oInitialOrientation,
                           Float fAngularStepSize, Float fAngularBoxSideLength, 
                           COverlapFunction & oObjectiveFunction,
                           Int nMaxMCStep,
                           Int nMaxFailedRestarts,
                           Float fMaxConvergenceCost );



    //---------------------------------------------------------------------
    //  AdaptiveSamplingZeroTemp
    //
    //   Purpose:  Local orientation search is performed by adaptively sampling
    //             the orientation space.
    //
    //
    //   Parameters:
    //         fCostFnAngularResolution - The size of the cost function peaks in
    //         the orientation space.  As a rule of thumb, it's the radius of
    //         the cost function that's above the "background."
    //---------------------------------------------------------------------
    std::pair<SMatrix3x3, Float>
    AdaptiveSamplingZeroTemp( const SMatrix3x3 & oInitialOrientation,
                              Float fCostFnAngularResolution,
                              Float fSearchRegionAngularSideLength,
                              COverlapFunction & oObjectiveFunction,
                              Int nMaxMCStep,
                              Int NumStratifiedSamples,
                              Float fConvergenceVariance, 
                              Float fMaxConvergenceCost );

  };
}



#endif
