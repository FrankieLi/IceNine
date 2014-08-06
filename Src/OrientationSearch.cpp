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
//  OrientationSearch.cpp
//
//  Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of orientation search and
//            optimization functions
//
////////////////////////////////////////////////////////////


#include "OrientationSearch.h"


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
  //  -- Given *UNIT vector* oRef, oVec, return a matrix that rotates oVec -> oRef
  //
  //--------------------------------------------------------------------------------------------------------
  SMatrix3x3 GetAlignmentRotation( const SVector3 &oVecDir, const SVector3 &oRefDir )
  {
    
    SVector3 oAxis  = Cross( oVecDir, oRefDir );
    Float    fAngle = acos( Dot( oVecDir, oRefDir ) );

    oAxis.Normalize();
    SMatrix3x3 oRes;
    oRes.BuildRotationAboutAxis( oAxis, fAngle );
    return oRes;
  }

}

//--------------------------------------------------------------------------------------------------------
//
//  namespace  OrientationOptimization
//
//--------------------------------------------------------------------------------------------------------
namespace OrientationOptimization
{

  //---------------------------------------------------------------------
  //  GetRandomVariable
  //
  //  Return a random variable that's between [fMin and fMax), with uniform
  //  distribution.
  //---------------------------------------------------------------------
  Float COrientationMC::GetRandomVariable( Float fMin, Float fMax )
  {
    Float fScale = fMax - fMin;
    return fScale * oRandomReal() + fMin;
  }
  
  //--------------------------------------------------------------------------------------------------------
  //  Public:   ZeroTemperatureOptimization
  //--------------------------------------------------------------------------------------------------------
  std::pair<SMatrix3x3, Float>
  COrientationMC::ZeroTemperatureOptimization( const SMatrix3x3 & oInitialOrientation,
                                               Float fAngularRadius, 
                                               COverlapFunction & oObjectiveFunction,
                                               Int nMaxMCStep )
  {
    SMatrix3x3 oOptimalState = oInitialOrientation;
    COverlapFunction::ValueType fMinimizedCost = oObjectiveFunction( oInitialOrientation );
    //-----------------
    //  Converting fRadius -- an angular
    //  distance to a radius in R^3 of the
    //  interpolation parameter space
    //-----------------
    Float fRadius = tan( fAngularRadius ) / sqrt( 12.0 ); // sqrt(12) = 2 * sqrt(3)
    
    for ( Int i = 0; i < nMaxMCStep; i ++ )
    {
      Float fX =  GetRandomVariable( -fRadius, fRadius );
      Float fY =  GetRandomVariable( -fRadius, fRadius );
      Float fZ =  GetRandomVariable( -fRadius, fRadius );
      SQuaternion q = oUniformGridGen.GetNearIdentityPoint( fX, fY, fZ );
      SMatrix3x3 oDelta = q.GetRotationMatrix3x3();
      SMatrix3x3 oCurrentState = oDelta * oOptimalState;
      
      COverlapFunction::ValueType fCurrentCost = oObjectiveFunction( oCurrentState );
      if( fCurrentCost <  fMinimizedCost )
      {
        oOptimalState = oCurrentState;
        fMinimizedCost = fCurrentCost;
      }
    }
    
    // Named return value optimization
    std::pair<SMatrix3x3, Float> oRes = std::make_pair( oOptimalState, fMinimizedCost ); 
    return oRes;
  }

  //--------------------------------------------------------------------------------------------------------
  //  Public:   ZeroTemperatureOptimization
  //--------------------------------------------------------------------------------------------------------
  boost::tuple<SMatrix3x3, Float, Float>
  COrientationMC::ZeroTemperatureOptimizationWithVariance( const SMatrix3x3 & oInitialOrientation,
                                                           Float fAngularRadius, 
                                                           COverlapFunction & oObjectiveFunction,
                                                           Int nMaxMCStep )
  {
    
    SMatrix3x3 oOptimalState = oInitialOrientation;
    COverlapFunction::ValueType fMinimizedCost = oObjectiveFunction( oInitialOrientation );

    Float Average           = fMinimizedCost;   // Running variance calculation due to Knuth, via Wikipedia
    int   NumSampledPoints  = 1;
    Float Variance          = 0;
    
    //-----------------
    //  Converting fRadius -- an angular
    //  distance to a radius in R^3 of the
    //  interpolation parameter space
    //-----------------
    Float fRadius = tan( fAngularRadius ) / sqrt( 12.0 ); // sqrt(12) = 2 * sqrt(3)


    for ( Int i = 0; i < nMaxMCStep; i ++ )
    {
      Float fX =  GetRandomVariable( -fRadius, fRadius );
      Float fY =  GetRandomVariable( -fRadius, fRadius );
      Float fZ =  GetRandomVariable( -fRadius, fRadius );
      SQuaternion q = oUniformGridGen.GetNearIdentityPoint( fX, fY, fZ );
      SMatrix3x3 oDelta = q.GetRotationMatrix3x3();
      SMatrix3x3 oCurrentState = oDelta * oOptimalState;
      
      COverlapFunction::ValueType fCurrentCost = oObjectiveFunction( oCurrentState );
      if( fCurrentCost <  fMinimizedCost )
      {
        oOptimalState = oCurrentState;
        fMinimizedCost = fCurrentCost;
      }

      //----------------
      //  variance calculation
      //  - calculate teh variance of the sampled region
      //----------------
      NumSampledPoints ++;
      Float Delta = fCurrentCost - Average;
      Average = Average + Delta / Float( NumSampledPoints );
      Variance += Delta * ( fCurrentCost - Average );
    }

    if( NumSampledPoints > 1 )
      Variance = Variance / Float( NumSampledPoints - 1);
    else
      Variance = -1.0;
    
    // Named return value optimization
    boost::tuple<SMatrix3x3, Float, Float> oRes = boost::make_tuple( oOptimalState, fMinimizedCost, Variance ); 
    return oRes;
  }

  //--------------------------------------------------------------------------------------------------------
  //  Public:   AdaptiveSamplingZeroTemp
  //
  //   This is essetially depth first search, non-recursive
  //
  //--------------------------------------------------------------------------------------------------------
  std::pair<SMatrix3x3, Float>
  COrientationMC::AdaptiveSamplingZeroTemp( const SMatrix3x3 & oInitialOrientation,
                                            Float fCostFnAngularResolution,
                                            Float fSearchRegionAngularSideLength,
                                            COverlapFunction & oObjectiveFunction,
                                            Int nMaxMCStep,
                                            Int NumMaxStratifiedSamples,
                                            Float fConvergenceVariance,
                                            Float fMaxConvergenceCost )
  {
    SMatrix3x3 oGlobalOptimalState = oInitialOrientation;
    SMatrix3x3 oCurrentState =  oInitialOrientation;
    COverlapFunction::ValueType fGlobalMinCost = oObjectiveFunction( oInitialOrientation );
    
    Float InitialSubregionRadius  = tan( fSearchRegionAngularSideLength ) / sqrt( 48.0 ); // sqrt(48) = 4 * sqrt(3)   ( random start radius )
    
    Int nTotalStepTaken = 0;
    Int NumGlobalPointsTaken = 0;
    Float fCurrentVariance;
    Float SubregionRadius = InitialSubregionRadius;
    
    while( nTotalStepTaken < nMaxMCStep )
    {
      Int NumSubregionSteps = Int( pow( ceil(SubregionRadius / fCostFnAngularResolution), 2.7 ) ); // Number of steps it takes to adequately
                                                                                                 // sample each subregions without missing
                                                                                                 // a cost function (giving it less more steps)

      NumSubregionSteps = std::max( NumSubregionSteps, 10 ); // at least 20 steps
      Float fCurrentCost;
      SMatrix3x3 oTmpResOrient;
      
      boost::tie( oTmpResOrient, fCurrentCost, fCurrentVariance )
        =  ZeroTemperatureOptimizationWithVariance( oCurrentState, SubregionRadius, 
                                                    oObjectiveFunction, NumSubregionSteps );
      if( fCurrentVariance > fConvergenceVariance )
        nMaxMCStep += NumSubregionSteps;
      
      nTotalStepTaken += NumSubregionSteps;
      if( fCurrentCost >=  fGlobalMinCost )           // restart with new position
      {
        //----------
        // Backtrack/Reset search parameters
        //----------
        SubregionRadius = std::min( static_cast<Float>( 2.0 * SubregionRadius), fSearchRegionAngularSideLength );
        NumGlobalPointsTaken ++;
        
        Float fX =  GetRandomVariable( -SubregionRadius, SubregionRadius );
        Float fY =  GetRandomVariable( -SubregionRadius, SubregionRadius );
        Float fZ =  GetRandomVariable( -SubregionRadius, SubregionRadius );
        SQuaternion q = oUniformGridGen.GetNearIdentityPoint( fX, fY, fZ );
        SMatrix3x3 oDelta = q.GetRotationMatrix3x3();
        oCurrentState = oDelta * oInitialOrientation;
      }
      else                                   // continuation of search, with narrowing of steps size
      {
        fGlobalMinCost      = fCurrentCost;
        oGlobalOptimalState = oTmpResOrient;
        oCurrentState       = oTmpResOrient;
        SubregionRadius     *= Float( 0.5 );    // reduce step size
      }
      
      if( (fGlobalMinCost < fMaxConvergenceCost)  && fabs( fCurrentVariance ) < fConvergenceVariance )    // convergence  -- move criterion to user defined
        break;
    }
    
        
    std::cout << RADIAN_TO_DEGREE( LatticeSymmetry::GetMisorientation( LatticeSymmetry::CCubicSymmetry::Get(), oGlobalOptimalState, oInitialOrientation )  )
	      << "\t|R_0| " << RADIAN_TO_DEGREE( InitialSubregionRadius )
	      << "\t|Step | " << nTotalStepTaken 
	      << "\t|GlbPts| " << NumGlobalPointsTaken
	      << "\t|R_c| " << RADIAN_TO_DEGREE( SubregionRadius )
	      << "\t|MCS| " << nMaxMCStep 
	      << "\t|Var| " << fCurrentVariance 
	      << "\t|C_i| " << oObjectiveFunction( oInitialOrientation )
	      << "\t|C_f|" << fGlobalMinCost 
	      << "\t[ " << RadianToDegree( oInitialOrientation.GetEulerAngles() ) << "]"
	      << std::endl; 
    
    // Named return value optimization
    std::pair<SMatrix3x3, Float> oRes = std::make_pair( oGlobalOptimalState, fGlobalMinCost ); 
    return oRes;
  }

  
  //--------------------------------------------------------------------------------------------------------
  //  Public:   RandomRestartZeroTemp
  //
  //  Parameter:  fDriftDistance is the side of the box of the local search space.
  //
  //--------------------------------------------------------------------------------------------------------
  std::pair<SMatrix3x3, Float>
  COrientationMC::RandomRestartZeroTemp( const SMatrix3x3 & oInitialOrientation,
                                         Float fAngularStepSize, Float fAngularBoxSideLength, 
                                         COverlapFunction & oObjectiveFunction,
                                         Int nMaxMCStep,
                                         Int nMaxFailedRestarts,
                                         Float fMaxConvergenceCost )
  {
    SMatrix3x3 oGlobalOptimalState = oInitialOrientation;
    SMatrix3x3 oCurrentState =  oInitialOrientation;
    COverlapFunction::ValueType fGlobalMinCost = oObjectiveFunction( oInitialOrientation );
    
    
    Float fCurAngularStepSize = fAngularStepSize;
    
    Int nTotalStepTaken = 0;
    Int nSuccessiveRestarts = 0;
    Int nMinErgodicSteps = 2 * pow( fAngularBoxSideLength / fCurAngularStepSize, 3 );  // Number of steps required to be ergodic
    Int nOptimizationSteps = nMinErgodicSteps;
    
    while(  nTotalStepTaken < nMaxMCStep )
    {
     
      nOptimizationSteps = std::min( nOptimizationSteps, ( nMaxMCStep - nTotalStepTaken ) );
      nOptimizationSteps = std::max( nOptimizationSteps, 0 );
      
      Float fCurrentCost;
      SMatrix3x3 oTmpResOrient;
      boost::tie( oTmpResOrient, fCurrentCost ) =  ZeroTemperatureOptimization( oCurrentState, fCurAngularStepSize, 
                                                                                oObjectiveFunction, nOptimizationSteps );
      nTotalStepTaken += nOptimizationSteps;
      
      if( fCurrentCost >=  fGlobalMinCost )           // restart with new position
      {
        Float fRadius     = tan( fAngularBoxSideLength ) / sqrt( 48.0 ); // sqrt(48) = 4 * sqrt(3)   ( random start radius )
        Float fX =  GetRandomVariable( -fRadius, fRadius );
        Float fY =  GetRandomVariable( -fRadius, fRadius );
        Float fZ =  GetRandomVariable( -fRadius, fRadius );
        SQuaternion q = oUniformGridGen.GetNearIdentityPoint( fX, fY, fZ );
        SMatrix3x3 oDelta = q.GetRotationMatrix3x3();
        oCurrentState = oDelta * oInitialOrientation;

        //----------
        // Reset search parameters
        //----------
        fCurAngularStepSize = fAngularStepSize;  
        nSuccessiveRestarts ++;
        nOptimizationSteps = nMinErgodicSteps;
      }
      else                                   // continuation of search, with narrowing of steps size
      {
        nSuccessiveRestarts = 0;
        fGlobalMinCost      = fCurrentCost;
        oGlobalOptimalState = oTmpResOrient;
        
        oCurrentState        = oTmpResOrient;
        fCurAngularStepSize *= Float( 0.5 );    // reduce step size
      }

      if( fGlobalMinCost < fMaxConvergenceCost )    // convergence  -- move criterion to user defined
        break;                                      // cheating -- should really allow convergence instead
      
      if( nSuccessiveRestarts > nMaxFailedRestarts )   // search failed change this to user defined
        break;
    }
    
    // Named return value optimization
    std::pair<SMatrix3x3, Float> oRes = std::make_pair( oGlobalOptimalState, fGlobalMinCost ); 
    return oRes;
  }


    
}
