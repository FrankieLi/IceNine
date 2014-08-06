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
//----------------------------------------------------------
//
//  ContinuousSearch.h
//
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Implementation of continuous search through templates.
//
//
//----------------------------------------------------------


#ifndef CONTINUOUSE_SEARCH_H
#define CONTINUOUSE_SEARCH_H
#include "OrientationSearch.h"
#include "SearchDetails.h"
#include "OverlapInfo.h"
#include "SimulationData.h"
#include "3dMath.h"
#include <limits>
#include "boost/tuple/tuple.hpp"

namespace OrientationSearch
{

  namespace Details
  {
    //--------------------------------------------------
    //  CalculateSearchParameter
    //--------------------------------------------------
    template< class SSearchParameter >
    std::pair<Float, Float> CalculateSearchParameter( const SSearchParameter & s )
    {
      Float BoxWidth = s.fLocalGridRadius / pow( Float(2), s.nLocalResolution );
      Float Radius   = BoxWidth * s.fMCRadiusScaleFactor;
      return std::make_pair( BoxWidth, Radius );
    }
  }  // Details
  
  //----------------------------------------------------------
  //  Continuous local search process
  //   Optimizes the orientation of a sample point assuming
  //   that the cost function has a locally continuous (or at least
  //   not ridiculously difficult) landscape.
  //
  //
  //   Cost function in this case converts SOverlapInfo -> float [0, 1].
  //   Cost of 0 being the best.
  //----------------------------------------------------------
  template< class CostFnT >
  class ContinuousLocalOptimizationProcess
  {
  public:
    typedef typename CostFnT::SamplePointT SamplePointT;
  private:
    
    CostFnT  CostFunction;    // Cost function dictates how the overlap is calculated
    OrientationOptimization::COrientationMC oMCOptimizer;   // need singleton
    
    //----------------
    //  Search Parameters
    //----------------
    Float fSearchBoxWidth;
    Float fRadius;
    Int   nMaxMCStep;
    Float fMaxConvergenceCost;
    Int   nSuccessiveRestarts;

    ContinuousLocalOptimizationProcess();  // this is not legal
  public:
    typedef SearchDetails::SCandidate   SCandidate;
    typedef SamplePointT ResultType;

    ContinuousLocalOptimizationProcess( CostFnT CostFn_,
                                        Float _SearchBoxWidth, Float _Radius,
                                        Int _MaxMCStep, Float _MaxConvergenceCost,
                                        Int _SuccessiveRestarts )
      : CostFunction( CostFn_ ),
        fSearchBoxWidth( _SearchBoxWidth ),
        fRadius        ( _Radius ),
        nMaxMCStep     ( _MaxMCStep ),
        fMaxConvergenceCost( _MaxConvergenceCost ),
        nSuccessiveRestarts( _SuccessiveRestarts )
    {}
    
    template< class SearchParam >
    ContinuousLocalOptimizationProcess( CostFnT CostFn_,
                                        const SearchParam & LocalSearchParam )
      : CostFunction( CostFn_ )
    {
      
      boost::tie( fSearchBoxWidth, fRadius )
        = Details::CalculateSearchParameter( LocalSearchParam );
      nMaxMCStep          = LocalSearchParam.nMaxMCSteps;
      fMaxConvergenceCost = LocalSearchParam.fMaxConvergenceCost;
      nSuccessiveRestarts = LocalSearchParam.nSuccessiveRestarts;
    }
    
    void Reinitialize( Float _SearchBoxWidth, Float _Radius,
                       Int _MaxMCStep, Float _MaxConvergenceCost,
                       Int _SuccessiveRestarts )
    {
      fSearchBoxWidth     = _SearchBoxWidth;
      fRadius             = _Radius;
      nMaxMCStep          = _MaxMCStep;
      fMaxConvergenceCost = _MaxConvergenceCost;
      nSuccessiveRestarts = _SuccessiveRestarts;
    }

    //-----------------------------------------------------
    //  Need to change this to remove voxel dependency
    //-----------------------------------------------------
    void operator() ( ResultType & oRes, const SamplePointT & v,
                      const SCandidate & oCandidate )
    {
      SMatrix3x3 oNewOrient;
      Float fCurrentCost;
      boost::tie( oNewOrient, fCurrentCost)
        = oMCOptimizer.RandomRestartZeroTemp( oCandidate.oOrientation,
                                              fRadius, fSearchBoxWidth,
                                              CostFunction, nMaxMCStep,
                                              nSuccessiveRestarts,
                                              fMaxConvergenceCost );
      
      oRes = v;
      oRes.oOrientMatrix = oNewOrient;
      oRes.fCost = fCurrentCost;
    }

   
    //-----------------------------------------------------
    //  Need to change this to remove voxel dependency
    //-----------------------------------------------------
    ResultType VarianceMinimizingOptimization( const SamplePointT & v,
                                               Float ConvergenceVariance )
    {
      SMatrix3x3 oNewOrient;
      Float fCurrentCost;
      boost::tie( oNewOrient, fCurrentCost)
        = oMCOptimizer.AdaptiveSamplingZeroTemp( v.oOrientMatrix,
                                                 DEGREE_TO_RADIAN( 0.5 ), fSearchBoxWidth,
                                                 CostFunction, nMaxMCStep,
                                                 nSuccessiveRestarts,
                                                 ConvergenceVariance,
                                                 fMaxConvergenceCost );
      ResultType oRes = v;
      oRes.oOrientMatrix = oNewOrient;
      oRes.fCost = fCurrentCost;
      return oRes;
    }

    
  };

  //----------------------------------------------------------
  //  Adapater function
  //----------------------------------------------------------
  template< class OverlapCostFn >
  class MCCostAdapter : public OrientationOptimization::COverlapFunction 
  {
  public:
    typedef OverlapCostFn VoxelCostFn;
    typedef typename VoxelCostFn::SamplePointT SamplePointT;
  private:
    OverlapCostFn CostFn;
    CSample & oSample;
    SamplePointT oVoxel;
    const XDMSimulation::CSimulationRange & oRangeToIndexMap;
    const DetectorListT & oDetectorList;
    const CSimulationData & oExpData;
    typedef CostFunctions::SOverlapInfo SOverlapInfo;

  public:

    MCCostAdapter( OverlapCostFn CostFn_, CSample & s_, const SamplePointT & v_,
                   const XDMSimulation::CSimulationRange & oRangeMap_,
                   const DetectorListT & _DetList,
                   const CSimulationData & oData_ )
      : CostFn( CostFn_ ), oSample( s_ ), oVoxel( v_ ),
        oRangeToIndexMap( oRangeMap_ ),
        oDetectorList( _DetList ),
        oExpData( oData_ ) {}

    //------------------------------------------
    //  Actual wrapper
    //------------------------------------------
    Float operator()( const SMatrix3x3 & oNewOrient )
    {
      oVoxel.oOrientMatrix = oNewOrient;
      SOverlapInfo oOverlapInfo
        = CostFn( oVoxel, oSample, oDetectorList,
                  oRangeToIndexMap, oExpData );

      if( oOverlapInfo.nPixelOnDetector <= 0 )
        return 1;
      else
        return Float( 1 ) - oOverlapInfo.fQuality;
    }

    Float operator()( const SamplePointT & v )
    {
      RUNTIME_ASSERT( 0, "Call to operator()( SVoxel ) Unexpected\n ");
      return std::numeric_limits<Float>::signaling_NaN();
    }
  };

  //------------------------------------------
  //  Define convergence as meeting a specific hit ratio
  //------------------------------------------
  template< class OverlapCostFn >
  class HitRatioConvergenceFn
  {
  public:
    typedef typename OverlapCostFn::SamplePointT SamplePointT;
  private:
    OverlapCostFn  CostFn;
    CSample        & oSample;
    const DetectorListT   & oDetList;
    const XDMSimulation::CSimulationRange & oRangeToIndexMap;
    const CSimulationData & oExpData;
    
    typedef CostFunctions::SOverlapInfo SOverlapInfo;

    Float fMinConvergenceHitRatio;
  public:

    HitRatioConvergenceFn( OverlapCostFn CostFn_, CSample & s_,
                           const XDMSimulation::CSimulationRange & oRangeMap_,
                           const DetectorListT & DetList_,
                           const CSimulationData & oData_,
                           Float fConvergenceHitRatio_ )
      : CostFn( CostFn_ ), oSample( s_ ),
        oDetList( DetList_ ),
        oRangeToIndexMap( oRangeMap_ ),
        oExpData( oData_ ), fMinConvergenceHitRatio( fConvergenceHitRatio_ ) {}
    
    //------------------------------------------
    //  Actual wrapper
    //------------------------------------------
    bool operator()( const SamplePointT & v )
    {
      SOverlapInfo oOverlapInfo
        = CostFn( v, oSample, oDetList, oRangeToIndexMap, oExpData );

      if( oOverlapInfo.nPixelOnDetector <= 0 )
        return false;
      else
        return ( Float( oOverlapInfo.nPixelOverlap ) / Float( oOverlapInfo.nPixelOnDetector )
                 >= fMinConvergenceHitRatio ) ; 
    }
  };
  
  //----------------------------------------------------------
  //  LocalOptimization
  //
  //  TODO:  Consider deprecating this.  This template
  //         really isn't doing all that much, and can be
  //         replaced by a functor applicator function (i.e., FOR_EACH) 
  //----------------------------------------------------------
  template< class ObjectiveFnT >
  class LocalOptimizationFn
  {
  public:
    typedef typename ObjectiveFnT::SamplePointT SamplePointT;
    typedef OrientationSearch::ContinuousLocalOptimizationProcess< ObjectiveFnT > LocalMCSearch;
    typedef ObjectiveFnT CostFunction;
  private:

    typedef SearchDetails::SCandidate       SCandidate;
    typedef SearchDetails::SSearchParameter SSearchParameter;
    
  public:
   
    LocalOptimizationFn(){}

    //--------------------------------------------------
    //  FindOptimal
    //--------------------------------------------------
    template< class ConvergenceFn, class CandidateIter > 
    std::pair<SamplePointT, bool> FindOptimal( const SamplePointT & oVoxel,
                                               CandidateIter pCur, CandidateIter pEnd,
                                               ConvergenceFn ConvergenceCheck,
                                               ObjectiveFnT & ObjectiveFn, 
                                               const SSearchParameter & oSearchParam )
    {
      LocalMCSearch      LocalSearch( ObjectiveFn, oSearchParam );
      bool   bConverged = false;
      SamplePointT oRes = oVoxel;
      Float fBestCost = 1000;
      int nBest = 0;
      int nCur = 0;
      for( ; pCur != pEnd; ++ pCur )
      {
        SamplePointT oCurrent;
        LocalSearch( oCurrent, oVoxel, *pCur );
        if( fBestCost > oCurrent.fCost )
        {
          nBest = nCur;
          fBestCost = oCurrent.fCost;
          oRes      = oCurrent;
          bConverged = ConvergenceCheck( oCurrent ); 
          if( bConverged )
            break;
        }
	//	std::cout << nCur << " --- FinalCost " << oCurrent.fCost << " " << RadianToDegree( oCurrent.oOrientMatrix.GetEulerAngles() ) <<  std::endl;
        nCur++;
      }
      //      std::cout << "Best found at " << nBest << std::endl;
      return std::make_pair( oRes, bConverged );
    }

    //--------------------------------------------------
    //  CalculateCost -
    //  Purpose:  Calculate cost for each of the candidates
    //
    //--------------------------------------------------
    template< class CandidateIter > 
    void CalculateCost( const SamplePointT & oVoxel,
                        CandidateIter pCur, CandidateIter pEnd,
                        ObjectiveFnT & ObjectiveFn, 
                        const SSearchParameter oSearchParam )
    {
      LocalMCSearch      LocalSearch( ObjectiveFn, oSearchParam );
      SamplePointT oRes = oVoxel;
      for( ; pCur != pEnd; ++ pCur )
      {
        LocalSearch( oRes, oVoxel, *pCur );
        pCur->oOrientation = oRes.oOrientMatrix;
        pCur->fCost        = oRes.fCost;
      }
    }

    //--------------------------------------------------
    //  Optimize
    //    - Optimize the orientation of the voxel with respect
    //      to the ObjectiveFn.
    //
    //      ConvergenceVariance should really be in SearchParameter - a hack currently
    //
    //--------------------------------------------------
    SamplePointT Optimize( const SamplePointT & Voxel,
                           ObjectiveFnT & ObjectiveFn, 
                           SSearchParameter oSearchParam,
                           Float ConvergenceVariance )  
    {
      oSearchParam.fMaxConvergenceCost = 0;   // this is for final optimization
      LocalMCSearch      LocalSearch( ObjectiveFn, oSearchParam );
      return LocalSearch.VarianceMinimizingOptimization( Voxel, ConvergenceVariance );
    }

  };
} //  namespace Orientation Search

#endif
