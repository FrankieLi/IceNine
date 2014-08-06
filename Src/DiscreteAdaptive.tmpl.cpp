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

#include <ctime>

namespace Reconstruction
{
  
  //----------------------------------------
  //  RunDiscreteSearch
  //
  //  Action:  Perform discrete orientation search and return a
  //           set of candidates
  //
  //  Need to specify sample symmetry 
  //----------------------------------------
  template< class SamplePointT >
  vector< typename DiscreteRefinement<SamplePointT>::SCandidate>
  DiscreteRefinement<SamplePointT>::RunDiscreteSearch( const SamplePointT & oVoxel,
                                                       vector<SMatrix3x3>::const_iterator pFirst,
                                                       vector<SMatrix3x3>::const_iterator pEnd,
                                                       const vector<SMatrix3x3> & oLocalGrid,
                                                       int nQMax,
                                                       LocalVoxelCostFn & CompleteVoxelCostFn,
                                                       Float fAngularRadius  )
  {
    typedef typename ReconstructionSetup::EtaAngularFilter PeakFilterT;
    typedef typename GlobalSearchCostFunctions::VoxelOverlapFn VoxelCostFn;

    //------------------
    //  Last parameter is the radius of the pixel based overlap cost function
    const int nPixelRadius = 3;
    GlobalSearchCostFunctions GlobalSearchCostFns( oSetup.EtaThresholdFilter(), oSimulator, nPixelRadius );
    VoxelCostFn  &  VoxelCostFunction = GlobalSearchCostFns.VoxelCostFn(); 

    AcceptAnyPeakOverlap CandidateAcceptFn;
    OrientationSearch::DiscreteSearchFn<VoxelCostFn, AcceptAnyPeakOverlap>
      DiscreteSearch( oSimulator,
                      oSetup.SampleGeometry(),
                      oSetup.ExperimentalSetup().GetDetectorList(),
                      oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                      oSetup.Data(), oLocalGrid, oVoxel.nPhase,
                      VoxelCostFunction, CandidateAcceptFn );
    
//     vector<SCandidate> OldCandidates
//       = DiscreteSearch.GetCandidates( oVoxel,
//                                       pFirst, pEnd,
//                                       nQMax );
    
    AdaptedLocalVoxelCostFn WrappedCostFn( CompleteVoxelCostFn,
                                           oSetup.SampleGeometry(),
                                           oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                           oSetup.ExperimentalSetup().GetDetectorList(),
                                           oSetup.Data() );
    LatticeSymmetry::CSymmetry * pSym = oSetup.SampleGeometry().GetSampleSymmetry( 1 );
    vector< vector< SCandidate > > oCandidatesByFZ
      = DiscreteSearch.GetSpacedCandidates( oVoxel,
                                            pFirst, pEnd,
                                            nQMax, WrappedCostFn, fAngularRadius,
                                            *pSym );
    
    vector< SCandidate > oSearchPointCandidates;
    for( Size_Type i = 0; i < oCandidatesByFZ.size(); ++ i)
      oSearchPointCandidates.insert( oSearchPointCandidates.end(), oCandidatesByFZ[i].begin(),
                                     oCandidatesByFZ[i].end() );

    std::cout << " || Num candiates " << oSearchPointCandidates.size() << std::endl;
    
    return oSearchPointCandidates;
  }
  
  //----------------------------------------
  //  ReconstructSample
  //----------------------------------------
  template< class SamplePointT >
  std::pair< SamplePointT, Int>
  DiscreteRefinement<SamplePointT>::ReconstructVoxel( const SamplePointT & oInputVoxel )
  {
    //-------------------------------------------------
    //  Build a cost function algorithm (type)
    typedef typename LocalSearchCostFunctions::VoxelOverlapFn                 VoxelCostFn;
    typedef typename OrientationSearch::MCCostAdapter< VoxelCostFn >          MCCostFn;
    typedef typename OrientationSearch::HitRatioConvergenceFn< VoxelCostFn >  HitRatioBasedConvergence;

    //  Instatiate the compuation of the cost function
    PeakFilterT Filter = oSetup.EtaThresholdFilter();
    LocalSearchCostFunctions LocalSearchCostFns( Filter, oSimulator );
    VoxelCostFn    & VoxelCostFunction = LocalSearchCostFns.VoxelCostFn();
    OrientationSearch::LocalOptimizationFn< MCCostFn >  LocalOptimization;
    
    SSearchParameter oTrialSampleSearch;
    SSearchParameter oCompleteLocalSearch;
    SearchDetails::Utilities::InitializeLocalSearchParameter( oTrialSampleSearch, InputParam );
    SearchDetails::Utilities::InitializeLocalSearchParameter( oCompleteLocalSearch, InputParam );
    oTrialSampleSearch.nMaxMCSteps         = 10;  // WARNING -- hard coded
    oTrialSampleSearch.nSuccessiveRestarts = 5;
    
    
    MCCostFn ObjectiveFn( VoxelCostFunction, oSetup.SampleGeometry(),
                          oInputVoxel,
                          oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                          oSetup.ExperimentalSetup().GetDetectorList(), oSetup.Data() );
    typedef typename vector<SMatrix3x3>::const_iterator SMatIter;
    SMatIter pFirst = oSetup.DiscreteSearchPoints().begin();
    SMatIter pEnd = oSetup.DiscreteSearchPoints().end();
    Int nQMax = 5 + InputParam.nMinLocalResolution;   // HACK need to be at least 5 to be informative...  make this a parameter
    vector<SMatrix3x3> oInterestingPoints;
    Float fDiameter = InputParam.fLocalOrientationGridRadius;

    
    vector<SCandidate> oSearchPointCandidates;
    int nHighResElements    = oSetup.InputParameters().nMaxDiscreteCandidates;
    for( int nLevel = 0; nLevel <= InputParam.nMaxLocalResolution; nLevel ++ )
    {
      vector<SMatrix3x3> oLocalGrid;
      OrientationSearch::Utilities::GenerateLocalGrid( oLocalGrid, fDiameter, 0, 1 );   // diameter is really diagonal
      
      oSearchPointCandidates  = RunDiscreteSearch( oInputVoxel, pFirst, pEnd,
                                                   oLocalGrid,  nQMax,
                                                   VoxelCostFunction, fDiameter  );      
      
      nQMax ++;
      for( Size_Type i = 0; i < oSearchPointCandidates.size(); ++ i)
      {
        SamplePointT oTestPoint = oInputVoxel;
        oTestPoint.oOrientMatrix = oSearchPointCandidates[i].oOrientation;
        SOverlapInfo oInfo = VoxelCostFunction( oTestPoint,  oSetup.SampleGeometry(),
                                                oSetup.ExperimentalSetup().GetDetectorList(),
                                                oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                                oSetup.Data() );
        oSearchPointCandidates[i].fCost = float(1) - CostFunctions::Utilities::GetConfidence( oInfo ); 
      }
      
      int nMaxPartialElements = oSearchPointCandidates.size();
      std::sort( oSearchPointCandidates.begin(), oSearchPointCandidates.end() );
      if( oSearchPointCandidates.size() > 0 )
      {
        int nPartialElements = std::min( int( oSearchPointCandidates.size() ), nMaxPartialElements );
        
        oTrialSampleSearch.fLocalGridRadius = ( fDiameter / 3.0 ) ;  // this radius comes from the spacing given by RunDiscretSecret
                                                                     // r = diameter/2, and there are 8 points in the volume equally
                                                                     // spaced, so each volume is around 4.  To make things overlap,
                                                                     // we fudge it to 3.5
        
        oTrialSampleSearch.fLocalGridRadius 
          = std::max( oTrialSampleSearch.fLocalGridRadius, Float( DEGREE_TO_RADIAN( 0.2 ) ) ); // Lower bound by 0.2 Degrees, cost fn width
        
        

        //--------------------------------------
        //
        //   Change this from CalculateCost, which optimizes on the cost, to something that
        //   optimizes for variance may help improve the performance/reliability.  At the
        //   vary minimum, it might be useful to flat out throw away all sample region with
        //   variance below a certain threshold.
        //--------------------------------------
        LocalOptimization.CalculateCost( oInputVoxel,
                                         oSearchPointCandidates.begin(),
                                         oSearchPointCandidates.begin() + nPartialElements,
                                         ObjectiveFn, 
                                         oTrialSampleSearch );
	
        fDiameter /= 1.5;  // fudge factor at this point
                           // respace candidate after small local optimization
        std::sort( oSearchPointCandidates.begin(), oSearchPointCandidates.end() );
        nHighResElements = int( std::floor( oSearchPointCandidates.size() / 4 ) ); // fudge factor
        Int nElements    = std::min( Int( oSearchPointCandidates.size() ), nHighResElements );
        oInterestingPoints.clear();
        for( int i = 0; i < nElements; i ++ )
          oInterestingPoints.push_back( oSearchPointCandidates[i].oOrientation );
        pFirst = oInterestingPoints.begin();
        pEnd = oInterestingPoints.end();
      }
    }

    std::cout << " |  Before final optimization " << std::endl;
    oCompleteLocalSearch.fLocalGridRadius = std::max( Float( fDiameter / 3.0 ), Float( DEGREE_TO_RADIAN( 0.2 ) ) ); // lower bound
    
    int nLocalOptCandidates = std::min( static_cast<Size_Type>( oSetup.InputParameters().nMaxDiscreteCandidates ),
                                        oSearchPointCandidates.size() );

    std::cout << " nLocalOptCandidates " << nLocalOptCandidates << std::endl;
    HitRatioBasedConvergence ConvergenceChecker( VoxelCostFunction, oSetup.SampleGeometry(),
                                                 oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                                 oSetup.ExperimentalSetup().GetDetectorList(),
                                                 oSetup.Data(), 1.0 ); // convergance ratio set to 1.0
    bool bConverged = false;
    SamplePointT Result;
    boost::tie( Result, bConverged )
      =      LocalOptimization.FindOptimal( oInputVoxel,
                                            oSearchPointCandidates.begin(),
                                            oSearchPointCandidates.begin() + nLocalOptCandidates,
                                            ConvergenceChecker,
                                            ObjectiveFn, 
                                            oCompleteLocalSearch );

    std::cout << "  After Final Optimization " << std::endl;
    
    Result = LocalOptimization.Optimize( Result, ObjectiveFn, oCompleteLocalSearch, 0.02 * 0.02 );  // second parameter is variance, which is set to 0.02^2
    std::cout << "  After Variance Optimization " << std::endl;
    SOverlapInfo oInfo = VoxelCostFunction( Result,  oSetup.SampleGeometry(),
                                            oSetup.ExperimentalSetup().GetDetectorList(),
                                            oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                            oSetup.Data() );

    std::cout << "  After Voxel Cost calculation " << std::endl;
    Result.fConfidence        = CostFunctions::Utilities::GetConfidence( oInfo );
    Result.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio( oInfo );
    Result.fCost              = static_cast<Float>( 1 ) - oInfo.fQuality;
    std::cout << "[ \t"
              << Result.fConfidence << "\t "
              << Result.fPixelOverlapRatio << " \t"
              << Float(1) - Result.fCost << " ]" << std::endl;
    
    return std::make_pair( Result,
                           CostFunctions::Utilities::ClassifyConvergence( Result, oCompleteLocalSearch ) );
  }

  //----------------------------------------
  // EvaluateOverlapInfo
  //----------------------------------------
  template< class SamplePointT >
  typename DiscreteRefinement<SamplePointT>::SOverlapInfo
  DiscreteRefinement<SamplePointT>::EvaluateOverlapInfo( const SamplePointT & oVoxel )
  {
    //-------------------------------------------------
    //  Build a cost function algorithm (type)
    typedef typename LocalSearchCostFunctions::VoxelOverlapFn                 VoxelCostFn;
    typedef typename OrientationSearch::MCCostAdapter< VoxelCostFn >          MCCostFn;
    typedef typename OrientationSearch::HitRatioConvergenceFn< VoxelCostFn >  HitRatioBasedConvergence;
    
    //  Instatiate the compuation of the cost function
    PeakFilterT Filter = oSetup.EtaThresholdFilter();
    LocalSearchCostFunctions LocalSearchCostFns( Filter, oSimulator );
    VoxelCostFn    & VoxelCostFunction = LocalSearchCostFns.VoxelCostFn();
    
    SOverlapInfo oInfo = VoxelCostFunction( oVoxel,  oSetup.SampleGeometry(),
                                            oSetup.ExperimentalSetup().GetDetectorList(),
                                            oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                            oSetup.Data() );
    return oInfo;
  }
  
  //----------------------------------------
  //  ReconstructSample
  //----------------------------------------
  template< class SamplePointT >
  typename DiscreteRefinement<SamplePointT>::SOverlapInfo
  DiscreteRefinement<SamplePointT>::LocalOptimization( SamplePointT & oInputVoxel )
  {
    typedef typename LocalSearchCostFunctions::VoxelOverlapFn                 VoxelCostFn;
    typedef typename OrientationSearch::MCCostAdapter< VoxelCostFn >          MCCostFn;
    typedef typename OrientationSearch::HitRatioConvergenceFn< VoxelCostFn >  HitRatioBasedConvergence;

    PeakFilterT Filter = oSetup.EtaThresholdFilter();
    LocalSearchCostFunctions LocalSearchCostFns( Filter, oSimulator );
    VoxelCostFn    & VoxelCostFunction = LocalSearchCostFns.VoxelCostFn();
    
    MCCostFn ObjectiveFn( VoxelCostFunction, oSetup.SampleGeometry(),
                          oInputVoxel,
                          oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                          oSetup.ExperimentalSetup().GetDetectorList(), oSetup.Data() );
    
    SSearchParameter oLocalSearchParam;
    SearchDetails::Utilities::InitializeLocalSearchParameter( oLocalSearchParam, InputParam );

    //    oLocalSearchParam.fLocalGridRadius = DEGREE_TO_RADIAN( 0.75 ); 
    typename OrientationSearch::LocalOptimizationFn< MCCostFn >::LocalMCSearch LocalSearch( ObjectiveFn, oLocalSearchParam );
    
    SCandidate oStartPoint;
    //oStartPoint.oOrientation = oInputVoxel.oOrientMatrix;
    SamplePointT oResult = LocalSearch.VarianceMinimizingOptimization( oInputVoxel, 0.02 * 0.02 ); // second parameter is variance 

    SOverlapInfo oInfo = VoxelCostFunction( oResult,  oSetup.SampleGeometry(),
                                            oSetup.ExperimentalSetup().GetDetectorList(),
                                            oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                            oSetup.Data() );
    oInputVoxel.fConfidence        = CostFunctions::Utilities::GetConfidence( oInfo );
    oInputVoxel.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio( oInfo );
    oInputVoxel.fCost              = Float(1) - oInfo.fQuality;
    oInputVoxel.oOrientMatrix      = oResult.oOrientMatrix;
    
    return oInfo;
  }

  
}
