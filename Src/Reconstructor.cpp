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
//  Reconstructor.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//   Purpose:  Interface and implementations of simple reconstructors
//
////////////////////////////////////////////////////////////


#include "Reconstructor.h"


namespace Reconstruction
{

  //--------------------------------------------------
  //  DEBUG_OutputVoxelCost()
  //
  //  Output the cost of this voxel as a function of orientation
  //--------------------------------------------------
  void BasicVoxelReconstructor::DEBUG_OutputVoxelCost( const SVoxel & TestVoxel )
  {
    typedef LocalSearchCostFunctions::VoxelOverlapFn                 VoxelCostFn;
    typedef OrientationSearch::MCCostAdapter< VoxelCostFn >          MCCostFn;
    typedef OrientationSearch::HitRatioConvergenceFn< VoxelCostFn >  HitRatioBasedConvergence;
    
    //  Instatiate the compuation of the cost function
    PeakFilterT Filter = oSetup.EtaThresholdFilter();
    LocalSearchCostFunctions LocalSearchCostFns( Filter, oSimulator );
    VoxelCostFn    & VoxelCostFunction = LocalSearchCostFns.VoxelCostFn();

    std::ofstream SearchPointCost("SearchPointCost.txt");

    vector<SMatrix3x3> SamplePoints;
    int nLine = 0;
    typedef vector<SMatrix3x3>::const_iterator MatrixIter;
    for( MatrixIter MatIt =  oSetup.DiscreteSearchPoints().begin();
         MatIt != oSetup.DiscreteSearchPoints().end(); ++ MatIt )
    {
      for( int n = 0; n <= 2; n ++ )
      {
        vector<SMatrix3x3> & CurrentLevel = AllLevelLocalGrids[n];
        for( int m = 0; m < CurrentLevel.size(); m ++ )
        {
          SVoxel TestPoint = TestVoxel;
          TestPoint.oOrientMatrix = CurrentLevel[m] * (*MatIt);
          SOverlapInfo oInfo = VoxelCostFunction( TestPoint,  oSetup.SampleGeometry(),
                                                  oSetup.ExperimentalSetup().GetDetectorList(),
                                                  oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                                  oSetup.Data() );
          
          SearchPointCost << CostFunctions::Utilities::GetConfidence( oInfo ) << " "
                          <<  GeneralLib::RadianToDegree( TestPoint.oOrientMatrix.GetEulerAngles() ) << std::endl;
          nLine ++;
        }
      }
    }
    std::cout << nLine << std::endl;
    for( int n = 0; n <= 3; n ++ )
    {
      vector<SMatrix3x3> & CurrentLevel = AllLevelLocalGrids[n];
      for( int m = 0; m < CurrentLevel.size(); m ++ )
      {
        SVoxel TestPoint = TestVoxel;
        TestPoint.oOrientMatrix = CurrentLevel[m] * TestVoxel.oOrientMatrix;
        SOverlapInfo oInfo = VoxelCostFunction( TestPoint,  oSetup.SampleGeometry(),
                                                oSetup.ExperimentalSetup().GetDetectorList(),
                                                oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                                oSetup.Data() );
        
        SearchPointCost << CostFunctions::Utilities::GetConfidence( oInfo ) << " "
                        <<  GeneralLib::RadianToDegree( TestPoint.oOrientMatrix.GetEulerAngles() ) << std::endl;
        
      }
    }
    
    
    std::cout << "DEBUG COMPLETED, exit without running reconstruction" << std::endl;
    exit(0);
    //-------------------------------------------------------
    
  }

  
  //----------------------------------------
  //  RunDiscreteSearch
  //
  //  Action:  Perform discrete orientation search and return a
  //           set of candidates
  //----------------------------------------
  vector<BasicVoxelReconstructor::SCandidate>
  BasicVoxelReconstructor::RunDiscreteSearch( const SVoxel & oVoxel,
                                              const vector<SMatrix3x3> & oLocalGrid )
  {
    //---------------------------------------
    //  DEBUG
    //   DEBUG_OutputVoxelCost( oVoxel );
    //---------------------------------------
    
    typedef ReconstructionSetup::EtaAngularFilter PeakFilterT;
    typedef GlobalSearchCostFunctions::VoxelOverlapFn VoxelCostFn;

    //------------------
    //  Last parameter is the radius of the pixel based overlap cost function
    const int nPixelRadius = 1;
    GlobalSearchCostFunctions GlobalSearchCostFns( oSetup.EtaThresholdFilter(), oSimulator, nPixelRadius );
        
    VoxelCostFn  &  VoxelCostFunction = GlobalSearchCostFns.VoxelCostFn(); 
    
    AcceptAnyPeakOverlap CandidateAcceptFn;
    const Int nQMaxDiscrete  = 5;
    OrientationSearch::DiscreteSearchFn<VoxelCostFn, AcceptAnyPeakOverlap>
      DiscreteSearch( oSimulator,
                      oSetup.SampleGeometry(),
                      oSetup.ExperimentalSetup().GetDetectorList(),
                      oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                      oSetup.Data(), oLocalGrid, oVoxel.nPhase,
                      VoxelCostFunction, CandidateAcceptFn );

    vector<SCandidate> oSearchPointCandidates
      = DiscreteSearch.GetCandidates( oVoxel,
                                      oSetup.DiscreteSearchPoints().begin(),
                                      oSetup.DiscreteSearchPoints().end(),
                                      nQMaxDiscrete );
    //-------------------------------------------------------
    //  For testing
    
//    std::ofstream SearchPointCost("SearchPointCost.txt");
//    for( int i = 0; i < oSearchPointCandidates.size(); i ++ )
//    {
//      SearchPointCost << oSearchPointCandidates[i].fCost << " " <<  GeneralLib::RadianToDegree(oSearchPointCandidates[i].oOrientation.GetEulerAngles()) << std::endl;
//    }
//    std::cout << "DEBUG COMPLETED, exit without running reconstruction" << std::endl;
//    exit(0);
    //-------------------------------------------------------
    return oSearchPointCandidates;
  }
  
  //----------------------------------------
  //  ReconstructSample
  //
  //  Basic exhaustive search with iterative deepening.
  //
  //----------------------------------------
  std::pair<SVoxel, Int>  BasicVoxelReconstructor::ReconstructVoxel( const SVoxel & oInputVoxel )
  {
    //-------------------------------------------------
    //  Build a cost function algorithm (type)
    typedef LocalSearchCostFunctions::VoxelOverlapFn                 VoxelCostFn;
    typedef OrientationSearch::MCCostAdapter< VoxelCostFn >          MCCostFn;
    typedef OrientationSearch::HitRatioConvergenceFn< VoxelCostFn >  HitRatioBasedConvergence;

    //  Instatiate the compuation of the cost function
    PeakFilterT Filter = oSetup.EtaThresholdFilter();
    LocalSearchCostFunctions LocalSearchCostFns( Filter, oSimulator );
    VoxelCostFn    & VoxelCostFunction = LocalSearchCostFns.VoxelCostFn();
    OrientationSearch::LocalOptimizationFn< MCCostFn >  LocalOptimization;
    
    SSearchParameter oTrialSampleSearch;
    SSearchParameter oCompleteLocalSearch;
    SearchDetails::Utilities::InitializeLocalSearchParameter( oTrialSampleSearch, InputParam );
    SearchDetails::Utilities::InitializeLocalSearchParameter( oCompleteLocalSearch, InputParam );
    oTrialSampleSearch.nMaxMCSteps         = 20;  // WARNING -- hard coded
    oTrialSampleSearch.nSuccessiveRestarts = 0;
    
    HitRatioBasedConvergence ConvergenceChecker( VoxelCostFunction, oSetup.SampleGeometry(),
                                                 oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                                 oSetup.ExperimentalSetup().GetDetectorList(),
                                                 oSetup.Data(), InputParam.fMaxDeepeningHitRatio );
    
    MCCostFn ObjectiveFn( VoxelCostFunction, oSetup.SampleGeometry(),
                          oInputVoxel,
                          oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                          oSetup.ExperimentalSetup().GetDetectorList(), oSetup.Data() );
    SVoxel Result;
    for( int nLevel = 0; nLevel <= InputParam.nMaxLocalResolution; nLevel ++ )
    {
      vector<SCandidate> oSearchPointCandidates  = RunDiscreteSearch( oInputVoxel,
                                                                      AllLevelLocalGrids[ nLevel ] );
      LocalOptimization.CalculateCost( oInputVoxel,
                                       oSearchPointCandidates.begin(),
                                       oSearchPointCandidates.end(),
                                       ObjectiveFn, 
                                       oTrialSampleSearch );
      std::sort( oSearchPointCandidates.begin(), oSearchPointCandidates.end() );
      Int nElements = std::min( Int( oSearchPointCandidates.size() ), InputParam.nMaxDiscreteCandidates );
      
      bool bConverged = false;
      std::tie( Result, bConverged )
        =      LocalOptimization.FindOptimal( oInputVoxel,
                                              oSearchPointCandidates.begin(),
                                              oSearchPointCandidates.begin() + nElements,
                                              ConvergenceChecker,
                                              ObjectiveFn, 
                                              oCompleteLocalSearch );
      if( bConverged )
        break;
    }
    SOverlapInfo oInfo = VoxelCostFunction( Result,  oSetup.SampleGeometry(),
                                            oSetup.ExperimentalSetup().GetDetectorList(),
                                            oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                            oSetup.Data() );
    
    Result.fConfidence        = CostFunctions::Utilities::GetConfidence( oInfo );
    Result.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio( oInfo );
    std::cout << "[ " << Result.fConfidence << " " << Result.fPixelOverlapRatio << " " << Float(1) - Result.fCost << " ]" << std::endl;
    return std::make_pair( Result,
                           CostFunctions::Utilities::ClassifyConvergence( Result, oCompleteLocalSearch ) );
  }


  //----------------------------------------
  //  ReconstructSample
  //----------------------------------------
  BasicVoxelReconstructor::SOverlapInfo
  BasicVoxelReconstructor::LocalOptimization( SVoxel & oInputVoxel )
  {
    typedef LocalSearchCostFunctions::VoxelOverlapFn                 VoxelCostFn;
    typedef OrientationSearch::MCCostAdapter< VoxelCostFn >          MCCostFn;
    typedef OrientationSearch::HitRatioConvergenceFn< VoxelCostFn >  HitRatioBasedConvergence;

    PeakFilterT Filter = oSetup.EtaThresholdFilter();
    LocalSearchCostFunctions LocalSearchCostFns( Filter, oSimulator );
    VoxelCostFn    & VoxelCostFunction = LocalSearchCostFns.VoxelCostFn();

    MCCostFn ObjectiveFn( VoxelCostFunction, oSetup.SampleGeometry(),
                          oInputVoxel,
                          oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                          oSetup.ExperimentalSetup().GetDetectorList(), oSetup.Data() );
    
    SSearchParameter oLocalSearchParam;
    SearchDetails::Utilities::InitializeLocalSearchParameter( oLocalSearchParam, InputParam );
    OrientationSearch::LocalOptimizationFn< MCCostFn >::LocalMCSearch LocalSearch( ObjectiveFn, oLocalSearchParam  );
        
    HitRatioBasedConvergence ConvergenceChecker( VoxelCostFunction, oSetup.SampleGeometry(),
                                                 oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                                 oSetup.ExperimentalSetup().GetDetectorList(),
                                                 oSetup.Data(), InputParam.fMaxDeepeningHitRatio );
    SCandidate oStartPoint;
    oStartPoint.oOrientation = oInputVoxel.oOrientMatrix;
    SVoxel oResult;
    LocalSearch( oResult, oInputVoxel, oStartPoint ); 

    SOverlapInfo oInfo = VoxelCostFunction( oResult,  oSetup.SampleGeometry(),
                                            oSetup.ExperimentalSetup().GetDetectorList(),
                                            oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                            oSetup.Data() );
    oResult.fConfidence        = CostFunctions::Utilities::GetConfidence( oInfo );
    oResult.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio( oInfo );
    oResult.fCost              = Float(1) - oInfo.fQuality;
    
    oInputVoxel = oResult;
    return oInfo;
  }
 
  //----------------------------------------
  // EvaluateOverlapInfo
  //----------------------------------------
  BasicVoxelReconstructor::SOverlapInfo BasicVoxelReconstructor::EvaluateOverlapInfo( const SVoxel & oVoxel )
  {
    //-------------------------------------------------
    //  Build a cost function algorithm (type)
    typedef LocalSearchCostFunctions::VoxelOverlapFn                 VoxelCostFn;
    typedef OrientationSearch::MCCostAdapter< VoxelCostFn >          MCCostFn;
    typedef OrientationSearch::HitRatioConvergenceFn< VoxelCostFn >  HitRatioBasedConvergence;
    
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
  
} // end namespace


