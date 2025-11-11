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
#include "DiscreteAdaptive.h"


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
  vector<DiscreteRefinement::SCandidate>
  DiscreteRefinement::RunDiscreteSearch( const SVoxel & oVoxel,
                                         vector<SMatrix3x3>::const_iterator pFirst,
                                         vector<SMatrix3x3>::const_iterator pEnd,
                                         const vector<SMatrix3x3> & oLocalGrid,
                                         int nQMax,
                                         LocalVoxelCostFn & CompleteVoxelCostFn,
                                         Float fAngularRadius  )
  {
    typedef ReconstructionSetup::EtaAngularFilter PeakFilterT;
    typedef GlobalSearchCostFunctions::VoxelOverlapFn VoxelCostFn;

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
    
    vector<SCandidate> OldCandidates
      = DiscreteSearch.GetCandidates( oVoxel,
                                      pFirst, pEnd,
                                      nQMax );

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

//       = DiscreteSearch.GetSpacedCandidates( oVoxel,
//                                             pFirst, pEnd,
//                                             nQMax, VoxelCostFunction, fAngularRadius,
//                                             *pSym );

    vector< SCandidate > oSearchPointCandidates;
    for( Size_Type i = 0; i < oCandidatesByFZ.size(); ++ i)
      oSearchPointCandidates.insert( oSearchPointCandidates.end(), oCandidatesByFZ[i].begin(),
                                     oCandidatesByFZ[i].end() );
    
    return oSearchPointCandidates;
  }
  
  //----------------------------------------
  //  ReconstructSample
  //----------------------------------------
  std::pair<SVoxel, Int>
  DiscreteRefinement::ReconstructVoxel( const SVoxel & oInputVoxel )
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
    oTrialSampleSearch.nMaxMCSteps         = 10;  // WARNING -- hard coded
    oTrialSampleSearch.nSuccessiveRestarts = 5;
    
    HitRatioBasedConvergence ConvergenceChecker( VoxelCostFunction, oSetup.SampleGeometry(),
                                                 oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                                 oSetup.ExperimentalSetup().GetDetectorList(),
                                                 oSetup.Data(), 1 );
    
    MCCostFn ObjectiveFn( VoxelCostFunction, oSetup.SampleGeometry(),
                          oInputVoxel,
                          oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                          oSetup.ExperimentalSetup().GetDetectorList(), oSetup.Data() );
    typedef vector<SMatrix3x3>::const_iterator SMatIter;
    SMatIter pFirst = oSetup.DiscreteSearchPoints().begin();
    SMatIter pEnd = oSetup.DiscreteSearchPoints().end();
    Int nQMax = 4;
    vector<SMatrix3x3> oInterestingPoints;
    Float fRadius = InputParam.fLocalOrientationGridRadius;
    
    vector<SCandidate> oSearchPointCandidates;
    int nHighResElements    = 100;
    for( int nLevel = 0; nLevel <= InputParam.nMaxLocalResolution; nLevel ++ )
    {
      vector<SMatrix3x3> oLocalGrid;
      OrientationSearch::Utilities::GenerateLocalGrid( oLocalGrid, fRadius, 0, 1 );
      
      oSearchPointCandidates  = RunDiscreteSearch( oInputVoxel, pFirst, pEnd,
                                                   oLocalGrid,  nQMax,
                                                   VoxelCostFunction, fRadius  );
      nQMax ++;
      for( Size_Type i = 0; i < oSearchPointCandidates.size(); ++ i)
      {
        SVoxel oTestPoint = oInputVoxel;
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
        fRadius /= 1.5;  // fudge factor at this point
        oTrialSampleSearch.fLocalGridRadius = fRadius;
        LocalOptimization.CalculateCost( oInputVoxel,
                                         oSearchPointCandidates.begin(),
                                         oSearchPointCandidates.begin() + nPartialElements,
                                         ObjectiveFn, 
                                         oTrialSampleSearch );

        // respace candidate after small local optimization
        
        std::sort( oSearchPointCandidates.begin(), oSearchPointCandidates.end() );
        Int nElements = std::min( Int( oSearchPointCandidates.size() ), nHighResElements );
        
        oInterestingPoints.clear();
        for( int i = 0; i < nElements; i ++ )
          oInterestingPoints.push_back( oSearchPointCandidates[i].oOrientation );
        pFirst = oInterestingPoints.begin();
        pEnd = oInterestingPoints.end();
        nHighResElements  /= 1;
      }
    }

    oCompleteLocalSearch.fLocalGridRadius = fRadius;
    int nLocalOptCandidates = std::min( static_cast<Size_Type>( 80 ),
                                        oSearchPointCandidates.size() );   // -------- WARNING hard code                                       
    bool bConverged = false;
    SVoxel Result;
    std::tie( Result, bConverged )
      =      LocalOptimization.FindOptimal( oInputVoxel,
                                            oSearchPointCandidates.begin(),
                                            oSearchPointCandidates.begin() + nLocalOptCandidates,
                                            ConvergenceChecker,
                                            ObjectiveFn, 
                                            oCompleteLocalSearch );
    SOverlapInfo oInfo = VoxelCostFunction( Result,  oSetup.SampleGeometry(),
                                            oSetup.ExperimentalSetup().GetDetectorList(),
                                            oSetup.ExperimentalSetup().GetRangeToIndexMap(),
                                            oSetup.Data() );

        
    Result.fConfidence        = CostFunctions::Utilities::GetConfidence( oInfo );
    Result.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio( oInfo );
    //std::cout << "[ \t" << Result.fConfidence << "\t " << Result.fPixelOverlapRatio << " \t" << Float(1) - Result.fCost << " ]" << std::endl;

    return std::make_pair( Result,
                           CostFunctions::Utilities::ClassifyConvergence( Result, oCompleteLocalSearch ) );
  }

  //----------------------------------------
  // EvaluateOverlapInfo
  //----------------------------------------
  DiscreteRefinement::SOverlapInfo DiscreteRefinement::EvaluateOverlapInfo( const SVoxel & oVoxel )
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
  
  //----------------------------------------
  //  ReconstructSample
  //----------------------------------------
  DiscreteRefinement::SOverlapInfo
  DiscreteRefinement::LocalOptimization( SVoxel & oInputVoxel )
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
    OrientationSearch::LocalOptimizationFn< MCCostFn >::LocalMCSearch LocalSearch( ObjectiveFn, oLocalSearchParam );
    
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
    oInputVoxel.fConfidence        = CostFunctions::Utilities::GetConfidence( o// Info );
//     oInputVoxel.fPixelOverlapRatio = CostFunctions::Utilities::GetHitRatio( oInfo );
//     oInputVoxel.fCost              = Float(1) - oInfo.fQuality;    
//     return oInfo;
//   }


  
// }
