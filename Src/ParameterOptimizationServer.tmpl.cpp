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


namespace ParallelReconstructor
{

  namespace ParameterOptimization
  {

    //----------------------------------
    //  LargeVolumeOptimization
    //
    //  This is VERY slow -- use this to home into
    //  closer answers
    //
    //----------------------------------
    template< class SamplePointT, class SamplePointGrid >
    vector< SLargeScaleOpt >
    ParameterOptimizationServer< SamplePointT, SamplePointGrid >
    ::GetLargeVolumeCandidates( const SEnergyOpt            &oStartEnergyLoc,
                                const vector<SDetParamMsg>  &vDetParams,
                                const vector<SStepSizeInfo> & vOptMaxDeviation )
    {
      typedef GeometricOptimizationBase<SamplePointT> Base;
      VoxelQueue.RandomizeReset();
      vector<SamplePointT> oLargeSearchList;

      VoxelQueue.Get( LocalSetup.InputParameters().nParamMCGlobalSearchElements, std::back_inserter( oLargeSearchList ) );
      vector<SLargeScaleOpt>        vCandidates;
      vector<vector<SDetParamMsg> > oDetParamList( nProcessingElements );
      vector<Float>                 oEnergyList  ( nProcessingElements );

      GET_LOG( osLogFile ) << "Large Scale Volume Candidate Search " << std::endl;
      GET_LOG( osLogFile ) << "Beam Energy " << oStartEnergyLoc.fBeamEnergy << std::endl;

      for( Size_Type i = 0; i < vOptMaxDeviation.size(); i ++ )
      {
        GET_LOG( osLogFile ) << "---------------------" << std::endl;
        GET_LOG( osLogFile ) << "Euler Angles Steps:         " << vOptMaxDeviation[i].oEulerSteps << std::endl;
        GET_LOG( osLogFile ) << "vOptMaxDeviation[i].fBeamCenterJ : " << vOptMaxDeviation[i].fBeamCenterJ << std::endl;
        GET_LOG( osLogFile ) << "vOptMaxDeviation[i].fBeamCenterK : " << vOptMaxDeviation[i].fBeamCenterK << std::endl;
        GET_LOG( osLogFile ) << "vOptMaxDeviation[i].xPos :         " << vOptMaxDeviation[i].oDetectorPos.m_fX << std::endl;
        GET_LOG( osLogFile ) << "---------------------" << std::endl;
      }


      for( Int nClientID = 1; nClientID < nProcessingElements; nClientID ++ )
      {
        vector<SDetParamMsg> vNewDetectorLoc;
        SEnergyOpt oNewEnergyLoc;
        if ( nClientID > 1 ) // have one of them start from the origin
        {
          vNewDetectorLoc           = Base::RandomMoveDet( vDetParams, vOptMaxDeviation );
          oNewEnergyLoc.fBeamEnergy = oStartEnergyLoc.fBeamEnergy
                                    + oRandomReal( -oStartEnergyLoc.fEnergyStep, oStartEnergyLoc.fEnergyStep );
        }
        else
        {
          vNewDetectorLoc           = vDetParams;
          oNewEnergyLoc.fBeamEnergy = oStartEnergyLoc.fBeamEnergy;
        }

        oEnergyList  [ nClientID ] = oNewEnergyLoc.fBeamEnergy;
        oDetParamList[ nClientID ] = vNewDetectorLoc;

        Base::Comm.SendCommand( 0, nClientID, XDMParallel::SET_EXP_PARAM );
        Base::SendExpParameters( nClientID, oNewEnergyLoc, vNewDetectorLoc );
        GET_LOG( osLogFile ) << " Sending " << oLargeSearchList.size() << " voxels to client " << nClientID << std::endl;
        Base::Comm.SendCommand( 0, nClientID, XDMParallel::FIT_MC_LIST );
        Base::Comm.SendWorkUnitList( nClientID, oLargeSearchList );
      }

      // listen for result
      Int nClientsLeft = nProcessingElements - 1;
      vector<SLargeScaleOpt> oCandidateList;
      while( nClientsLeft > 0 )
      {
        Int nCommand, nClientID;
        Base::Comm.RecvCommand( &nClientID, &nCommand );
        RUNTIME_ASSERT( nCommand == XDMParallel::REPORT_MC_LIST, "Server ERROR!  Wrong comamnd recv'd \n");
        vector< SParamOptMsg<SamplePointT> > vOptResults;
        Base::Comm.RecvWorkUnitList( nClientID, vOptResults );

        // calculate cost
        Float fCost   = 0;
        Int nFitted   = 0;
        Int nUnfitted = 0;
        for( Size_Type i = 0; i < vOptResults.size(); i ++ )
        {
          if( vOptResults[i].bConverged )
          {
            fCost += Float(1) - vOptResults[i].oOverlapInfo.fQuality;
            nFitted ++;
          }
          else
          {
            fCost += Float(1);
            nUnfitted ++;
          }
        }
        GET_LOG( osLogFile ) << "Client " << nClientID << " Fitted "
                             << nFitted << " Unfitted " << nUnfitted << " cost = " << fCost << std::endl;
        if( fCost < static_cast<Float>( oLargeSearchList.size() ) )
        {
          SLargeScaleOpt oNewPoint;
          oNewPoint.fCost      = fCost / static_cast<Float>( oLargeSearchList.size() );
          oNewPoint.fEnergy    = oEnergyList[ nClientID ];
          oNewPoint.vDetParams = oDetParamList[ nClientID ];
          vCandidates.push_back( oNewPoint );
        }
        nClientsLeft --;
      }
      GET_LOG( osLogFile ) << " Finished Large Scale Optimization: Num Candidates = " << vCandidates.size() << std::endl;
      return vCandidates;
    }



    //---------------------------------------------------------------------------
    //  GetNConvergedElements
    //---------------------------------------------------------------------------
    template<  class SamplePointT, class SamplePointGrid >
    std::pair< std::vector< SParamOptMsg< SamplePointT > >, bool >
    ParameterOptimizationServer< SamplePointT, SamplePointGrid >
    ::GetNConvergedElements( )
    {
      typedef GeometricOptimizationBase<SamplePointT> Base;
      const Int nMaxElements = nProcessingElements * LocalSetup.InputParameters().nOptNumElementPerPE * 3;
      if( VoxelQueue.Size() > nMaxElements )
      {
        GET_LOG( osLogFile ) << "Number of elements from structure: " << VoxelQueue.Size() << std::endl;
        GET_LOG( osLogFile ) << "   exceeds the maximum number of elements to be fitted: " << nMaxElements << std::endl;
        GET_LOG( osLogFile ) << "Limiting to " <<   nMaxElements << " elements " << std::endl;
        VoxelQueue.RandomizedSetMaxElements( nMaxElements );
      }

      std::queue<int> WaitQueue;
      for( int i = 1; i <= (nProcessingElements - 1); i ++ )  // start from 1, since 0 is Server
        WaitQueue.push( i );
      typedef XDMParallel::MultiElementDistribution<SamplePointT>       MultiElementDistributor;
      MultiElementDistributor SingleVoxelDistributor(1);
      Utilities::WorkUnitDistribution( VoxelQueue, SingleVoxelDistributor, WaitQueue, Base::Comm, osLogFile, 0,  XDMParallel::FIT_MC );
      GET_LOG( osLogFile ) << "Finished with Initial Work Unit Distribution " << std::endl;

      typedef SParamOptMsg< SamplePointT > ParamOptMsg;
      vector<ParamOptMsg> oConvergedSamplePoints;
      Int nElementsToFit  = LocalSetup.InputParameters().nOptNumElementPerPE;   // parameters?
      int NumClients      = nProcessingElements -1;
      Int NumElementsUsed = NumClients;
      while( ! VoxelQueue.Empty() || WaitQueue.size() < NumClients
             || Int( oConvergedSamplePoints.size() ) >= nElementsToFit )
      {
        Int                  nClientPE;
        Int                  nCommand;
        vector<ParamOptMsg> TmpResult;
        ParamOptMsg          oOptResult;
        NumElementsUsed ++;
        Base::Comm.RecvCommand( &nClientPE, &nCommand );
        RUNTIME_ASSERT( nCommand == XDMParallel::REPORT_MC_LIST, "\nDriver::Server: Unknown command from client recv \n" );
        Base::Comm.RecvWorkUnitList( nClientPE, TmpResult );
        RUNTIME_ASSERT( TmpResult.size() == 1, "[This error check shall be removed after debugging] Error: Size mismatch in PMC\n " );
        oOptResult = TmpResult[0];
        if( oOptResult.bConverged )
          oConvergedSamplePoints.push_back( oOptResult );

        if ( oConvergedSamplePoints.size() < nElementsToFit && (NumElementsUsed < nMaxElements) ) // if there's still voxels left
          Utilities::WorkUnitDistribution( VoxelQueue, SingleVoxelDistributor, WaitQueue, Base::Comm, osLogFile, 0,  XDMParallel::FIT_MC );
        else                                     // tell client that we're done
          Base::Comm.SendCommand( nMyID, nClientPE, XDMParallel::WAIT );
      }

      if( oConvergedSamplePoints.size() > nElementsToFit )
      {
        GET_LOG( osLogFile ) << "Needed " << nElementsToFit << " Ended with "
                             <<   oConvergedSamplePoints.size()  << std::endl;
        std::sort( oConvergedSamplePoints.begin(), oConvergedSamplePoints.end()  ); // sort by quality, ascending
        int nExtraPoints = oConvergedSamplePoints.size() - nElementsToFit;
        oConvergedSamplePoints.erase( oConvergedSamplePoints.begin(),
                                      oConvergedSamplePoints.begin() + nExtraPoints ); // keep only the last nElementsToFit
        GET_LOG( osLogFile ) << "List has been trimed to the best " << oConvergedSamplePoints.size() << " elements " << std::endl;
      }
      return  std::make_pair( oConvergedSamplePoints, oConvergedSamplePoints.size() >= nElementsToFit );
    }

    //----------------------------------
    //
    //  LocalParameterOptimization
    //
    //----------------------------------
    template<  class SamplePointT, class SamplePointGrid >
    boost::tuple< std::vector< SDetParamMsg  >, Float, SEnergyOpt >
    ParameterOptimizationServer< SamplePointT, SamplePointGrid >
    ::LocalParamOptimization( const SEnergyOpt            & oEnergyLoc,
                              const vector<SDetParamMsg>  & vDetParams,
                              const vector< SParamOptMsg<SamplePointT> >  & vSamplePoints,
                              const vector<SStepSizeInfo> & vClientStepSizeInfo,
                              const vector<SStepSizeInfo> & vOptMaxDeviation )
    {
      typedef GeometricOptimizationBase<SamplePointT> Base;
      SMonteCarloParam oMCParam;
      oMCParam.fTemperature     = LocalSetup.InputParameters().fParameterMCTemperature;
      oMCParam.nCoolingSteps    = LocalSetup.InputParameters().nNumParamOptSteps * LocalSetup.InputParameters().fCoolingFraction;
      oMCParam.fDTemperature    = oMCParam.fTemperature / oMCParam.nCoolingSteps;
      oMCParam.nThermalizeSteps = LocalSetup.InputParameters().nNumParamOptSteps * LocalSetup.InputParameters().fThermalizeFraction;
      oMCParam.nMaxIterations   = LocalSetup.InputParameters().nNumParamOptSteps;
      oMCParam.eSearchMethod    = static_cast<SO3SearchMethod>( LocalSetup.InputParameters().nOrientationSearchMethod);

      Base::Comm.SendCommand  ( nMyID, 1, nProcessingElements - 1, XDMParallel::EVAL_OVERLAP );
      Base::Comm.BcastSend    ( nMyID,    oMCParam            );
      Base::Comm.BcastSendList( nMyID,    vSamplePoints       );
      Base::Comm.BcastSendList( nMyID,    vClientStepSizeInfo );

      GET_LOG( osLogFile ) << "Begin Local Optimization " << std::endl;
      //-----------------------------------
      //  Send to All  (making this structured would help)
      //-----------------------------------
      for ( Int nClientID = 1; nClientID < nProcessingElements; nClientID ++ )
      {
        vector< SDetParamMsg > vDetectorStartLoc;
        SEnergyOpt oEnergyStartLoc;
        oEnergyStartLoc.fEnergyStep = oEnergyLoc.fEnergyStep;
        if ( nClientID > 1 ) // have one of them start from the origin
        {
          vDetectorStartLoc           = Base::RandomMoveDet( vDetParams, vOptMaxDeviation );
          oEnergyStartLoc.fBeamEnergy = oEnergyLoc.fBeamEnergy
                                      + oRandomReal( -oEnergyLoc.fEnergyStep,
                                                      oEnergyLoc.fEnergyStep );
        }
        else
        {
          vDetectorStartLoc           = vDetParams;
          oEnergyStartLoc.fBeamEnergy = oEnergyLoc.fBeamEnergy;
        }
        Base::SendExpParameters( nClientID, oEnergyStartLoc, vDetectorStartLoc );
      }

      vector<SDetParamMsg>   vBestParam;
      SEnergyOpt             oBestEnergyParam;
      Float fBestCost        = 2;
      Int   nBestID          = nMyID;
      Int   nReported        = 0;
      while( nReported < nProcessingElements - 1 )
      {
        Int nClientID;
        vector< SDetParamMsg > vReturnedList;
        Float fCost;
        Int nCommand;
        SEnergyOpt oNewEnergyParam;
        Base::Comm.RecvCommand( &nClientID, &nCommand );
        Base::Comm.RecvWorkUnit( nClientID, &fCost );
        Base::Comm.RecvWorkUnit( nClientID, &oNewEnergyParam );
        Base::Comm.RecvWorkUnitList( nClientID, vReturnedList );
        if( fBestCost > fCost )
        {
          fBestCost         = fCost;
          nBestID           = nClientID;
          vBestParam        = vReturnedList;
          oBestEnergyParam  = oNewEnergyParam;
        }
        GET_LOG( osLogFile ) << "Optimal from: " << nClientID << " " << fCost <<  std::endl;
        nReported ++;
      }
      return boost::make_tuple( vBestParam, fBestCost, oBestEnergyParam );
    }


    //---------------------------------------------------------------------------
    //  RunParameterOpt
    //
    //  A copy of the voxel list on purpose
    //---------------------------------------------------------------------------
    template< class SamplePointT, class SamplePointGrid >
    void ParameterOptimizationServer<  SamplePointT, SamplePointGrid >::Process( )
    {
      typedef GeometricOptimizationBase<SamplePointT> Base;
      //      typedef ParameterOptimizationServer< Reconstructor, SamplePointT, SamplePointGrid > Self;
      RUNTIME_ASSERT( LocalSetup.InputParameters().fThermalizeFraction >= 0, "ParallelReconstruction:  fThermalizeFraction < 0" );
      RUNTIME_ASSERT( LocalSetup.InputParameters().fCoolingFraction >= 0,    "ParallelReconstruction:  fCoolingFraction < 0" );
      RUNTIME_ASSERT( ( LocalSetup.InputParameters().fThermalizeFraction + LocalSetup.InputParameters().fCoolingFraction <= 1),
                      "ParallelReconstruction: fThermalizeFraction + fCoolingFraction > 1 ");

      vector<SStepSizeInfo>   vGlobalMaxDeviation = LocalSetup.ExperimentalSetup().GetOptimizationInfo();
      vector<SDetParamMsg>      vCurrentDetParams = LocalSetup.ExperimentalSetup().GetExperimentalParameters();
      const vector<SStepSizeInfo> & vMinStepSizes = LocalSetup.ExperimentalSetup().GetDetectorSensitivity();
      //-----------------------------------
      // reduce step size by the cube root of the number or processors
      // (This is the approximate volume taken by each processor)
      //-----------------------------------
      Float fScale = Float ( 1 )
                   / ( pow( (nProcessingElements - 1)  * LocalSetup.InputParameters().nNumParamOptSteps,
                            Float( 1 ) / Float ( 3 ) ) );

      GET_LOG( osLogFile ) << "Scaling Factor for parameter MC: " << fScale << " " << nProcessingElements << std::endl;

      vector<SStepSizeInfo> vClientStepSizeInfo = vGlobalMaxDeviation;
      Base::ScaleParamOptMsg( vClientStepSizeInfo, vMinStepSizes, fScale );

      Int                  nIter           = 0;
      Float                fGlobalBestCost = 2;

      const Float fReductionScale = std::min( LocalSetup.InputParameters().fSearchVolReductionFactor /
                                              pow( Float( nProcessingElements -1 ), Float(1) / Float(3) ),
                                              Float( 0.8 ) );  // 4 is the fudge factor

      vector< SParamOptMsg<SamplePointT> > vSamplePoints;
      VoxelQueue.RandomizeReset();
      SEnergyOpt oEnergyLoc;
      oEnergyLoc.fBeamEnergy = LocalSetup.ExperimentalSetup().GetBeamEnergy();
      oEnergyLoc.fEnergyStep = LocalSetup.InputParameters().BeamEnergyWidth / Float(2) * fScale;

      SEnergyOpt            oBestEnergyLoc     = oEnergyLoc;
      vector<SDetParamMsg>   vBestParams        = vCurrentDetParams;

      Int                   nLocalRestarts     = 0;
      Int                   nLargeStepRestarts = 0;
      bool                  bFitNewVoxels      = true;

      while ( nIter < LocalSetup.InputParameters().nParameterRefinements )
      {
        Bool bFitSuccess = false;
        if( bFitNewVoxels )
        {
          boost::tie(vSamplePoints, bFitSuccess) = GetNConvergedElements( );
          VoxelQueue.RandomizeReset();
          bFitNewVoxels = false;
        }
        Bool bNeedRestart = true;
        if( bFitSuccess || (!bFitNewVoxels) )
        {
          // Run LocalParamOptimization
          SEnergyOpt            OptimizedEnergyLoc;
          vector<SDetParamMsg>   OptimizedParams;
          Float                 fCurrentCost = 2;
          boost::tie( OptimizedParams, fCurrentCost, OptimizedEnergyLoc ) =
            LocalParamOptimization( oEnergyLoc, vCurrentDetParams,
                                    vSamplePoints, vClientStepSizeInfo,
                                    vGlobalMaxDeviation );

          GET_LOG( osLogFile ) << "[Best Cost, New Cost] " << fGlobalBestCost << " " << fCurrentCost << std::endl;
          if( fCurrentCost < fGlobalBestCost )   // if something's found, refine
          {
            fGlobalBestCost       = fCurrentCost;
            vCurrentDetParams     = OptimizedParams;
            vBestParams           = OptimizedParams;
            oBestEnergyLoc        = OptimizedEnergyLoc;
            oEnergyLoc            = OptimizedEnergyLoc;
            GET_LOG( osLogFile ) << "Decided to refine, reduced by: " << fReductionScale << std::endl;
            oEnergyLoc.fEnergyStep *= fReductionScale;
            Base::ScaleParamOptMsg( vGlobalMaxDeviation, vMinStepSizes, fReductionScale );
            Base::ScaleParamOptMsg( vClientStepSizeInfo, vMinStepSizes, fReductionScale );
            bNeedRestart = false;
          }
        }

        if( bNeedRestart )   // is this even useful?
        {
          bFitNewVoxels = true;
          VoxelQueue.RandomizeReset();
          if(  nLocalRestarts < LocalSetup.InputParameters().nMaxParamMCLocalRestarts )
          {
            vCurrentDetParams      = Base::RandomMoveDet( vCurrentDetParams, vClientStepSizeInfo );
            oEnergyLoc.fBeamEnergy = oEnergyLoc.fBeamEnergy
                                   + oRandomReal( -oEnergyLoc.fEnergyStep, oEnergyLoc.fEnergyStep );
            SetClientParameters( oEnergyLoc, vCurrentDetParams );
            nLocalRestarts ++;
          }
          else if ( nLargeStepRestarts < LocalSetup.InputParameters().nMaxParamMCGlobalRestarts )
          {
            nLocalRestarts       = 0;
            Bool bCandidateFound = false;
            while( nLargeStepRestarts < LocalSetup.InputParameters().nMaxParamMCGlobalRestarts && ! bCandidateFound )
            {
              nLargeStepRestarts ++;
              vector< SLargeScaleOpt > oCandidateList =
                GetLargeVolumeCandidates( oEnergyLoc, vCurrentDetParams, vGlobalMaxDeviation );
              if( oCandidateList.size() > 0 )
              {
                bCandidateFound = true;
                std::sort( oCandidateList.begin(), oCandidateList.end() );  // sort ascending order by cost
                vCurrentDetParams      = oCandidateList[0].vDetParams;
                oEnergyLoc.fBeamEnergy = oCandidateList[0].fEnergy;
                SetClientParameters( oEnergyLoc, vCurrentDetParams );
              }
            }
          }
        }
        nIter ++;
      }  //----------------------------------------

      //-----------------------------------
      // save mic file
      //-----------------------------------

      RUNTIME_ASSERT( 0, "Need to rewrite the part that saves to mic");

      VoxelQueue.ClearSolution();
      for( Size_Type i = 0; i < vSamplePoints.size(); i ++ )
        VoxelQueue.Push( vSamplePoints[i].oVoxel );

      //      Base::WriteFitResult( VoxelQueue.Solution(), ".opt");
      Base::SetExperimentalParameters( oBestEnergyLoc.fBeamEnergy, vBestParams );

      GET_LOG( osLogFile ) << "Best Energy: " << oBestEnergyLoc.fBeamEnergy << std::endl;
      Base::Comm.SendCommand( nMyID, 1, nProcessingElements -1, XDMParallel::PROCESS_DONE );
      //      VoxelQueue.ClearSolution();


    }

  } // end namespace ParameterOptimization

}
