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
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL LAB BE LIABLE FOR
// ANY
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
namespace ParallelReconstructor {

namespace ParameterOptimization {
//---------------------------------------------------------------------------
//  ParameterMC
//---------------------------------------------------------------------------
template <class SamplePointT, class Reconstructor, class SamplePointGrid>
void ParameterOptimizationClient<SamplePointT, Reconstructor,
                                 SamplePointGrid>::ParameterMC() {
  typedef ParameterOptimizationClient<SamplePointT, Reconstructor,
                                      SamplePointGrid>
      Self;
  typedef GeometricOptimizationBase<SamplePointT> Base;

  SMonteCarloParam oMCParam;
  typedef SParamOptMsg<SamplePointT> ParamOptMsg;

  vector<ParamOptMsg> vSamplePoints;
  vector<SStepSizeInfo> vStepSizeInfo;
  Comm.BcastRecv(nServerPE, &oMCParam);
  Comm.BcastRecvList(nServerPE, vSamplePoints);
  Comm.BcastRecvList(nServerPE, vStepSizeInfo);
  nParamSO3SearchMethod = oMCParam.eSearchMethod;

  SEnergyOpt oEnergyParam;
  vector<SDetParamMsg> vDetectorParams;
  RecvExpParameters(oEnergyParam, vDetectorParams);
  Float fBestEnergy = oEnergyParam.fBeamEnergy;
  Float fEnergyStep = oEnergyParam.fEnergyStep / Float(2);
  Base::SetExperimentalParameters(fBestEnergy, vDetectorParams);

  //----------------
  //  Evaluate starting point
  //----------------
  for (Size_Type i = 0; i < vSamplePoints.size(); i++)
    vSamplePoints[i].oOverlapInfo =
        pReconstructor->LocalOptimization(vSamplePoints[i].oVoxel);

  Float fBestCost = EvaluateCost(vSamplePoints);

  vector<Float> MaxVisitedL(vDetectorParams.size(),
                            std::numeric_limits<Float>::min());
  vector<Float> MinVisitedL(vDetectorParams.size(),
                            std::numeric_limits<Float>::max());
  vector<Float> MaxVisitedJ(vDetectorParams.size(),
                            std::numeric_limits<Float>::min());
  vector<Float> MinVisitedJ(vDetectorParams.size(),
                            std::numeric_limits<Float>::max());
  vector<Float> MaxVisitedK(vDetectorParams.size(),
                            std::numeric_limits<Float>::min());
  vector<Float> MinVisitedK(vDetectorParams.size(),
                            std::numeric_limits<Float>::max());
  std::pair<Float, Float> CostRange;
  CostRange.first = fBestCost;
  CostRange.second = fBestCost;

  for (Int nIter = 0; nIter < oMCParam.nMaxIterations; nIter++) {
    // perturb parameters
    vector<SDetParamMsg> vNewDetParams =
        Base::RandomMoveDet(vDetectorParams, vStepSizeInfo);
    Float fNewEnergy = fBestEnergy + oRandomReal(-fEnergyStep, fEnergyStep);

    //----------------------------------------------------------
    //  range logging
    for (Int n = 0; n < vDetectorParams.size(); ++n) {
      MaxVisitedL[n] =
          std::max(vDetectorParams[n].oPosition.m_fX, MaxVisitedL[n]);
      MinVisitedL[n] =
          std::min(vDetectorParams[n].oPosition.m_fX, MinVisitedL[n]);
      MaxVisitedJ[n] =
          std::max(vDetectorParams[n].fBeamCenterJ, MaxVisitedJ[n]);
      MinVisitedJ[n] =
          std::min(vDetectorParams[n].fBeamCenterJ, MinVisitedJ[n]);
      MaxVisitedK[n] =
          std::max(vDetectorParams[n].fBeamCenterK, MaxVisitedK[n]);
      MinVisitedK[n] =
          std::min(vDetectorParams[n].fBeamCenterK, MinVisitedK[n]);
    }
    //----------------------------------------------------------

    Base::SetExperimentalParameters(fNewEnergy, vNewDetParams);

    for (Size_Type i = 0; i < vSamplePoints.size(); i++)
      vSamplePoints[i].oOverlapInfo =
          pReconstructor->LocalOptimization(vSamplePoints[i].oVoxel);

    Float fCost = EvaluateCost(vSamplePoints);
    Float fRand = oRandomReal();

    CostRange.first = std::min(fCost, CostRange.first);
    CostRange.second = std::max(fCost, CostRange.second);
    if (fCost < fBestCost) {
      vDetectorParams = vNewDetParams;
      fBestCost = fCost;
      fBestEnergy = fNewEnergy;
    } else if (oMCParam.fTemperature > 0 &&
               exp((fBestCost - fCost) / oMCParam.fTemperature) > fRand) {
      vDetectorParams = vNewDetParams;
      fBestCost = fCost;
      fBestEnergy = fNewEnergy;
    }
    if (nIter > oMCParam.nThermalizeSteps && oMCParam.fTemperature >= 0)
      oMCParam.fTemperature -=
          oMCParam.fDTemperature;  // start cooling after thermaling
  }

  GET_LOG(osLogFile) << "  Best Cost: " << fBestCost << " Total Range [ "
                     << CostRange.first << ", " << CostRange.second << "]"
                     << std::endl;
  GET_LOG(osLogFile) << "  L, j, k Distance Range Traversed: " << std::endl;
  for (Int n = 0; n < vDetectorParams.size(); ++n) {
    GET_LOG(osLogFile) << " Det " << n << " --------------------- "
                       << std::endl;
    GET_LOG(osLogFile) << " [ " << MinVisitedL[n] << ", " << MaxVisitedL[n]
                       << " ]  Delta = " << MaxVisitedL[n] - MinVisitedL[n]
                       << std::endl;
    GET_LOG(osLogFile) << " [ " << MinVisitedJ[n] << ", " << MaxVisitedJ[n]
                       << " ]  Delta = " << MaxVisitedJ[n] - MinVisitedJ[n]
                       << std::endl;
    GET_LOG(osLogFile) << " [ " << MinVisitedK[n] << ", " << MaxVisitedK[n]
                       << " ]  Delta = " << MaxVisitedK[n] - MinVisitedK[n]
                       << std::endl;
  }

  oEnergyParam.fBeamEnergy = fBestEnergy;
  Comm.SendCommand(nMyID, nServerPE, XDMParallel::REPORT_MC);
  Comm.SendWorkUnit(nServerPE, fBestCost);
  Comm.SendWorkUnit(nServerPE, oEnergyParam);
  Comm.SendWorkUnitList(nServerPE, vDetectorParams);
}

//----------------------------------------------------
//  FitMCElementList
//  Fit a list of elements
//----------------------------------------------------
template <class SamplePointT, class Reconstructor, class SamplePointGrid>
void ParameterOptimizationClient<SamplePointT, Reconstructor,
                                 SamplePointGrid>::FitMCElements() {
  typedef ParameterOptimizationClient<SamplePointT, Reconstructor,
                                      SamplePointGrid>
      Self;
  typedef GeometricOptimizationBase<SamplePointT> Base;

  typedef SParamOptMsg<SamplePointT> ParamOptMsg;
  vector<SamplePointT> oWorkUnitList;
  GET_LOG(osLogFile) << "[Client] Recving Work Unit List:" << std::endl;
  Comm.RecvWorkUnitList(nServerPE, oWorkUnitList);
  RUNTIME_ASSERT(
      nServerPE == 0,
      "Client recving messages from units other than ROOT!  STOP \n");
  vector<ParamOptMsg> vOptResult(oWorkUnitList.size());

  GET_LOG(osLogFile) << "[Client] Rec'd Work Unit List:" << oWorkUnitList.size()
                     << " units " << std::endl;

  for (Size_Type i = 0; i < oWorkUnitList.size(); i++) {
    SamplePointT Result;
    Int nResultCode;
    boost::tie(Result, nResultCode) =
        pReconstructor->ReconstructVoxel(oWorkUnitList[i]);
    if (nResultCode == SearchDetails::CONVERGED ||
        nResultCode == SearchDetails::PARTIAL)
      vOptResult[i].bConverged = true;
    else
      vOptResult[i].bConverged = false;
    vOptResult[i].oVoxel = Result;
    vOptResult[i].oOverlapInfo = pReconstructor->EvaluateOverlapInfo(Result);
  }

  GET_LOG(osLogFile) << "[Client] Finished fitting:" << oWorkUnitList.size()
                     << " units " << std::endl;

  Comm.SendCommand(nMyID, nServerPE, XDMParallel::REPORT_MC_LIST);
  Comm.SendWorkUnitList(nServerPE, vOptResult);
  GET_LOG(osLogFile) << " Finished MC fitting element list " << std::endl;
}

//----------------------------------------------------
//  ParameterOptimizationClient::Process
//----------------------------------------------------
template <class SamplePointT, class Reconstructor, class SamplePointGrid>
void ParameterOptimizationClient<SamplePointT, Reconstructor,
                                 SamplePointGrid>::Process() {
  typedef ParameterOptimizationClient<SamplePointT, Reconstructor,
                                      SamplePointGrid>
      Self;
  typedef GeometricOptimizationBase<SamplePointT> Base;

  GET_LOG(osLogFile)
      << "[ParameterOptimization Client]: Beginning Parameter Optimization"
      << std::endl;
  Int nCommand;
  time_t oStartTime, oStopTime;
  do {
    int nSourcePE;
    Comm.RecvCommand(&nSourcePE, &nCommand);
    switch (nCommand) {
      case XDMParallel::FIT_MC: {
        GET_LOG(osLogFile) << "FIT_MC " << std::endl;
        time(&oStartTime);
        FitMCElements();
        time(&oStopTime);
        double oTimeDiff = difftime(oStopTime, oStartTime);
        GET_LOG(osLogFile) << "Time elapsed for MC Fit (sec): " << oTimeDiff
                           << std::endl;
      } break;

      case XDMParallel::EVAL_OVERLAP:
        ParameterMC();
        break;

      case XDMParallel::FIT_MC_LIST:
        RUNTIME_ASSERT(0, "FIT_MC_LIST is no longer available \n");
        exit(0);
        break;

      case XDMParallel::SET_EXP_PARAM: {
        GET_LOG(osLogFile) << "SET_EXP_PARAM " << std::endl;
        SEnergyOpt oEnergyOpt;
        vector<SDetParamMsg> oDetParams;
        RecvExpParameters(oEnergyOpt, oDetParams);
        GET_LOG(osLogFile) << "Request for changes in parameter " << std::endl;
        GET_LOG(osLogFile) << "oEnergyOpt.fBeamEnergy "
                           << oEnergyOpt.fBeamEnergy << std::endl;
        for (Size_Type nDet = 0; nDet < oDetParams.size(); nDet++) {
          GET_LOG(osLogFile) << "--------------------------------" << std::endl;
          GET_LOG(osLogFile)
              << " New Euler     "
              << RadianToDegree(oDetParams[nDet].oOrientation.GetEulerAngles())
              << std::endl;
          GET_LOG(osLogFile) << " New Det Pos   "
                             << oDetParams[nDet].oPosition.m_fX << std::endl;
          GET_LOG(osLogFile) << " Beam Center J "
                             << oDetParams[nDet].fBeamCenterJ << std::endl;
          GET_LOG(osLogFile) << " Beam Center K "
                             << oDetParams[nDet].fBeamCenterK << std::endl;
          GET_LOG(osLogFile) << "--------------------------------" << std::endl;
        }
        Base::SetExperimentalParameters(oEnergyOpt.fBeamEnergy, oDetParams);
      } break;
      case XDMParallel::WAIT:
        GET_LOG(osLogFile) << "Client waiting for server " << std::endl;
        break;

      case XDMParallel::PROCESS_DONE:
        GET_LOG(osLogFile) << "PROCESS_DONE " << std::endl;
        break;

      default:
        RUNTIME_ASSERT(0,
                       "\n Parameter Optimization: Recv-ed unknown command\n");
        break;
    };
  } while (nCommand != XDMParallel::PROCESS_DONE);
}

// template< class SamplePointT, class R, class G >
// void LazyBFSClient<SamplePointT, R, G>::Initialize( )
template <class SamplePointT, class Reconstructor, class SamplePointGrid>
void ParameterOptimizationClient<SamplePointT, Reconstructor,
                                 SamplePointGrid>::Initialize() {
  typedef
      typename LazyBFSClient<SamplePointT, Reconstructor, SamplePointGrid>::Mic
          Mic;
  boost::shared_ptr<Mic> pMic =
      boost::dynamic_pointer_cast<Mic>(this->LocalSetup.ReconstructionRegion());
  VoxelQueue.Initialize(*pMic, LocalSetup.MinSideLength());

  std::cout << "Client:  Finished with initialization" << std::endl;
  // broadcast
  // int nCommand;
  // Comm.BcastRecvCommand(0, &nCommand);
  // RUNTIME_ASSERT(nCommand == BEGIN_BFS_FIT,
  //                "Unexpected command recv'd in LazyBFSClient,
  //                Initialize\n");
  this->Simulator.Initialize(LocalSetup.ExperimentalSetup());
  this->pReconstructor =
      ReconstructorPtr(new Reconstructor(Simulator, LocalSetup));
}

}  // namespace ParameterOptimization

}  // namespace ParallelReconstructor
