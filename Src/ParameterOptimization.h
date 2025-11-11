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
//  ParameterOptimization.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of the server side of the paralle IceNine
// 
//
//
////////////////////////////////////////////////////////////


#ifndef _XDM_PARAMETER_OPTIMIZATION_H_
#define _XDM_PARAMETER_OPTIMIZATION_H_

#include "ParallelReconstructor.h"
#include "DiscreteAdaptive.h"
#include <memory>
#include <queue>
#include "XDMCommCore.h"
#include <limits>
#include "SimulationData.h"
#include <ctime>
#include "XDMParallel.h"
namespace ParallelReconstructor
{

  namespace ParameterOptimization
  {

    typedef CXDMDetectorFactory::SDetParameters SDetParamMsg;
    typedef CXDMExperimentSetup::SStepSizeInfo  SStepSizeInfo;
    
    //------------------------------------------
    //  This replaced CommCore's ParallelDetails
    // 
    //------------------------------------------
    template< class SamplePointT >
    struct SParamOptMsg
    {
      SamplePointT oVoxel;
      CostFunctions::SOverlapInfo oOverlapInfo;
      Bool bConverged;
      bool operator< ( const SParamOptMsg<SamplePointT> & oRHS ) const
      {
        return oOverlapInfo.fQuality < oRHS.oOverlapInfo.fQuality;
      }
    };
    
    struct SMonteCarloParam
    {
      Float fTemperature;
      Float fDTemperature;
      Int   nCoolingSteps;
      Int   nThermalizeSteps;
      Int   nMaxIterations;
      SO3SearchMethod eSearchMethod;
    };
    
    //----------------------------------
    //  SEnergyOpt
    //----------------------------------
    struct SEnergyOpt
    {
      Float fBeamEnergy;
      Float fEnergyStep;
    };
    
    //----------------------------------
    //  LargeScaleOpt
    //----------------------------------
    struct SLargeScaleOpt
    {
      Float        fCost;
      Float        fEnergy;
      vector< SDetParamMsg > vDetParams;
      
      bool operator< ( const SLargeScaleOpt & oRHS ) const
      {
        return fCost < oRHS.fCost;
      }
    };
    
    //-------------------------------------------------------------------
    //  Core structure of Parallel IceNine
    //
    //  Might be able to remove template of SamplePointT
    //-------------------------------------------------------------------
    template< class SamplePointT > 
    class GeometricOptimizationBase
    {
    public:
      typedef Reconstruction::ReconstructionSetup ReconstructionSetup;
      typedef typename CXDMDetectorFactory::SDetParameters SDetParamMsg;
      typedef typename CXDMExperimentSetup::SStepSizeInfo  SStepSizeInfo;
      typedef vector< SDetParamMsg >                 DetParamMsgList;
      typedef vector< SStepSizeInfo >                StepSizeList;
      

      XDMParallel::XDMCommunicator Comm;      
      LogStream          & osLogFile;
      mutable CUniformRandomReal   oRandomReal;
      ReconstructionSetup        & LocalSetup;
      SO3SearchMethod              nParamSO3SearchMethod;

    public:
      //---------------------------------------------------------------------
      //  SendExpParameters
      //---------------------------------------------------------------------
      void SendExpParameters( Int nClientPE,
                              const SEnergyOpt & oEnergyParam,
                              const vector<SDetParamMsg> &vDetectorParams )
      {
        Comm.SendWorkUnit    ( nClientPE, oEnergyParam  );
        Comm.SendWorkUnitList( nClientPE, vDetectorParams );
      }

      //---------------------------------------------------------------------
      //  ScaleParamOptMsg
      //---------------------------------------------------------------------
      void ScaleParamOptMsg( SStepSizeInfo & oMsg,
                             const SStepSizeInfo & oMinStepSize,
                             Float fScale ) const;
      
      //---------------------------------------------------------------------
      //  ScaleParamOptMsg
      //---------------------------------------------------------------------
      void ScaleParamOptMsg( vector< SStepSizeInfo >       & vMsg,
                             const vector< SStepSizeInfo > & vMinStepSize,
                             Float fScale ) const
      {
        for( Size_Type i = 0; i < vMsg.size(); i ++ ) 
          ScaleParamOptMsg( vMsg[i], vMinStepSize[i], fScale );
      }
      
      //---------------------------------------------------------------------------
      //  RandomMoveOneDet
      //---------------------------------------------------------------------------
      virtual SDetParamMsg RandomMoveDet( const SDetParamMsg  & oCurParam, 
                                          const SStepSizeInfo & oOptStepSize ) const;
  
      
      //---------------------------------------------------------------------------
      //  DecoupledRandomMoveDet
      //
      //  Randomly move detector based on the step size and reconstructions parameters
      //---------------------------------------------------------------------------
      DetParamMsgList
      DecoupledRandomMoveDet( const DetParamMsgList  & oParamList,
                              const vector<SStepSizeInfo> & oOptStepSizeList )
      {
        DetParamMsgList oNewParamList;
        oNewParamList.resize( oParamList.size() );
        for ( Size_Type nDet = 0; nDet < oParamList.size(); nDet ++ )
          oNewParamList[ nDet ] =  RandomMoveDet( oParamList[ nDet ], oOptStepSizeList[ nDet ] );
        return oNewParamList;
      }

      //---------------------------------------------------------------------------
      //  ConstrainedRandomMoveDet
      //  -- moving detector based on constraint set  -- only the first detector is
      //  unconstrained.
      //---------------------------------------------------------------------------
      DetParamMsgList
      ConstrainedRandomMoveDet( const DetParamMsgList & oParamList,
                                const vector<SStepSizeInfo> & oOptStepSizeList ) const;

  
      //----------------------------------
      //  General RandomMove 
      //----------------------------------
      DetParamMsgList RandomMoveDet( const DetParamMsgList  & oParamList,
                                     const vector<SStepSizeInfo> & oOptStepSizeList )
      {
        if ( LocalSetup.InputParameters().bConstrainedParamMC )
          return ConstrainedRandomMoveDet( oParamList, oOptStepSizeList );
        else
          return DecoupledRandomMoveDet( oParamList, oOptStepSizeList );
      }

    public:
      GeometricOptimizationBase( ReconstructionSetup & oReconstructionSetup,
                                 LogStream & osLogFile_ )
        : osLogFile( osLogFile_ ),
          LocalSetup( oReconstructionSetup ) {}
      
      //----------------------------------
      // SetExperimentalParameters 
      //----------------------------------
//       void SetExperimentalParameters( const Float fNewEnergy,
//                                       const DetParamMsgList &vNewDetParams )
//       {
//         for( Size_Type nDetID = 0; nDetID < vNewDetParams.size(); nDetID ++ )
//           LocalSetup.ExperimentalSetup().SetExperimentalParameters( vNewDetParams[ nDetID ], nDetID );
//         LocalSetup.ExperimentalSetup().SetBeamEnergy( fNewEnergy );
//       }
        
      //----------------------------------
      // SetExperimentalParameters
      //   ParameterList is usually vector<SDetParamMsg>
      //   But SDetParamMsg could take SamplePointT to be
      //   other than SVoxel.
      //----------------------------------
      template< class ParameterList >
      void SetExperimentalParameters( const Float fNewEnergy,
                                      const ParameterList & vNewDetParams )
      {
        for( Size_Type nDetID = 0; nDetID < vNewDetParams.size(); nDetID ++ )
          LocalSetup.ExperimentalSetup().SetExperimentalParameters( vNewDetParams[ nDetID ], nDetID );
        LocalSetup.ExperimentalSetup().SetBeamEnergy( fNewEnergy );
      }

      virtual void Process() = 0;  // Just to make this an abstract class
    };  // end GeometricOptimizationBase

    //================================================================================================

    
    //--------------------------------------------------------------------------
    //  ParameterOptimization
    //
    //
    //--------------------------------------------------------------------------
    template< class SamplePointT = SVoxel,
              //              class Reconstructor = Reconstruction::DiscreteRefinement<SamplePointT>,
              class SamplePointGrid = MicAnalysis::CMicGrid  >
    class ParameterOptimizationServer :
      public ParallelReconstructorServer,
      public GeometricOptimizationBase< SamplePointT >
     
    {

    public:
      typedef ReconstructionStrategies::RandomizedQueue< SamplePointT,SamplePointGrid > ReconstructionStrategy;
      typedef XDMParallel::SingleElementDistribution                                             SingleElementDistributor;
      typedef MicFile<SamplePointT>                                                              Mic;
      typedef std::shared_ptr<MicIOBase>                                                       MicPtr;
      typedef vector< SDetParamMsg >                  DetParamMsgList;
      typedef GeometricOptimizationBase< SamplePointT >  OptimizationBase;
      
      
    private:
      
      void SetClientParameters( const SEnergyOpt & oEnergyLoc,
                                const vector<SDetParamMsg> & vCurrentDetParams ) 
      {
        for( Int nClientID = 1; nClientID < nProcessingElements; nClientID ++ )  
        {
          OptimizationBase::Comm.SendCommand( 0, nClientID, XDMParallel::SET_EXP_PARAM );
          OptimizationBase::SendExpParameters( nClientID, oEnergyLoc, vCurrentDetParams );
        }
      }
      
    protected:
    
      int nMyID;
      int nProcessingElements;
      ReconstructionSetup & LocalSetup;  // Sever doesn't do anything to this other than to send
      LogStream & osLogFile;
      ReconstructionStrategy VoxelQueue;
      //      ReconstructionStrategy ResultQueue;
      mutable CUniformRandomReal   oRandomReal;
      
  
      //---------------------------------------------------------------------------
      //  GetNConvergedElements
      //  Purpose:  Distribute work unit to all clients, and try to find
      //            n voxels that are fitted, specified by the parameter
      //            optimization.
      //---------------------------------------------------------------------------
      std::pair< std::vector< SParamOptMsg< SamplePointT > >, bool > GetNConvergedElements( );

      //----------------------------------
      //
      // LocalParamOptization
      //
      // Algorithm:  Given M fitted voxels from the origin,
      //             each of the N processing elements will
      //             perturb the experimental parameter, then
      //             locally optimize the orientation of each
      //             voxel.  Note that the orientation can only
      //             drift so far, and therefore bad initial
      //             initial orientation from the voxels fitted
      //             at the origin will lead to unfittable starting
      //             point.
      //
      //  TODO:      Send along a list of candidate list so
      //             that the local perturbation can go a bit
      //             further and becomes just a bit more robust.
      //----------------------------------
      std::tuple< DetParamMsgList, Float, SEnergyOpt >
      LocalParamOptimization( const SEnergyOpt            & oEnergyLoc,
                              const vector<SDetParamMsg>  & vDetParams,
                              const vector< SParamOptMsg<SamplePointT> >  & vSamplePoints,
                              const vector<SStepSizeInfo> & vClientStepSizeInfo,
                              const vector<SStepSizeInfo> & vOptMaxDeviation );

      //----------------------------------
      //  GetLargeVolumeCandidates
      //
      //  This is VERY slow -- use this to home into
      //  closer answers
      //
      //----------------------------------
      vector< SLargeScaleOpt >
      GetLargeVolumeCandidates( const SEnergyOpt            &oStartEnergyLoc,
                                const vector<SDetParamMsg>  &vDetParams,
                                const vector<SStepSizeInfo> &vOptMaxDeviation );

    public:
    
      ParameterOptimizationServer( ReconstructionSetup & oReconstructionSetup,
                                   int nID, int nPE,
                                   LogStream & os )
        : GeometricOptimizationBase< SamplePointT >( oReconstructionSetup, os ),
          nMyID( nID ), nProcessingElements( nPE ),
          LocalSetup( oReconstructionSetup ),
          osLogFile( os ) {}
      
      //--------------------------------------------------
      //  This is the interface layer
      //--------------------------------------------------
      void Initialize( )
      {
        typedef typename LazyBFSServer<SamplePointT, SamplePointGrid>::Mic Mic;
        std::shared_ptr< Mic > pMic= std::dynamic_pointer_cast<Mic>( LocalSetup.ReconstructionRegion() );
        
        VoxelQueue.Initialize( *pMic,
                               LocalSetup.MinSideLength(),
                               LocalSetup.InputParameters().SampleCenter,
                               LocalSetup.InputParameters().SampleRadius );
        // broadcast start
        OptimizationBase::Comm.BcastSendCommand( 0, XDMParallel::OPT_PARAM );
      }
    
      //------------------------------------------
      //  Initialize with explicit set of points that we'd want to
      //  do optimization with.
      //------------------------------------------
      void Initialize( MicPtr StartPoint );
    
      void Process();

      //------------------------------------------
      //  TODO:  return iterator pair
      //------------------------------------------
      MicPtr Solution()
      {
        MicPtr Result = MicPtr( new Mic() );
        *Result = VoxelQueue.Solution();
        return Result;
      }
    };
  

    //===============================================================================


    //--------------------------------------------------------------------------
    //  ParameterOptimizationClient
    //   This should be ready for testing.
    //--------------------------------------------------------------------------
    template< class SamplePointT = SVoxel,
              class Reconstructor = Reconstruction::DiscreteRefinement<SamplePointT>,
              class SamplePointGrid = MicAnalysis::CMicGrid  >
    class ParameterOptimizationClient
      : public ParallelReconstructorClient,
        public GeometricOptimizationBase< SamplePointT >
    {
    public:
      typedef std::shared_ptr< Reconstructor >            ReconstructorPtr;
      typedef ReconstructionStrategies
      ::BreadthFirstStrategy<SamplePointT, SamplePointGrid> ReconstructionStrategy;  // probably don't really need a queue
      
      typedef typename ReconstructionStrategy::SamplePointPtr  SamplePointPtr;
      typedef MicFile<SamplePointT>                            Mic;
      typedef GeometricOptimizationBase< SamplePointT >  OptimizationBase;
      typedef SParamOptMsg<SamplePointT>                 ParamOptMsg;
    private:

      SO3SearchMethod              nParamSO3SearchMethod;
      XDMParallel::XDMCommunicator Comm;      
      LogStream          & osLogFile;
      mutable CUniformRandomReal   oRandomReal;
      ReconstructionSetup        & LocalSetup;

      //      ReconstructionSetup   & LocalSetup;
      CSimulation            Simulator;
      ReconstructionStrategy VoxelQueue;
      ReconstructorPtr       pReconstructor;
      
      Int nMyID;
      const int nServerPE;

      vector<SamplePointPtr> Fit( const SamplePointT & oCenter );
      
      //----------------------------------
      //  Helper functions
      //----------------------------------
    
      //----------------------------------
      //  RecvExpParameters
      //----------------------------------
      void RecvExpParameters( SEnergyOpt & oEnergyParam,
                              vector<SDetParamMsg> &vDetectorParams )
      {
        OptimizationBase::Comm.RecvWorkUnit    ( nServerPE, &oEnergyParam  );
        OptimizationBase::Comm.RecvWorkUnitList( nServerPE, vDetectorParams );
      }
      
      //----------------------------------
      //  EvaluteCost
      //----------------------------------
      Float EvaluateCost( const vector<ParamOptMsg> & vParamInfo )
      {
        if( vParamInfo.size() == 0 )
          return Float(1);
        Float fCost = 0; 
        for( Size_Type i = 0; i < vParamInfo.size(); i++ )
        {
          Float fHitRatio   = CostFunctions::Utilities::GetHitRatio( vParamInfo[i].oOverlapInfo );
          Float fConfidence = CostFunctions::Utilities::GetConfidence( vParamInfo[i].oOverlapInfo );;
          fCost += ( 1.0 - vParamInfo[i].oOverlapInfo.fQuality );
          GET_LOG( osLogFile ) << "Element " << i << " Quality " << vParamInfo[i].oOverlapInfo.fQuality
                               << " Confidence " << fConfidence
                               << " HitRatio " << fHitRatio <<  std::endl;
        }
        GET_LOG( osLogFile ) << "Averaged Cost: " << fCost/ Float( vParamInfo.size() ) << std::endl;
        return fCost / Float( vParamInfo.size() );
      }
    
      void FitMCElements();
      void ParameterMC();
    
    public:
     ParameterOptimizationClient(ReconstructionSetup &oReconstructionSetup,
                                 int nID, LogStream &os)
         : GeometricOptimizationBase<SamplePointT>(oReconstructionSetup, os),
           LocalSetup(oReconstructionSetup),
           nMyID(nID),
           nServerPE(0),
           osLogFile(os) {}  // server PE is always 0

     void Initialize();

     //----------------------------------------------------
     //  Process
     //----------------------------------------------------
     void Process();
      
    };

  }
}
#include "ParameterOptimization.tmpl.cpp"
#include "ParameterOptimizationClient.tmpl.cpp"
#include "ParameterOptimizationServer.tmpl.cpp"
#endif
