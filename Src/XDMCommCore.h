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
//  XDMParallel.h
//  Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:  Implementation of basic communication wrapper via MPI.
//
////////////////////////////////////////////////////////////

#ifndef _XDM_COMM_CORE_H_
#define _XDM_COMM_CORE_H_

#include <sstream>
#include "Serializer.h"
#include "MicIO.h"
#include <cstdio>
#include <ctime>
#include <random>
#include <tuple>
#include <memory>
#include "Sampling.h"
#include "Symmetry.h"
#include "XDMParallel.h"
#include "Utilities.h"
#include "Reconstructor.h"
#include "ReconstructionSetup.h"
#include <queue>

namespace XDMParallel
{

  namespace core
  {
    //----------------------------------
    //  Experimental Parameter for one detector
    //
    //----------------------------------
    typedef CXDMDetectorFactory::SDetParameters SDetParamMsg;
    typedef CXDMExperimentSetup::SStepSizeInfo  SStepSizeInfo;
    typedef std::map< int, SVoxel > ClientVoxelMapT;
    struct SParamOptMsg
    {
      SVoxel oVoxel;
      CostFunctions::SOverlapInfo oOverlapInfo;
      Bool bConverged;
      
      bool operator< ( const SParamOptMsg & oRHS ) const
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
      vector<SDetParamMsg> vDetParams;
      
      bool operator< ( const SLargeScaleOpt & oRHS ) const
      {
        return fCost < oRHS.fCost;
      }
    };
    
    //-------------------------------------------------------------------
    //  Core structure of Parallel IceNine
    //-------------------------------------------------------------------
    class CParallelDetails
    {
    public:
      typedef Reconstruction::ReconstructionSetup ReconstructionSetup;
      
    protected:
      XDMCommunicator Comm;      
      LogStream osLogFile;
      // CMic      oReconstructedMic;
      double    fTimeBetweenSave;
      Int       nNumParamOptSteps;
      Int       nNumOptElementPerPE;
      Float     fParameterMCTemperature;
      Float     fThermalizeFraction;
      Float     fCoolingFraction;
      Int       nParameterRefinements;
      Float     fEnergyDeviation;
      
      Bool          bConstrainedParamMC;
      Float         fDetDistDeviation;    
      Float         fDetOrientDeviation;     // in radians

      vector< SStepSizeInfo >  oOptimizationConstrains;
      SVector3                 oDetOrientEulerDeviation;
      vector<Float> vDetectorShifts;   // distances between L1 and L2, L2, and L3, and so on  
      
      Int   nMaxParamMCLocalRestarts;
      Int   nMaxParamMCGlobalRestarts;
      Int   nParamMCGlobalSearchElements;
      Float fSearchVolReductionFactor;

      SVector3 oSampleCenter;   // Maybe put these guys in sample?  Need to change CSample?
      Float    fSampleRadius; 
      Float    fMinAccelerationThreshold;
      
      vector<Int>     oCommandList;
      SO3SearchMethod nParamSO3SearchMethod;

      
      CUniformRandomReal oRandomReal;

      ReconstructionSetup oSetup;
      
      //---------------------------------------------------------------------
      //  ScaleParamOptMsg
      //---------------------------------------------------------------------
      void ScaleParamOptMsg( SStepSizeInfo & oMsg,
                             const SStepSizeInfo & oMinStepSize,
                             Float fScale )
      {
        for( Int i = 0; i < 3; i ++ )
        {
          oMsg.oEulerSteps[i]  = std::max( oMsg.oEulerSteps[i]  * fScale, oMinStepSize.oEulerSteps[i]  );
          oMsg.oDetectorPos[i] = std::max( oMsg.oDetectorPos[i] * fScale, oMinStepSize.oDetectorPos[i] );
        }
        oMsg.fBeamCenterJ   = std::max( oMsg.fBeamCenterJ * fScale, oMinStepSize.fBeamCenterJ );
        oMsg.fBeamCenterK   = std::max( oMsg.fBeamCenterK * fScale, oMinStepSize.fBeamCenterK );
        oMsg.fPixelHeight   = std::max( oMsg.fPixelHeight * fScale, oMinStepSize.fPixelHeight );
        oMsg.fPixelWidth    = std::max( oMsg.fPixelWidth  * fScale, oMinStepSize.fPixelWidth  );
        oMsg.fAngularRadius = std::max( oMsg.fAngularRadius * fScale, oMinStepSize.fAngularRadius );
      }
      
      //---------------------------------------------------------------------
      //  ScaleParamOptMsg
      //---------------------------------------------------------------------
      void ScaleParamOptMsg( vector< SStepSizeInfo >       & vMsg,
                             const vector< SStepSizeInfo > & vMinStepSize,
                             Float fScale )
      {
        for( Size_Type i = 0; i < vMsg.size(); i ++ ) 
          ScaleParamOptMsg( vMsg[i], vMinStepSize[i], fScale );
      }
      
      //---------------------------------------------------------------------------
      //  RandomMoveOneDet
      //---------------------------------------------------------------------------
      SDetParamMsg RandomMoveDet( const SDetParamMsg  & oCurParam, 
                                  const SStepSizeInfo & oOptStepSize )
      {
        SMatrix3x3 oDelta;
        Float fRadius = -1; 
        if( nParamSO3SearchMethod == eUniformQuaternion )
        {
          const Float fAngularRadius = oOptStepSize.fAngularRadius;
          UniformGrid::CQuaternionGrid oUniformGridGen;
          fRadius = tan( fAngularRadius ) / sqrt( 48.0 ); // sqrt(48) = 4 * sqrt(3)   ( random start radius )
          Float fX =  oRandomReal( -fRadius, fRadius );
          Float fY =  oRandomReal( -fRadius, fRadius );
          Float fZ =  oRandomReal( -fRadius, fRadius );
          SQuaternion q = oUniformGridGen.GetNearIdentityPoint( fX, fY, fZ );
          oDelta = q.GetRotationMatrix3x3();
        }
        else if ( nParamSO3SearchMethod == eConstraintedEuler )
        {
          Float fX =  oRandomReal( -oOptStepSize.oEulerSteps.m_fX, oOptStepSize.oEulerSteps.m_fX  );
          Float fY =  oRandomReal( -oOptStepSize.oEulerSteps.m_fY, oOptStepSize.oEulerSteps.m_fY );
          Float fZ =  oRandomReal( -oOptStepSize.oEulerSteps.m_fZ, oOptStepSize.oEulerSteps.m_fZ );
          oDelta.BuildActiveSmallRotation( fX / Float( 2 ), fY / Float( 2 ), fZ / Float( 2 ) );
        }
        else
        {
          RUNTIME_ASSERT( 0, "Unknown orientation search method in RandomMove\n");
        }
        SDetParamMsg oNewParameters = oCurParam;
        //
        // WARNING!!!  Only the x-direction is perturbed in the position
        //             of the detector!!
        //             The combination of beam center and detector positions
        //             are redundant.  The *natural* coordinates from the experiment
        //             however, is the beam center + detector distance.  This needs
        //             to resovled somehow.
        oNewParameters.oOrientation    = oDelta * oNewParameters.oOrientation;
        oNewParameters.oPosition.m_fX += oRandomReal( -oOptStepSize.oDetectorPos.m_fX / Float(2), oOptStepSize.oDetectorPos.m_fX / Float(2) );
        oNewParameters.fBeamCenterJ   += oRandomReal( -oOptStepSize.fBeamCenterJ / Float(2), oOptStepSize.fBeamCenterJ / Float(2) );
        oNewParameters.fBeamCenterK   += oRandomReal( -oOptStepSize.fBeamCenterK / Float(2), oOptStepSize.fBeamCenterK / Float(2) );
        oNewParameters.fPixelHeight   += oRandomReal( -oOptStepSize.fPixelHeight / Float(2), oOptStepSize.fPixelHeight / Float(2) );
        oNewParameters.fPixelWidth    += oRandomReal( -oOptStepSize.fPixelWidth  / Float(2), oOptStepSize.fPixelWidth / Float(2) );
    

        return oNewParameters;
      }

      
      //---------------------------------------------------------------------------
      //  DecoupledRandomMoveDet
      //
      //  Randomly move detector based on the step size and reconstructions parameters
      //---------------------------------------------------------------------------
      vector< SDetParamMsg >
      DecoupledRandomMoveDet( const vector<SDetParamMsg>  & oParamList,
                              const vector<SStepSizeInfo> & oOptStepSizeList )
      {
        vector<SDetParamMsg> oNewParamList;
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
      vector< SDetParamMsg >
      ConstrainedRandomMoveDet( const vector<SDetParamMsg>  & oParamList,
                                const vector<SStepSizeInfo> & oOptStepSizeList )
      {
        vector<SDetParamMsg> oNewParamList;
        oNewParamList.resize( oParamList.size() );

        std::stringstream ss;
        ss << "ConstrainedRandomMoveDet:  Number of L distances does not match the number of difference in Ls \n"
           << " vDetectorShifts.size() = " <<  vDetectorShifts.size()
           << " oParamList.size() = " <<  oParamList.size() << std::endl;
        
        RUNTIME_ASSERT( vDetectorShifts.size() == oParamList.size() -1,
                        ss.str() );

        for ( Size_Type nDet = 0; nDet < oParamList.size(); nDet ++ )
        {
          SStepSizeInfo oCurStepSizeInfo = oOptStepSizeList[nDet];
          SDetParamMsg oNewParam = oParamList[nDet];
          if( nDet != 0 )
          {
            Float        fMinAngle = std::min( fDetOrientDeviation, oCurStepSizeInfo.fAngularRadius );
            Float        fMinStep  = std::min( fDetDistDeviation,   oCurStepSizeInfo.oDetectorPos.m_fX );
            SVector3     oMinEulerStep;
            oMinEulerStep.m_fX = std::min( oDetOrientEulerDeviation.m_fX, oCurStepSizeInfo.oEulerSteps.m_fX );
            oMinEulerStep.m_fY = std::min( oDetOrientEulerDeviation.m_fY, oCurStepSizeInfo.oEulerSteps.m_fY );
            oMinEulerStep.m_fZ = std::min( oDetOrientEulerDeviation.m_fZ, oCurStepSizeInfo.oEulerSteps.m_fZ );

            oCurStepSizeInfo.oEulerSteps    = oMinEulerStep;
            oCurStepSizeInfo.fAngularRadius = fMinAngle;
            oCurStepSizeInfo.oDetectorPos   = SVector3( fMinStep, 0, 0 );
            
            oCurStepSizeInfo.fBeamCenterJ = std::min( oOptimizationConstrains[ nDet -  1].fBeamCenterJ,
                                                      oCurStepSizeInfo.fBeamCenterJ ); 
            oCurStepSizeInfo.fBeamCenterK = std::min( oOptimizationConstrains[ nDet - 1 ].fBeamCenterK,
                                                      oCurStepSizeInfo.fBeamCenterK );	    
            oNewParam.oPosition.m_fX        = oNewParamList[ nDet - 1].oPosition.m_fX
                                            + vDetectorShifts[ nDet - 1 ];
            oNewParam.oOrientation          = oNewParamList[ nDet - 1 ].oOrientation;
            oNewParam.fBeamCenterJ          = oNewParamList[ nDet - 1 ].fBeamCenterJ;
            oNewParam.fBeamCenterK          = oNewParamList[ nDet - 1 ].fBeamCenterK;
          }
      
          oNewParamList[ nDet ] = RandomMoveDet( oNewParam, oCurStepSizeInfo );
        }
        return oNewParamList;
      }
  
      //----------------------------------
      //  General RandomMove 
      //----------------------------------
      vector< SDetParamMsg > RandomMoveDet( const vector<SDetParamMsg>  & oParamList,
                                            const vector<SStepSizeInfo> & oOptStepSizeList )
      {
        if ( bConstrainedParamMC )
          return ConstrainedRandomMoveDet( oParamList, oOptStepSizeList );
        else
          return DecoupledRandomMoveDet( oParamList, oOptStepSizeList );
      }

      //----------------------------------
      // SpecifyOptions
      //----------------------------------
      void SpecifyOptions( const CConfigFile & oConfigFile )
      {
        fTimeBetweenSave         = static_cast<double>( oConfigFile.nSecondsBetweenSave );
        nNumParamOptSteps        = oConfigFile.nNumParamOptSteps;
        nNumOptElementPerPE      = oConfigFile.nOptNumElementPerPE;
        fParameterMCTemperature  = oConfigFile.fParameterMCTemperature;
        fThermalizeFraction      = oConfigFile.fThermalizeFraction;
        fCoolingFraction         = oConfigFile.fCoolingFraction;
    
        nParamSO3SearchMethod    = static_cast<SO3SearchMethod>( oConfigFile.nOrientationSearchMethod );
        oCommandList             = oConfigFile.vCommands;

    
        nParameterRefinements    = oConfigFile.nParameterRefinements;
        fEnergyDeviation         = oConfigFile.BeamEnergyWidth;
        
        bConstrainedParamMC      = oConfigFile.bConstrainedParamMC;
        fDetDistDeviation        = oConfigFile.fDetDistDeviation;
    
        fDetOrientDeviation      = oConfigFile.fDetOrientDeviationSO3;     // in radians
        oOptimizationConstrains  = oSetup.ExperimentalSetup().GetOptimizationConstrains();
    
        oDetOrientEulerDeviation = oConfigFile.oDetOrientDeviationEuler;
        vDetectorShifts              = oConfigFile.vDetDistSpacing; 

        nMaxParamMCLocalRestarts     = oConfigFile.nMaxParamMCLocalRestarts;
        nMaxParamMCGlobalRestarts    = oConfigFile.nMaxParamMCGlobalRestarts;
        nParamMCGlobalSearchElements = oConfigFile.nParamMCGlobalSearchElements;
        fSearchVolReductionFactor    = oConfigFile.fSearchVolReductionFactor;
        fMinAccelerationThreshold    = oConfigFile.fMinAccelerationThreshold;
    
        if( bConstrainedParamMC )
          RUNTIME_ASSERT( Int( vDetectorShifts.size() ) == oConfigFile.nDetectors - 1,
                          " ERROR:  Constrained Parameter MC requires N -1 detector shifts specification for N detectors.\n");
      }
      
      CParallelDetails( ){}
    public:

      CParallelDetails( const CConfigFile & oConfigFile )
      {  }
      
      //----------------------------------
      // SetExperimentalParameters 
      //----------------------------------
      void SetExperimentalParameters( const Float fNewEnergy,
                                      const vector<SDetParamMsg> &vNewDetParams )
      {
        for( Size_Type nDetID = 0; nDetID < vNewDetParams.size(); nDetID ++ )
          oSetup.ExperimentalSetup().SetExperimentalParameters( vNewDetParams[ nDetID ], nDetID );
        oSetup.ExperimentalSetup().SetBeamEnergy( fNewEnergy );
      }

      //------------------------------------------------------------------------
      //  Serializer for MPI
      //------------------------------------------------------------------------
      Bool Save   ( CSerializer   & oSerialBuf ) const
      {
        Bool bSuccess;
        bSuccess = oSetup.Save( oSerialBuf );
        return bSuccess;
      }
      
      //------------------------------------------------------------------------
      // Restore
      //------------------------------------------------------------------------
      Bool Restore( CDeserializer & oSerialBuf )
      {
        Bool bSuccess;
        bSuccess = oSetup.Restore( oSerialBuf );
        return bSuccess;
      }
  
      
    }; // ParallelDetails
  } // namespace core


  
  //---------------------------------------------------
  //  SingleElementDistribution
  //---------------------------------------------------
  struct SingleElementDistribution
  {
    template< class Strategy >
    int operator()( XDMCommunicator & Comm, Int nCurrentPE,
                    Strategy & VoxelQueue ) const
    {
      if( ! VoxelQueue.Empty() )   // throw exception otherwise
      {
        Comm.SendWorkUnit( nCurrentPE, VoxelQueue.First() );
        VoxelQueue.Pop();
        return 1;
      }
      return 0;
    }
  };

  //---------------------------------------------------
  //  MultiElementDistribution
  //  Distribute a maximum of nElements
  //---------------------------------------------------
  template <class SamplePointT >
  class MultiElementDistribution
  {
  private:
    MultiElementDistribution();
    int nElements;
  public:
    MultiElementDistribution( int n_ ): nElements( n_ ) {}
    template< class Strategy >
    int operator()( XDMCommunicator & Comm, Int nCurrentPE,
                    Strategy & VoxelQueue ) const
    {
      vector<SamplePointT> WorkUnitList;
      for( int n = 0; n < nElements && ! VoxelQueue.Empty(); n ++ )
      {
        WorkUnitList.push_back( VoxelQueue.First() );
        VoxelQueue.Pop();
      }
      Comm.SendWorkUnitList( nCurrentPE, WorkUnitList );
      return WorkUnitList.size();
    }

    void SetMaxElements( int n ) { nElements = n; }
  };

}

#endif
