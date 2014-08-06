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
  //===============================================================================
  namespace ParameterOptimization
  {
  //---------------------------------------------------------------------
  //  ScaleParamOptMsg
  //---------------------------------------------------------------------
  template< class SamplePointT > 
  void  GeometricOptimizationBase<SamplePointT>
  ::ScaleParamOptMsg( GeometricOptimizationBase<SamplePointT>::SStepSizeInfo & oMsg,
                      const GeometricOptimizationBase<SamplePointT>::SStepSizeInfo & oMinStepSize,
                      Float fScale ) const
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
    //  RandomMoveDet
    //---------------------------------------------------------------------
  template< class SamplePointT > 
  SDetParamMsg GeometricOptimizationBase<SamplePointT>::
  RandomMoveDet( const GeometricOptimizationBase<SamplePointT>::SDetParamMsg  & oCurParam, 
                 const GeometricOptimizationBase<SamplePointT>::SStepSizeInfo & oOptStepSize ) const
  {
    SMatrix3x3 oDelta;
    Float fRadius = -1;
    SO3SearchMethod nParamSO3SearchMethod = static_cast<SO3SearchMethod>( LocalSetup.InputParameters().nOrientationSearchMethod ); 
    if( nParamSO3SearchMethod  == eUniformQuaternion )
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
    typename GeometricOptimizationBase<SamplePointT>::SDetParamMsg oNewParameters = oCurParam;
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
  //  ConstrainedRandomMoveDet
  //  -- moving detector based on constraint set  -- only the first detector is
  //  unconstrained.
  //---------------------------------------------------------------------------
  template< class S >
  typename GeometricOptimizationBase<S>::DetParamMsgList
  GeometricOptimizationBase<S>::
  ConstrainedRandomMoveDet( const GeometricOptimizationBase<S>::DetParamMsgList  & oParamList,
                            const GeometricOptimizationBase<S>::StepSizeList     & oOptStepSizeList ) const
  {
    DetParamMsgList oNewParamList;
    oNewParamList.resize( oParamList.size() );

    vector<Float> vDetectorShifts = LocalSetup.InputParameters().vDetDistSpacing; 

    std::stringstream ss;
    ss << "ConstrainedRandomMoveDet:  Number of L distances does not match the number of difference in Ls \n"
       << " vDetectorShifts.size() = " <<  vDetectorShifts.size()
       << " oParamList.size() = " <<  oParamList.size() << std::endl;
        
    RUNTIME_ASSERT( vDetectorShifts.size() == oParamList.size() -1,
                    ss.str() );

    for ( Size_Type nDet = 0; nDet < oParamList.size(); nDet ++ )
    {
      typename GeometricOptimizationBase<S>::SStepSizeInfo oCurStepSizeInfo = oOptStepSizeList[nDet];
      typename GeometricOptimizationBase<S>::SDetParamMsg  oNewParam = oParamList[nDet];
      if( nDet != 0 )
      {
        Float        fMinAngle = std::min( LocalSetup.InputParameters().fDetOrientDeviationSO3,
                                           oCurStepSizeInfo.fAngularRadius );
        
        Float        fMinStep  = std::min( LocalSetup.InputParameters().fDetDistDeviation,
                                           oCurStepSizeInfo.oDetectorPos.m_fX );
        SVector3     oMinEulerStep;
        SVector3 oDetOrientEulerDeviation = LocalSetup.InputParameters().oDetOrientDeviationEuler;
        oMinEulerStep.m_fX = std::min( oDetOrientEulerDeviation.m_fX, oCurStepSizeInfo.oEulerSteps.m_fX );
        oMinEulerStep.m_fY = std::min( oDetOrientEulerDeviation.m_fY, oCurStepSizeInfo.oEulerSteps.m_fY );
        oMinEulerStep.m_fZ = std::min( oDetOrientEulerDeviation.m_fZ, oCurStepSizeInfo.oEulerSteps.m_fZ );

        oCurStepSizeInfo.oEulerSteps    = oMinEulerStep;
        oCurStepSizeInfo.fAngularRadius = fMinAngle;
        oCurStepSizeInfo.oDetectorPos   = SVector3( fMinStep, 0, 0 );

        vector< SStepSizeInfo >  oOptimizationConstrains  = LocalSetup.ExperimentalSetup().GetOptimizationConstrains();

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
  
  //===============================================================================

  }// end namespace ParameterOptimization
}
