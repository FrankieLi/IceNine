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
#include "ExperimentSetup.h"





//-----------------------------------------------------------------------
//
//  Private:
//  
//  GetMaxQ()   
//
//  Return an estimate on the maximum magnitude of scattering vector given our experimental geometry
//
//
//-----------------------------------------------------------------------
Float CXDMExperimentSetup::GetMaxQ( const CDetector & oDetector, const CSample &oSample ) const
{

  //------------------------------
  //
  // By our convention, coordinate origin is either top left or
  // top right corner.  i.e., It is the point that's furthest away
  // from the beam.  MaxQ is defined by how far "up"
  // coordinate origin 
  //
  //
  //------------------------------
  
  SVector3 oKOutMaxDir = oDetector.GetDetectorCoordinateOrigin();    // detector coordinate origin in lab frame
  oKOutMaxDir += oSample.GetLocation( );
  oKOutMaxDir.Normalize();
  SVector3 oQMax = GetReciprocalVector( oKOutMaxDir );
  return oQMax.GetLength();
}

//-----------------------------------------------------------------------
//
//   GetSampleSymmetry
//
//
//-----------------------------------------------------------------------
const CSymmetry* CXDMExperimentSetup::GetSampleSymmetry()    const
{
  switch( oExpConfigFile.SampleSymmetry )
  {
    case LatticeSymmetry::eCubic:
      return & ( LatticeSymmetry::CCubicSymmetry::Get() );
    case LatticeSymmetry::eHexagonal:
      return & ( LatticeSymmetry::CHexagonalSymmetry::Get() );
      break;
    case LatticeSymmetry::eTetragonal:
      return & ( LatticeSymmetry::CTetragonalSymmetry::Get() );
    default:
      stringstream ss;
      ss << "Symmetry " << oExpConfigFile.SampleSymmetry << " is not recognized." << std::endl;
      RUNTIME_ASSERT( 0, ss.str() );
      return NULL;   // this will cause a segfault - this is supposed to be a fatal error
  }
}

//-----------------------------------------------------------------------
//
//  Private:
//  CXDMExperiment::SetDetectionLimit
//  
//  Set limits on reflection (scattering) vector based on experimental geometry  
//
//
//  TODO:  Have an intensity dependent one as well
//
//
//-----------------------------------------------------------------------
void CXDMExperimentSetup::SetDetectionLimit( CUnitCell & oResCellStruct,
                                             const CDetector & oDetector, const CSample &oSample ) const
{

  const vector<SVector3> & oRecpVecs = oResCellStruct.GetReciprocalVectorList();
  RUNTIME_ASSERT( oRecpVecs.size() == 3, 
                  "[CUnitCell::GetReflectionVectorList]: Recipricol Vector not initialied\n");
  Float fMaxQ = GetMaxQ( oDetector, oSample );


  // THIS ONLY WORKS FOR CUBIC --
  std::cout << "WARNING -- QMax is set by Cubic method still!!!!! " << std::endl;
  fMaxQ = min( fMaxQ, oExpConfigFile.fMaxQ );  // take the minimum of user input and computed Q
  Int nMaxH = (Int) ( fMaxQ / oRecpVecs[0].GetLength() );
  Int nMaxK = (Int) ( fMaxQ / oRecpVecs[1].GetLength() );
  Int nMaxL = (Int) ( fMaxQ / oRecpVecs[2].GetLength() ); 

  oResCellStruct.SetReflectionVectorLimits( nMaxH, nMaxK, nMaxL, fMaxQ, 
					    oExpConfigFile.fMinAmplitudeFraction );

}

//-----------------------------------------------------------------------
//  SetConfigFile
//-----------------------------------------------------------------------
void CExperimentSetup::SetConfigFile( const CConfigFile & oConfigFile )
{
  oExpConfigFile = oConfigFile;
  
  oBeamDirection.m_fX = oConfigFile.BeamDirection[0];
  oBeamDirection.m_fY = oConfigFile.BeamDirection[1];
  oBeamDirection.m_fZ = oConfigFile.BeamDirection[2];
  fBeamDeflectionChiLaue = oConfigFile.BeamDeflectionChiLaue;
  
  fMinAcceptedIntensityFraction =  oConfigFile.fMinAmplitudeFraction;
  
  fBeamEnergy = oConfigFile.BeamEnergy;
  fBeamEnergyWidth = oConfigFile.BeamEnergyWidth;
  fEtaLimit = oConfigFile.fEtaLimit;
  
  bInitialized = true;
}

//-----------------------------------------------------------------------
//
//  Public:
//  CExperimentSetup::CExperimentSetup
//  
//  Input:  oConfigFile  -  a properly initialized ExpConfigFile
//
//-----------------------------------------------------------------------
CExperimentSetup::CExperimentSetup( const CConfigFile & oConfigFile ): bInitialized(true)
{
  SetConfigFile( oConfigFile );
}

//-----------------------------------------------------------------------
//
//  ReadRotationInterval
//
//-----------------------------------------------------------------------
Size_Type CXDMExperimentSetup::ReadRotationInterval( Size_Type nNumDetector )
{    
  string oOmegaFile( oExpConfigFile.sRotationRangeFilename );
  
  RUNTIME_ASSERT( InitFileIO::ReadRotationIntervalFiles( vOmegaRangeList, vFileRangeList, nNumDetector, oOmegaFile ),
                  "CXDMExperimentSetup: ERROR:  Simulation failed!! Failed to read omega file\n" );
  
  stringstream ssErr;
  ssErr << "ERROR:  Config file specified  "
        <<  oExpConfigFile.nDetectors
        << " detectors, but "
        << vFileRangeList.size()
        << " file ranges (implying " << vFileRangeList.size() << " detectors) specified in omega file" << std::endl;
  RUNTIME_ASSERT( static_cast<int>( vFileRangeList.size() ) == oExpConfigFile.nDetectors,
                  ssErr.str() );
  
  // will change later.  currently using existing omega file to make sure that everything works.
  for ( Size_Type i = 0; i < vOmegaRangeList.size(); i ++ )
  {
    if ( vOmegaRangeList[i].fHigh < vOmegaRangeList[i].fLow )
    {
      DEBUG_ALERT( 0,  "WARNING:  High and low range are in reversed order!  Swapping them!" );
      std::swap( vOmegaRangeList[i].fHigh, vOmegaRangeList[i].fLow );
    }
    stringstream ss;
    ss << "Simulating ranges [" << vOmegaRangeList[i].fLow
       << ", " << vOmegaRangeList[i].fHigh << "]" << std::endl;
    DEBUG_ALERT( 0, ss.str() );
  }

  RUNTIME_ASSERT( vOmegaRangeList.size() > 0,
                  "ERROR:  Simulation failed!! Omega Range not specified\n");

  //
  //  Note:  We are currently doing multiple wedge structure by artificially
  //         inflating the omega list, and fill most of them with empties.
  Size_Type nLastIndex = vOmegaRangeList.size() - 1;
  Float fHigh;
  Float fLow;
  
  if ( vOmegaRangeList[ 0 ].fHigh > vOmegaRangeList[ nLastIndex ].fHigh ) // flipped
  {
    fHigh =  vOmegaRangeList[ nLastIndex ].fLow;   // flipped on purpose
    fLow  =  vOmegaRangeList[ 0 ].fHigh;
  }
  else
  {
    fHigh =  vOmegaRangeList[ nLastIndex ].fHigh;  // normal
    fLow  =  vOmegaRangeList[ 0 ].fLow;    
  }

  Float fWidth = vOmegaRangeList[0].fHigh - vOmegaRangeList[0].fLow;
  oRangeToIndexMap.Set( fLow, fHigh, fWidth, vOmegaRangeList ); 
  
  return vOmegaRangeList.size();
}

//-----------------------------------------------------------------------
//
//  ReadDetectorInfo
//
//-----------------------------------------------------------------------
Size_Type CXDMExperimentSetup::ReadDetectorInfo( )
{
  InitFileIO::CDetectorFile oDetParser;
  bool bSuccess = oDetParser.Parse( vDetectorList, oExpConfigFile.sDetectorFilename );
  RUNTIME_ASSERT( bSuccess, ( "ERROR:  Simulation failed!! Failed to read detector file\n" + oExpConfigFile.sDetectorFilename).c_str() );
  RUNTIME_ASSERT( vDetectorList.size() > 0,
                  "ERROR:  Simulation failed!! No detector specified\n" );
  return vDetectorList.size();
}

//-----------------------------------------------------------------------
//  ReadStepSizeFile
//-----------------------------------------------------------------------
vector<CXDMExperimentSetup::SStepSizeInfo>
CXDMExperimentSetup::ReadStepSizeFile( const string & sFilename )
{
  vector<SStepSizeInfo> vStepSizeList;
  vector< InitFileIO::CDetectorInfo > vOptInfo;
  InitFileIO::CDetectorFile oDetParser;
  RUNTIME_ASSERT( oDetParser.ParseOptimizationInfo( vOptInfo, sFilename ),
                  "ERROR:  Simulation failed!! Failed to read optimization file\n"  );
  RUNTIME_ASSERT( vOptInfo.size() > 0,
                  "ERROR:  Simulation failed!! No optimization specified\n" );
  
  vStepSizeList.resize( vOptInfo.size() );
  for( Size_Type i = 0; i < vOptInfo.size(); i ++ )
  {
    vStepSizeList[i].oEulerSteps  = vOptInfo[i].vLabFrameOrientation;
    vStepSizeList[i].oDetectorPos = vOptInfo[i].oLabFrameLocation;
    vStepSizeList[i].fBeamCenterJ = vOptInfo[i].fBeamCenterJ;
    vStepSizeList[i].fBeamCenterK = vOptInfo[i].fBeamCenterK;
    vStepSizeList[i].fPixelHeight = vOptInfo[i].fPixelHeight;
    vStepSizeList[i].fPixelWidth  = vOptInfo[i].fPixelWidth;
        
    SMatrix3x3 oEulerStep;
    oEulerStep.BuildActiveEulerMatrix( vStepSizeList[i].oEulerSteps.m_fX,
                                       vStepSizeList[i].oEulerSteps.m_fY,
                                       vStepSizeList[i].oEulerSteps.m_fZ );
    SMatrix3x3 oIdentity;
    oIdentity.SetIdentity();
    //-----------------
    //  This is *DETECTOR* symmetry not sample symmetry.
    //-----------------
    vStepSizeList[i].fAngularRadius
      = GetMisorientation( LatticeSymmetry::CCubicSymmetry::Get(), oIdentity, oEulerStep );
  }
  
  return vStepSizeList;
}

//-----------------------------------------------------------------------
//
//  Initialize Sample
//
//  Note that detector information is needed to properly set the detection limit
//  of the reciprical vectors.  (perhaps changing this would be nice?)
//
//-----------------------------------------------------------------------
void CXDMExperimentSetup::InitializeSample( CSample & oSample, const CDetector & oDetector )
{
  
  oSample.Rotate( oExpConfigFile.SampleOrientation[0],
                  oExpConfigFile.SampleOrientation[1],
                  oExpConfigFile.SampleOrientation[2] );
	
  oSample.Translate( SVector3(oExpConfigFile.SampleLocation[0], 
                              oExpConfigFile.SampleLocation[1], 
                              oExpConfigFile.SampleLocation[2]) );

  // Input sample and structure information
  if ( ! oSample.LoadSample( oExpConfigFile.SampleFilename, oExpConfigFile.nMicGridType ) )
  {
    exit(0);  // Replace with exceptions
  }

  //////////
  //
  // TODO:  Generalize to allow multiple structure
  //
  //

  // the zeroth element is the 'empty' element (do not generate scattering)
  CUnitCell oEmptyStructure;
  oEmptyStructure.InitializeCoordinateSystem();
  oEmptyStructure.SetReflectionVectorLimits(0, 0, 0, 0, 0);
  oSample.AddCrystalStructure( oEmptyStructure );
  
  CUnitCell oCellStructure;
  if ( !InitFileIO::ReadCrystalStructureFile( oCellStructure, oExpConfigFile.StructureFilename ) )
    exit( 0 );

//   oCellStructure.WriteReflectionVectorList( "TestCryStructure.txt", 0  );
//   RUNTIME_ASSERT( 0, " DEBUGGING RUN ONLY - quit at ExperimentSetup.cpp:288\n");
  // Set limits on oSample's scattering vector list
  SetDetectionLimit( oCellStructure, oDetector, oSample );

  const CSymmetry* pSym = GetSampleSymmetry();
  oCellStructure.SetUniqueReflectionVectorList( *pSym );
  oSample.AddCrystalStructure( oCellStructure );
  oSample.SetSampleSymmetry( oExpConfigFile.SampleSymmetry );
}

//-----------------------------------------------------------------------
//
//  InitializeExperiment
//
//  TODO:  Generalize to volume in the future
//-----------------------------------------------------------------------
void CXDMExperimentSetup::InitializeExperiment( )
{
  Size_Type nNumDetector = ReadDetectorInfo( );
  ReadRotationInterval( nNumDetector );
  vOptimizationInfo     = ReadStepSizeFile( oExpConfigFile.sOptimizationFilename );
  vDetectionSensitivity = ReadStepSizeFile( oExpConfigFile.sDetectionLimitFilename );
  vOptimizationConstrains = ReadStepSizeFile( oExpConfigFile.sOptimizationConstrainFilename );
  
  //--------  Error handling ---------------------------------------------
  stringstream ssErr;
  ssErr << "ERROR:  Config file specified  "
        <<  oExpConfigFile.nDetectors
        << " but only "
        << vDetectorList.size()
        << " detectors specified in detector file" << std::endl;
  RUNTIME_ASSERT( static_cast<int>( vDetectorList.size() ) == oExpConfigFile.nDetectors,
                  ssErr.str() );

  ssErr.str("");
  ssErr << "ERROR:  Config file specified  "
        <<  oExpConfigFile.nDetectors
        << " but only "
        << oExpConfigFile.vDetDistSpacing.size()
        << " detector spacing specified.  ( Requirement: nSpacing = nDet -1 )" << std::endl;
  RUNTIME_ASSERT( int( oExpConfigFile.nDetectors - 1 )
                  == int( oExpConfigFile.vDetDistSpacing.size() ),
                  ssErr.str() );
  
  for( Size_Type i = 0; i < oExpConfigFile.vDetDistSpacing.size(); i ++ )
  {
    ssErr.str("");
    ssErr << "Spacing between detector " << i << " "
          <<  i + 1 << " is not specified " << std::endl; 
    RUNTIME_ASSERT( oExpConfigFile.vDetDistSpacing[i] >= 0,
                    ssErr.str() );
  }
  
  RUNTIME_ASSERT( vOptimizationInfo.size() == vDetectorList.size(),
                  "ERROR: Need to specify stepsize for each detector distance in optimization file\n" ); 
 
  RUNTIME_ASSERT( vDetectionSensitivity.size() == vDetectorList.size(),
                  "ERROR: Need to specify stepsize for each detector distance in sensitivity file\n" );
  
  RUNTIME_ASSERT( int( vOptimizationConstrains.size() ) == int( vDetectorList.size() ) -1,
                  "Error:  Number of constrains must equal to NumDetectors -1 in constrain file " );
}

//-----------------------------------------------------------------------
//
//  GetNextSample (questionable if this is needed)
//
//-----------------------------------------------------------------------
const string & CXDMExperimentSetup::GetNextSample() const
{
  DEBUG_ASSERT(bInitialized, "[CXDMExperimentSetup::GetNextSample]: Experiment not initialized\n");
  return oExpConfigFile.SampleFilename;
}


/////////////////////////////////////////////////////////////////////
// 
//  A C C E S S O R S
//
/////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------
//  Public
//  CXDMExperimentSetup::GetBeamEnergy
//
//-----------------------------------------------------------------------
Float CXDMExperimentSetup::GetBeamEnergy() const
{
  DEBUG_ASSERT(bInitialized, "[CXDMExperimentSetup::GetBeamEnergy]: Experiment not initialized\n");
  return fBeamEnergy;
}

//-----------------------------------------------------------------------
//  Public
//  CXDMExperimentSetup::GetEtaLimit
//
//-----------------------------------------------------------------------
Float CXDMExperimentSetup::GetEtaLimit() const
{
  DEBUG_ASSERT(bInitialized, "[CXDMExperimentSetup::GetEtaLimit]: Experiment not initialized\n");
  return fEtaLimit;
}

//-----------------------------------------------------------------------
//  Public:
//  CXDMExperimentSetup::GetBeamEnergyWidth
//
//-----------------------------------------------------------------------
Float CXDMExperimentSetup::GetBeamEnergyWidth() const
{
  DEBUG_ASSERT(bInitialized, "[CXDMExperimentSetup::GetBeamEnergy]: Experiment not initialized\n");
  return fBeamEnergyWidth;
}

//-----------------------------------------------------------------------
//  Public:
//  CXDMExperimentSetup::GetBeamEnergyWidth
//
//-----------------------------------------------------------------------
const SVector3 & CXDMExperimentSetup::GetXrayBeamDirection() const
{
  DEBUG_ASSERT(bInitialized, "[CXDMExperimentSetup::GetXrayBeamDirection]: Experiment not initialized\n");
  return oBeamDirection;
}

//-----------------------------------------------------------------------
//   Public:
//   CXDMExperimentSetup::GetMinAcceptedAmplitude
//
//-----------------------------------------------------------------------
Float CXDMExperimentSetup::GetMinAcceptedIntensityFraction() const
{
  DEBUG_ASSERT(bInitialized, "[CXDMExperimentSetup::GetMinAcceptedIntensity]: Experiment not initialized\n");
  return fMinAcceptedIntensityFraction;
}

//-----------------------------------------------------------------------
//   Public:
//   CXDMExperimentSetup::GetBeamDeflectionChiLaue
//
//-----------------------------------------------------------------------
Float CXDMExperimentSetup::GetBeamDeflectionChiLaue() const
{
  DEBUG_ASSERT(bInitialized, "[CXDMExperimentSetup::GetBeamDeflectionChiLaue]: Experiment not initialized\n");
  return fBeamDeflectionChiLaue;
}

//-----------------------------------------------------------------------
//  GetReciprocalVector
//
//  Purpose:  Given oKOutDir, output direction of the momentum vector,
//            a reciprocal lattice vector G that produced this output
//            direction that is compatible with the experiment is returned.
//            It is assumed that scattering is completely elastic.
//-----------------------------------------------------------------------
SVector3 CXDMExperimentSetup::GetReciprocalVector( const SVector3 & oKOutDir ) const
{
  Float    fWaveNumber = fBeamEnergy * PhysicalConstants::keV_over_hbar_c_in_ang;
  SVector3 oKIn        = this->oBeamDirection * fWaveNumber;         // Incoming momentum vector
  SVector3 oKOut       = oKOutDir * fWaveNumber;                     // Elastic scattering => k_in = k_out
  
  SVector3 oRecipVec   = oKOut - oKIn;
  
  return oRecipVec;
}

//-----------------------------------------------------------------------
//  public:  Save
//
//  Action:  Save all of simulation related variables into CSerializer.
//           Since the order of this matters to the Deserializer, these two
//           functions must be identical.
//
//  Precondition:  This object has already been initialized.  Otherwise, bogus
//                 state information will be saved.
//
//-----------------------------------------------------------------------
Bool CXDMExperimentSetup::Save( CSerializer   & oSerialBuf ) const
{
  
  Bool bSuccess = true;
  bSuccess = oExpConfigFile.Save( oSerialBuf );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactObj( bInitialized );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactObj( oBeamDirection );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactObj( fMinAcceptedIntensityFraction );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactObj( fBeamEnergy );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactObj( fBeamEnergyWidth );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactObj( fBeamDeflectionChiLaue );
  
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactVector( vOmegaRangeList   );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactVector( vFileRangeList    );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactVector( vDetectorList     );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactVector( vOptimizationInfo );
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactVector( vOptimizationConstrains );
  
  bSuccess = bSuccess &&  oSerialBuf.InsertCompactVector( vDetectionSensitivity );
  bSuccess = bSuccess &&  oSerialBuf.InsertComplexObj( oRangeToIndexMap );
  
 return bSuccess;
}

//-----------------------------------------------------------------------
//  public: Restore
//-----------------------------------------------------------------------
Bool CXDMExperimentSetup::Restore( CDeserializer & oSerialBuf )
{
  Bool bSuccess = true;
  bSuccess = oExpConfigFile.Restore( oSerialBuf );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & bInitialized );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & oBeamDirection );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fMinAcceptedIntensityFraction );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fBeamEnergy );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fBeamEnergyWidth );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fBeamDeflectionChiLaue );
  
  // clear everything existing
  vOmegaRangeList.clear();
  vFileRangeList.clear();
  vDetectorList.clear();
  vOptimizationInfo.clear();
  vOptimizationConstrains.clear();
  vDetectionSensitivity.clear();
  
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vOmegaRangeList );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vFileRangeList  );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vDetectorList   );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vOptimizationInfo );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vOptimizationConstrains );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vDetectionSensitivity );
  bSuccess = bSuccess && oSerialBuf.GetComplexObj( & oRangeToIndexMap );
  
  return bSuccess;
}

//--------------------------------------------------------------------
//  SetExperimentalParameters
//
//  Purpose:  Used to change experimental parameters after initialization.
//
//  Note that the position of the parameter in the vector corrosponds to
//  the detector number.
//
//--------------------------------------------------------------------
void CXDMExperimentSetup::SetExperimentalParameters( const CXDMDetectorFactory::SDetParameters & oDetParam,
                                                     Size_Type nDetectorID )
{
  RUNTIME_ASSERT( nDetectorID < vDetectorList.size(), "ERROR:  Unknown Detector ID\n" );

  CXDMDetectorFactory::ModifyImageParameters( vDetectorList[ nDetectorID ],
                                              oDetParam.oPosition,
                                              oDetParam.fBeamCenterJ, oDetParam.fBeamCenterK,
                                              oDetParam.fPixelWidth, oDetParam.fPixelHeight,
                                              oDetParam.oOrientation );
}

//--------------------------------------------------------------------
//  GetExperimentalParameters
//--------------------------------------------------------------------
vector< CXDMDetectorFactory::SDetParameters >
CXDMExperimentSetup::GetExperimentalParameters( ) const
{
  vector< CXDMDetectorFactory::SDetParameters > vDetectorParamList;

  for( Size_Type i = 0; i < vDetectorList.size(); i ++ )
  {
    const CDetector & oCurDetector = vDetectorList[i];
    vDetectorParamList.push_back( CXDMDetectorFactory::GetImageParameters( oCurDetector ) );
  }
  
  return vDetectorParamList;
}

//--------------------------------------------------------------------
//  WriteDetectorFile
//  Purpose:  Output the detector file to the ostream.
//
//--------------------------------------------------------------------
void CXDMExperimentSetup::WriteDetectorFile( ostream & os ) const
{
  for( Size_Type i = 0; i < vDetectorList.size(); i ++ )
  {
    InitFileIO::CDetectorInfo oDetInfo = CXDMDetectorFactory::GetImageInfo( vDetectorList[i] );
    os << "{" << std::endl;
    oDetInfo.Print( os );
    os << "}" << std::endl << std::endl;
  }
}
