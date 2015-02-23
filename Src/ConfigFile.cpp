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
///////////////////////////////////////////////////////
//
//  ConfigFile.cpp
//
//  Implementation of a parser for configuration files.
//  A simple tokenizer is implemented.
//
//
//
//  TODO:   Jun 22, 2009.  Change all atoi and so forth to
//          lexical_cast, with try-catch around it.  This will
//          produce a much more robust error catching code for
//          the configuration file.
//
///////////////////////////////////////////////////////

#include "ConfigFile.h"

using std::cout;
using std::cerr;
using std::endl;

//-------------------------------------------------
//
//  Private
//  CConfigFile::Parse - simple parser from a character buffer
//
//-------------------------------------------------
bool CConfigFile::Parse( const string & sBuf )
{
  const char * pSO3SearchTags[] = { "ConstrainedEuler", "UniformSO3" };
  const char * pGridTypeSearchTag[] = { "Triangular", "Square" };
  
  RUNTIME_ASSERT( sizeof(pSO3SearchTags)/sizeof(char*)  == eNumSO3SearchTypes,
                  "[ConfigFile] ERROR: Number of keywords for SO3 Search methods is not up to date \n  " );
  
  RUNTIME_ASSERT( sizeof(pGridTypeSearchTag)/sizeof(char*)  == eNumGridTypes,
                  "[ConfigFile] ERROR: Number of keywords for NumGridTypesnot up to date \n  " );
  const char* sKeywordStrings[] = {
    // comment
    "#",  

    // infile information
    "InfileBasename", "InfileExtension", "FileNumStart",
    "FileNumEnd", "InFileType", "InfileSerialLength",
		
    // Outfile information
    "OutfileBasename", "OutfileExtension",
    "OutfileSerialLength", "OutStructureBasename",

    // Beam Information
    "BeamEnergy", "BeamEnergyWidth", "BeamHeight", 
    "BeamDirection",
    "BeamDeflectionChiLaue",
        
    // Detector Information
    "DetectorFilename", "EtaLimit",

    // Sample Information
    "SampleLocation", "SampleRadius", "SampleCenter", "SampleOrientation", "SampleFilename",
    "StructureFilename", "FundamentalZoneFilename", "SampleSymmetry",

    "MaxInitSideLength",
    "MinSideLength",
    
    // Simulation information
    "MinAmplitudeFraction", "MaxQ",

    // Initialization files
    "RotationRangeFilename",

    // Search algorithm parameters
    "LocalOrientationGridRadius",
    "MinLocalResolution",
    "MaxLocalResolution",
    "MaxDiscreteCandidates",
    "MaxAcceptedCost",
    "MaxConvergenceCost",
    "MaxDeepeningHitRatio",
    "MaxMCSteps",
    "MCRadiusScaleFactor",
    "SuccessiveRestarts",

    //  Acceleration parameters
    "MinAccelerationThreshold",
    
    // Saving options
    "SecondsBetweenSave",

    // Parameter Optimization
    "OptimizationFilename",
    "OptimizationConstrainFilename",
    "DetectionLimitFilename",
    "NumParameterOptimizationSteps",
    "NumElementToOptimizePerPE",
    "ParameterMCInitTemperature",
    "OrientationSearchMethod",
    "CoolingFraction",
    "ThermalizeFraction",
    "ParameterRefinements",

    "NumDetectors",
    "DetectorSpacing",
    "DetectorSpacingDeviation",
    "DetectorOrientationDeviationInEuler",
    "DetectorOrientationDeviationInSO3",
    
    "ParamMCMaxLocalRestarts",
    "ParamMCMaxGlobalRestarts",
    "ParamMCNumGlobalSearchElements",
    "ConstrainedOptimization",
    "SearchVolumeReductionFactor",

    //  Scattering Vector Based Initial Guesses
    "ConsistencyError",
    "BraggFilterTolerance",
    
    //  Fit type
    "DEBUG_FirstSearchType",
    "RunParameterOptimization",
    "RunReconstruction",
    "RunAdpReconstruction",
    "SelectBoundaryVoxels",
    "IntensityDecomposition",
    "LazyBFS",
    "LazyStrain",
    "LocalOrientationOptimization",
    "DEBUG_LastSearchType",
    //  end search types
    //-----------------------------
    
    "BCPeakDetectorOffset",

    //------------------------
    //  Strain related options
    //------------------------
    "EnableStrain",
    "StrainOptConfigFilename",
    
    //------------------------
    //  Partial Result Loading
    //------------------------
    "PartialResultFilename",
    "PartialResultAcceptanceConfidence",

    //------------------------
    "GridType"
  };

  enum EKeyword{
    // comment
    eComment, 
		
    // Infile information
    eInfileBasename, eInFileExtension, eFileNumStart,
    eFileNumEnd, eInFileType, eInFileSerialLength,
		
    // Outfile information
    eOutfileBasename, eOutfileExtension,
    eOutfileSerialLength, eOutStructureBasename,
    
    // Beam Information
    eBeamEnergy, eBeamEnergyWidth, eBeamHeight,
    eBeamDirection,
    eBeamDeflectionChiLaue,
    
    // Detector Information
    eDetectorFilename, eEtaLimit,

    // Sample Information
    eSampleLocation, eSampleRadius, eSampleCenter, eSampleOrientation, 
    eSampleFilename, 

    eStructureFilename, eFundamentalZoneFilename, eSampleSymmetry,

    eMaxInitSideLength,
    eMinSideLength,

    // Simulation information
    eMinAmplitudeFraction,
    eMaxQ,

    // Initialization files
    eRotationRangeFilename,
    
    // Search algorithm parameters
    eLocalOrientationGridRadius,
    eMinLocalResolution,
    eMaxLocalResolution,
    eMaxDiscreteCandidates,
    eMaxAcceptedCost,
    eMaxConvergenceCost,
    eMaxDeepeningHitRatio,
    eMaxMCSteps,
    eMCRadiusScaleFactor,
    eSuccessiveRestarts,

    eMinAccelerationThreshold,
    
    // file saving options
    eSecondsBetweenSave,

    //  Optimization parameters
    eOptimizationFilename,
    eOptimizationConstrainFilename,
    eDetectionLimitFilename,
    eNumParameterOptimizationSteps,
    eNumElementToOptimizePerPE,
    eParameterMCInitTemperature,
    eOrientationSearchMethod,

    eCoolingFraction,
    eThermalizeFraction,
    eParameterRefinements,

    eNumDetectors,
    eDetectorSpacing,
    eDetectorSpacingDeviation,
    eDetectorOrientationDeviationInEuler,
    eDetectorOrientationDeviationInSO3,
    
    eParamMCMaxLocalRestarts,
    eParamMCMaxGlobalRestarts,
    eParamMCNumGlobalSearchElements,
    eConstrainedOptimization,
    eSearchVolumeReductionFactor,

    eConsistencyError,
    eBraggFilterTolerance,
    
    //  search types
    eDEBUG_FirstSearchType,
    eRunParameterOptimization,
    eRunReconstruction,
    eRunAdpReconstruction,
    eSelectBoundaryVoxels,
    
    eIntensityDecomposition,
    eLazyBFS,
    eLazyStrain,
    eLocalOrientationOptimization,

    eDEBUG_LastSearchType,
    //  end search type

    eBCPeakDetectorOffset,

    eEnableStrain,
    eStrainOptConfigFilename,
    ePartialResultFilename,
    ePartialResultAcceptanceConfidence,

    eGridType,
    
    // error check
    eNumKeywords
  };
	
  const char* sFileTypes[] = {"bin", "tif", "ascii"};

  bool *vInitializationCheck;
  bool *vRequirementCheck;
  vInitializationCheck = new bool[eNumKeywords];   // a check list to see which variable is not initialized
  vRequirementCheck    = new bool[eNumKeywords];   // a check list to see if requirement is met

  for ( Int i = 0; i < eNumKeywords; i ++ )
  {
    vInitializationCheck[i] = false;
    vRequirementCheck[i]    = false;
  }
  //------------------------------------------------------------------
  //
  //  OPTIONAL KEYWORDS 
  //  check list must be set to true for optional keywords
  //
  // a list of variables that are optional
  vInitializationCheck[eComment]                  = true;
  vInitializationCheck[eDEBUG_FirstSearchType]    = true;
  vInitializationCheck[eDEBUG_LastSearchType]     = true;
  vInitializationCheck[eRunParameterOptimization] = true;
  vInitializationCheck[eRunReconstruction]        = true;
  vInitializationCheck[eRunAdpReconstruction]     = true;
  vInitializationCheck[eSelectBoundaryVoxels]     = true;
  vInitializationCheck[eIntensityDecomposition]   = true;
  vInitializationCheck[eLazyBFS]                  = true;
  vInitializationCheck[eLazyStrain]               = true;
  vInitializationCheck[eLocalOrientationOptimization] = true;
  vInitializationCheck[ePartialResultFilename] = true;
  vInitializationCheck[ePartialResultAcceptanceConfidence] = true;
  vInitializationCheck[eEnableStrain] = true;
  vInitializationCheck[eStrainOptConfigFilename] = true;
  vInitializationCheck[eGridType]                = true;   // Default to Triangular
  nMicGridType = eTriangular;                     //  Default grid type
  //------------------------------------------------------------------

  RUNTIME_ASSERT( sizeof(sKeywordStrings)/sizeof(char*) == eNumKeywords,
                  "[ConfigFile] ERROR:  Keyword array is not up to date" );

  RUNTIME_ASSERT( sizeof(sFileTypes)/sizeof(char*) == eNumFileTypes, 
                  "[ConfigFile] ERROR:  File type Keyword array is not up to date" );
	
  vector<vector<string> > vsTokens;
  GeneralLib::Tokenize( vsTokens, sBuf, " \t\n");

  for(Size_Type i = 0; i < vsTokens.size(); i ++)
  {
    Size_Type iFirstToken;

    if ( vsTokens[i].size() == 0)
    {
      iFirstToken = eComment;
    }
    else if ( vsTokens[i][0].find_first_of( sKeywordStrings[ eComment ] ) == 0) // if first token is comment
    {
      iFirstToken = eComment;
    }
    else
    {
      // Identify the keyword in the beginning of the line (i.e., vertex? texture?
      for(iFirstToken = 0; iFirstToken < eNumKeywords; iFirstToken ++)
      {
        if(strcmp(sKeywordStrings[iFirstToken], vsTokens[i][0].c_str()) == 0)
          break;
      }
      CONFIG_DEBUG( std::cout << i + 1 << " " << sKeywordStrings[iFirstToken] << " "
                    << vsTokens[i][0].c_str() << " Token " << iFirstToken << endl); 
    }

    //----------------
    //  Error handling
    //----------------
    if( (iFirstToken < eDEBUG_FirstSearchType 
	|| iFirstToken > eDEBUG_LastSearchType ) && iFirstToken != eComment )
    {
      if( vsTokens[i].size() <= 1 )
      {
	std::cerr << "Keyword " << vsTokens[i][0] 
		  << " requires at least one argument. Please check the config file at around line: "
		  << i << endl;
	RUNTIME_ASSERT( false, " Config File has a typo somewhere!" );
      }
    }
   
    switch(iFirstToken)   // look at first token of each line
    {
      // Comment
      case eComment:  
        break;
			
        ////
        //  Infile Information
        //
        ////

      case eInfileBasename:
        InFileBasename = vsTokens[i][1];
        break;
      case eFileNumStart:
        FileStartNum = atoi(vsTokens[i][1].c_str());
        break;
      case eFileNumEnd:
        FileEndNum = atoi(vsTokens[i][1].c_str());
        break;
      case eInFileType:
        InFileType = eNumFileTypes;
        for(UInt j = 0; j < eNumFileTypes; j++){
          if (sFileTypes[j] == vsTokens[i][1])
            InFileType = (EFileType)j;
        }
        if (InFileType == eNumFileTypes){
          cerr << "[ConfigFile] Error:  Unrecognized file type" << endl;
          exit(0);
        }else{
          break;			
        }
      case eInFileSerialLength:
        InFileSerialLength = atoi( vsTokens[i][1].c_str() );
        break;
        ////
        //  Outfile Information
        ////
      case eOutfileBasename:
        OutFileBasename = vsTokens[i][1];
        break;
      case eOutfileExtension:
        OutFileExt = vsTokens[i][1];
        break;
      case eOutStructureBasename:
        OutStructureBasename = vsTokens[i][1];
        break;
      case eOutfileSerialLength:
        OutFileSerialLength = atoi( vsTokens[i][1].c_str() );
        break;
      case eInFileExtension:
        InFileExt = vsTokens[i][1];
        break;

        ////
        //  Beam Information
        ////
      case eBeamEnergy:
        BeamEnergy = atof(vsTokens[i][1].c_str());
        break;
      case eBeamEnergyWidth:
        BeamEnergyWidth = atof(vsTokens[i][1].c_str());
        break;
      case eBeamHeight:
        BeamHeight = atof(vsTokens[i][1].c_str());
        break;
      case eBeamDirection:
        BeamDirection = InitFileIO::ExtractVector( vsTokens[i], i );
        break;
      case eBeamDeflectionChiLaue:
        BeamDeflectionChiLaue = atof( vsTokens[i][1].c_str() );
        break;        
        ////
        //  Detector Information
        ////
      case eDetectorFilename:
        sDetectorFilename = vsTokens[i][1];
        break;
      case eEtaLimit:
        fEtaLimit = DEGREE_TO_RADIAN( atof( vsTokens[i][1].c_str() ) );
        break;
        ////
        //  Sample Position
        ////
      case eSampleLocation:
        SampleLocation = InitFileIO::ExtractVector( vsTokens[i], i );
        break;
      case eSampleRadius:
        SampleRadius = InitFileIO::ExtractReal( vsTokens[i], i );
        break;
      case eSampleCenter:
        SampleCenter = InitFileIO::ExtractVector( vsTokens[i], i );
        break;
      case eSampleOrientation:
        SampleOrientation = InitFileIO::ExtractVector( vsTokens[i], i );
        break;			
      case eSampleFilename:
        SampleFilename = vsTokens[i][1];
        break;
      case eStructureFilename:
        StructureFilename = vsTokens[i][1];
        break;
      case eFundamentalZoneFilename:
        FundamentalZoneFilename = vsTokens[i][1];
        break;
      case eSampleSymmetry:
        SampleSymmetry = LatticeSymmetry::ParseSymmetry( vsTokens[i][1] );
        if( SampleSymmetry == LatticeSymmetry::eNumSymmetry )
        {
          cerr << "[ConfigFile] Error: Line " << i  << "  Unknown symmetry " << endl;
          exit(0);
        }
        break;
      case eMaxInitSideLength:
        fMaxInitSideLength = atof( vsTokens[i][1].c_str() );
        break;
      case eMinSideLength:
        fMinSideLength = atof( vsTokens[i][1].c_str() );
        break;
      case eMinAmplitudeFraction:
        fMinAmplitudeFraction = atof( vsTokens[i][1].c_str() );
        if( fMinAmplitudeFraction > 1 || fMinAmplitudeFraction < 0){
          cerr << "[ConfigFile] Error: Line " << i  << "  fMinAmplitude must be between [0, 1]" << endl;
          exit(0);
        }
        
        break;
      case eMaxQ:
        fMaxQ = atof( vsTokens[i][1].c_str() );
        if( fMaxQ < 0){
          cerr << "[ConfigFile] Error: Line " << i  << "  MaxQ must be > 0" << endl;
          exit(0);
        }
        break;
      case eRotationRangeFilename:
        sRotationRangeFilename = vsTokens[i][1];
        break;

        //-----------------------------
        // Search algorithm parameters
        //-----------------------------
      case eLocalOrientationGridRadius:
        fLocalOrientationGridRadius = DEGREE_TO_RADIAN( atof( vsTokens[i][1].c_str() ) );
        break;
      case eMinLocalResolution:
        nMinLocalResolution = atoi( vsTokens[i][1].c_str() );
        break;
      case eMaxLocalResolution:
        nMaxLocalResolution = atoi( vsTokens[i][1].c_str() );
        break;
      case eMaxDiscreteCandidates:
        nMaxDiscreteCandidates = atoi( vsTokens[i][1].c_str() );
        break;
      case eMaxAcceptedCost:
        fMaxAcceptedCost = atof( vsTokens[i][1].c_str() );
        break;
      case eMaxConvergenceCost:
        fMaxConvergenceCost = atof( vsTokens[i][1].c_str() );
        break;
      case eMaxDeepeningHitRatio:
        fMaxDeepeningHitRatio = atof( vsTokens[i][1].c_str() );
        break;
      case eMaxMCSteps:
        nMaxMCSteps = atoi( vsTokens[i][1].c_str() );
        break;
      case eMCRadiusScaleFactor:
        fMCRadiusScaleFactor = atof( vsTokens[i][1].c_str() );
        break;
      case eSuccessiveRestarts:
        nSuccessiveRestarts =  atoi( vsTokens[i][1].c_str() );
        break;
      case eMinAccelerationThreshold:
        fMinAccelerationThreshold =  atof( vsTokens[i][1].c_str() );
        break;
      case eSecondsBetweenSave:
        nSecondsBetweenSave =  atoi( vsTokens[i][1].c_str() );
        break;
      case eOptimizationFilename:
        sOptimizationFilename =  vsTokens[i][1];
        break;
      case eOptimizationConstrainFilename:
        sOptimizationConstrainFilename = vsTokens[i][1];
        break;
      case eDetectionLimitFilename:
        sDetectionLimitFilename =  vsTokens[i][1];
        break;
      case eNumParameterOptimizationSteps:
        nNumParamOptSteps = atoi( vsTokens[i][1].c_str() );
        break;
      case eNumElementToOptimizePerPE:
        nOptNumElementPerPE = atoi( vsTokens[i][1].c_str() );
        break;
      case eParameterMCInitTemperature:
        fParameterMCTemperature = atof( vsTokens[i][1].c_str() );
        break;
      case eOrientationSearchMethod:

        for( Int nTagID = 0; nTagID < eNumSO3SearchTypes; nTagID ++ )
          if( strcmp( pSO3SearchTags[ nTagID ], vsTokens[i][1].c_str() ) == 0 )
            nOrientationSearchMethod = nTagID;

        if( nOrientationSearchMethod == eNumSO3SearchTypes )
        {
          cerr << "Unexpected search method for orientation search of detector parameter at line: \n" << i << endl; 
          exit( 0 );
          return false;
        }
        break;
      case eGridType:

        for( Int nTagID = 0; nTagID < eNumGridTypes; nTagID ++ )
          if( strcmp( pGridTypeSearchTag[ nTagID ], vsTokens[i][1].c_str() ) == 0 )
            nMicGridType= nTagID;

        if( nMicGridType ==  eNumGridTypes)
        {
          cerr << "Unexpected grid type detected at line: \n" << i << endl; 
          exit( 0 );
          return false;
        }
        std::cout << "Mic grid type selected: " << pGridTypeSearchTag[ nMicGridType ] << std::endl;
        break;        
      case eCoolingFraction:
        fCoolingFraction = atof( vsTokens[i][1].c_str() );
        break;
      case eThermalizeFraction:
        fThermalizeFraction = atof( vsTokens[i][1].c_str() );
        break;
      case eParameterRefinements:
        nParameterRefinements = atoi( vsTokens[i][1].c_str() );
        break;

      case eNumDetectors:
        nDetectors =  InitFileIO::ExtractInt( vsTokens[i], i );
        RUNTIME_ASSERT( nDetectors > 0, "Number of detectors defined must be greater than 0\n");
        vDetDistSpacing.resize( nDetectors -1 );
        std::fill( vDetDistSpacing.begin(), vDetDistSpacing.end(), -1 );
        break;
      case eDetectorSpacing:
        {
          RUNTIME_ASSERT( nDetectors > 0,
                          " Number of detectors is not yet defined.  It needs to be defined before this line\n");
          Int nCurDet = atoi( vsTokens[i][1].c_str() );
          RUNTIME_ASSERT( nCurDet < int( nDetectors ) -1,
                          "Either not enough detectors specified, or there is an attempt to define spacing between a non-existing detector.\n Note: We count from 0.\n");
          vDetDistSpacing[ nCurDet ] = atof( vsTokens[i][2].c_str() );
        }
        break;
      case eDetectorSpacingDeviation:
        fDetDistDeviation = InitFileIO::ExtractReal( vsTokens[i], i );
        break;
      case eDetectorOrientationDeviationInEuler:
        oDetOrientDeviationEuler = DegreeToRadian( InitFileIO::ExtractVector( vsTokens[i], i ) );
        break;
      case eDetectorOrientationDeviationInSO3:
        fDetOrientDeviationSO3 = DEGREE_TO_RADIAN( InitFileIO::ExtractReal( vsTokens[i], i ) );
        break;
      case eParamMCMaxLocalRestarts:
        nMaxParamMCLocalRestarts = InitFileIO::ExtractInt( vsTokens[i], i );
        break;
      case eParamMCMaxGlobalRestarts:
        nMaxParamMCGlobalRestarts = InitFileIO::ExtractInt( vsTokens[i], i );
        break;
      case eParamMCNumGlobalSearchElements:
        nParamMCGlobalSearchElements = InitFileIO::ExtractInt( vsTokens[i], i );
        break;
      case eConstrainedOptimization:
        {
          Int nChoice  = InitFileIO::ExtractInt( vsTokens[i], i );
          RUNTIME_ASSERT( nChoice == 1 || nChoice == 0,
                          "Error:  Unknown choice for ConstrainedOptimization");
          if( nChoice == 1 )
            bConstrainedParamMC = true;
          else
            bConstrainedParamMC = false;
          break;
        }
      case eSearchVolumeReductionFactor:
        fSearchVolReductionFactor = InitFileIO::ExtractReal( vsTokens[i], i );
        break;

      case eConsistencyError:
        fConsistencyError = DEGREE_TO_RADIAN( InitFileIO::ExtractReal( vsTokens[i], i ) );
        break;
      case eBraggFilterTolerance:
        fBraggFilterTol   = DEGREE_TO_RADIAN( InitFileIO::ExtractReal( vsTokens[i], i ) );
        break;
      case eRunParameterOptimization:
        vCommands.push_back( XDMParallel::OPT_PARAM );
        break;
      case eRunReconstruction:
        vCommands.push_back( XDMParallel::FIT );
        break;
      case eRunAdpReconstruction:
        vCommands.push_back( XDMParallel::FIT_ADP );
        break;
      case eSelectBoundaryVoxels:
        RUNTIME_ASSERT( vsTokens[i].size() == 4,
                        "Usage: SelectBoundaryVoxels MaxCost AngleThreshold(degrees) Radius (mm)\n" );

        fMaxCost           = atof( vsTokens[i][1].c_str() );
        fMaxMisorientation = DEGREE_TO_RADIAN( atof( vsTokens[i][2].c_str() ) );
        fMaxRadius         = atof( vsTokens[i][3].c_str() );

        vCommands.push_back( XDMParallel::USE_BND );
        break;

      case eIntensityDecomposition:
        vCommands.push_back( XDMParallel::FIT_INT_DECOMP );
      case eLazyBFS:
        vCommands.push_back( XDMParallel::FIT_LBFS );
        break;
      case eLocalOrientationOptimization:
        vCommands.push_back( XDMParallel::FIT_LOCAL_OPTIMIZE );
        break;
      case eLazyStrain:
        vCommands.push_back( XDMParallel::FIT_STRAIN_LBFS );
        break;        
      case eBCPeakDetectorOffset:
        nBCPeakDetectorOffset = atoi( vsTokens[i][1].c_str() );
        break;

      case eEnableStrain:
        {
          RUNTIME_ASSERT( vsTokens[i].size() == 2,
                          "Usage:  EnableStrain  [ 1 | 0]\n" );
          if( atoi( vsTokens[i][1].c_str() ) == 1 )
            bStrainEnabled = true;
          else
            bStrainEnabled = false;
          break;  
        }
      case eStrainOptConfigFilename:
        {
          bStrainEnabled = true;  
          RUNTIME_ASSERT( vsTokens[i].size() == 2,
                          "Usage: Missing Strain Optimization Config Filename \n" );
          string StringOptFilename = vsTokens[i][1];
          RUNTIME_ASSERT( StrainOptSetup.Read( StringOptFilename ), "Failed to read StrainOptConfigFile" );;
          break;
        }
      case ePartialResultFilename:
        RUNTIME_ASSERT( vsTokens[i].size() == 2,
                        "Usage: Missing path for PartialResultFilename \n" );
        bUsePartialResult     = true;
        PartialResultFilename = vsTokens[i][1];
        break;

      case ePartialResultAcceptanceConfidence: 
        RUNTIME_ASSERT( vsTokens[i].size() == 2,
                        "Usage: Missing min accptance restart confidence \n" );
        bUsePartialResult     = true;
        std::cout << "Partial Result Acceptance " << vsTokens[i][0] << " " << vsTokens[i][1] << std::endl;
        fPartialResultAcceptanceConf = InitFileIO::ExtractReal( vsTokens[i], i );
        break;
        
      default:
        {
          cerr << "[ConfigFile] Error: syntax not recognized:  Line " << i  << endl;
          exit(0);
          return false;
        }
        // insert error stuff here
        
    }// end switch
    
    // set check off this particular token
    vRequirementCheck[ iFirstToken ]    = true;
    vInitializationCheck[ iFirstToken ] = true;

    
    if( bUsePartialResult &&
        ( !    vInitializationCheck[ePartialResultFilename] 
          || ! vInitializationCheck[ePartialResultAcceptanceConfidence]
          )
        ) 
    {
      cerr << "[ConfigFile] Error!" << std::endl
           << "Restart with partial results requires specification of both filename and minimum acceptable confidance"
           <<  std::endl;
      exit(0);
    }
    
  }

  //-------------------------------------------------
  //  follow requirements
  RUNTIME_ASSERT( bStrainEnabled == vRequirementCheck[ eStrainOptConfigFilename ],
                  "If strain is enabled, StrainOpConfigFilename must be specified\n");
  //-------------------------------------------------
  bool bMissingTokens = false;
  for ( Int i = 0; i < eNumKeywords; i ++ )
    if( ! vInitializationCheck[ i ] )
    {
      cerr << "[ConfigFile] Missing variable: \'" << sKeywordStrings[i] << "\' not optional"  << endl;
      bMissingTokens = true;
    }
  delete [] vInitializationCheck;
  delete [] vRequirementCheck;

  if ( bMissingTokens )
    return false;
  
  return true;

}

///////////////////////////
//
//  Public
//  CConfigFile::InputConfigParameters - User funtion to initialize the config file
//  
///////////////////////////
bool CConfigFile::InputConfigParameters( const string &filename )
{
  Size_Type  nBufferSize = 0;
  char *pBuffer = InitFileIO::ReadFileToBuf( nBufferSize, filename);
    
  if( nBufferSize<= 0 || !pBuffer )
  {
    std::cerr << "ERROR: Unable to read config file" << endl; 
    return false;
  }
	
  if( !Parse( string( pBuffer, nBufferSize ) ) )
  {
    std::cerr << "ERROR: Unable to parse config file" << endl;
    return false;
  }
  
  if( pBuffer )
    delete[] pBuffer;
  
  return true;
}

///////////////////////////
//
//  Public
//  CConfigFile::PrintFile
//  
///////////////////////////
void CConfigFile::PrintFile()
{
  cout << "[PrintFile]: Printing ConfigFile" << endl;
  cout << "FileStartNum " << FileStartNum << endl
       << " FileEndNum " << FileEndNum << endl
       << " OutFileExt " << OutFileExt << endl
       << " InFileExt " << InFileExt << endl
       << " OutFileBasename " << OutFileBasename << endl
       << " InFileBasename " << InFileBasename << endl

	// Beam Information
       << " BeamEnergy " << BeamEnergy << endl
       << " BeamEnergyWidth " << BeamEnergyWidth << endl
       << " BeamHeight " << BeamHeight << endl

	// Detector Information;
       << "Detector Filename " << sDetectorFilename << endl

    // Sample Information
       << " SampleLocation " << SampleLocation[0] << " "
       << SampleLocation[1] << " " << SampleLocation[2] << endl
       << " SampleOrientation " << RADIAN_TO_DEGREE( SampleOrientation[0] )
       << " " << RADIAN_TO_DEGREE( SampleOrientation[1] ) << " " 
       << RADIAN_TO_DEGREE( SampleOrientation[2] ) << endl
       << " " << SampleFilename << endl;
	
}

//----------------------------------------------------------------
//  Save
//
//  no strings or files related stuff are saved
//
//----------------------------------------------------------------
Bool CConfigFile::Save  ( CSerializer & oSerialBuf ) const
{
  bool bSuccess;
  // Beam Information
  bSuccess = oSerialBuf.InsertCompactObj( BeamEnergy );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( BeamEnergyWidth );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( BeamHeight );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( BeamDirection );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( BeamDeflectionChiLaue );
  
    
  //  Search information
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nNumParamOptSteps );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nOptNumElementPerPE );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fParameterMCTemperature );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nOrientationSearchMethod );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fCoolingFraction );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fThermalizeFraction );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nParameterRefinements );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nDetectors );
  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( vDetDistSpacing );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fDetDistDeviation );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( oDetOrientDeviationEuler );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fDetOrientDeviationSO3 );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nMaxParamMCLocalRestarts );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nMaxParamMCGlobalRestarts );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nParamMCGlobalSearchElements );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( bConstrainedParamMC );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fSearchVolReductionFactor );

  //  Scattering Vector Based Initial Guesses
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fConsistencyError );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fBraggFilterTol );
  
  // Sample Information
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( SampleLocation );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( SampleRadius );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( SampleCenter );
  
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( SampleOrientation );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( SampleSymmetry );
  
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxInitSideLength );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMinSideLength );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMinAmplitudeFraction );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxQ );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fEtaLimit );
 
  // Search algorithm parameters
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fLocalOrientationGridRadius );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nMinLocalResolution );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nMaxLocalResolution );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nMaxDiscreteCandidates );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxAcceptedCost );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxConvergenceCost );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxDeepeningHitRatio );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nMaxMCSteps );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMCRadiusScaleFactor );
  
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nSuccessiveRestarts );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMinAccelerationThreshold );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nSecondsBetweenSave );

  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( vCommands );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxCost );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxMisorientation );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxRadius );

  //------------------ Strain related objets;
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( bStrainEnabled );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( StrainOptSetup );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( bUsePartialResult );


  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nMicGridType );

  //----------------- Restart information
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fPartialResultAcceptanceConf );
  return bSuccess;
}

//----------------------------------------------------------------
//  Restore
//
//
//
//----------------------------------------------------------------
Bool CConfigFile::Restore( CDeserializer & oSerialBuf )
{
  bool bSuccess;
  // Beam Information
  bSuccess = oSerialBuf.GetCompactObj( &BeamEnergy );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &BeamEnergyWidth );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &BeamHeight );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &BeamDirection );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &BeamDeflectionChiLaue );
  
    
  //  Search information
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nNumParamOptSteps );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nOptNumElementPerPE );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fParameterMCTemperature );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nOrientationSearchMethod );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fCoolingFraction );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fThermalizeFraction );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nParameterRefinements );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nDetectors );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vDetDistSpacing );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fDetDistDeviation );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &oDetOrientDeviationEuler );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fDetOrientDeviationSO3 );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nMaxParamMCLocalRestarts );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nMaxParamMCGlobalRestarts );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nParamMCGlobalSearchElements );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &bConstrainedParamMC );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fSearchVolReductionFactor );

  //  Scattering Vector Based Initial Guesses
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fConsistencyError );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fBraggFilterTol );
  
  // Sample Information
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &SampleLocation );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &SampleRadius );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &SampleCenter );
  
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &SampleOrientation );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &SampleSymmetry );
  
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxInitSideLength );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMinSideLength );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMinAmplitudeFraction );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxQ );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fEtaLimit );
 
  // Search algorithm parameters
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fLocalOrientationGridRadius );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nMinLocalResolution );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nMaxLocalResolution );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nMaxDiscreteCandidates );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxAcceptedCost );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxConvergenceCost );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxDeepeningHitRatio );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nMaxMCSteps );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMCRadiusScaleFactor );
  
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nSuccessiveRestarts );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMinAccelerationThreshold );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nSecondsBetweenSave );

  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vCommands );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxCost );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxMisorientation );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fMaxRadius );

  //------------------ Strain related objets;
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &bStrainEnabled );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &StrainOptSetup );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & bUsePartialResult );
  
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &nMicGridType );


  //----------------- Restart information
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &fPartialResultAcceptanceConf );
  return bSuccess;
}

    
