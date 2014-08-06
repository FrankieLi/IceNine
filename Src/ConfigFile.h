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
//  ConfigFile.h
//
//  Classes:  CConfigFile
//
//  Purpose:  A container for initialization information
//            for simulations.  Contains a parser to read
//            simple text based configuration files.
//
//
//  This file needs some restructuring
//
///////////////////////////////////////////////////////



#ifndef ICENINE_CONFIG_H
#define ICENINE_CONFIG_H


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include "Debug.h"
#include "Types.h"
#include "3dMath.h"
#include "Error.h"
#include "Parser.h"
#include "XDMParallel.h"
#include "InitFilesIO.h"
#include "MicIO.h"
#include "Symmetry.h"
#include "Serializer.h"
#include "OptimizationConfig.h"

using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;

//-------------------------
//  SO3SearchMethod
//
//  Tags for possible methods of the
//  search algorithm.
//-------------------------
enum SO3SearchMethod{
  eConstraintedEuler, eUniformQuaternion,
  eNumSO3SearchTypes
};

enum EFileType{
  eBin, eTif, eASCII, eNumFileTypes
};


//------------------------------------------
//
//   Class CConfigFile
//
//------------------------------------------
class CConfigFile
{
  
private:

  bool Parse( const string & sBuf);	

public:
  
  UInt FileStartNum;
  UInt FileEndNum;
  string InFileExt;
  EFileType InFileType;
  string InFileBasename;
  Int    InFileSerialLength;

  string OutFileExt;
  string OutFileBasename;
  Int    OutFileSerialLength;

  string OutStructureBasename;
	
  // Beam Information
  Float BeamEnergy;
  Float BeamEnergyWidth;
  Float BeamHeight;
  SVector3 BeamDirection;
  Float BeamDeflectionChiLaue;  // deflection angle, Chi Laue

  // Detector Information;
  string sDetectorFilename;


  //------------------------------------------
  // Parameter Optimization
  //
  //  TODO:  clean this up in the future -- put them into a seperate initialization file
  //
  string sOptimizationFilename;
  string sOptimizationConstrainFilename;
  string sDetectionLimitFilename;
  
  Int    nNumParamOptSteps;
  Int    nOptNumElementPerPE;
  Float  fParameterMCTemperature;
  Int    nOrientationSearchMethod;
  Float  fCoolingFraction;      // how much time spent on cooling (linear)
  Float  fThermalizeFraction;   // how much time spent on Thermalizing
  Int    nParameterRefinements;  // how many times to refine a parameter search

  Int            nDetectors;     // this has to occur before everything below!
  vector<Float>  vDetDistSpacing;
  Float          fDetDistDeviation;
  SVector3       oDetOrientDeviationEuler;
  Float          fDetOrientDeviationSO3;

  Int   nMaxParamMCLocalRestarts;
  Int   nMaxParamMCGlobalRestarts;
  Int   nParamMCGlobalSearchElements;
  Bool  bConstrainedParamMC;
  Float fSearchVolReductionFactor;
  //------------------------------------------

  //------------------------------------------
  //  Scattering Vector Based Initial Guesses

  Float fConsistencyError; // tolerance on what is considered as a consistent
                           // rotational matrix
  
  Float fBraggFilterTol;   // tolerance on the omega interval for bragg acception
  //
  //
  //------------------------------------------
  
  // Sample Information
  SVector3  SampleLocation;
  Float     SampleRadius;
  SVector3  SampleCenter;
  
  SVector3  SampleOrientation;
  string SampleFilename;
  string StructureFilename;
  string FundamentalZoneFilename;
  ESymmetryT SampleSymmetry;
  
  Float fMaxInitSideLength;   // maxmimum sidelength allowed initially
  Float fMinSideLength;          // minimum sidelength allowed
  
  
  // Simulation parameters
  Float fMinAmplitudeFraction;  // max fraction of amplitude to be taken 
  Float fMaxQ;
  Float fEtaLimit;              // symmetric limit for eta

  // Initialization File Information
  string sRotationRangeFilename;  // a.k.a:  omega files
  
  
  // Search algorithm parameters
  Float fLocalOrientationGridRadius;
  Int   nMinLocalResolution;
  Int   nMaxLocalResolution;

  Int   nMaxDiscreteCandidates; // Maximum number of candidates after the initial discrete search
                                  // Clearly, with better cost function, this can become a smaller number.
  Float fMaxAcceptedCost;
  Float fMaxConvergenceCost;
  Float fMaxDeepeningHitRatio;   // Hit ratio below which iterative deepening will continue
                                 // to refine.
  Float nMaxMCSteps;
  Float fMCRadiusScaleFactor;    //  A fudge factor used to scale the radius of
                                 //  the Monte Carlo step.  This is to be taken out
                                 //  once further testing tells us more about the
                                 //  search algorithm.
  Int   nSuccessiveRestarts;
  
  Float fMinAccelerationThreshold;   //  Threshold (currently hit ration) before acceleration kicks in
        
  // File Saving options
  Int   nSecondsBetweenSave; 
  
  // Pipeline scheduling
  vector<Int>  vCommands;        // TODO:  Generalization of this will require a second file for pipelines only
    
  vector< vector<string> > TokenizedConfigFile;
  
  // Backward compatibility related functions  (everything starts with BC)
  //  (To keep compatible with XDMMPI)
  //
  Int nBCPeakDetectorOffset; 

  //------------- For boundary selection
  Float fMaxCost;
  Float fMaxMisorientation;
  Float fMaxRadius;
  //------------- end For boundary selection

  Bool bStrainEnabled;
  
  OptimizationConfig::StrainConfigFile StrainOptSetup;
  
  //------------- For Partial Results
  bool   bUsePartialResult;
  float  fPartialResultAcceptanceConf;
  string PartialResultFilename;
  //------------- For Partial Results


  

  Int    nMicGridType;
  

  
  
  //-------------------------------------------------
  //  CConfigFile default constructor
  //-------------------------------------------------
  CConfigFile():FileStartNum(0), FileEndNum(0),
                InFileExt(""), InFileBasename(""), OutFileExt(""), OutFileBasename(""),
                OutStructureBasename(""),
                BeamEnergy(0), BeamEnergyWidth(0), BeamHeight(0),
                nDetectors( -1 ), bUsePartialResult( false )  {}

  //-------------------------------------------------
  //  InputConfigParameters
  //-------------------------------------------------
  bool InputConfigParameters(const string &filename);
  
  //-------------------------------------------------
  //  PrintFile
  //-------------------------------------------------
  void PrintFile(  );


  //-------------------------------------------------
  //  Save and restore - this is limited to the
  //  values but not the buffers, strings, or file related
  //  parameters.
  //
  //
  //-------------------------------------------------
  Bool Save  ( CSerializer   & oSerialBuf ) const;
  Bool Restore( CDeserializer & oSerialBuf );
    
};



#endif
