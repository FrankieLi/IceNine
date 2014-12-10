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
////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//   ExperimentSetup.h
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:  *Container class* for experimental parameters.  Specifically to hide from the user the
//             particular files input output, parsing, and error handling.  It differs from Simulation.h
//             in that no computation is done here.
//
//             Currently CExperimentSetup is the abstract class that defines the minimum requirement for
//             for an x-ray experiment.  CXDMExpermentSetup contains specific parameters and initialization
//             routines used for the 3DXDM, or now HEDM experiment.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


 
#ifndef _EXPERIMENT_SETUP_H_
#define _EXPERIMENT_SETUP_H_

#include <string>
#include "Sample.h"
#include "Detector.h"
#include "Debug.h"
#include "SimulationData.h"
#include "DetectorFile.h"
#include "Serializer.h"
#include "Symmetry.h"

#define MAX_QUADTREE_DEPTH 10  // maximum tree depth

using namespace XDMSimulation;

class CExperimentSetup
{
protected:

  CConfigFile oExpConfigFile;
  bool bInitialized;

  // Experiment setup and parameters
  SVector3 oBeamDirection;
  Float fMinAcceptedIntensityFraction;
  Float fBeamEnergy;
  Float fBeamEnergyWidth;
  Float fBeamDeflectionChiLaue;
  Float fEtaLimit;     // absolute value of limit of angle range from azimuthal
  
public:
  
  //---------------------------------
  //  Default Constructor:  Use with extreme caution!
  //  This must be used with a delayed "SetConfigFile"
  //---------------------------------
  CExperimentSetup(){};
  CExperimentSetup( const CConfigFile & oConfigFile );

  //---------------------------------
  //  SetConfigFile
  //---------------------------------
  void SetConfigFile( const CConfigFile & oConfigFile );
  
  virtual ~CExperimentSetup() {};
  
  //--------------------------------------------------------------------
  //  PositionExperiment:  
  //
  //  Action: Correctly position the sample and image at the proper locations
  //  and orientations.  Also initialize oSample and oDetImage to the correct
  //  parameters with oExpConfigFile
  //--------------------------------------------------------------------
  virtual void InitializeExperiment(CSample & oSample, CDetector & vDetImage) = 0;
  
  //--------------------------------------------------------------------
  //  GetNextSample:  (should write this as iterator) 
  //
  //  Action: Iteratores through the list of samples
  //--------------------------------------------------------------------	
  virtual const string & GetNextSample() const = 0;

  //--------------------------------------------------------------------
  //  Accessors
  //--------------------------------------------------------------------
  virtual Float GetBeamEnergy() const = 0;
  virtual Float GetBeamEnergyWidth() const = 0;
  virtual const SVector3 & GetXrayBeamDirection() const = 0;
  virtual Float GetMinAcceptedIntensityFraction() const = 0;
  virtual Float GetEtaLimit() const = 0;
  
  //--------------------------------------------------------------------
  //  GetChiLaue  - returns deflection angle of x-ray beam
  //--------------------------------------------------------------------
  virtual Float GetBeamDeflectionChiLaue() const = 0;
  
  //--------------------------------------------------------------------
  //  Save and Restore
  //
  //  Purpose:  Pack all of the useful data from this class
  //            into a *char array for transport through MPI.
  //--------------------------------------------------------------------
  virtual Bool Save  ( CSerializer   & oSerialBuf ) const = 0;
  virtual Bool Restore( CDeserializer & oSerialBuf ) = 0;
  
};

//--------------------------------------------------------------------
// CXDMExperimentSetup  A builder class that implements CExperimentSetup
// The initialization will follow the convention and setup of the XDM setup
//--------------------------------------------------------------------
class CXDMExperimentSetup : public CExperimentSetup
{

public:
  struct SStepSizeInfo
  {
    SVector3 oEulerSteps;
    SVector3 oDetectorPos;
    Float fBeamCenterJ;
    Float fBeamCenterK;
    Float fPixelHeight;
    Float fPixelWidth;
    Float fAngularRadius; // angular search radius
  };
  
private:
 
  void SetDetectionLimit( CUnitCell & oCellStruct,
                          const CDetector & oDetImage,
                          const CSample &oSample ) const;
  
  //--------------------------------------------------------------------
  //  I/O
  //--------------------------------------------------------------------
  Size_Type ReadRotationInterval( Size_Type nNumDetectors );
  Size_Type ReadDetectorInfo( );
  vector<SStepSizeInfo> ReadStepSizeFile ( const string & sFilename );
    
  //  for each vOmegaRangeList, there exist a filename (number)
  vector< SRange >     vOmegaRangeList; //  perhaps merge into a map?
                                        //  This seems to be a bit ugly
  
  vector< SIntRange >  vFileRangeList;  // simply a set of ranges of files for
                                        // different detector distances, move to
                                        // a detector distance structure
  
  vector<CDetector>     vDetectorList;            // a vector of detectors
  vector<SStepSizeInfo> vOptimizationInfo;        // A vector of optimization step sizes
  vector<SStepSizeInfo> vOptimizationConstrains;  // A vector of constrains between detectors
  vector<SStepSizeInfo> vDetectionSensitivity;  // A vector that describes the sensitivity of the detector geometry
  CSimulationRange      oRangeToIndexMap;   // an object that is able to map from an angle
                                            // to the proper file number and range index

public: 

  CXDMExperimentSetup(){};
  CXDMExperimentSetup( const CConfigFile & oConfigFile ): CExperimentSetup( oConfigFile ) {};
  ~CXDMExperimentSetup() {};
  
  void InitializeExperiment( );
  void InitializeExperiment(CSample & oSample, CDetector & vDetImage)
  {
    DEBUG_ASSERT(0, "NOT IMPLEMENTED:  InitializeExperiment(CSample & oSample, CDetector & vDetImage) \n");
  }

  void InitializeSample( CSample & oSample, const CDetector & oDetector );
  
  //---------------------------------
  // Geometry functions
  //
  // TODO:  Generalize belowed function to take only two points and determine maxQ
  //---------------------------------
  Float GetMaxQ( const CDetector & oDetImage, const CSample &oSample ) const;
  
  //---------------------------------
  // Accessors
  //---------------------------------
  Float GetBeamEnergy()                   const;
  Float GetEtaLimit()                     const;
  Float GetBeamEnergyWidth()              const;
  const SVector3 & GetXrayBeamDirection() const;
  Float GetMinAcceptedIntensityFraction() const;
  Float GetBeamDeflectionChiLaue()        const;
  const string & GetNextSample()          const;
  const CSymmetry* GetSampleSymmetry()    const;
  
  const vector< SRange > &    GetOmegaRangeList()       const { return vOmegaRangeList    ; } 
  const vector< SIntRange > & GetFileRangeList()        const { return vFileRangeList     ; } 
  const vector< CDetector > & GetDetectorList()         const { return vDetectorList      ; } 
  const vector< SStepSizeInfo > & GetOptimizationInfo() const { return vOptimizationInfo  ; }
  const vector< SStepSizeInfo > & GetOptimizationConstrains() const { return vOptimizationConstrains  ; }
  const vector< SStepSizeInfo > & GetDetectorSensitivity()    const { return vDetectionSensitivity ; } 
  const CConfigFile &         GetInputParameters() const { return oExpConfigFile; }
  
  const CSimulationRange &    GetRangeToIndexMap() const { return oRangeToIndexMap; } 


  //---------------------------------------------
  //  GetDetector - Direct access to a pointer to detector
  // 
  //---------------------------------------------
  CDetector* GetDetector( int n ) 
  {
    RUNTIME_ASSERT( n < vDetectorList.size(), 
		    " GetDetector: Tried to access detector n > vDetectorList.size()\n ");
    return &( vDetectorList[ n ] );
  }

  //----------------------------------
  //  Backward compatibility related accessors
  //----------------------------------
  Int GetBCDetectorOffset() { return oExpConfigFile.nBCPeakDetectorOffset; }
  
  //----------------------------------
  // Utilities
  //----------------------------------

  //-----------------
  //  GetReciprocalVector
  //
  //  Purpose:  Given oKOutDir, output direction of the momentum vector,
  //            a reciprocal lattice vector G that produced this output
  //            direction that is compatible with the experiment is returned.
  //            It is assumed that scattering is completely elastic.
  //
  //  Input:    oKOutDir is a *unit vector* in the direction of the output momentum.
  //
  //-----------------
  SVector3 GetReciprocalVector( const SVector3 & oKOutDir )    const;

  //----------------------------------
  //  Serializer for MPI
  //----------------------------------
  Bool Save  ( CSerializer   & oSerialBuf ) const;
  Bool Restore( CDeserializer & oSerialBuf );

  //----------------------------------
  // Mutator
  //----------------------------------

  //--------------------------------------------------------------------
  //  SetExperimentalParameters
  //
  //  Purpose:  Used to change experimental parameters after initialization.
  //
  //  NOTE:
  //  1.  nDetectorID ranges from [0, NumDetectors)
  //  2.  The detector limit will be reset.  i.e., MaxQ and so forth will be
  //      recalculated.
  //
  //--------------------------------------------------------------------
  void SetExperimentalParameters( const CXDMDetectorFactory::SDetParameters & oNewDetParam, Size_Type nDetectorID );

  //--------------------------------------------------------------------
  //  SetBeamEnergy
  //--------------------------------------------------------------------
  void SetBeamEnergy( Float fEnergy ){ fBeamEnergy = fEnergy; }
  
  //--------------------------------------------------------------------
  //  GetExperimentalParameters
  //  Purpose:  Return a set of experimental parameters that are modifiable
  //            A vector is returned corrpsonding to the detector list is returned
  //--------------------------------------------------------------------
  vector< CXDMDetectorFactory::SDetParameters > GetExperimentalParameters( ) const;
  
  //--------------------------------------------------------------------
  //  WriteDetectorFile
  //  Purpose:  Output the detector file to the ostream.
  //
  //--------------------------------------------------------------------
  void WriteDetectorFile( ostream & os ) const;
  
};
#endif
