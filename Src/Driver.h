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
//  Driver.h
//
//  Author:  Frankie Li (sfli@cmu.edu)
//   
//  Purpose:  Driver program for the entire suite of simulation tools
//            This is the top level of the program, and therefore there
//            should be no computation directly in this object. The driver  
//            should also be the location where exceptions are caught
//            and handled.  Stacktrace should be printed in this object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef _DRIVER_H_
#define _DRIVER_H_

#include "DetectorCalibration.h"
#include "ConfigFile.h"
#include "Debug.h"
#include <string>
#include <vector>
#include <fstream>

// For simulations
#include "3dMath.h"
#include "Sample.h"
#include "Detector.h"
#include "Types.h"
#include "Debug.h"
#include <string>
#include "BBox.h"
#include "ExperimentSetup.h"
#include "PhysicalConstants.h"

#include "ImageData.h"
#include "Quaternion.h"

#include "ForwardSimulation.h"
#include "SerialReconstruction.h"

// #include "ReconstructionAnalysis.h"
#include "ParallelDriver.h"
#include "PaintGrid.h"

#include "GrainReconstruction.h"
#include "DetectorCalibration.h"
//#include "UnitTest.h"

using namespace GeneralLib;
using std::ofstream;
using std::ifstream;
using std::string;
using std::vector;

using std::cout;

////////////////////////////////////////////
//
//  SimulateDetectorImage
//
//
////////////////////////////////////////////
int XDMDriver( int argc, char* args[] );
void SimulateDetectorImage( const string &sConfigFilename );
void ReconstructSample( const string & sConfigFilename );
void ReconstructGrain( const string & sConfigFilename );
void ReoptimizeSample( const string & sConfigFilename );
void CostCalculation  ( const string & sConfigFilename );

void DetectorOptimization  ( const string & sConfigFilename );
void PaintGrid( const string & sConfigFilename,
                const string & sPaintGridFilename,
                int nThreads );


//----------------------------------
//   ParallelReconstruction
//
//   The basic setup for parallel operations.
//----------------------------------
int ParallelReconstruction( int argc, char* args[],
                            const string & sConfigFilename );



#endif
