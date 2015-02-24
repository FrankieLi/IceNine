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
//  Driver.cpp
//
//  Author:  Frankie Li (sfli@cmu.edu)
//   
//  Purpose:  Implementation of Driver.h
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Driver.h"


//---------------------------------------------------------------------------
//
//  Driver
//
//  A rudimentary driver program
//
//  TODO:  Switch to Boost command shell library
//
//---------------------------------------------------------------------------
int XDMDriver( int argc, char* args[] )
{
  if( argc < 3 )
  {
    cerr << "USAGE: " << args[0]
         << " [ (s)imulate | (r)econstruct | (p)arallel reconstruction ]  <Config File> " << endl;
    return -1;
  }

  if ( args[1][0] == 's')
  {
    SimulateDetectorImage( args[2] );
  }
  else if ( args[1][0] == 'r')
  {
    if( argc == 3 )
    {
      ReconstructSample( args[2] );
    }
    if( argc == 4 )  // running test
    {
      if( args[4] == NULL )     
      {
	std::cerr 
	<< "Serial reconstruction optional parameters: [(P)arameter Optimization | (G)rain Reconstruction) " 
	<< std::endl;
      }
      if( args[4][0] == 'P' )
      {
      }
      if( args[4][0] == 'G' )
      {
	ReconstructGrain( args[2] );
      }
    }
  }
  else if ( args[1][0] == 'p')
  {
    ParallelReconstruction( argc, args, args[2] );
  }
  else if ( args[1][0] == 'G')
  {
    if( argc != 4  || args[3] == NULL )
      std::cerr << "USAGE: " << args[0]
                << "G<Num Threads> <Config File> <PaintGrid Config>" << std::endl;
    
    int nThreads = atoi( ++args[1] );
    std::cout << "Running nThreads " << nThreads << std::endl;
    PaintGrid( args[2], args[3], nThreads );
  }
  else
  {
    cerr << "Undefined program option: "  << args[1][0] << endl;
    return -1;
  }

  return 0;
}

//---------------------------------------------------------------------------
//
//  SimulateDetectorImage
//
//
//---------------------------------------------------------------------------
void SimulateDetectorImage( const string &sConfigFilename )
{
  CConfigFile oConfigFile;
  if( !oConfigFile.InputConfigParameters( sConfigFilename ) )
  {
    std::cerr << "Experiment setup failed: Unable to read file: " << sConfigFilename << std::endl;
    exit(0);
  }

  CXDMForwardSimulation oSimulation( oConfigFile );
  if( oConfigFile.bStrainEnabled )
  {
    oSimulation.SimulateDetectorImagesOptimized( true );
  }
  else
  {
    oSimulation.SimulateDetectorImagesOptimized();
  }
}

//---------------------------------------------------------------------------
//
//  ReconstructSample
//
//---------------------------------------------------------------------------
void ReconstructSample( const string & sConfigFilename )
{
  CConfigFile oConfigFile;
  if( !oConfigFile.InputConfigParameters( sConfigFilename ) )
  {
    std::cerr << "Experiment setup failed: Unable to read file: " << sConfigFilename << std::endl;
    exit(0);
  }
  Reconstruction::SerialReconstruction oReconstruction( oConfigFile );
  oReconstruction.ReconstructSample();
}

//---------------------------------------------------------------------------
//
//  PaintGrid
//
//---------------------------------------------------------------------------
void PaintGrid( const string & sConfigFilename,
                const string & sPaintGridFilename,
                int nThreads )
{
  CConfigFile oConfigFile;
  if( !oConfigFile.InputConfigParameters( sConfigFilename ) )
  {
    std::cerr << "Experiment setup failed: Unable to read file: " << sConfigFilename << std::endl;
    exit(0);
  }
  
  CPaintGrid S03GridPainter( oConfigFile, nThreads );
  S03GridPainter.Run( sPaintGridFilename );
}

//---------------------------------------------------------------------------
//   ParallelReconstruction
//
//--------------------------------------------------------------------------- 
int ParallelReconstruction( int argc, char* args[],
                            const string & sConfigFilename )
{
  XDMParallel::CParallelProcessing oPReconstructor;
  oPReconstructor.Process( argc, args, sConfigFilename );
  
  return 0;
}

//---------------------------------------------------------------------------
//  ReconstructGrain
//---------------------------------------------------------------------------
void ReconstructGrain( const string & sConfigFilename )
{
  CConfigFile oConfigFile;
  if( !oConfigFile.InputConfigParameters( sConfigFilename ) )
  {
    std::cerr << "Experiment setup failed: Unable to read file: " << sConfigFilename << std::endl;
    exit(0);
  }
  Reconstruction::GrainReconstruction oReconstruction( oConfigFile );

  oReconstruction.ReconstructRandomGrain(  );


}
