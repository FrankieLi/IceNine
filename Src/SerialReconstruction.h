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
//  SerialReconstruction.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//   Purpose:  Simple serial reconstruction
//
////////////////////////////////////////////////////////////

#ifndef _SERIAL_RECONSTRUCTION_H
#define _SERIAL_RECONSTRUCTION_H

#include "Reconstructor.h"
#include "DiscreteAdaptive.h"
#include <ctime>
#include <boost/shared_ptr.hpp>

namespace Reconstruction
{
   //----------------------------------------
  //
  //  SerialReconstruction
  //
  //  This is the root of all serial reconstruction calls, at least
  //  for the simple version.
  //
  //----------------------------------------
  class SerialReconstruction
  {
  private:
    typedef ReconstructionStrategies::UniformSubidivisonSequentialSelection ReconstructionStrategy;
    typedef SearchDetails::SCandidate       SCandidate;
    typedef SearchDetails::SSearchParameter SSearchParameter;
    
    ReconstructionSetup     oSetup;
    ReconstructionStrategy  VoxelQueue;   // acts as a queue
    CSimulation oSimulator;
    SerialReconstruction();
  public:

    SerialReconstruction( const CConfigFile & oConfigFile )
    {
      RUNTIME_ASSERT( oConfigFile.nMicGridType == eTriangular, "Non-triangular grid not implemented for serial reconstruction\n" );
      oSetup.InitializeWithDataFiles( oConfigFile );
      boost::shared_ptr<CMic> pMic= boost::dynamic_pointer_cast<CMic>( oSetup.ReconstructionRegion());
      VoxelQueue.Initialize( *pMic,
                             oSetup.MinSideLength() );
      oSimulator.Initialize( oSetup.ExperimentalSetup() );
    }
    
    //----------------------------------------
    // ReconstructSample
    //----------------------------------------
    void ReconstructSample()
    {
      BasicVoxelReconstructor     Reconstructor( oSimulator, oSetup );
      DiscreteRefinement<SVoxel>  AdpReconstructor( oSimulator, oSetup );
      std::cout << "Num Elements in Queue " << VoxelQueue.Size() << std::endl;
      while( VoxelQueue.Size() != 0 )
      {
        SVoxel Result;
        Int nCode;
        time_t oStartTime, oStopTime;
        time( & oStartTime );
        boost::tie( Result, nCode ) = Reconstructor.ReconstructVoxel( VoxelQueue.First() );
        time( & oStopTime );
        double oTimeDiff = difftime( oStopTime, oStartTime );
        std::cout << "Normal Reconstruction " << oTimeDiff << " sec " << std::endl;
        
        time( & oStartTime );
        boost::tie( Result, nCode ) = AdpReconstructor.ReconstructVoxel( VoxelQueue.First() );
        time( & oStopTime );
        oTimeDiff = difftime( oStopTime, oStartTime );
        std::cout << "Adaptive Reconstruction " << oTimeDiff << " sec " << std::endl;

        VoxelQueue.Pop();
        VoxelQueue.Push( Result );
      }  
    }

    //----------------------------------------
    // ReconstructedMic
    //----------------------------------------
    CMic GetReconstructedMic() const
    {
      return VoxelQueue.Solution();
    }
  };
}


#endif
