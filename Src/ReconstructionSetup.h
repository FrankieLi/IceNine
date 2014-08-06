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
//  ReconstructionSetup.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//   Purpose:  This is the data portion of a simple reconstruction.
//             ReconstructionSetup (probably should be renamed) is
//             where all initializations take place.
//
//
////////////////////////////////////////////////////////////


#ifndef RECONSTRUCTION_SETUP_H
#define RECONSTRUCTION_SETUP_H

#include <boost/function.hpp>
#include "boost/multi_array.hpp"
#include "Quaternion.h"
#include "Serializer.h"
#include "Symmetry.h"
#include <boost/tuple/tuple.hpp>
#include <limits>
#include "OverlapInfo.h"
#include "SimulationData.h"
#include "ExperimentSetup.h"
#include "PeakFilters.h"

namespace Reconstruction
{

  class ReconstructionSetup
  {
  public:
    typedef CSimulationData::ImageIteratorT  ImageIteratorT;
    typedef CSample::MicIOPtrT               MicIOPtrT;

  private:
    //------------------------------------------------------------------------
    //
    //  Private: ReadExperimentalData
    //
    //  Experimental data is read based on information from oExpSetup
    //
    //------------------------------------------------------------------------
    Bool ReadExperimentalData( XDMSimulation::CSimulationData & oExpData, bool bServerMode = false );
    
  public:
    
    //------------------------------------------------------------------------
    //
    //  A reconstruction must first be initialized, either with data files 
    //  or a restore using a Deserializer.  This function must be called
    //  before any other function!
    //
    //
    //  bServerMode is used by server - Enabling server mode disables the building
    //  of the sparse matrix and the search tree, as the server will never use them.
    //------------------------------------------------------------------------
    virtual void InitializeWithDataFiles( const CConfigFile & oConfigFile, bool bServerMode = false );

    
    CXDMExperimentSetup  oExpSetup;
    vector<SMatrix3x3>   oFZOrientationList;     // remove this - this is decided by specific reconstruction methods
    CSimulationData      oExpData;
    CSample              oReconstructionLayer;   //  TODO:  Move voxel information
                                                  //  away from the rest of sample info?

    MicIOPtrT            pPartialResult;
    typedef HEDM::XDMSimpleEtaFilter EtaAngularFilter; 
    typedef HEDM::XDMEtaAcceptFn     EtaIntensityFilter;
    EtaAngularFilter   FPeakFilter;         //  need ability to choose different filters
    EtaIntensityFilter FPeakIntensityAccept;

    //-------------------------------
    //  Accessor
    //-------------------------------
    const CXDMExperimentSetup & ExperimentalSetup()    const { return oExpSetup;   }
    const EtaAngularFilter    & EtaThresholdFilter()   const { return FPeakFilter; }
    const CConfigFile         & InputParameters()      const { return oExpSetup.GetInputParameters(); }
    const vector<SMatrix3x3>  & DiscreteSearchPoints() const { return oFZOrientationList; }
    Float MinSideLength() const  { return InputParameters().fMinSideLength; }


    //----------- mutator --------------
    CXDMExperimentSetup & ExperimentalSetup()  { return oExpSetup;   }
    CSimulationData & Data()                   { return oExpData; }
    
    //-------------------------------
    //  These guys need to become const
    //-------------------------------
    CSample &         SampleGeometry()         { return oReconstructionLayer; }
    const CSimulationData & Data()       const { return oExpData; }
    
    const MicIOPtrT ReconstructionRegion() const { return oReconstructionLayer.GetMic(); }
    MicIOPtrT       ReconstructionRegion()       { return oReconstructionLayer.GetMic(); }
    
    const MicIOPtrT PartialResult() const { return pPartialResult; }
    MicIOPtrT       PartialResult()       { return pPartialResult; }
    
    Bool Save   ( CSerializer   & oSerialBuf ) const;
    Bool Restore( CDeserializer & oSerialBuf ) ;

    
    //-------------------------------
    // For debugging use
    //-------------------------------
    template< typename ImageMapT >
    void SetSimulationData( const ImageMapT & ImageMap )
    {

      for( int i = 0; i < ImageMap.shape()[0]; i ++ )
      {
        for( int j = 0; j < ImageMap.shape()[1]; j ++ )
        {
          oExpData.mImageMap[i][j].Overwrite( ImageMap[i][j] );
        }
      }
    }

    

  };
    
}


#endif
