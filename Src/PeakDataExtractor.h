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
//  PeakDataExtractor
//  Author:      Frankie Li
//  e-mail:      sfli@cmu.edu
//  Purpose:     Implementation of cross correlation metric for diffraction
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _PEAK_EXTRACTOR_H
#define _PEAK_EXTRACTOR_H

#include "Simulation.h"
#include "SimulationData.h"
#include "DetectorFile.h"
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include "Raster.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "PeakFilters.h"
#include <boost/shared_ptr.hpp>



namespace IceNineAnalysis
{
////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  class ExtracPeaks
//
//        Given a collection of voxels, extract the set of peaks overlapping in the detector images.
//////////////////////////////////////////////////////////////////////////////////////////////////////////
class PeakDataExtractor
{

public:
  class OutputRaster :
    public CRaster< IntensityT, boost::numeric::ublas::matrix<IntensityT> >
  {
  public:
    //--------------------------
    //   Add Triangle
    //--------------------------
    void AddTriangle( const Point &v0, const Point & v1, const Point & v2, Float fIntensity )
    {
      RasterizeTriangle( v0, v1, v2, fIntensity );
    }
  };
  typedef boost::multi_array<OutputRaster, 2>     ImageMap;
  
  
  
private:

  //------------------------------------------------------
  //
  //  SimulatePeaks
  //  Input:   oCurrentLayer, oOuputImage, range 
  //  Output:  oDetector
  //
  //  TODO:   Also, fix the non-const of oCurrentLayer
  //
  //------------------------------------------------------
  void SimulatePeaks( ImageMap & oSimData, const DetectorListT & vDetectorList,
                      CSample oCurrentLayer,  const CSimulationRange &oRangeToIndexMap );

  boost::shared_ptr<CSimulation> pSimulator;        // make this a reference/pointers soon

  boost::shared_ptr<ReconstructionSetup> pSetup; // should probably be parameter instead of global variable in this scope

public:
  
  PeakDataExtractor( boost::shared_ptr<ReconstructionSetup> _pSetup,
			  boost::shared_ptr<CSimulation> _pSimulator )
    :pSetup( _pSetup), pSimulator( _pSimulator ) {}
  
  
  //------------------------------------------------------
  //
  //------------------------------------------------------
  void GeneratePeaks( ImageMap & oSimData, const DetectorListT & vDetectorList,
		      CSample oCurrentLayer, const CSimulationRange &oRangeToIndexMap )

    
  
  

};


}
#endif
