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
//  PaintGrid.h
//  Author:      Frankie Li
//  e-mail:      sfli@cmu.edu/li31@llnl.gov
//  Purpose:     Implementation of J. V. Bernier's Paint Grid algorithm
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _PAINT_GRID_H
#define _PAINT_GRID_H

#include "Reconstructor.h"
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
#include "InitFilesIO.h"
#include "BBox.h"
#include <sstream>

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  class CPaintGrid
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
class CPaintGrid : public CSimulation
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
  
  typedef boost::multi_array<OutputRaster, 2>  ImageMap;
  typedef boost::multi_array<Float, 2>         EtaOmegaMapT;
  typedef PBRMath::BBox2D                      BBoxT;
  typedef std::vector<EtaOmegaMapT>            EtaOmega_2ThetaMapT;


  typedef boost::multi_array<Float, 3>         IntensityMapT;
  
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
  std::vector<Float> Paint( const DetectorListT & vDetectorList,
                            const vector<SMatrix3x3> & FZList,
                            CSample oCurrentLayer );
  
  //------------------------------------------------------
  //  DirectPaint - directly result in "intensities"
  //------------------------------------------------------
  void DirectPaint( const DetectorListT & vDetectorList,
                    const vector<SMatrix3x3> & FZList,
                    CSample oCurrentLayer );


  //------------------------------------------------------
  //  DilateIntensityMap
  //  Dilation of Two-theta rings - strictly in the two theta,
  //  in the radial direction.
  //------------------------------------------------------
  void DilateIntensityMap( Float fOmegaWidth, Float fEtaWidth, Float f2ThetaWidth );


  //------------------------------------------------------
  //  BBoxTwoEtaBBox
  //------------------------------------------------------
  //BBoxT GetEtaTwoThetaBoundingBox( const BBoxT & PeakPixelBBox ) const;
  void PixelToEtaTwoTheta( Float & Eta, Float & TwoTheta,
                           Float  x, Float y,
                           const CDetector & Detector,
                           const SVector3 & Center) const;
  
  //------------------------------------------------------
  //  ConstructEtaOmegaMap
  //
  //   - Construct the eta-omega map by iterating through all detector images
  //
  //------------------------------------------------------
  void ConstructEtaOmegaMap( const DetectorListT & vDetectorList,
                             CSample oCurrentLayer );
  
  
  //------------------------------------------------------
  // ReadOrientationTestPoints
  //------------------------------------------------------
  void ReadOrientationTestPoints( const std::string & OrientationFilename,
                                  std::vector<SMatrix3x3 > & FZList );

  //------------------------------------------------------
  // InitializeIntensityMap
  //------------------------------------------------------
  
  void InitializeIntensityMap( Float EtaBinSize,
                               Float f2ThetaBinSize, Float f2ThetaPeakWidth,
                               const vector<CRecpVector> & oRecipVectors ); // need to be replaced
  
  //------------------------------------------------------
  // ConstructTwoThetaRange
  //------------------------------------------------------
  std::vector<SRange> ConstructTwoThetaRanges( Float ThetaMapWidth,
                                               Float ThetaWidth,
                                               const vector<CRecpVector> & oRecipVectors );
  
  int NumThreads;
  IntensityMapT IntensityMap;   // Mapping (2theta, omega, eta) -> Intensity 
  
  CSimulationRange    OmegaRangeToIndexMap;
  CSimulationRange    EtaRangeToIndexMap;
  CSimulationRange    TwoThetaRangeToIndexMap;
  

  Reconstruction::ReconstructionSetup     oSetup;
  
  CSimulation         oSimulator;        // make this a reference/pointers soon
  CConfigFile         oSetupFile;    
  


public:
  
  CPaintGrid( const CConfigFile & oConfigFile, int nThreads );

  bool Run( const std::string & PaintGridConfig );

};
#endif
