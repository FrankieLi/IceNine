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
//  Simulation.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//   Purpose:  This is the container for major computations for the forward modelling simulation.
//             Note that parameters are input first via initialization using an external file.  Note
//             also that this class is specific to CXDMSimulation, meaning that the simulation is specific
//             to ONE experiment.  Future interests in simulating other systems will mean a replacement of
//             the simulation class.  (Note that the geometry class and ExperimentSetup class may still be
//             kept in tact.)
//
//
//             The intent for this class is to be general.  Therefore optimization techniques are not used
//             unless they abstraction maybe preserved.   The justification is that this is a class built to
//             test algorithms.  Therefore each component must work independently without some hidden global
//             variable.  Transparacny and abstraction are taken precedence.  Another version of this code
//             may be implemented with speed optimization, and one could easily call it CXDMFastSimulation.
//             It would probably be a good idea to write an abstract class for general simulation to keep
//             everything organized.
//
////////////////////////////////////////////////////////////
#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "3dMath.h"
#include "Sample.h"
#include "Detector.h"
#include "DetectorFile.h"
#include "Types.h"
#include "Debug.h"
#include <string>
#include "BBox.h"
#include "ExperimentSetup.h"
#include "PhysicalConstants.h"
#include "ConfigFile.h"
#include "ImageData.h"
#include "SimulationData.h"
#include "boost/multi_array.hpp"
#include "ProgramLimits.h"
#include <tuple>

#include "DiffractionCore.h"    // All inline functions are here
#include "PeakFilters.h"

using namespace GeneralLib;
using namespace XDMSimulation;

using namespace DiffractionCore;
#ifndef NUM_VERTICES
#define NUM_VERTICES 3
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 
//  class CXDMSimulationState
//
//  A class used to document the state of the simulation
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
class CXDMSimulationState
{
private:
  string sLogFilename;
  ofstream oLogFile;

public:

  ~CXDMSimulationState();
  
  Float fOmega;                              // simulating rotation angle
  vector <Point> vCurrentPixelList;          // The list of pixels generated from the scattering
  SVoxel oVoxel;                             // scattering "origin."  The voxel being processed
  Int    nVoxelID;                           // the ID of the current voxel
  SVector3 oScatteringVector;                // the scattering vector that produced this peak
  Int    nRecipVecID;                        // the ID of the current reciprocal vector
  
  bool InitializeLogFile( string sFilename );
  void PrintToLog();
  void Finalize() { oLogFile.close(); };
  void Clear() { vCurrentPixelList.clear(); } ;   // may need a better method
  
};



typedef vector<CDetector> DetectorListT;
typedef boost::multi_array< Point, 2 > PointMapT;


typedef boost::multi_array< Point,     3 > ProjectedPixelMapT;
typedef boost::multi_array< Size_Type, 2 > PixelCountMapT;
typedef boost::multi_array< Float,     2 > PixelIntensityMapT;


typedef const ProjectedPixelMapT::const_array_view<2>::type ConstProjectedPixelView;
typedef const PixelCountMapT    ::const_array_view<1>::type ConstPixelCountView;
typedef const PixelIntensityMapT::const_array_view<1>::type ConstPixelIntensityView;

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 
//  class CSimulation
//
//  Base class of all simulation.  This contains all of the necessary functions to actually run
//  simulations of all types.
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
class CSimulation   // to be parameterized by acceptance function?
{

public:
  
  typedef CSimulationData::ImageMapT::array_view<1>::type ImageListT;
  typedef CSimulationData::ImageMapT::const_array_view<1>::type Const_ImageListT;
  
protected:  // was private
  
  //  CXDMSimulationState oSimulationState;   // Used to save the state of the simulation

  SVector3 oBeamDirection;
  Float    fBeamEnergy;
  Float    fBeamDeflectionChiLaue;
  Bool bSimulationInitialized;


  struct FTrivialFilter
  {
    Float fIntensity;
    FTrivialFilter(): fIntensity( 0 ) {}
    std::pair<bool, Float> operator() ( const SVector3 & d ) const
    { return std::make_pair( true, fIntensity ); }
  };
  
  //-------------------------------------------------------------------------------------
  //    H E L P E R   F U N C T I O N S  ( for class specialization )
  //-------------------------------------------------------------------------------------

  FTrivialFilter TrivialAcceptFilter;
  
public:
  
  //------------------------------------------------------------------------------
  //  GetProjectedVertices
  //  Purpose:  Project Vertices based on given geometry.  Note that certain projections
  //            will 
  //------------------------------------------------------------------------------
  template < class ProjectedPixelMap, class PixelCountMap,
             class PixelIntensityMap, class VertexIter,
             class PeakFilterFn >
  void GetProjectedVertices( ProjectedPixelMap & vProjectedPixels,
                             PixelCountMap     & vNumPixelHit,
                             PixelIntensityMap & vPixelIntensity,
                             const DetectorListT & oDetectorList, const CSample &oSample,
                             VertexIter pStart, VertexIter pEnd,  const SVector3 & oNormal,
                             PeakFilterFn FPeakAcceptFilter ) const;
               
  //-------------------------------------------------------------------------------------
  //  ProjectVoxel - Fast version that takes in a view of the Simulation data ImageMap
  //  Note that oOutImage must correspond to the ordering of oDetectorList
  //
  //
  //  TODO:  To be come specialization of the templated version soon.
  //-------------------------------------------------------------------------------------
  bool ProjectVoxel( ImageListT & oOutImage, const DetectorListT & oDetectorList,
                     const CSample &oSample, const SVoxel &oVoxel,
                     const SVector3 &oNormal, Float fIntensity );

  //-------------------------------------------------------------------------------------
  //  ProjectVoxel
  //  Purpose:  Templated, generalized ProjectVoxel.  A functional is taken as parameter
  //            to filter vertices.
  //-------------------------------------------------------------------------------------
  template < class OutputImageListT, class PeakFilterFn >
  bool ProjectVoxel( OutputImageListT & oOutImage, const DetectorListT & oDetectorList,
                     const CSample &oSample, const SVoxel &oVoxel,
                     const SVector3 &oNormal, PeakFilterFn FPeakFilter );
  
  //-------------------------------------------------------------------------------------
  //
  //  Explicit Functions - These are functions that are slow but work correctly
  //
  //-------------------------------------------------------------------------------------
  bool ProjectVoxel( CImageData & oOutputImage, const CDetector &oDetector,
                     const CSample &oSample, const SVoxel &voxel,
                     const SVector3 &oNormal, Float fIntensity);
  
  //-------------------------------------------------------------------------------------
  //
  //  GetScatteringOmega
  //
  //
  //  Input:  oScatteringVec - scattering vector of the sample in the sample frame,
  //          i.e., omega = 0
  //          Note scattering magitude is added for speed optimization
  //
  //  Ouput:  Possibly two solutions for which the sample will result in detectable
  //          scatting peaks on the CCD.  Returns true if Bragg scatting may be detected,
  //          false otherwise.
  //
  //
  //-------------------------------------------------------------------------------------
  bool GetScatteringOmegas(Float &fOmega1, Float &fOmega2,
                           const SVector3 &oScatteringVec, const Float &fScatteringVecMag ) const;

  //-------------------------------------------------------------------------------------
  //
  //  GetObservablePeaks
  //
  //  Purpose:  Generate a vector of peak info list, specified by the current
  //            simulation setup.  To speed things up, a preallocated vector is passed
  //            in to retrive the results.
  //
  //  Comment:  a reserve will be placed on oResultPeakInfo to reduce memory allocation
  //            related slow down.  The cost is approximately <1% in runtime when compiled
  //            under optimized mode. (It's practically neglegible when compiled with the
  //            intel compiler. (This can be optimized
  //            by a one time allocation of oResultPeakInfo as a vector in the callee.)
  //-------------------------------------------------------------------------------------
  void GetObservablePeaks( vector<SPeakInfo> & oResultPeakInfo,
                           const SMatrix3x3 &oVoxelOrient,
                           const vector<CRecpVector> & oRecipVectors,
                           const DetectorListT & vDetectorList ) const;

  //-------------------------------------------------------------------------------------
  //  Default constructor -- use with caution!!!
  //-------------------------------------------------------------------------------------
  CSimulation():bSimulationInitialized( false ){};


  //-------------------------------------------------------------------------------------
  //  ProjectVertex
  //
  //  To be deprecated
  //-------------------------------------------------------------------------------------
  bool ProjectVertex( Point &p, const CDetector &oDetector, 
                      const CSample &oSample, const SVector3 &oVectex,
                      const SVector3 &oNormal ) const;

  //-------------------------------------------------------------------------------------
  //  ProjectVertex -- Generalized form
  //
  //  FPeakFilter is a function of the form 
  //  FPeakFilter: SVector3 -> bool
  //-------------------------------------------------------------------------------------
  template < class PeakFilterFn >
  bool ProjectVertex( Point &p, const CDetector &oDetector, 
                      const CSample &oSample, const SVector3 &oVectex,
                      const SVector3 &oNormal, PeakFilterFn FPeakFilter ) const;
public:

  CSimulation( const CXDMExperimentSetup & oExpParam ) { Initialize( oExpParam );  }
  //-------------------------------------------------------------------------------------
  //  THIS IS FOR GDB!!!  
  //
  //  TODO:  Remove this upon release!
  //-------------------------------------------------------------------------------------
  SMatrix3x3 DEBUG_T;
  SVector3   DEBUG_V;
  mutable Int DEBUG_I;

  void Initialize( const CXDMExperimentSetup & oExpParam )
  {
    bSimulationInitialized = true;
    oBeamDirection         = oExpParam.GetXrayBeamDirection();
    fBeamEnergy            = oExpParam.GetBeamEnergy();
    fBeamDeflectionChiLaue = oExpParam.GetBeamDeflectionChiLaue();
  }
  
  Float GetBeamEnergy() const { return fBeamEnergy; }
};

#include "Simulation.tmpl.cpp"

#endif
