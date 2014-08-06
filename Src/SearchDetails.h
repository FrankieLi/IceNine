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
//------------------------------------------------------------------------
//   SearchDetails
//
//------------------------------------------------------------------------

#ifndef SEARCH_DETAILS_H
#define SEARCH_DETAILS_H
#include "3dMath.h"
#include "Types.h"
#include <vector>
#include "ConfigFile.h"
namespace SearchDetails
{
  using namespace GeneralLib;
  //------------------------------------------------------------------------
  //  SScatteringVectorEntry
  //
  //  This is the "scattering vector" in the lab frame which is used to
  //  calculate the elastic X-ray diffraction.
  //------------------------------------------------------------------------
  struct SScatteringVectorEntry
  {
    SVector3   oReciprocalVecDir;
    Float      fReciprocalVecMag; 
    Float      vObservedOmega[2];
    Bool       bObserved;
    SMatrix3x3 vOmegaRotationMatrix[2];
    SVector3   vLabRecpVecDir[2];  
  };

  
  //------------------------------------------------------------------------
  //  AcceptanceCriteron
  //  Knobs for acceptance - for testing and debugging purpose
  //
  //------------------------------------------------------------------------
  struct SAcceptCriterion
  {
    Float fConfidence; // nPeakOverlap/nPeakOnDetector
    Int nPeakOverlap;
    Float fOverlapRatio;
    Int nPixelOverlap;
    Int nReciprocalVectorHit;
  };
  
  typedef std::vector<SScatteringVectorEntry>  ScatteringVectorListT;


  //------------------------------------------------------------------------
  //  A candidate in the searh
  //------------------------------------------------------------------------
  struct SCandidate
  {
    SMatrix3x3 oOrientation; // change to quaternion?
    Float fCost;             // initial cost of the system
    
    SCandidate() {};
    SCandidate( const SMatrix3x3 & oMat, Float f )
      : oOrientation( oMat ), fCost( f ) {};
    
    //---------------------------
    //  less than
    //---------------------------
    Bool operator< ( const SCandidate & oRHS) const
    {
      return fCost < oRHS.fCost; 
    };
  };

  //  typedef std::pair< SSearchPoint, Int > CandidatePair;
  
  //---------------------------
  //  SSearchParameter
  //
  //  Private tuning parameters for
  //  tuning of the search.
  //---------------------------
  struct SSearchParameter
  {
    //------------------------------------------------------------------------------
    //  Fittign function     Cost function is ranged betwee [0, 1)
    //------------------------------------------------------------------------------
    Float fMaxAcceptedCost;      // Maximum cost where a fitting voxel is considered good
    Float fMaxConvergenceCost;   // Maximum cost where an element is considered converged
    Float fMaxDeepeningHitRatio; // Maximum hit ratio below which iterative deepening will
    // kick in.
    //------------------------------------------------------------------------------
    //  Randomize Best First Search
    //------------------------------------------------------------------------------
    Float fMaxBestFirstCost;    // Maximum cost to trigger best first search
      
    Int   nLocalResolution;
    Float fLocalGridRadius;     // L_infinity radius of local grid, i.e., its diameter

    Int   nMaxMCSteps;          // Maximum allowed monte carlo steps for orientation optimization
      
    Float fMCRadiusScaleFactor; //  A fudge factor used to scale the radius of
    //  the Monte Carlo step.  This is to be taken out
    //  once further testing tells us more about the
    //  search algorithm.
      
    Int   nSuccessiveRestarts;  //  A fudge factor.  We really should be calculating this
    //  based on search space volume information.
  };
  
  enum ReconstructionStatus{
    CONVERGED, PARTIAL, NON_ACCEPTED, NONE
  };



  namespace Utilities
  {
    void InitializeLocalSearchParameter( SSearchParameter & oSearchParam,
                                         const CConfigFile oSetupFile );
  }


  
} // end namespace SearchDetails

#endif
