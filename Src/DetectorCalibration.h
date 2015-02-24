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
//
//   DetectorCalibration.h
//   Author:   Frankie Li  (li31@llnl.gov)
//
//   Purpose:  New, refactor detector calibration class.
//             
//
//
////////////////////////////////////////////////////////////
#ifndef DETECTOR_CALIBRATION_H_
#define DETECTOR_CALIBRATION_H_


#include "Reconstructor.h"
#include "DiscreteAdaptive.h"
#include <ctime>
#include <boost/shared_ptr.hpp>
#include "GrainReconstruction.h"
#include <Eigen/Dense>

namespace Calibration
{
  
  
  //-------------------------------------------
  //  class DetectorCalibration
  //        Given a region of interest, auto-magically figure out
  //        the geometric parameters of the detector based on
  //        some previously tested heuristics. 
  //-------------------------------------------
  class DetectorCalibration
  {

  private:
    typedef SVoxel SamplePointT;
    

    typedef Eigen::Matrix<Float, Eigen::Dynamic, 1> VectorN;
    typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> MatrixNM;
    
    
    //----------------------------------------------
    //    SVector3   oPosition;
    //    Float      fBeamCenterJ;
    //    Float      fBeamCenterK;
    //    SMatrix3x3 oOrientation
    //----------------------------------------------
    VectorN GetDetParam( const CDetector & Det )
    {
      VectorN v;
      
      return v;
    }

    //----------------------------------------------
    //----------------------------------------------
    void SetDetParam( CDetector & Det, const VectorN & v )
    {

    }

    
    //----------------------------------------------
    //  GetAggregateCost
    //----------------------------------------------

    //----------------------------------------------
    //  LocalOptimization
    //----------------------------------------------

  public:
    
    //----------------------------------------------
    //  
    //
    //
    //----------------------------------------------


    //----------------------------------------------
    
    
  };
}


#endif
