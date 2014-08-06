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
#ifndef _DETECTOR_FILE_H_
#define _DETECTOR_FILE_H_



#include <string>
#include "Detector.h"
#include "Parser.h"
#include <vector>
#include "BBox.h"


// forward declerationx
class CDetector;

namespace InitFileIO
{
  //--------------------------
  //
  //  A compact form of the detector
  //  parameters.  i.e., only the most
  //  essential of the geometry.   (This is also acting as a factory)
  //
  //  NOTE:     This sturcture is to be kept as a P.O.D.
  //--------------------------
  class CDetectorInfo
  {
  public:
    SVector3 vJUnitVector;
    SVector3 vKUnitVector;
    
    Float fBeamCenterJ;
    Float fBeamCenterK;
    
    SVector3   oLabFrameLocation;
    SMatrix3x3 oLabFrameOrientMatrix;
    SVector3   vLabFrameOrientation;   // Temporary hack so that the parameter optimization works
    
    Int nNumJPixels;
    Int nNumKPixels;
    
    Float fPixelHeight;
    Float fPixelWidth;
    
    Bool ParseDetectorBlock( const string & sBuf );
    CDetector GetDetector() const;
    void Print( ostream &os ) const;
  };
  
  //--------------------------
  //  CDetectorFile
  //--------------------------
  class CDetectorFile
  {
  private:
    
    //-----------------------------------------
    //  FindDetectorBlock
    //
    //  Find the number of { }  blocks.
    //  An error will be thrown when { isn't followed
    //  with an }, or vice versa
    //-----------------------------------------
    Bool FindDetectorBlock( vector<string> &vBlocks, const string & sBuf );
    
  public:
    Bool Parse( vector<CDetector> &oDetectorList, const string & oFilename );

    Bool ParseOptimizationInfo( vector<CDetectorInfo> &oDetStepSizeInfo, const string & oFilename );
  };

}

#endif
