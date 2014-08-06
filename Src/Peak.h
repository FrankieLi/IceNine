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
//   Peak.h
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:  Provide basic functionalities for peak manipulation, such as finding
//             max and min of a peak.
//
//
////////////////////////////////////////////////////////////

#ifndef CPEAK_H_
#define CPEAK_H_

#include <vector>
#include <algorithm>
#include <ostream>
#include "Types.h"
#include "Pixel.h"
#include "Serializer.h"
#include "BBox.h"
using std::ostream;
using std::vector;
using std::sort;
using std::endl;

using namespace GeneralLib;


//-------------------------------------------------------------------------------------
//
//  CLASS CPeak:  Base class of the container type Peak.
//
//-------------------------------------------------------------------------------------
class CPeak
{
private:
	       
  Float dMaxIntensity;
  Float dMinIntensity;
  bool PixelsSorted;
	
public:
  vector<Pixel> vPixelList;
  CPeak():dMaxIntensity(-1),dMinIntensity(-1) ,PixelsSorted(false), vPixelList() {};
  
  //-------------------------
  //  function: SetPixelThreshold(fraction)
  //            remove all peak below fraction of peak height
  //-------------------------
  void SetPixelThreshold(Float dIntensity);
  void FindRange();
  Float GetMaxIntensity() const;
  Float GetMinIntensity() const;
	
  void GetCenterOfIntensity(Float &row, Float &col);
	
  UInt NumPixel() const {return vPixelList.size();}
  void WritePeaks(ostream &os);
  void WritePeaksWithID(ostream &os, Int id);

  //-------------------------
  //  Save and Restore
  //-------------------------
  bool Save   ( CSerializer & oSerialBuf ) const;
  bool Restore( CDeserializer & oSerialBuf );
  
};


//-------------------------------------------------------------------------------------
//  PeakLookup
//-------------------------------------------------------------------------------------
class PeakLookup
{
public:
  vector<CPeak> vPeaks;
  
  inline void operator() (const Point &p, const CPeak &peak){
    vPeaks.push_back(peak);
  }
};


//-------------------------------------------------------------------------------------
//
//
//
//                      Class  CDetectorPeak
//
//
//
//-------------------------------------------------------------------------------------
class CDetectorPeak: public CPeak
{
private:
  bool bInitialized;
public:
  BBox2D oBBox;
  void CalculateBoundingBox();
  
  //-------------------------------------------------------------------------
  //
  //  CDetectorPeak::GetBoundingBox
  //
  //  Calculate the bounding box of this peak
  //
  //-------------------------------------------------------------------------
  inline BBox2D GetBoundingBox() const
  {
    return oBBox;
  }

  
  inline void GetBoundingBox(BBox2D & oBox ) const
  {
    oBox = oBBox;
  }

  inline void AddPointToBBox( const Pixel & p )
  {
    using namespace PBRMath;
    Point pTest( p.x, p.y );
    oBBox = Union( oBBox, pTest );
  }

  //-------------------------
  //  Save and Restore
  //-------------------------
  bool Save   ( CSerializer & oSerialBuf ) const;
  bool Restore( CDeserializer & oSerialBuf );
  
};


#endif
