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
//   Peak.cpp
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:  Implementation of Peak.h
//
//
////////////////////////////////////////////////////////////
#include "Peak.h"




//-------------------------------------------------------------------------
//
//  CPeak::FindRange
//            Find both max and min pixel intensity
//            of the pixel list
//
//-------------------------------------------------------------------------
void CPeak::FindRange()
{
	
	
  if( vPixelList.size() == 1){
    dMaxIntensity = vPixelList[0].fIntensity;
    dMinIntensity = vPixelList[0].fIntensity;
    return;
  }
	
  dMaxIntensity = MIN_FLOAT;
  dMinIntensity = MAX_FLOAT;


  for(UInt i = 0; i < vPixelList.size() -1; i += 2)
  {
    if( vPixelList[i].fIntensity > vPixelList[i+1].fIntensity )
    {
      if(vPixelList[i].fIntensity > dMaxIntensity)
        dMaxIntensity = vPixelList[i].fIntensity;
			
      if(vPixelList[i+1].fIntensity < dMinIntensity)
        dMinIntensity = vPixelList[i+1].fIntensity;
		
    }
    else
    {
			
      if(vPixelList[i+1].fIntensity > dMaxIntensity)
        dMaxIntensity = vPixelList[i+1].fIntensity;
			
      if(vPixelList[i].fIntensity < dMinIntensity)
        dMinIntensity = vPixelList[i].fIntensity;
    }
	
  }
  
}

//-------------------------------------------------------------------------
// function: SetPixelThreshold
//
//           Set threshold of intensity 
//           to be considered as part
//           of the peak
//-------------------------------------------------------------------------
void CPeak::SetPixelThreshold(Float dIntensity)
{
	
  if(dMaxIntensity < 0)  // intensity should really be positive
  {
    FindRange();
  }
	

  // sort in assending order
  if ( !PixelsSorted )
  {
    sort(vPixelList.begin(), vPixelList.end());
  }
	
  //
  // Find first element that meets the threshold
  // criteria
  //

  UInt i;
  for(i = 0; i < vPixelList.size(); i ++ )
  {
    if(vPixelList[i].fIntensity > dIntensity)
      break;
    
  }
	
  if( i > 0 )
  {
    vPixelList.erase(vPixelList.begin(), vPixelList.begin() + i - 1);
  }
	
	
}

//-------------------------------------------------------------------------
// GetMaxIntensity
//-------------------------------------------------------------------------
Float CPeak::GetMaxIntensity() const 
{
  return dMaxIntensity;
}

//-------------------------------------------------------------------------
//  GetMinIntensity
//-------------------------------------------------------------------------
Float CPeak::GetMinIntensity() const 
{
  return dMinIntensity;
}

//-------------------------------------------------------------------------
// GetCenterOfIntensity
//
// Note - a better averaging scheme may be needed, as we lose
// precision when suming a large number, like Intensity.
//
//-------------------------------------------------------------------------
void CPeak::GetCenterOfIntensity(Float &x, Float &y)
{
  Float xI = 0;
  Float yI = 0;
  Float intensitySum = 0;
	
  for( UInt i = 0; i < vPixelList.size(); i ++ )
  {
    xI += vPixelList[i].x * vPixelList[i].fIntensity;
    yI += vPixelList[i].y * vPixelList[i].fIntensity;
    intensitySum += vPixelList[i].fIntensity;
  }
  
  x = xI / intensitySum;
  y = yI / intensitySum;
	
}

//-------------------------------------------------------------------------
//  WritePeaks
//-------------------------------------------------------------------------
void CPeak::WritePeaks(ostream &os)
{
  for( UInt i = 0; i < vPixelList.size(); i++ )
  {
    os << vPixelList[i].x << " " 
       << vPixelList[i].y << " " 
       << vPixelList[i].fIntensity << endl;
  }
}

//-------------------------------------------------------------------------
//  WritePeaksWithID
//-------------------------------------------------------------------------
void CPeak::WritePeaksWithID(ostream &os, Int id)
{
  for( UInt i = 0; i < vPixelList.size(); i++ )
  {
    os << vPixelList[i].x << " " 
       << vPixelList[i].y << " " 
       << vPixelList[i].fIntensity << " "
       << id << endl;
  }
}


//-------------------------------------------------------------------------
//  Save
//-------------------------------------------------------------------------
bool CPeak::Save   ( CSerializer & oSerialBuf ) const
{
  bool bSuccess;
  bSuccess = oSerialBuf.InsertCompactObj( dMaxIntensity );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( dMinIntensity );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( PixelsSorted );
  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( vPixelList );
  DEBUG_ASSERT( bSuccess, "CPeak::Save failed\n");
  return bSuccess;
}

//-------------------------------------------------------------------------
//  Restore
//-------------------------------------------------------------------------
bool CPeak::Restore( CDeserializer & oSerialBuf )
{
  bool bSuccess;
  vPixelList.clear();
  
  bSuccess = oSerialBuf.GetCompactObj( &dMaxIntensity );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &dMinIntensity );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &PixelsSorted );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( vPixelList );
  DEBUG_ASSERT( bSuccess, "CPeak::Restore failed\n");
  return bSuccess;
}


///////////////////////////////////////////////////////////////////////////
//
//
//
//                      Class  CDetectorPeak
//
//
//
///////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
//
//  CDetectorPeak::GetBoundingBox
//
//  Calculate the bounding box of this peak
//
//-------------------------------------------------------------------------
void CDetectorPeak::CalculateBoundingBox()
{
  Int maxX = 0;
  Int maxY = 0;
  Int minX = MAX_INT;
  Int minY = MAX_INT;

  for(Size_Type i = 0; i < vPixelList.size(); i ++)
  {
    if( vPixelList[i].x > maxX )
      maxX = vPixelList[i].x;
    if(vPixelList[i].y > maxY)
      maxY = vPixelList[i].y;

    if(vPixelList[i].x < minX)
      minX = vPixelList[i].x;
    if(vPixelList[i].y < minY)
      minY = vPixelList[i].y;
		
  }
  oBBox.pMax.x = (Float)maxX + EPSILON;
  oBBox.pMax.y = (Float)maxY + EPSILON;
  
  oBBox.pMin.x = (Float)minX - EPSILON;
  oBBox.pMin.y = (Float)minY - EPSILON;
  
}

//--------------------------------------------------
//  Save 
//--------------------------------------------------
bool CDetectorPeak::Save   ( CSerializer & oSerialBuf ) const
{
  bool bSuccess;
  bSuccess = CPeak::Save( oSerialBuf );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( bInitialized );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( oBBox );

  DEBUG_ASSERT( bSuccess, "CDetectorPeak::Save failed\n");
  return bSuccess;
}

//--------------------------------------------------
//  Restore
//--------------------------------------------------
bool CDetectorPeak::Restore( CDeserializer & oSerialBuf )
{
  bool bSuccess;
  bSuccess = CPeak::Restore( oSerialBuf );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &bInitialized );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &oBBox );
  
  DEBUG_ASSERT( bSuccess, "CDetectorPeak::Restore failed\n");
  return bSuccess;
}
