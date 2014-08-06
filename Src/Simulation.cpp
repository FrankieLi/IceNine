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
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  File:   Simulation.cpp
//
//  Author: Frankie Li 
//  email:  sfli@cmu.edu
//
//  Purpose:  Implementation of basic simulation functions.  (Basic being
//            those that will be used by almost all types of X-ray
//            sychrotron simulation.)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Simulation.h"

//------------------------------------------------------------------------------
//
//  Private:  GetScatteringOmegas
//
//
//  Comment:  This function really just answers the question, "At what omega is the
//  bragg condition satistfied for a specific scattering vector.  This does NOT say
//  anything about the diffracted beam hitting the detector, and hence detectable.
//
//  Converted from Bob's program.
//
//  Hand optimized after profiling shows that 50% of the time is spent in this function
//
//------------------------------------------------------------------------------
bool CSimulation::GetScatteringOmegas( Float &fOmegaRes1, Float &fOmegaRes2,
                                       const SVector3 &oScatteringVec,
                                       const Float &fScatteringVecMag ) const
{
  ///
  ///  NOTE:  reciprical vectors are measured in angstrom                               
  ///  k  = 2pi/lambda  = E/ (h-bar c)                                                  
  ///
  
  Float fWavenumber =  PhysicalConstants:: keV_over_hbar_c_in_ang * fBeamEnergy;
  
  Float fSinTheta = fScatteringVecMag / ( (Float)2.0 * fWavenumber);   // Bragg angle
  Float fCosChi = oScatteringVec.m_fZ / fScatteringVecMag;             // Tilt angle of G relative to z-axis
  Float fSinChi = sqrt( (Float)1.0 - fCosChi * fCosChi );

  Float fSinChiLaue = sin( fBeamDeflectionChiLaue );     // ! Tilt angle of k_i (+ means up)
  Float fCosChiLaue = cos( fBeamDeflectionChiLaue );

  Float fNumerator = fSinTheta + fCosChi * fSinChiLaue;
  Float fDenom     = fSinChi * fCosChiLaue;

  if( fabs( fNumerator ) <= fabs( fDenom ) )
  {
    // [-pi:pi]: angle to bring G to nominal position along +y-axis
    Float fDeltaOmega0 = atan2( oScatteringVec.m_fX, oScatteringVec.m_fY);
    
    //  [0:pi/2] since arg >0: phi goes from above to Bragg angle
    Float fDeltaOmega_b1 = asin( fNumerator / fDenom );
    
    Float fDeltaOmega_b2 = PI -  fDeltaOmega_b1;
    
    fOmegaRes1 = fDeltaOmega_b1 + fDeltaOmega0;  // oScatteringVec.m_fY > 0
    fOmegaRes2 = fDeltaOmega_b2 + fDeltaOmega0;  // oScatteringVec.m_fY < 0
    
    if ( fOmegaRes1 > PI )          // range really doesn't matter
      fOmegaRes1 -=  2.f * PI;

    if ( fOmegaRes1 < -PI)
      fOmegaRes1 +=  2.f * PI;

    if ( fOmegaRes2 > PI)
      fOmegaRes2 -= 2.f * PI;
        
    if ( fOmegaRes2 < -PI)
      fOmegaRes2 += 2.f * PI;
    
    return true;
  }
  else
  {
    fOmegaRes1 = fOmegaRes2 = 0;     // too close to rotation axis to be illumination
    return false;
  }
  
}

//-------------------------------------------------------------------------------------
//
//  GetObserablePeaks
//
//  Purpose:  Generate a vector of peak info list, specified by the current
//            simulation setup.  To speed things up, a preallocated vector is passed
//            in to retrive the results.
//
//-------------------------------------------------------------------------------------
void CSimulation::GetObservablePeaks( vector<SPeakInfo> &vPeakInfo,
                                      const SMatrix3x3 &oVoxelOrientation,
                                      const vector<CRecpVector> & oRecipVectors,
                                      const DetectorListT & vDetectorList ) const
{
  vPeakInfo.reserve( oRecipVectors.size() * 2 );
  for( Size_Type nRecipIndex =0; nRecipIndex < oRecipVectors.size();  nRecipIndex++ )   
  {
    SVector3 oScatteringVec  = oRecipVectors[nRecipIndex].v; 	
    oScatteringVec = oVoxelOrientation * oScatteringVec;

    Float fOmegas[2];
    bool  bPeakObservable = GetScatteringOmegas( fOmegas[0], fOmegas[1], 
                                                 oScatteringVec,
                                                 oRecipVectors[nRecipIndex].fMag );   // This is VERY expensive
    if( bPeakObservable )
    {
      oScatteringVec.Normalize();
      vPeakInfo.push_back( SPeakInfo( bPeakObservable, oScatteringVec, fOmegas[0], oRecipVectors[nRecipIndex].fMag ) );
      vPeakInfo.push_back( SPeakInfo( bPeakObservable, oScatteringVec, fOmegas[1], oRecipVectors[nRecipIndex].fMag ) );
    }
  }
}


//------------------------------------------------------------------------------
//
//  ProjectVoxel  -- The fast version
//
//  Precondition:  oDetectorList and oOutImageList have the same ordering:  i.e., 
//                 oOutImageList[n] <---- corresponds --->  oDetectorList[n]
//------------------------------------------------------------------------------
bool CSimulation::ProjectVoxel( ImageListT & oOutImageList, const DetectorListT & oDetectorList,
                                const CSample &oSample, const SVoxel &oVoxel,
                                const SVector3 &oNormal, Float fIntensity )
{
  TrivialAcceptFilter.fIntensity = fIntensity;
  return ProjectVoxel( oOutImageList, oDetectorList, oSample, oVoxel, oNormal, TrivialAcceptFilter );
}

//------------------------------------------------------------------------------
//
//
//
//  CSimulation::ProjectVertex
//
//
//  Input:  oVertex  - input vertex
//          oNormal  - normal of the surface
//          oDetector - image plane of the projection
//
//
//  Output:  p  - point will be filled with the location ImageCol, ImageRow, which
//                is a point on the image plane, oDetector
//  return false if projection was not successful, i.e., point does not lie on the plane; true otherwise
//
//------------------------------------------------------------------------------
bool CSimulation::ProjectVertex( Point &p, const CDetector &oDetector, 
                                 const CSample &oSample, const SVector3 &oVertex,
                                 const SVector3 &oNormal) const
{
  CRay oReflectedRay = DiffractionCore::GetReflectedRay( oSample, oVertex, oNormal, oBeamDirection );
  Bool bIntersected  = DiffractionCore::GetIlluminatedPixel( p, oDetector, oReflectedRay );
  return bIntersected;
}

//------------------------------------------------------------------------------
//
//
// Public:
// CSimulation::ProjectVoxel
//
// Input:  oSample, oCurrnetVoxel, oNormal
// Output:  oDetector - with pixels filled where voxel is projected
//
// Purpose:  determine location of the detector that is lit by the x-ray diffraction from the voxel
//
//------------------------------------------------------------------------------
bool CSimulation::ProjectVoxel( CImageData & oOutputImage, const CDetector &oDetector,
                                const CSample &oSample, const SVoxel &oVoxel,
                                const SVector3 &oNormal, Float fIntensity )
{
  Point ProjectedPixels[3];	
  Int nNumPoints = 0;

  // for each vertex of the voxel
  for(UInt j = 0; j < 3; j ++)
  {  
    if( ProjectVertex( ProjectedPixels[ nNumPoints ], oDetector, oSample, oVoxel.pVertex[j], oNormal) )
      nNumPoints ++;
    else
      return false;
  }
	
  // should add more than just pixels
  oOutputImage.AddTriangle( ProjectedPixels[0], ProjectedPixels[1], ProjectedPixels[2], fIntensity );
  return true;
}

//------------------------------------------------------------------------
//
//  CXDMSimulationState
//
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//  Public:  CXDMSimulationState::~CXDMSimulationState
//------------------------------------------------------------------------
CXDMSimulationState::~CXDMSimulationState()
{
  if( oLogFile.is_open() )
    oLogFile.close();
}

//------------------------------------------------------------------------
//
//  Public:  CXDMSimulationState::InitializeLogFile
//
//  Initialize the log file of the simulation state
//  
//------------------------------------------------------------------------
bool CXDMSimulationState::InitializeLogFile( string sFilename )
{
  if ( oLogFile.is_open() )
    return false;
  
  oLogFile.open( sFilename.c_str() );
  return ( oLogFile.good() );
}

//------------------------------------------------------------------------
//
//  Public:  CXDMSimulationState::Print
//
//  Print (dump) the current state of the simulation into the log file
//
//  TODO:  Figure out a better format of the log file?
//  
//------------------------------------------------------------------------
void CXDMSimulationState::PrintToLog()
{
  //  for ( Size_Type i = 0; i < vCurrentPixelList.size(); i ++ )

  if ( vCurrentPixelList.size() > 0 ) 
    oLogFile << nVoxelID << " "
             << nRecipVecID << " "
             << fOmega   << " "
             << oScatteringVector.m_fX << " " << oScatteringVector.m_fY
             << " " << oScatteringVector.m_fZ
             << vCurrentPixelList[0].x << " "
             << vCurrentPixelList[0].y << endl;
}




