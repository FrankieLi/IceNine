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
//  Simulation.tmpl.cpp
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//   Purpose:  Implementation of templated functions from Simulation.h
//
//
////////////////////////////////////////////////////////////



//------------------------------------------------------------------------------
//
//  ProjectVoxel  -- The fast AND generalized version
//
//  Precondition:  oDetectorList and oOutImageList have the same ordering:  i.e., 
//                 oOutImageList[n] <---- corresponds --->  oDetectorList[n]
//
//  Parameter:     FPeakFilter processes the filter
//
//  TODO:    Possible generalization of using PeakFilter for vertices.
//
//
//------------------------------------------------------------------------------
template <class OutputImageListT, class PeakFilterFn>
bool CSimulation::ProjectVoxel( OutputImageListT & oOutImageList,
                                const DetectorListT & oDetectorList,
                                const CSample &oSample, const SVoxel &oVoxel,
                                const SVector3 &oNormal,
                                PeakFilterFn FPeakAcceptFilter ) 
{
 
  DEBUG_ASSERT( oOutImageList.size() == oDetectorList.size(),
                  "ERROR:  CSimulation::ProjectVoxel: num of images != num detectors\n" );

  DEBUG_ASSERT( MAX_NUM_DETECTORS >  oDetectorList.size(),
                "CSimulation::ProjectVoxel:  Error:  MAX_NUM_DETECTORS <= oDectorList.size() \n ");

  //-----------------------------------------------------------------
  // TODO:  Remove these into a "buffer area" in the CSimulation space
  //-----------------------------------------------------------------
  Point        vProjectedPixels[ MAX_NUM_DETECTORS ][ NUM_VERTICES ];
  Size_Type    vNumPixelHit    [ MAX_NUM_DETECTORS ];
  Float        vPixelIntensity [ MAX_NUM_DETECTORS ];
  
  GetProjectedVertices( vProjectedPixels, vNumPixelHit, vPixelIntensity, oDetectorList,
                        oSample, &oVoxel.pVertex[0], &oVoxel.pVertex[0] + NUM_VERTICES,
                        oNormal, FPeakAcceptFilter );
  bool bDetectorHit =  false;   
  for ( Size_Type nDet = 0; nDet < oDetectorList.size(); nDet ++ )
  {
    if ( vNumPixelHit[ nDet ] == NUM_VERTICES )
    {
      for( Size_Type i = 0; i < NUM_VERTICES; i ++ )
        if ( oDetectorList[ nDet ].InRange( vProjectedPixels[ nDet ][i].x, vProjectedPixels[ nDet ][i].y ) )
          bDetectorHit = true;
      
      // should add more than just pixels  ( note intensity value is current flat )
      oOutImageList[ nDet ].AddTriangle( vProjectedPixels[ nDet ][ 0 ],
                                         vProjectedPixels[ nDet ][ 1 ],
                                         vProjectedPixels[ nDet ][ 2 ],
                                         vPixelIntensity [ nDet ] );
    }
    
  }
  return bDetectorHit; 
}
  
//------------------------------------------------------------------------------
//  GetProjectedVertices
//  Purpose:  Project Vertices based on given geometry.  Note that certain projections
//            will
//
//  Note that intensity calculation must be done here -- have to refactor or something.
//  NumPixelHit indicates how many pixel hitting a specific detector for this one peak.
//
//------------------------------------------------------------------------------
template < class ProjectedPixelMap, class PixelCountMap,
           class PixelIntensityMap, class VertexIter,
           class PeakFilterFn >
void CSimulation::GetProjectedVertices( ProjectedPixelMap & vProjectedPixels,
                                        PixelCountMap     & vNumPixelHit,
                                        PixelIntensityMap & vPixelIntensity,
                                        const DetectorListT & oDetectorList, const CSample &oSample,
                                        VertexIter pStart, VertexIter pEnd,  const SVector3 & oNormal,
                                        PeakFilterFn FPeakAcceptFilter ) const
{
  
  for ( Size_Type i = 0; i < oDetectorList.size(); i ++)
    vNumPixelHit[i] = 0;
  
  SVector3 oReflectedRayDir = DiffractionCore::GetReflectedRayDir( oSample, oNormal, oBeamDirection );   

  bool bAccept;
  Float fIntensity;
  boost::tie( bAccept, fIntensity ) = FPeakAcceptFilter( oReflectedRayDir ); 
  if( ! bAccept )
    return;
  
  for ( VertexIter pCur = pStart; pCur != pEnd; ++ pCur ) // for each vertex of the voxel
  {
    CRay oReflectedRay = DiffractionCore::BuildReflectedRay( oSample, *pCur, oReflectedRayDir );
    for ( Size_Type nDet = 0; nDet < oDetectorList.size(); nDet ++ )
    {
      Int nCurrentPoint = vNumPixelHit[ nDet ];
      Bool bIntersected = DiffractionCore::GetIlluminatedPixel( vProjectedPixels[ nDet ][ nCurrentPoint ],
                                                                oDetectorList[ nDet ] , oReflectedRay );
      if ( bIntersected )
      {
        vNumPixelHit[ nDet ] ++;
        vPixelIntensity[ nDet ] = fIntensity;
      }
      else
      {
        // short circuit -- if one pixel doesn't interset the detector plane
        // (planes are not clipped), than no pixel should intersect the
        // detector plane  [ WARNING -- this is NOT general -- move to CXDMReconstructSample!!!!!]
        // THIS ASSUMES THAT DETECTOR PLANES ARE PARALLEL!!!  
        
        return;
      }
    }  
  }
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
template < class PeakFilterFn >
bool CSimulation::ProjectVertex( Point &p, const CDetector &oDetector, 
                                 const CSample &oSample, const SVector3 &oVertex,
                                 const SVector3 &oNormal, PeakFilterFn FPeakFilter ) const
{
  CRay oReflectedRay = DiffractionCore::GetReflectedRay( oSample, oVertex, oNormal, oBeamDirection );
  Bool bAccept = FPeakFilter( oReflectedRay.GetDirection() );
  if( bAccept )
  {
    Bool bIntersected = DiffractionCore::GetIlluminatedPixel( p, oDetector, oReflectedRay );
    return bIntersected;
  }
  return false;
}
