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
//  DiffractionCore.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//  Purpose:  This contains all of utilities needed for the simulations
//            of diffractions
//
////////////////////////////////////////////////////////////

#ifndef DIFFRACTION_CORE_H
#define DIFFRACTION_CORE_H

#include "3dMath.h"
#include "Sample.h"
#include "Detector.h"
#include "Voxel.h"
#include "SearchDetails.h"
namespace DiffractionCore
{
  using namespace GeneralLib;

  //--------------------------------------
  // SPeakInfo
  //--------------------------------------
  struct SPeakInfo
  {
    Bool     bObservable;
    SVector3 oScatteringDir;
    Float    fOmega;
    Float    fRecipVecMag;
    SPeakInfo( Bool b,
               const SVector3 & v,
               Float f, Float m ):
      bObservable( b ),
      oScatteringDir( v ),
      fOmega( f ),
      fRecipVecMag( m ) {}
  };

  
  //-----------------------------------------------------------------------------------
  // I N L I N E   F U N C T I O N S
  //-----------------------------------------------------------------------------------
  
  //-------------------------------------------------------------------------------------
  //
  //  GetReflectionVector
  //  Input:   Incident vector and the normal vector
  //  Output:  The reflected vector (based on geometric optics)
  //
  //-------------------------------------------------------------------------------------
  inline SVector3 GetReflectionVector( const SVector3 &R_in, 
                                       const SVector3 &normal )
  {
    // What we want is:
    // R_out = R_in - 2 * Dot(R_in, normal) * normal;
    
    SVector3 R_out;
    R_out = R_in - (Float) (2.0) * Dot( R_in, normal ) * normal;
    return R_out;
  }

  //------------------------------------------------------------------------------
  //
  //  GetReflectedRayDir 
  //
  //  Returns the direction (i.e., unit vector) of the reflected ray given a normal.
  //  Note that the oNormal is in the reference frame of oSample, and the reflected ray
  //  will be in the reference frame of the lab.
  //
  //  The purpose of this function is to allow decoupling of reflection vector
  //  calculation, which is redundantly called for each of the vertices of each voxel.
  //  Since each of the vertex is assumed to have the same normal vector (constant interpolation),
  //  there is no reason to recalculate the same reflection 3 times.
  //
  //------------------------------------------------------------------------------
  inline SVector3 GetReflectedRayDir( const CSample &oSample, const SVector3& oNormal,
                                      const SVector3 & oBeamDir )
  {
    SVector3 oLabNormal = oSample.ToLabFrame( oNormal ) ;
    SVector3 oOutDir = GetReflectionVector( oBeamDir, oLabNormal );
    return oOutDir;
  }

  //------------------------------------------------------------------------------
  //
  //  BuildReflectedRay
  //
  //  Returns a reflected ray given the sample, vertex, and the reflection direction.
  //  Note that this simply means that a ray is drawn from the vertex in the lab frame
  //  along the direction of the reflected ray.  Note that the reflected ray will be
  //  in the reference frame of the lab.
  //
  //------------------------------------------------------------------------------
  inline CRay BuildReflectedRay( const CSample &oSample, const SVector3 &oVertex,
                                 const SVector3 &oRefDir)
  {
    CRay oReflectedRay;                                  // construct a Ray
    SVector3 oStartPoint = oSample.ToLabFrame( oVertex );    // convert starting point to lab frame
    oReflectedRay.SetDirection( oRefDir );
    oReflectedRay.SetStart( oStartPoint );
    return oReflectedRay;
  }

  //------------------------------------------------------------------------------
  //  Private:  GetReflectedRay
  //            A helper function that calculates the reflected vector given a
  //            sample and the vertex.
  //            A Reflected ray from originated from the oVertex will be computed.
  //  NOTE:     Both oVertex and oNormal will be rotated into the sample frame.
  //
  //------------------------------------------------------------------------------
  inline CRay GetReflectedRay( const CSample &oSample, const SVector3 &oVertex,
                               const SVector3 &oNormal, const SVector3 & oBeamDir )
  {  
    SVector3 oRefDir   = GetReflectedRayDir( oSample, oNormal, oBeamDir );
    CRay oReflectedRay = BuildReflectedRay ( oSample, oVertex, oRefDir);
    return oReflectedRay;
  }
  

  //------------------------------------------------------------------------------
  // Private:   GetIlluiminatedPixel
  //            A helper function that returns the point where the detector is lit
  //            by an incident ray, oIncidentRay.  (Incident upon the detector.)
  //            Note that intersection will be calculated in this function.
  //------------------------------------------------------------------------------
  inline Bool GetIlluminatedPixel( Point & p, const CDetector &oDetector,
                                   const CRay & oIncidentRay )
  {
    Float fT;
    if( oDetector.Intersects( oIncidentRay, fT ) )      
    {
      Float ImageRow, ImageCol;
      SVector3 oIntersectLoc = oIncidentRay.Evaluate( fT );   // Get location of intersection in lab frame
      oDetector.LabToPixel( ImageRow, ImageCol, oIntersectLoc );    
      p.x = ImageCol;
      p.y = ImageRow;
      return true;
    }
    return false;
  }
  
} // namespace Physics

  //--------------------------------------------------
  //  HEDM specific calculations
  //
  //  Purpose:  Shortcuts and approximations used in HEDM
  //
  //  W A R N I N G
  //  - Some/most of these functions are very coordinate
  //    specific.  Unless you're implementing some kind of
  //    new algorithm for HEDM, do NOT use this.
  //--------------------------------------------------
namespace HEDM
{
  //------------------------------------------------------------------------
  // InterpolateScatteringOmegas
  //
  //
  // WARNING:  This is a coordinate system dependent function!!!
  //
  //
  // WARNING:  This function operates under the assumption that fGLabX is the x-coordinate
  //           of *SOME* G vector in the lab frame.
  //
  // IF one chooses to have RecpVec to be a G vector in the lab frame, then
  // oTransformation = Omega * Delta * Omega^-1
  //
  // OTHERWISE - if G vector is in the sample frame, then oTransformation = Omega * Delta
  //
  // Bottom line is that we have to preserve the following equation:
  //
  //   GL' = dOmega * Omega * Delta * Omega^-1 GL 
  //
  //
  //
  // We are assuming that we're doing small rotation of the reciprocal vector.
  // This means that the new reciprocal vector, GL' will have the following properties:
  //
  //  dot( GL', K_in ) = dot( GL, K_in ) 
  //  GL' = dOmega * Omega * Delta * Omega^-1 GL 
  //
  //  where dOmega is an infinitesimal rotation matrix around the z axis.
  //
  //  TODO:  How to interpolate omegas that wasn't there to begin with?
  //
  // WARNING:  This is an interpolation function, so use it as such!
  //
  // Error of this function grows rapidly as fTop/fBottom -> 0.5
  // This is essentially an inverse tangent function, so be wary of the sign.
  // ( Note that there's a sign correction section commented out )
  //------------------------------------------------------------------------
  inline bool InterpolateScatteringOmegas(Float & fOmegaRes,
                                          const SVector3 & oRecpVec,
                                          const SMatrix3x3 & sMat,
                                          Float fGLabX )
  {
    Float fTop =   sMat.m[0][0] * oRecpVec.m_fX - fGLabX
      + sMat.m[0][1] * oRecpVec.m_fY
      + sMat.m[0][2] * oRecpVec.m_fZ;
      
    Float fBottom = sMat.m[1][0] * oRecpVec.m_fX
      + sMat.m[1][1] * oRecpVec.m_fY
      + sMat.m[1][2] * oRecpVec.m_fZ;
      
    fOmegaRes = fTop / fBottom;
      
    if( fabs( fOmegaRes ) > 1 ) 
      return false;
      
    return true;
  }
  
  //--------------------------------------------------------------------------------------------------------
  //  PrecomputeScatteringVector
  //
  //  Purpose:  A short cut to computing all of the scattering vectors for each of the
  //            orientational sample point.  This method locally interpolates 
  //
  //--------------------------------------------------------------------------------------------------------
  template< typename DetectorListT, typename SimRangeToIndexMapT,
            typename ReciprocalVecIter, typename SimulatorT >
  inline void PrecomputeScatteringVector( SearchDetails::ScatteringVectorListT  & ResultScatteringVectors,
                                          const SMatrix3x3          & oFZRotMatrix,
                                          const DetectorListT       & vDetectorList,
                                          const SimRangeToIndexMapT & oRangeToIndexMap,
                                          ReciprocalVecIter           pRecipLatticeVector,
                                          ReciprocalVecIter           pEnd,
                                          SimulatorT                  Simulator )
  {
    ResultScatteringVectors.reserve( pEnd - pRecipLatticeVector );

    // DEBUG  (SPECIFIC TO THIS EXPERIMENT)
    SVector3 oRotAxis( 0, 0, 1 );
    // END DEBUG

    for( Size_Type n = 0 ; pRecipLatticeVector != pEnd;  ++pRecipLatticeVector, ++ n )   
    {
      SVector3 oRecipVec  = pRecipLatticeVector->v; 	
      oRecipVec.Transform( oFZRotMatrix );
        
      Float fOmegaRes[2];
      bool  bPeakObservable = Simulator.GetScatteringOmegas( fOmegaRes[0], fOmegaRes[1],
                                                             oRecipVec,
                                                             pRecipLatticeVector->fMag);   // This is VERY expensive
  
      ResultScatteringVectors[n].bObserved = bPeakObservable;
      if( bPeakObservable )
      {
        SearchDetails::SScatteringVectorEntry oRes;
        
        SVector3 oRecipVecDir = oRecipVec;
        oRes.fReciprocalVecMag = oRecipVec.GetLength();
        oRecipVecDir.Normalize();
        oRes.oReciprocalVecDir = oRecipVecDir;
        oRes.vObservedOmega[0] = fOmegaRes[0];
        oRes.vObservedOmega[1] = fOmegaRes[1];
          
        oRes.vOmegaRotationMatrix[0].BuildActiveEulerMatrix( fOmegaRes[0], 0, 0 );
        oRes.vOmegaRotationMatrix[1].BuildActiveEulerMatrix( fOmegaRes[1], 0, 0 );

        oRes.vLabRecpVecDir[0] = oRes.vOmegaRotationMatrix[0] * oRecipVecDir;
        oRes.vLabRecpVecDir[1] = oRes.vOmegaRotationMatrix[1] * oRecipVecDir;
        ResultScatteringVectors.push_back( oRes );
      }
    }
  }
    
  //-------------------------------------------------------------------------------------
  //
  //  InterpolateObservablePeaks
  //
  //  Purpose:  Generate a vector of peak info list, specified by the current
  //            simulation setup.  To speed things up, a preallocated vector is passed
  //            in to retrive the results.
  //  Comment:  a reserve will be placed on oResultPeakInfo to reduce memory allocation
  //            related slow down.  The cost is approximately <1% in runtime when compiled
  //            under optimized mode. (It's practically neglegible when compiled with the
  //            intel compiler. (This can be optimized
  //            by a one time allocation of oResultPeakInfo as a vector in the callee.)
  //
  //-------------------------------------------------------------------------------------
  template< typename PeakListT>
  inline void InterpolateObservablePeaks( PeakListT & oResultPeakInfo,
                                   const SearchDetails::ScatteringVectorListT & CandidateLatticeVectors,
                                   const SMatrix3x3 & oSmallRotation )
  {
    oResultPeakInfo.reserve( CandidateLatticeVectors.size() * 2 );
    for( Size_Type n = 0; n <  CandidateLatticeVectors.size(); n ++  )   
    {
      SVector3 oSampleRecpVecDir  = CandidateLatticeVectors[ n ].oReciprocalVecDir;
      oSampleRecpVecDir = oSmallRotation * oSampleRecpVecDir;
        
      for( Int i = 0; i < 2; i ++ )
      {
        // run interpolation with G_sample, or reciprocal lattice vector in sample frame
        const SMatrix3x3 & oZRotation = CandidateLatticeVectors[ n ].vOmegaRotationMatrix[i];   // use reference - no modify needed
          
        Float fDeltaOmega;
        SVector3 vLabRecpVecDir = CandidateLatticeVectors[ n ].vLabRecpVecDir[i];   // using S = Omega with Sample Frame vector (can be saved)
        bool bPeakObservable = InterpolateScatteringOmegas( fDeltaOmega, oSampleRecpVecDir,
                                                            oZRotation, vLabRecpVecDir.m_fX );  // (very cheap)
        Float fOmega = CandidateLatticeVectors[ n ].vObservedOmega[i];
        Float fOmegaNew = fOmega + fDeltaOmega;
          
        if( bPeakObservable )
          oResultPeakInfo.push_back( DiffractionCore::SPeakInfo( bPeakObservable, oSampleRecpVecDir,
                                                                 fOmegaNew, CandidateLatticeVectors[ n ].fReciprocalVecMag ) );
      }
    }
  }
    
} // ----------------------------- Namespace HEDM -----------------------------



#endif
