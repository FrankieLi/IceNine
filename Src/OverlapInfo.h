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
//  File:    OverlapInfo.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Data structure of overlap counts and cost plus the implementation
//           of cost functions.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef OVERLAP_INFO_H_
#define OVERLAP_INFO_H_

#include "Simulation.h"
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/function.hpp>
#include "boost/multi_array.hpp"
#include "Quaternion.h"
#include "OrientationSearch.h"
#include "Sampling.h"
#include "Serializer.h"
#include "Symmetry.h"
#include <boost/tuple/tuple.hpp>
#include <limits>
#include "SearchDetails.h"
#include "PeakFilters.h"
namespace CostFunctions
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //  struct SOverlapInfo
  //
  //  some of this logging info should be disabled in high performance run
  //
  //  TODO:  Make this derivable -- specialization for each class.
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  struct SOverlapInfo
  {
    Int nPixelOverlap;          // number of pixels overlapping experimental data
    Int nPixelOnDetector;       // Number of simulated pixels lit on the detector
  
    Int nPeakOverlap;           // number of Bragg peaks *qualified* as overlapping data 
    Int nPeakOnDetector;        // number of Bragg peaks *qualified* as being on the detector
 
    //    Int nIndexablePeaks;        // for PathTraceSearch
    Int nDetectorsOverlap;      // number of detectors hit by this peak   (doesn't get incremented)

    Int nDetectorsHit;   // for other uses
  
    Float fQuality;             // Quality = 1 - cost;
    Float fBestPeakHit;         // DEBUG tmp added to debug
    
    Float nPeaksCounted;
  
    //----------------------
    //  Initialize all variables
    //----------------------
    void Initialize()
    {
      nPixelOverlap     = 0;
      nPixelOnDetector  = 0;
      nPeakOverlap      = 0;
      nPeakOnDetector   = 0;
      fQuality          = 0;
      nDetectorsOverlap = 0;

      nDetectorsHit     = 0;
    
      nPeaksCounted     = 0;
      //      fMeanAngularError = 0;
    }

    //----------------------
    //  An increment/ update function
    //----------------------
    void UpdateCounts( const SOverlapInfo & oRHS )
    {
      nPeakOverlap     += oRHS.nPeakOverlap;
      nPixelOverlap    += oRHS.nPixelOverlap;
      nPixelOnDetector += oRHS.nPixelOnDetector;
      nPeakOnDetector  += oRHS.nPeakOnDetector; 
      nPeaksCounted    +=  oRHS.nPeaksCounted;
    };
  
    //----------------------
    //  UpdateQuality
    //  nNumPoints is the total number of points before this point
    //----------------------
    void UpdateQuality( const SOverlapInfo & oRHS, Int nNumDetectors, Int &nNumPoints )
    {
      Float fCurQuality;
      if( oRHS.nPixelOnDetector > 0  )
      {
        if ( oRHS.nDetectorsOverlap > 0 )
        {
          Float fPixelOverlapRatio = Float( oRHS.nPixelOverlap )     / Float( oRHS.nPixelOnDetector );
          Float fDetRatio          = Float( oRHS.nDetectorsOverlap ) / Float ( nNumDetectors );
          fCurQuality              = fPixelOverlapRatio * fDetRatio;
        }
        else
        {
          fCurQuality = 0;
        }
        fQuality = fQuality +  ( fCurQuality - fQuality  ) / Float( nNumPoints + 1 );
        nNumPoints ++;
      }
    }


    //----------------------
    //  UpdateSmoothedQuality
    //  nNumPoints is the total number of points before this point
    //----------------------
    void UpdateQuality( const SOverlapInfo & oRHS, Int nNumDetectors,
                        Float fWeight, Float & TotalWeight )
    {
      Float fCurQuality;
      if( oRHS.nPixelOnDetector > 0  )
      {
        if ( oRHS.nDetectorsOverlap > 0 )
        {
          Float fPixelOverlapRatio = Float( oRHS.nPixelOverlap )     / Float( oRHS.nPixelOnDetector );
          Float fDetRatio          = Float( oRHS.nDetectorsOverlap ) / Float ( nNumDetectors );
          fCurQuality              = fPixelOverlapRatio * fDetRatio;
        }
        else
        {
          fCurQuality = 0;
        }
        fQuality = fQuality +  fWeight * ( fCurQuality - fQuality  ) / Float( TotalWeight + fWeight );
        TotalWeight += fWeight;
      }
    }

    
  };
  

  //---------------------------------------------------------------------
  //  Utilities
  //---------------------------------------------------------------------
  namespace Utilities
  {
    //-------------------------------------------------------------------------------------
    //  InitializeFlags
    //-------------------------------------------------------------------------------------
    void InitializeFlags( Bool bFlags[ MAX_NUM_DETECTORS ], Int nSize, Bool bInitValue );

    //-------------------------------------------------------------------------------------
    //  CountQualifiedPeaks
    //-------------------------------------------------------------------------------------
    void CountQualifiedPeaks( SOverlapInfo & oOverlapInfo, Bool bDetectorLit[ MAX_NUM_DETECTORS ],
                              Bool bSpotOverlap[ MAX_NUM_DETECTORS ], Size_Type nDetectors );
  } // end namespace Utilites


  
  //---------------------------------------------------------------------
  //  SampleVertexOverlapCounter
  //
  //  This object counts the overlap of diffraction spot generated by
  //  a set of sample vertices (indicating sample points, vertices of a
  //  triangle, or whate have you) in a set of metric defined
  //---------------------------------------------------------------------
  template< class PeakAcceptFn, class DetectorOverlapFn, class SimulatorFn >
  class XDMSampleVertexOverlapCounter
  {
  private:
    PeakAcceptFn      FPeakIntensityAccept;
    DetectorOverlapFn DetOverlapCounter;
    Float             fWaveNumber;
    const SimulatorFn & oSimulator;
    //-------------------------------------------------------------------------------------
    //  GenerateProjectedPixels
    //-------------------------------------------------------------------------------------
    void GenerateProjectedPixels( ProjectedPixelMapT & oProjectedPixels,
                                  PixelCountMapT     & oNumPixelHit,
                                  PixelIntensityMapT & oPixelIntensity,
                                  const vector<SPeakInfo> & vPeakInfo,
                                  CSample & oCurLayer,
                                  const vector<SVector3> & oVertexList,
                                  const DetectorListT & vDetectorList,
                                  const CSimulationRange &oRangeToIndexMap );
  public:
              
    XDMSampleVertexOverlapCounter( PeakAcceptFn      PeakFn,
                                   DetectorOverlapFn DetFn,
                                   const SimulatorFn   & oSim ):
      FPeakIntensityAccept( PeakFn ), DetOverlapCounter( DetFn ), oSimulator( oSim ),
      fWaveNumber( PhysicalConstants::keV_over_hbar_c_in_ang * oSim.GetBeamEnergy() ){}
    
    //------------------------------------------------------------------------
    //
    //  CalculateDiiffractionOverlap
    //
    //  Purpose:   Given a list of vertices, calculate the overlap produce
    //             by their diffraction at the specified diffraction peaks on
    //             the given detectors.
    //------------------------------------------------------------------------
    SOverlapInfo CalculateDiffractionOverlap( const vector<SVector3> &oVertexList,
                                              const vector<SPeakInfo> &vPeakInfo,
                                              CSample & oCurrentLayer,
                                              const DetectorListT & vDetectorList,
                                              const CSimulationRange &oRangeToIndexMap,
                                              const CSimulationData & oExpData );

    //------------------------------------------------------------------------
    //
    //  CalculateDiffractionOverlap
    //
    //  Purpose:   Given a list of vertices, calculate the overlap produce
    //             by their diffraction at the specified diffraction peaks on
    //             the given detectors.
    //------------------------------------------------------------------------
    SOverlapInfo CalculateDiffractionOverlapSmoothed( const vector<SVector3> &oVertexList,
                                                      const vector<SPeakInfo> &vPeakInfo,
                                                      CSample & oCurrentLayer,
                                                      const DetectorListT & vDetectorList,
                                                      const CSimulationRange &oRangeToIndexMap,
                                                      const CSimulationData & oExpData,
                                                      int DeltaOmega );
    
    //------------------------------------------------------------------------
    //  Accessor:
    //     Simulator
    //------------------------------------------------------------------------
    const SimulatorFn & Simulator() const { return oSimulator; }
  };
  
  //-------------------------------------------------------------------------------------
  //
  //  DetectorOverlapCounter
  //
  //  Purpose:  Get the approximate number of voxel overlapping between the simulation
  //            and the experiment.  To get this approximation, the center of voxel
  //            is projected onto the detector.
  //
  //  Return:   SOverlapInfo with the peak counts and peak generated fiels filled.
  //
  //
  //
  //  OverlapCounterFn  - must implement the function:
  //
  //  !!! OverlapCounterFn must NOT have any side effects
  //
  //  template<ImageT, PixelIter>
  //  void operator()( SOverlapInfo & oOverlapInfo, Bool & bDetectorLit, Bool & bSpotOverlap,
  //                   const ImageT & oImage, PixelIter pFirst, PixelIter pEnd )
  //                            
  //
  //-------------------------------------------------------------------------------------
  template< typename OverlapCounterFn >
  class DetectorOverlapCounter
  {
  private:
    OverlapCounterFn  CountFn;
    DetectorOverlapCounter( );        
  public:
    DetectorOverlapCounter( OverlapCounterFn _Fn ):
      CountFn( _Fn ) {}
    template< typename ImageList,
              typename ProjectedPixelList,
              typename PixelCountList,
              typename PixelIntensityList >
    SOverlapInfo operator()( ImageList           & oExpImageList,
                             ProjectedPixelList  & oProjectedPixels,
                             PixelCountList      & oNumPixelHit,
                             PixelIntensityList  & oPixelIntensity,
                             Int                   nNumPixels,
                             const DetectorListT & oDetectorList )
    {
      Bool bDetectorLit[ MAX_NUM_DETECTORS ];   // TODO:  Make this into class global vector that is initialized once
      Bool bSpotOverlap[ MAX_NUM_DETECTORS ];
      Utilities::InitializeFlags( bDetectorLit, oDetectorList.size(), false );
      Utilities::InitializeFlags( bSpotOverlap, oDetectorList.size(), false );
      
      SOverlapInfo oOverlapInfo;
      oOverlapInfo.Initialize();
      for( Size_Type nDet = 0; nDet < oDetectorList.size(); nDet++ )
        if( oNumPixelHit[ nDet ]  > 0 )
          CountFn( oOverlapInfo, bDetectorLit[ nDet ], bSpotOverlap[ nDet ],
                   oExpImageList[ nDet ], & oProjectedPixels[ nDet ][ 0 ],
                   ( & oProjectedPixels[ nDet ][ 0 ] ) + nNumPixels );        
      
      Utilities::CountQualifiedPeaks( oOverlapInfo, bDetectorLit, bSpotOverlap, oDetectorList.size() );
      return oOverlapInfo;
    }
  };
  
  //-------------------------------------------------------------------------------------
  //  PixelOverlapCount
  //
  //  Calculate overlap using pixels only.
  //
  //-------------------------------------------------------------------------------------
  class PixelBasedPeakOverlapCounter
  {
  private:
    PixelBasedPeakOverlapCounter();  // disabled
  public:
    Int nRadius;
    PixelBasedPeakOverlapCounter( Int n_ ): nRadius( n_ ){}

    template< typename ImageT >
    inline bool ProcessPoint( bool &bDetectorLit, bool & bSpotOverlap,
                              const ImageT & oImage, Float x, Float y ) const
    {
      if( oImage.IsInBound( x, y ) )
      {
        bDetectorLit = true;
        if( oImage.IsBright( x, y ) )
        {
          bSpotOverlap    = true;
          return false;
        }
      }
      return true;
    }

    //--------------------------------------
    //  operator()
    //--------------------------------------
    template< typename ImageT, typename PixelIter>
    inline void operator()( SOverlapInfo & oOverlapInfo,
                            Bool & bDetectorLit, Bool & bSpotOverlap,
                            const ImageT & oImage,
                            PixelIter pFirst, PixelIter pEnd ) const
    {
      const Point & oDetPoint = *pFirst;
      bool bContinueSearch = true;
      if( nRadius == 0 )            // Not sure if this speeds anything up
      {
        ProcessPoint( bDetectorLit, bSpotOverlap, oImage,
                      oDetPoint.x, oDetPoint.y );
        return;
      }
      for( Int nDX = -nRadius; nDX <= nRadius && bContinueSearch; nDX++ )
        for( Int nDY = -nRadius; nDY <= nRadius && bContinueSearch; nDY ++ )
          bContinueSearch = ProcessPoint( bDetectorLit, bSpotOverlap, oImage,
                                          oDetPoint.x + nDX, oDetPoint.y + nDY );
    }
  };

  //-------------------------------------------------------------------------------------
  //  SearchTreeBasedPeakOverlapCounter
  //  Calculate overlap using search tree
  //
  //  NOTE - this CLEARLY REQUIRES that the ImageT to be i) searchable and ii) initialized
  //         with search structure.
  //
  //  Performance note:  There is NO reason at all to use this if the vertices formed is
  //                     around the size of a pixel (or 1-2 pixel radius) in the detector
  //                     space.  This is built for cases where projected shapes are
  //                     much larger than pixels, i.e., 10s of pixels.
  //-------------------------------------------------------------------------------------
  struct SearchTreeBasedPeakOverlapCounter
  {
    template< typename ImageT, typename PixelIter>
    inline void operator()( SOverlapInfo & oOverlapInfo,
                            Bool & bDetectorLit, Bool & bSpotOverlap,
                            const ImageT & oImage,
                            PixelIter pFirst, PixelIter pEnd ) const
    {
      const Point & p1 = *  pFirst;
      const Point & p2 = *( pFirst + 1 );
      const Point & p3 = *( pFirst + 2 );
      // bound check
      DEBUG_ASSERT( pEnd != &p1, "PeakOverlapCount Memory: Vertex does not exist \n" );
      DEBUG_ASSERT( pEnd != &p2, "PeakOverlapCount Memory: Vertex does not exist \n" );
      DEBUG_ASSERT( pEnd != &p3, "PeakOverlapCount Memory: Vertex does not exist \n" );
      
      Int nOverlap = oImage.HasOverlap( p1, p2, p3 );
      if ( nOverlap > 0 )
      {
        bDetectorLit = true;
        bSpotOverlap = true;
      }
      else
      {
        if( oImage.IsInBound( p1.x, p1.y )
            || oImage.IsInBound( p2.x, p2.y )
            || oImage.IsInBound( p3.x, p3.y ) )
        {
          bDetectorLit = true;
        }
      }
    }
  };

  // TrianglePixelOverlapCounter

  template< class T >
  struct ShapePixelOverlapCounter
  {
    public:
    
    ShapePixelOverlapCounter()
    {
      std::cout << "Error - wrong specialization called" << std::endl;
      exit(0);
    }
    
    template< typename ImageT, typename PixelIter>
      inline void operator()( SOverlapInfo & oOverlapInfo,
			      Bool & bDetectorLit, Bool & bSpotOverlap,
			      const ImageT & oImage,
			      PixelIter pFirst, PixelIter pEnd ) const
    {
      std::cout << "Incorrect specialization" << std::endl;
      exit(0);
    }
    // general, empty ShapePixelOverlapCounter
  };
  //---------------------------------------------------
  //  Count the number of pixels overlapping given the triangle
  //  is projected onto the image
  //---------------------------------------------------
  template<>
  struct ShapePixelOverlapCounter<SVoxel>
  {
  public:
    template< typename ImageT, typename PixelIter>
    inline void operator()( SOverlapInfo & oOverlapInfo,
                            Bool & bDetectorLit, Bool & bSpotOverlap,
                            const ImageT & oImage,
                            PixelIter pFirst, PixelIter pEnd ) const
    {
      const Point & p1 = *  pFirst;
      const Point & p2 = *( pFirst + 1 );
      const Point & p3 = *( pFirst + 2 );
      // bound check
      DEBUG_ASSERT( pEnd != &p1, "TrianglePixelOverlapCounter Memory: Vertex does not exist \n" );
      DEBUG_ASSERT( pEnd != &p2, "TrianglePixelOverlapCounter Memory: Vertex does not exist \n" );
      DEBUG_ASSERT( pEnd != &p3, "TrianglePixelOverlapCounter Memory: Vertex does not exist \n" );
      
      Int nPixelOverlap  = 0;
      Int nPixelProduced = 0;                     // need to figure out pixels produced anyways

      //----------------------------------------------
      //  Explaination of the use of const_cast<>
      //   The use of const_cast here is because oImage is
      //   really a const object that should not be changed.
      //   However, thanks to C++'s obscurity, the const-ness
      //   of this object makes it impossible to call non-const
      //   functions.  Were this limited to non-const functions
      //   that modifies some precomputation members, we could
      //   just use the keyword "mutable" in the class.  However,
      //   the problem comes from a function (Rastersize) that
      //   should switch from const to non-const function based
      //   on the function object passed in.  In another words,
      //   when one is rasterizing to "fill" the image, the
      //   object should not be const.  On the other hand, a
      //   rasterize to "check" should only require a const
      //   reference.  This is where const_cast comes in to
      //   remove the "const-ness" of the object without (logically)
      //   disrupting all callers of TrianglePixelOverlapCounter().
      //----------------------------------------------
      boost::tie( nPixelOverlap, nPixelProduced )
        = const_cast< ImageT &>( oImage ).GetTriangleOverlapProperty( p1, p2, p3 );
      
      oOverlapInfo.nPixelOverlap    += nPixelOverlap;
      oOverlapInfo.nPixelOnDetector += nPixelProduced;
      
      if( nPixelOverlap > 0 )
      {
        bDetectorLit = true;
        bSpotOverlap = true;
      }
      else
      {
        if(  oImage.IsInBound( p1.x, p1.y )
             || oImage.IsInBound( p2.x, p2.y )
             || oImage.IsInBound( p3.x, p3.y ) )
        {
          bDetectorLit = true;
        }
      }
    }
  };
  typedef ShapePixelOverlapCounter<SVoxel> TrianglePixelOverlapCounter;
  //---------------------------------------------------
  //  Count the number of pixels overlapping given the Square
  //  is projected onto the image
  //---------------------------------------------------
  template<>
  struct ShapePixelOverlapCounter<SquareVoxel>
  {
    template< typename ImageT, typename PixelIter>
    inline void operator()( SOverlapInfo & oOverlapInfo,
                            Bool & bDetectorLit, Bool & bSpotOverlap,
                            const ImageT & oImage,
                            PixelIter pFirst, PixelIter pEnd ) const
    {
      const Point & p1 = *  pFirst;
      const Point & p2 = *( pFirst + 1 );
      const Point & p3 = *( pFirst + 2 );
      const Point & p4 = *( pFirst + 3 );
      // bound check
      DEBUG_ASSERT( pEnd != &p1, "TrianglePixelOverlapCounter Memory: Vertex does not exist \n" );
      DEBUG_ASSERT( pEnd != &p2, "TrianglePixelOverlapCounter Memory: Vertex does not exist \n" );
      DEBUG_ASSERT( pEnd != &p3, "TrianglePixelOverlapCounter Memory: Vertex does not exist \n" );
      DEBUG_ASSERT( pEnd != &p4, "TrianglePixelOverlapCounter Memory: Vertex does not exist \n" );
      
      Int nPixelOverlap  = 0;
      Int nPixelProduced = 0;                     // need to figure out pixels produced anyways

      //----------------------------------------------
      //  Explaination of the use of const_cast<>
      //   The use of const_cast here is because oImage is
      //   really a const object that should not be changed.
      //   However, thanks to C++'s obscurity, the const-ness
      //   of this object makes it impossible to call non-const
      //   functions.  Were this limited to non-const functions
      //   that modifies some precomputation members, we could
      //   just use the keyword "mutable" in the class.  However,
      //   the problem comes from a function (Rastersize) that
      //   should switch from const to non-const function based
      //   on the function object passed in.  In another words,
      //   when one is rasterizing to "fill" the image, the
      //   object should not be const.  On the other hand, a
      //   rasterize to "check" should only require a const
      //   reference.  This is where const_cast comes in to
      //   remove the "const-ness" of the object without (logically)
      //   disrupting all callers of TrianglePixelOverlapCounter().
      //----------------------------------------------
      boost::tie( nPixelOverlap, nPixelProduced )
        = const_cast< ImageT &>( oImage ).GetSquareOverlapProperty( p1, p2, p3, p4 );
      
      oOverlapInfo.nPixelOverlap    += nPixelOverlap;
      oOverlapInfo.nPixelOnDetector += nPixelProduced;
      
      if( nPixelOverlap > 0 )
      {
        bDetectorLit = true;
        bSpotOverlap = true;
      }
      else
      {
        if(  oImage.IsInBound( p1.x, p1.y )
             || oImage.IsInBound( p2.x, p2.y )
             || oImage.IsInBound( p3.x, p3.y )
             || oImage.IsInBound( p4.x, p4.y ))
        {
          bDetectorLit = true;
        }
      }
    }

  };

  
}// end namespace CostFunction

#include "OverlapInfo.tmpl.cpp"
#endif
