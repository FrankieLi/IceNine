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
//  File:    OverlapInfo.tmpl.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Implementation of cost functions
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace CostFunctions
{

  namespace Utilities
  {
    //-------------------------------------------------------------------------------------
    //  InitializeFlags
    //-------------------------------------------------------------------------------------
    inline void InitializeFlags( Bool bFlags[ MAX_NUM_DETECTORS ], Int nSize, Bool bInitValue )
    {
      for( Int i = 0; i < nSize; i ++ )
        bFlags[ i ] = bInitValue;
    }

    //-------------------------------------------------------------------------------------
    //  CountQualifiedPeaks
    //
    //  TODO: Potentially change this function into a bitwise function instead of a bool array.  
    //  The consequence is that we won't need a for loop to check.  We can simply
    //  do bitwise operation to check for patterns.  The *only* valid patterns of
    //  valid peaks and valid overlap is:
    //
    //  1 1 1 1 1 ... 1
    //  1 1 1 1 1 ... 0
    //  1 1 1 1 0 ... 0
    //
    //  The number of zeros allowed should be a parameter.  (This sets the number of candidate
    //  to use for more detailed searches. )  (NOTE:  Profile this first.)
    //-------------------------------------------------------------------------------------
    inline void CountQualifiedPeaks( SOverlapInfo & oOverlapInfo,
                                     Bool bDetectorLit[ MAX_NUM_DETECTORS ],
                                     Bool bSpotOverlap[ MAX_NUM_DETECTORS ],
                                     Size_Type nDetectors )
    {
      Bool bValidSimPeak       = bDetectorLit[0];
      Bool bOnPreviousDetector = bDetectorLit[0];
      Int  nDetectorsOverlap   = 0;

      for( Size_Type i = 1; i < nDetectors && bValidSimPeak; i ++ )
      {
        if( bDetectorLit[i] )
        {
          bValidSimPeak = true;
        }
        else
        {
          if( bOnPreviousDetector )
          {
            bValidSimPeak = true;
            bOnPreviousDetector = false;
          }
          else
          {
            bValidSimPeak = false;
          }
        }
      }
  
      Bool bPreviousOverlap = bSpotOverlap[0];
      Bool bValidOverlap    = bSpotOverlap[0];

      if( bSpotOverlap[0] )
        nDetectorsOverlap ++;
  
      for( Size_Type i = 1; i < nDetectors && bValidOverlap; i ++ )
      {
        if( bSpotOverlap[i] )
        {
          nDetectorsOverlap ++;
          bValidOverlap = true;
        }
        else
        {
          if( bPreviousOverlap && ! bDetectorLit[i] )
          {
            bValidOverlap = true;
            bPreviousOverlap = false;
          }
          else
          {
            bValidOverlap = false;
          }
        }
      }
  
      if( bValidSimPeak )
        oOverlapInfo.nPeakOnDetector = 1;
      if( bValidOverlap )
        oOverlapInfo.nPeakOverlap    = 1;
  
      if( bValidSimPeak && bValidOverlap )
        oOverlapInfo.nDetectorsOverlap = nDetectorsOverlap;
    }
  }  // namespace Utilities


  //--------------------------------------------------------------------------------------------
  
  //------------------------------------------------------------------------
  //
  //  CalculateDiffractionOverlap
  //------------------------------------------------------------------------
  template< class PeakAcceptFn, class DetectorOverlapFn, class SimulatorFn >
  SOverlapInfo XDMSampleVertexOverlapCounter<PeakAcceptFn, DetectorOverlapFn, SimulatorFn>::
  CalculateDiffractionOverlap( const vector<SVector3> &oVertexList,
                               const vector<SPeakInfo> &vPeakInfo,
                               CSample & oCurrentLayer,
                               const DetectorListT & vDetectorList,
                               const CSimulationRange &oRangeToIndexMap,
                               const CSimulationData & oExpData )
  {
    const Int nNumPixels    = oVertexList.size();
    const Int nNumPeaks     = vPeakInfo.size();
    const Int nNumDetectors = vDetectorList.size();
    ProjectedPixelMapT oProjectedPixels( boost::extents[ nNumPeaks ][ nNumDetectors ][ nNumPixels ] );
    PixelCountMapT     oNumPixelHit    ( boost::extents[ nNumPeaks ][ nNumDetectors ] );
    PixelIntensityMapT oPixelIntensity ( boost::extents[ nNumPeaks ][ nNumDetectors ] );
  
    GenerateProjectedPixels( oProjectedPixels, oNumPixelHit, oPixelIntensity,
                             vPeakInfo, oCurrentLayer, oVertexList, vDetectorList, oRangeToIndexMap );
    SOverlapInfo oOverlapInfo;
    oOverlapInfo.Initialize();
    Int nPoints = 0;
    for( Size_Type nPeakNum = 0; nPeakNum < vPeakInfo.size(); nPeakNum ++ )
    {
      typedef boost::multi_array_types::index_range range;
      ConstProjectedPixelView  oPixelView
        = oProjectedPixels[ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ][ range( 0, nNumPixels )] ]; 
      ConstPixelCountView     oHitView
        = oNumPixelHit    [ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ] ];
      ConstPixelIntensityView oIntensityView
        = oPixelIntensity [ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ] ];

      Size_Type nOmegaIndex = oRangeToIndexMap( vPeakInfo[ nPeakNum ].fOmega );
      
      if ( nOmegaIndex != XDMSimulation::NoMatch )
      {
        typedef typename CSimulation::Const_ImageListT Const_ImageListT;
        Const_ImageListT oImList = oExpData.mImageMap[ boost::indices [ nOmegaIndex ][ range( 0, vDetectorList.size() ) ] ];

        SOverlapInfo oCurOverlapInfo =
          DetOverlapCounter( oImList, oPixelView, oHitView, oIntensityView, nNumPixels, vDetectorList );  
        oOverlapInfo.UpdateCounts ( oCurOverlapInfo );
        oOverlapInfo.UpdateQuality( oCurOverlapInfo, vDetectorList.size(), nPoints );   // nPoints is currently a hack
      }
    }
    return oOverlapInfo;
  }
  
  //------------------------------------------------------------------------
  //
  //  CalculateDiffractionOverlapSmoothed
  //------------------------------------------------------------------------
  template< class PeakAcceptFn, class DetectorOverlapFn, class SimulatorFn >
  SOverlapInfo XDMSampleVertexOverlapCounter<PeakAcceptFn, DetectorOverlapFn, SimulatorFn>::
  CalculateDiffractionOverlapSmoothed( const vector<SVector3> &oVertexList,
                                       const vector<SPeakInfo> &vPeakInfo,
                                       CSample & oCurrentLayer,
                                       const DetectorListT & vDetectorList,
                                       const CSimulationRange &oRangeToIndexMap,
                                       const CSimulationData & oExpData,
                                       int DeltaOmega )
  {
    const Int nNumPixels    = oVertexList.size();
    const Int nNumPeaks     = vPeakInfo.size();
    const Int nNumDetectors = vDetectorList.size();
    ProjectedPixelMapT oProjectedPixels( boost::extents[ nNumPeaks ][ nNumDetectors ][ nNumPixels ] );
    PixelCountMapT     oNumPixelHit    ( boost::extents[ nNumPeaks ][ nNumDetectors ] );
    PixelIntensityMapT oPixelIntensity ( boost::extents[ nNumPeaks ][ nNumDetectors ] );
  
    GenerateProjectedPixels( oProjectedPixels, oNumPixelHit, oPixelIntensity,
                             vPeakInfo, oCurrentLayer, oVertexList, vDetectorList, oRangeToIndexMap );
    SOverlapInfo oOverlapInfo;
    oOverlapInfo.Initialize();
    Float TotalWeight = 0;
    for( Size_Type nPeakNum = 0; nPeakNum < vPeakInfo.size(); nPeakNum ++ )
    {
      typedef boost::multi_array_types::index_range range;
      ConstProjectedPixelView  oPixelView
        = oProjectedPixels[ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ][ range( 0, nNumPixels )] ]; 
      ConstPixelCountView     oHitView
        = oNumPixelHit    [ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ] ];
      ConstPixelIntensityView oIntensityView
        = oPixelIntensity [ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ] ];

      Size_Type nOmegaIndex = oRangeToIndexMap( vPeakInfo[ nPeakNum ].fOmega );

      if ( nOmegaIndex != XDMSimulation::NoMatch )
      {
        for( int DeltaW = -DeltaOmega; DeltaW <= DeltaOmega; DeltaW ++ )
        {
          int nCurIndex = static_cast<int>( nOmegaIndex ) + DeltaW;
          if( oRangeToIndexMap.ValidIndex( nCurIndex ) )
          {
            typedef typename CSimulation::Const_ImageListT Const_ImageListT;
            Const_ImageListT oImList = oExpData.mImageMap[ boost::indices [ nCurIndex ][ range( 0, vDetectorList.size() ) ] ];
            
            SOverlapInfo oCurOverlapInfo =
              DetOverlapCounter( oImList, oPixelView, oHitView, oIntensityView, nNumPixels, vDetectorList );
            if( DeltaW == 0 )
              oOverlapInfo.UpdateCounts ( oCurOverlapInfo );

            if( DeltaW == 0 )
              oOverlapInfo.UpdateQuality( oCurOverlapInfo, vDetectorList.size(),
                                          1, TotalWeight );   // nPoints is currently a hack
            oOverlapInfo.UpdateQuality( oCurOverlapInfo, vDetectorList.size(),
                                        0.1,
                                        TotalWeight );   // nPoints is currently a hack
          }
        }
      }
    }
    return oOverlapInfo;
  }

  
  //------------------------------------------------------------------------
  //  GenerateProjectedPixels
  //
  //  Private
  //  Purpose:  A generalized function that generates projected pixels given
  //            a set of vertices.  The pixels are just projection of these
  //            vertices onto the detector based on ONE reciprocal lattice vector.
  //            (i.e., the normal, which is given in vPeakInfo)
  //
  //  Note:     Pixels generated are of the form:  
  //            [ nNumDetectors ][ nNumPixels ]
  //
  //            Also returned are oNumPixelHit, which is the number of pixel hit per detector,
  //            and oPixelIntensity, which is the intensity of the peak on the detector.
  //            Note that this function assumes a "flat shading model," in that
  //            the intensity is constant across one pixel and/or polygon delineated
  //            by the pixel.
  //
  //            Because this is a generalized function, it could take multiple vertex
  //            as an argument.  In that case, pixels must be specified in the counterclockwise
  //            sense.
  //
  //  TODO:     Generalize this to take FPeakIntensityAccept as a general function object.
  //
  //------------------------------------------------------------------------
  template< class PeakAcceptFn, class DetectorOverlapFn, class SimulatorFn >
  void XDMSampleVertexOverlapCounter<PeakAcceptFn, DetectorOverlapFn, SimulatorFn>::
  GenerateProjectedPixels( ProjectedPixelMapT & oProjectedPixels,
                           PixelCountMapT     & oNumPixelHit,
                           PixelIntensityMapT & oPixelIntensity,
                           const vector<SPeakInfo> & vPeakInfo,
                           CSample & oCurLayer,
                           const vector<SVector3> &oVertexList,
                           const DetectorListT & vDetectorList,
                           const CSimulationRange &oRangeToIndexMap )
  {
    const Int nNumPixels    = oVertexList.size();
    const Int nNumDetectors = vDetectorList.size();
    
    const SMatrix3x3 oCurOrientation = oCurLayer.GetOrientationMatrix();
    for( Size_Type nPeakNum = 0; nPeakNum < vPeakInfo.size(); nPeakNum ++ )
    {
      if( vPeakInfo[ nPeakNum ].bObservable )
      {
        Size_Type nOmegaIndex = oRangeToIndexMap( vPeakInfo[ nPeakNum ].fOmega );
        if ( nOmegaIndex != XDMSimulation::NoMatch )
        {
          oCurLayer.RotateZ( vPeakInfo[ nPeakNum ].fOmega );

          typedef boost::multi_array_types::index_range range;
          
          ProjectedPixelMapT::array_view<2>::type  oPixelView
            = oProjectedPixels[ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ][ range( 0, nNumPixels )] ]; 
          
          boost::multi_array< Size_Type, 3 >::array_view<1>::type  oHitView
            = oNumPixelHit    [ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ] ];
          
          boost::multi_array< Float,     2 >::array_view<1>::type  oIntensityView
            = oPixelIntensity[ boost::indices[ nPeakNum ][ range( 0, nNumDetectors ) ] ];
          
          const SVector3 & oNormal = vPeakInfo[ nPeakNum ].oScatteringDir;

          Float fSinTheta = vPeakInfo[ nPeakNum ].fRecipVecMag / ( Float(2) * fWaveNumber );
          FPeakIntensityAccept.SetSin2Theta( sin( Float( 2 ) * asin( fSinTheta ) ) );    // change this function object to save x,y points instead
          oSimulator.GetProjectedVertices( oPixelView, oHitView, oIntensityView,
                                           vDetectorList, oCurLayer,
                                           oVertexList.begin(), oVertexList.end(),
                                           oNormal, FPeakIntensityAccept );
          oCurLayer.SetOrientation( oCurOrientation );  // use this to reduce numerical errors
        }
      }
    }
  }
  
} // namespace CostFunctions
