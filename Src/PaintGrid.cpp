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
//--------------------------------------------------------------------------------------------------------
//  File:   ForwardSimulation.cpp
//
//
//--------------------------------------------------------------------------------------------------------

#include "PaintGrid.h"


//--------------------------------------------------------------------------------------------------------
//  Public:
//  CPaintGrid::CPaintGrid(const string & sConfigFilename)
//  Constructor
//
//  input:  sConfigFilename - contains file name of the config file with experimental setup and simulation
//                           geometry
// 
//--------------------------------------------------------------------------------------------------------
CPaintGrid::CPaintGrid( const CConfigFile & oConfigFile, int nThreads ) 
{
  NumThreads = nThreads;
  // omp_set_num_threads( NumThreads );

  RUNTIME_ASSERT(0, "omp_set_num_threads not working!");
  exit(0);
  oSetupFile = oConfigFile;                  // to be removed

  //--------------------------------------
  time_t oStartTime, oStopTime;
  time ( &oStartTime );
  //--------------------------------------
  std::cout << "Begin InitializeWithDataFiles " << std::endl;
  oSetup.InitializeWithDataFiles( oConfigFile );
  oSimulator.Initialize( oSetup.ExperimentalSetup() );

  time ( &oStopTime );
  std::cout << "Time elapsed for initialization (input) (sec): "
            << difftime( oStopTime, oStartTime ) << std::endl;
}

//--------------------------------------------------------------------------------------------------------
//
//  This is temporary
//
//
//--------------------------------------------------------------------------------------------------------
void CPaintGrid::ReadOrientationTestPoints( const std::string & OrientationFilename,
                                            std::vector<SMatrix3x3 > & FZList )
{

  //--------------------------------------
  time_t oStartTime, oStopTime;
  time ( &oStartTime );
  //--------------------------------------
  std::cout << "Begin Reading Orientation File " << std::endl;
  
  std::ifstream OrientationFile ;
  OrientationFile.open( OrientationFilename.c_str() );
  if( ! OrientationFile )
  {
    std::cout << "Cannot open orientation file" << std::endl;
    exit(0);
  }
  
  std::vector<GeneralLib::SVector3> RFList;
  do
  {
    GeneralLib::SVector3 Tmp;
    OrientationFile >> std::skipws >> Tmp.m_fX
                    >> std::skipws >> Tmp.m_fY
                    >> std::skipws >> Tmp.m_fZ;
    if( OrientationFile.good() )
      RFList.push_back( Tmp );
 
  }  
  while( OrientationFile.good() );
  time ( &oStopTime );
  std::cout << "Time elapsed for reading orientation grid (input) (sec): "
            << difftime( oStopTime, oStartTime ) << std::endl;
  time ( &oStartTime );
  FZList.resize( RFList.size() );
#pragma omp parallel for
#pragma omp nowait
  for( int i = 0; i < RFList.size(); i ++)
    {
      GeneralLib::SQuaternion q;

      RUNTIME_ASSERT(0, "PAINTGRID IS BROKEN RIGHT NOW\n");
      //     q.CreateFromRodrigues( RFList[i] );
      FZList[i] =  q.GetRotationMatrix3x3();
    }


  time ( &oStopTime );
  std::cout << "Time elapsed for conversion (input) (sec): "
            << difftime( oStopTime, oStartTime ) << std::endl;
}

//--------------------------------------------------------------------------------------------------------
// ConstructTwoThetaRange
//--------------------------------------------------------------------------------------------------------
std::vector<SRange> CPaintGrid::ConstructTwoThetaRanges( Float ThetaMapWidth, Float ThetaWidth,
                                                         const vector<CRecpVector> & oRecipVectors )
{
  Float HalfWidth = ThetaWidth / static_cast<Float>( 2 );
  std::vector<SRange> TwoThetaRangeList;
  Float fWavenumber =  PhysicalConstants:: keV_over_hbar_c_in_ang * oSetup.ExperimentalSetup().GetBeamEnergy();  

  std::vector<Float> TwoThetaList;
  for( int nRecipIndex = 0; nRecipIndex < oRecipVectors.size();  nRecipIndex++ )   
  {
    Float fSinTheta = oRecipVectors[nRecipIndex].fMag / (Float(2.0) * fWavenumber );
    Float f2Theta   =  Float( 2 ) * asin( fSinTheta );
    TwoThetaList.push_back( f2Theta );
  }
  
  // take care of overlaps

  std::sort( TwoThetaList.begin(), TwoThetaList.end() );

  std::vector<SRange> RegularizedTwoThetaRanges;
  Float fMax = *std::max_element( TwoThetaList.begin(), TwoThetaList.end() ) + HalfWidth;
  Float fMin = *std::min_element( TwoThetaList.begin(), TwoThetaList.end() ) - HalfWidth;

  int NumIntervals  = (fMax - fMin) / ThetaMapWidth;

  std::cout << "Num Intervals " << NumIntervals << std::endl;
  
  std::vector<int> Mask( NumIntervals, 0 );
  for( int i = 0; i < TwoThetaList.size(); i ++ )
  {

    Float fThetaLow  = TwoThetaList[i] - HalfWidth;
    Float fThetaHigh = TwoThetaList[i] + HalfWidth;
    
    
    int nStart = floor( (fThetaLow  - fMin) / ThetaMapWidth  );
    int nStop  = floor( (fThetaHigh - fMin) / ThetaMapWidth );
    for( int n = nStart; n <= nStop; n ++  )
    {
      Float fBinLow  = n       * ThetaMapWidth + fMin;
      Float fBinHigh = (n + 1) * ThetaMapWidth + fMin;
      if( ( fBinHigh >= fThetaLow )  && ( fBinLow < fThetaHigh ) )
        Mask[n] = 1;
    }
  }
  
  for( int i = 0; i < Mask.size(); i ++ )
  {
    if( Mask[i] > 0 )
    {
      Float fBinLow  = i       * ThetaMapWidth + fMin;
      Float fBinHigh = (i + 1) * ThetaMapWidth + fMin;
      RegularizedTwoThetaRanges.push_back(  SRange( fBinLow, fBinHigh ) );
      std::cout << fBinLow << " " << fBinHigh << std::endl;
    }
  }

  TwoThetaRangeToIndexMap.Set( fMin, fMax, ThetaMapWidth, RegularizedTwoThetaRanges ); 
  return RegularizedTwoThetaRanges;
}

//--------------------------------------------------------------------------------------------------------
//  InitializeIntensityMap
//--------------------------------------------------------------------------------------------------------
void CPaintGrid::InitializeIntensityMap( Float EtaBinSize,
                                         Float f2ThetaBinSize, Float f2ThetaPeakWidth,
                                         const vector<CRecpVector> & oRecipVectors )
{
  Float MinEta = -PI;
  Float MaxEta =  PI;
  int NumEtaIntervals = std::ceil( ( MaxEta - MinEta ) / EtaBinSize );
  EtaBinSize = ( MaxEta - MinEta ) / NumEtaIntervals;                 // Readjust bin size to fit.
  
  std::vector<SRange> TwoThetaRanges  = ConstructTwoThetaRanges( f2ThetaBinSize,
                                                                 f2ThetaPeakWidth,
                                                                 oRecipVectors );
  
  OmegaRangeToIndexMap = oSetup.ExperimentalSetup().GetRangeToIndexMap();
  
  IntensityMap.resize( boost::extents[ OmegaRangeToIndexMap.GetNumIntervals() ]
                       [ NumEtaIntervals  ]
                       [ TwoThetaRanges.size() ] );

  
  for( int i = 0; i < OmegaRangeToIndexMap.GetNumIntervals(); i ++ )
    for( int j = 0; j < NumEtaIntervals; j ++ )
      for( int k = 0; k < TwoThetaRanges.size(); k ++ )
        IntensityMap[i][j][k] = 0;
  
  
  std::vector<SRange> EtaRanges;
  for( int i = 0; i < NumEtaIntervals; i ++ )
    EtaRanges.push_back( SRange( static_cast<int>(i) * EtaBinSize + MinEta,
                                 static_cast<int>(i + 1) * EtaBinSize + MinEta ) );
  
  EtaRangeToIndexMap.Set( MinEta, MaxEta,
                          EtaBinSize, EtaRanges ); 
}

//--------------------------------------------------------------------------------------------------------
//  DialateTwoTheta
//--------------------------------------------------------------------------------------------------------
void CPaintGrid::DilateIntensityMap( Float fOmegaWidth, Float fEtaWidth, Float f2ThetaWidth )
{

  Float fOmegaHW  = fOmegaWidth / static_cast<Float>( 2.0 );
  Float fEtaHW    = fEtaWidth   / static_cast<Float>( 2.0 );
  Float f2ThetaHW = f2ThetaWidth / static_cast<Float>( 2.0 );

  int MaxNumDilation = 0;
  for( int i = 0; i < OmegaRangeToIndexMap.GetNumIntervals(); i ++ )
  {
    for( int j = 0; j < EtaRangeToIndexMap.GetNumIntervals(); j ++ )
    {
      for( int k = 0; k < TwoThetaRangeToIndexMap.GetNumIntervals(); k ++ )
      {
        if( IntensityMap[i][j][k] > 0 )
        {
          Float fOmegaCenter  = OmegaRangeToIndexMap.IndexToIntervalCenter   ( i );
          Float fEtaCenter    = EtaRangeToIndexMap.IndexToIntervalCenter     ( j );
          Float f2ThetaCenter = TwoThetaRangeToIndexMap.IndexToIntervalCenter( k );

          Float fOmegaLow = fOmegaCenter - fOmegaHW;
          Float fOmegaHi  = fOmegaCenter + fOmegaHW;

          Float fEtaLow = fEtaCenter - fEtaHW;
          Float fEtaHi  = fEtaCenter + fEtaHW;
          
          Float f2ThetaLow = f2ThetaCenter - f2ThetaHW;
          Float f2ThetaHi  = f2ThetaCenter + f2ThetaHW;
          
          typedef CSimulationRange::IndexIter IndexIter;
          std::pair<IndexIter, IndexIter> OmegaRange    = OmegaRangeToIndexMap   ( fOmegaLow,  fOmegaHi  );
          std::pair<IndexIter, IndexIter> EtaRange      = EtaRangeToIndexMap     ( fEtaLow,    fEtaHi    );
          std::pair<IndexIter, IndexIter> TwoThetaRange = TwoThetaRangeToIndexMap( f2ThetaLow, f2ThetaHi );
          int NumDilation = 0;
          for( IndexIter OmegaCur = OmegaRange.first; OmegaCur != OmegaRange.second; ++ OmegaCur )
          {
            for( IndexIter EtaCur = EtaRange.first; EtaCur != EtaRange.second; ++ EtaCur )
            {
              for( IndexIter ThetaCur = TwoThetaRange.first; ThetaCur != TwoThetaRange.second; ++ ThetaCur )
              {
                if( *OmegaCur != XDMSimulation::NoMatch
                    && *EtaCur != XDMSimulation::NoMatch
                    && *ThetaCur != XDMSimulation::NoMatch )
                {
                  IntensityMap[ *OmegaCur ][ *EtaCur ][ *ThetaCur ] ++;
                  NumDilation ++;
                }
              }
            }
          }
          MaxNumDilation = std::max( NumDilation, MaxNumDilation );
        }
      }
    }
  }
  std::cout << "Maximum number of pixels used for dilation: " << MaxNumDilation << std::endl;
}

//--------------------------------------------------------------------------------------------------------
//
//  CPaintGrid:: SimulateDetectorImagesOptimized
//
//--------------------------------------------------------------------------------------------------------
bool CPaintGrid::Run( const std::string & PaintGridConfig  )
{

  Float OmegaDilationWidth = DEGREE_TO_RADIAN( 1 );
  Float EtaDilationWidth   = DEGREE_TO_RADIAN( 1 );
  Float ThetaDilationWidth = DEGREE_TO_RADIAN( 1 );
  

  Float EtaBinSize   = DEGREE_TO_RADIAN( 0.1 );
  Float ThetaBinSize = DEGREE_TO_RADIAN( 0.1 );
  
  
  std::string OrientationGridFile;
  //  std::string OutputFilename;
  
  //-----------------------------------------------------------------
  //   WARNING:  This is NOT parsed!
  std::ifstream ConfigFile;
  ConfigFile.open( PaintGridConfig.c_str() );
  RUNTIME_ASSERT( ConfigFile.good(), " PaintGrid config file cannot be opened" );

  ConfigFile >> std::skipws >> OrientationGridFile;

  ConfigFile >> std::skipws >> OmegaDilationWidth
             >> std::skipws >> EtaDilationWidth
             >> std::skipws >> ThetaDilationWidth;
    

  ConfigFile >> std::skipws >> EtaBinSize
             >> std::skipws >> ThetaBinSize;

  //  ConfigFile >> std::skipws >> OutputFilename;
  
  std::cout << " WARNING - PaintGrid Config is not a parsed file!  Check output!" << std::endl;
  std::cout << " " << OrientationGridFile << std::endl;

  std::cout << " " << OmegaDilationWidth 
            << " " << EtaDilationWidth
            << " " << ThetaDilationWidth  << std::endl;
  
  
  std::cout << " " << EtaBinSize
            << " " << ThetaBinSize  << std::endl;
  //  std::cout << " " <<  OutputFilename << std::endl;
  
  OmegaDilationWidth = DEGREE_TO_RADIAN( OmegaDilationWidth );
  EtaDilationWidth   = DEGREE_TO_RADIAN( EtaDilationWidth );
  ThetaDilationWidth = DEGREE_TO_RADIAN( ThetaDilationWidth ); 

  EtaBinSize   = DEGREE_TO_RADIAN( EtaBinSize );
  ThetaBinSize = DEGREE_TO_RADIAN( ThetaBinSize ); 

  ConfigFile.close();

  
  //-----------------------------------------------------------------

  
  CSample oCurrentLayer;
  Initialize( oSetup.ExperimentalSetup() );
  
  vector<SMatrix3x3> FZList;
  ReadOrientationTestPoints( OrientationGridFile,  FZList );  //  WARNING:  Temporary hard code
  
  
  const vector<SRange>    & vOmegaRangeList  = oSetup.ExperimentalSetup().GetOmegaRangeList();
  const vector<CDetector> & oDetectorList    = oSetup.ExperimentalSetup().GetDetectorList();
  

  oSetup.ExperimentalSetup().InitializeSample( oCurrentLayer,  oDetectorList[ 0 ] );  // binding CurrentLayer to the first detector's limit

  std::shared_ptr< CMic > pMic= std::dynamic_pointer_cast<CMic>( oCurrentLayer.GetMic() );

  
  RUNTIME_ASSERT( pMic->VoxelListBegin() != pMic->VoxelListEnd(),
                  "ERROR:  Mic file is empty!");
  
  Int nCryStructIndex = pMic->VoxelListBegin()->nPhase;
  const vector<CUnitCell> & oCryStructList  = oCurrentLayer.GetStructureList();
  const vector<CRecpVector> & oRecipVectors = oCryStructList[ nCryStructIndex ].GetReflectionVectorList(); 

  InitializeIntensityMap( EtaBinSize,
                          ThetaBinSize,
                          ThetaDilationWidth,
                          oRecipVectors );
  //--------------------------------------
  time_t oStartTime, oStopTime;
  time ( &oStartTime );
  std::cout << "Begin Construct EtaOmega Map " << std::endl;
  ConstructEtaOmegaMap( oDetectorList,
                        oCurrentLayer );
  time ( &oStopTime );
  std::cout << "Time elapsed for ConstructEtaOmegaMap (sec): "
            << difftime( oStopTime, oStartTime ) << std::endl;

  time ( &oStartTime );
  DilateIntensityMap( OmegaDilationWidth,
                      EtaDilationWidth,
                      ThetaDilationWidth );
  time ( &oStopTime );
  std::cout << "Time elapsed for Dilation (sec): "
            << difftime( oStopTime, oStartTime ) << std::endl;
  
  std::cout << "Begin PaintGrid with " << FZList.size() << " points in the FZList " <<  std::endl;
  time ( &oStartTime );
  std::vector<Float> Intensities = Paint( oDetectorList, FZList, oCurrentLayer );
  time ( &oStopTime );
  std::cout << "Finished PainGrid " << std::endl;
  std::cout << "Time elapsed for PaintGrid (sec): "
            << difftime( oStopTime, oStartTime ) << std::endl;

  //-----------------------------------------------
  std::ofstream Output;
  std::stringstream ss;
  ss << oSetup.ExperimentalSetup().GetInputParameters().OutStructureBasename
     << ".PaintGridIntensity.txt";
  Output.open( ss.str().c_str() );
  for( int i = 0; i < Intensities.size(); i ++ )
    Output << Intensities[i] << std::endl;
  Output.close();
  //-----------------------------------------------
  
  return true;
}

//------------------------------------------------------
//  GetEtaTwoThetaBoundingBox
//------------------------------------------------------
void CPaintGrid::PixelToEtaTwoTheta( Float & Eta, Float & TwoTheta,
                                     Float x, Float y,
                                     const CDetector & Detector,
                                     const SVector3 & Center ) const
{
  GeneralLib::SVector3 ReflectedDir = Detector.PixelToLabCoordinate( x, y ) - Center;
  
  Eta      = atan2( ReflectedDir.m_fZ, -ReflectedDir.m_fY );
  TwoTheta = fabs( acos( ReflectedDir.m_fX / ReflectedDir.GetLength() ) );
}

//------------------------------------------------------
//  ConstructEtaOmegaMap
//
//   - Construct the eta-omega map by iterating through all detector images
//
//------------------------------------------------------
void CPaintGrid::ConstructEtaOmegaMap( const DetectorListT & vDetectorList,
                                       CSample oCurrentLayer )
{
  const CSimulationData & Data = oSetup.Data();

  int TotalNumPeaks = 0;
  int NumIntensityMapMarked = 0;
  
  RUNTIME_ASSERT( Data.nNumDetectors >= 1, "Need at least one detector!");
  std::cout << "Number of detector intervals: " << Data.nNumIntervals << std::endl;
  for( int i = 0; i < Data.nNumIntervals; i ++ )  // assumed to have single detector
  {
    const SMatrix3x3 oCurOrientation = oCurrentLayer.GetOrientationMatrix();
    SRange fTest = OmegaRangeToIndexMap.IndexToInterval( i );
    RUNTIME_ASSERT( fabs( OmegaRangeToIndexMap.IndexToIntervalCenter( i ) - (fTest.fHigh + fTest.fLow)/Float(2.0) ) < 0.0001, "Interval is broken!\n");
    oCurrentLayer.RotateZ( OmegaRangeToIndexMap.IndexToIntervalCenter( i ) );

    SVector3 Center = oCurrentLayer.ToLabFrame( SVector3( 0, 0, 0 ) );  // could change to arbitary center
      
    std::vector<CDetectorPeak> DetectorPeaks = Data.mImageMap[i][0].GetPeakList( );  // peaks are purged - BBox only
    TotalNumPeaks += DetectorPeaks.size();
    for( int j = 0; j < DetectorPeaks.size(); j ++ )
    {
      //       for( int k = 0; k < DetectorPeaks[j].vPixelList.size(); k ++ )  // can use BBox instead if we figure out the mapping is correct
      //       {
        Float Eta     [4];
        Float TwoTheta[4];
        BBoxT PixelBBox = DetectorPeaks[j].GetBoundingBox();

        PixelToEtaTwoTheta( Eta[0], TwoTheta[0], PixelBBox.pMin.x - 2.5, PixelBBox.pMin.y - 2.5, vDetectorList[0], Center );
        PixelToEtaTwoTheta( Eta[1], TwoTheta[1], PixelBBox.pMax.x + 2.5, PixelBBox.pMin.y - 2.5, vDetectorList[0], Center );
        PixelToEtaTwoTheta( Eta[2], TwoTheta[2], PixelBBox.pMin.x - 2.5, PixelBBox.pMax.y + 2.5, vDetectorList[0], Center );
        PixelToEtaTwoTheta( Eta[3], TwoTheta[3], PixelBBox.pMax.x + 2.5, PixelBBox.pMax.y + 2.5, vDetectorList[0], Center );
        
        Float MinEta = *( std::min_element( &Eta[0], &Eta[3] ) );
        Float MaxEta = *( std::max_element( &Eta[0], &Eta[3] ) );
        
        Float MinTwoTheta = *( std::min_element( &TwoTheta[0], &TwoTheta[3] ) );
        Float MaxTwoTheta = *( std::max_element( &TwoTheta[0], &TwoTheta[3] ) );
        
        typedef CSimulationRange::IndexIter IndexIter;
        std::pair<IndexIter, IndexIter> EtaRange      = EtaRangeToIndexMap     ( MinEta,      MaxEta      );
        std::pair<IndexIter, IndexIter> TwoThetaRange = TwoThetaRangeToIndexMap( MinTwoTheta, MaxTwoTheta );

        int nTried = 0;
        for( IndexIter EtaCur = EtaRange.first; EtaCur != EtaRange.second; ++ EtaCur )
        {
          for( IndexIter ThetaCur = TwoThetaRange.first; ThetaCur != TwoThetaRange.second; ++ ThetaCur )
          {
            
            if( *EtaCur != XDMSimulation::NoMatch && *ThetaCur != XDMSimulation::NoMatch )
            {
              nTried ++;
              IntensityMap[i][ *EtaCur ][ *ThetaCur ] ++;  // mark
              NumIntensityMapMarked ++;

            }
          }
        }
        
//         if( nTried <= 0 )
//         {
//           std::cout << "Failed: Angles " << MinEta << " " << MaxEta << " | " << MinTwoTheta << " " << MaxTwoTheta << std::endl;
//           std::cout << "Mask area has " << nTried << " points " << std::endl;
//           exit(0);
//         }
        // }
    }
    oCurrentLayer.SetOrientation( oCurOrientation );  // use this to reduce numerical errors
    
  }
  std::cout << "Total number of peaks inserted " << TotalNumPeaks
            << " non-zero histogram: " << NumIntensityMapMarked <<  std::endl;
}


//--------------------------------------------------------------------------------------------------------
//
//  CPaintGrid:: SimulatePeaks   -- Fast version
//
//  Note for reader:  GetReflectionVectorList is "lazy."  In another words, it doesn't recalculate unless
//  it needs to.  Therefore its complexity is really O(1).
//
//--------------------------------------------------------------------------------------------------------
std::vector<Float> CPaintGrid::Paint( const DetectorListT & vDetectorList,
                                      const vector<SMatrix3x3> & FZList,
                                      CSample oCurrentLayer  )
{
  
  std::cout << "Starting Paint " << std::endl; 
 
  Float fWavenumber =  PhysicalConstants:: keV_over_hbar_c_in_ang * oSetup.ExperimentalSetup().GetBeamEnergy();  // (Energy (keV)/( hbar c ) - in ang.
  // i.e., E = 64 => 64 keV/(hbar c) in final unit of ang
  const vector<CUnitCell> & oCryStructList = oCurrentLayer.GetStructureList();
  std::shared_ptr< CMic > pMic= std::dynamic_pointer_cast<CMic>( oCurrentLayer.GetMic() );
  
  std::cout << " Number of voxels " << pMic->GetNumVoxels() << std::endl;
  std::vector<int> IntersectDetectorCount( FZList.size(), 0);
  std::vector<int> EtaOmegaHitCount      ( FZList.size(), 0);

  std::ofstream DEBUG_PaintGridOutput( "DEBUG_PaintGrid.txt" );

  int NumIntersected = 0;
  for( vector<SVoxel>::const_iterator pCurVoxel = pMic->VoxelListBegin();
       pCurVoxel != pMic->VoxelListEnd(); pCurVoxel ++ )
  {    
    SVector3 oCenter = pCurVoxel->GetCenter();
    int nOrientCount    = 0;
    Int nCryStructIndex = pCurVoxel->nPhase;
    const vector<CRecpVector> & oRecipVectors = oCryStructList[ nCryStructIndex ].GetReflectionVectorList(); 

    const SMatrix3x3 oCurOrientation = oCurrentLayer.GetOrientationMatrix();    
        
#pragma omp parallel for private( oCurrentLayer )
#pragma omp nowait
    for( int n = 0; n < FZList.size(); n ++  )
    {
      if( omp_get_thread_num() == 0 && (n % 50000 == 0) )
        std::cout << n << "/" << FZList.size() << std::endl;
    
      CSample _ThreadSample = oCurrentLayer;
      SVoxel  v             = * pCurVoxel; 
      
      for( int nRecipIndex = 0; nRecipIndex < oRecipVectors.size();  nRecipIndex++ )   
      {
        v.oOrientMatrix = FZList[n];
        SVector3 oScatteringVec  = oRecipVectors[nRecipIndex].v; 	
        oScatteringVec.Transform( v.oOrientMatrix );    // g_hkl' = O * g_hkl
        
        Float fOmegaRes[2];
        //------------------------------------
        //  Magnitude is precomputed.  Rotation does not change magintude
        //------------------------------------
        bool  bPeakObservable = oSimulator.GetScatteringOmegas( fOmegaRes[0], fOmegaRes[1],
                                                                oScatteringVec, oRecipVectors[nRecipIndex].fMag );
        if( bPeakObservable )
        {
          oScatteringVec.Normalize();
          const SVector3 & oScatteringDir = oScatteringVec; 
          for ( int i = 0; i < 2; i ++ )   // using this for loop to enforce uniformity     
          {
            Size_Type nOmegaIndex = OmegaRangeToIndexMap( fOmegaRes[i] );
            if ( nOmegaIndex != XDMSimulation::NoMatch )
            {
              Float fSinTheta = oRecipVectors[nRecipIndex].fMag / ( Float(2.0) * fWavenumber );
              Float f2Theta   = Float( 2 ) * asin( fSinTheta );
              HEDM::XDMEtaAcceptFn FAcceptFn( - oSetup.ExperimentalSetup().GetEtaLimit(),
                                              oSetup.ExperimentalSetup().GetEtaLimit() );
              FAcceptFn.fSin2Theta = sin( f2Theta ) ;   
              _ThreadSample.RotateZ( fOmegaRes[i] );
              
              CRay oReflectedRay = DiffractionCore::GetReflectedRay( _ThreadSample, oCenter, oScatteringDir, oBeamDirection );
              SVector3 oReflectedRayDir = oReflectedRay.GetDirection();
              oReflectedRayDir.Normalize();
            
              Point p;
              Bool bIntersected = DiffractionCore::GetIlluminatedPixel( p, vDetectorList[0], oReflectedRay );
              
              if ( bIntersected && ( p.x >= 0 && p.x < vDetectorList[0].GetNumCols()
                                     && p.y >= 0 && p.y < vDetectorList[0].GetNumRows() )  )
              {
//                 {  //-------- DEBUG
//                   GeneralLib::SVector3 DEBUG_Reflected = vDetectorList[0].PixelToLabCoordinate( p.x, p.y );
//                   Float fT;
//                   vDetectorList[0].Intersects( oReflectedRay, fT);
//                   std::cout << " +  " << oReflectedRay.Evaluate( fT ) << std::endl;
//                   std::cout << " +  " << DEBUG_Reflected << std::endl;
//                   GeneralLib::SVector3 VecDiff = DEBUG_Reflected - oReflectedRay.Evaluate( fT );
//                   DEBUG_PaintGrid << VecDiff << std::endl;
//                 } //--------- DEBUG
                
                Float fEta = atan2( oReflectedRayDir.m_fZ, -oReflectedRayDir.m_fY  );   // convert to Joel's convention
                //                std::cout << fEta << " " << f2Theta << std::endl;
                PixelToEtaTwoTheta( fEta, f2Theta, p.x, p.y, vDetectorList[0], SVector3(0, 0,0) );
                //                std::cout << fEta << " " << f2Theta << std::endl;
                bool bAccept;
                Float fIntensity;
                std::tie( bAccept, fIntensity ) = FAcceptFn( oScatteringDir ); 
                if( bAccept )
                {
                  IntersectDetectorCount[n] ++;
                  
                  Size_Type EtaIdx      = EtaRangeToIndexMap( fEta );
                  Size_Type TwoThetaIdx = TwoThetaRangeToIndexMap( f2Theta );

                  if(EtaIdx == XDMSimulation::NoMatch)
                  {
                    std::cout << "EXCEPTION: Eta should never be out of range " << std::endl;
                    exit(0);
                  }
                  if( TwoThetaIdx == XDMSimulation::NoMatch)
                  {
                    std::cout << "EXCEPTION: 2 Theta out of range " << std::endl;
                  }
                  if( EtaIdx != XDMSimulation::NoMatch && TwoThetaIdx != XDMSimulation::NoMatch )
                    if( IntensityMap[ nOmegaIndex ][ EtaIdx ][ TwoThetaIdx ] > 0 )
                      EtaOmegaHitCount[n] ++;    
                  
                }
              }
              _ThreadSample.SetOrientation( oCurOrientation );  // use this to reduce numerical errors
            }
          }
        }
      }// end for each Reciprocal Vector
    }
  }

  int NumPaintedOverlap = 0;
  std::vector<Float> Completeness   ( FZList.size(), 0);
  //-----------
  //  TESTING
  for( int n = 0; n < FZList.size(); n ++  )
  {
    Completeness[n] = Float(EtaOmegaHitCount[n]) / Float(IntersectDetectorCount[n]);
    NumPaintedOverlap += EtaOmegaHitCount[n];
    DEBUG_PaintGridOutput
      << Float(EtaOmegaHitCount[n]) << " " <<  Float(IntersectDetectorCount[n])  << " "
      << Completeness[n] << std::endl;

  }
  std::cout << "FZList.size() " << FZList.size() << std::endl;
  std::cout << "NumPaintedOverlap " << NumPaintedOverlap << std::endl;

  DEBUG_PaintGridOutput.close();
  
  return Completeness;
}
