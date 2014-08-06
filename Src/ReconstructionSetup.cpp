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
//  ReconstructionSetup.cpp
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//   Purpose:  implementation of ReconstructionSetup
//
//
////////////////////////////////////////////////////////////
#include "ReconstructionSetup.h"

namespace Reconstruction
{
  //--------------------------------------------------------------------------------------------------------
  //  
  //  Private:  ReadExperimentalData
  //
  //--------------------------------------------------------------------------------------------------------
  Bool ReconstructionSetup::ReadExperimentalData( CSimulationData & oExpData, bool bServerMode )
  {
    const vector< SIntRange > &  vFileRangeList  = oExpSetup.GetFileRangeList();
    const vector< SRange    > &  vOmegaRangeList = oExpSetup.GetOmegaRangeList();
    const vector< CDetector > &  vDetectorList   = oExpSetup.GetDetectorList();
  
    RUNTIME_ASSERT( vFileRangeList.size() == vDetectorList.size() && vDetectorList.size() > 0,
                    "ERROR:  File range not specified for some detector(s) \n" );

    //
    //  TODO:  This needs to be more robust
    // 
    oExpData.Initialize( vOmegaRangeList.size(), vDetectorList.size(),
                         vDetectorList[0].GetNumCols(), vDetectorList[0].GetNumRows() );
  
    for ( Size_Type nDetector = 0; nDetector < vFileRangeList.size(); nDetector ++ )
    {
      Size_Type nNumFile = vFileRangeList[ nDetector ].nHigh - vFileRangeList[ nDetector ].nLow + 1;

      RUNTIME_ASSERT( nNumFile == vOmegaRangeList.size(),
                      "ERROR:  Number of files does not match specified number of omega inteval\n");

      //  Using -- openmp --  TODO - limit number of threads?
#pragma omp parallel for
      for ( Int nFileNum = vFileRangeList[nDetector].nLow;
            nFileNum <= vFileRangeList[nDetector].nHigh; nFileNum ++ )
      {
        Int nCurrentFileIndex = nFileNum - vFileRangeList[nDetector].nLow;
        stringstream ssFilename;
        ssFilename << InputParameters().InFileBasename
                   << InitFileIO::NumToSuffix( nFileNum, InputParameters().InFileSerialLength )
                   << "." << InputParameters().InFileExt << nDetector + oExpSetup.GetBCDetectorOffset();

        Bool bSuccess = false;
        switch( InputParameters().InFileType )
        {
          case eASCII:
            bSuccess= oExpData.mImageMap[ nCurrentFileIndex ][ nDetector ].ReadCXDMSimulationDataFile( ssFilename.str() );
            break;
          case eBin:
            bSuccess= oExpData.mImageMap[ nCurrentFileIndex ][ nDetector ].ReadCXDMSimulationUFFFile( ssFilename.str(), bServerMode );
            break;
          case eTif:
            RUNTIME_ASSERT( 0, "Error:  Tif input format not yet supported\n" );
            break;
          default:
            RUNTIME_ASSERT( 0, "Error:  Input File Type Not Recognized" );
        }
        DEBUG_ALERT( 0, ssFilename.str() + "\n" );
        RUNTIME_ASSERT( bSuccess, "Error!  Input data files cannot be read, please check file format and location\n" + ssFilename.str() );
      }
    }
    return true;
  }
  
  //------------------------------------------------------------------------
  //  InitializeWithDataFiles
  //
  //  Purpose:  This is the central point to read from all data/initialization
  //            file.
  //  Note:     This function or InitializeWithRestore must be called beore
  //            any reconstruction!
  //------------------------------------------------------------------------
  void ReconstructionSetup::InitializeWithDataFiles( const CConfigFile & oConfigFile, bool bServerMode )
  {
    oExpSetup.SetConfigFile( oConfigFile );
    std::cerr << "Before initialization " << std::endl;
    oExpSetup.InitializeExperiment( );   // read in  FZ files, omega files

    std::cerr << "Before reading experimental data " << std::endl;
    ReadExperimentalData( oExpData, bServerMode );    // read in experimental data
    std::cerr << "After reading experimental data " << std::endl;

    vector<SVector3> oFZEulerAngleList;
    bool bSuccess = InitFileIO::ReadFundamentalZoneFile( oFZEulerAngleList, InputParameters().FundamentalZoneFilename );
    RUNTIME_ASSERT( bSuccess, "InitializeWithDataFiles:  Fundamental Zone File Not Found!\n");
  
    for ( Size_Type i = 0; i < oFZEulerAngleList.size(); i ++)
    {
      SMatrix3x3 oOrient;
      oOrient.BuildActiveEulerMatrix( oFZEulerAngleList[i].m_fX,
                                      oFZEulerAngleList[i].m_fY,
                                      oFZEulerAngleList[i].m_fZ );

      oFZOrientationList.push_back( oOrient );
    }
  
  
    const vector< CDetector > &  vDetectorList = oExpSetup.GetDetectorList();
    oExpSetup.InitializeSample( oReconstructionLayer,
                                vDetectorList[ 0 ] );  // binding CurrentLayer to the first detector's limit
    //
    // TODO:  Initialization should become uniform (or group together somewhere)
    //
    FPeakFilter.fMinEta = - InputParameters().fEtaLimit;
    FPeakFilter.fMaxEta =   InputParameters().fEtaLimit;
    FPeakIntensityAccept.fMinEta = - InputParameters().fEtaLimit;
    FPeakIntensityAccept.fMaxEta =   InputParameters().fEtaLimit;

    pPartialResult = MicIOFactory::Create( InputParameters().nMicGridType );   // initialize partial result

    // Read partial result
    if( InputParameters().bUsePartialResult )
    {
      const string & filename =  InputParameters().PartialResultFilename;
      std::cout << "before reading partial result " << std::endl;
      //  select out the correct result
      RUNTIME_ASSERT( pPartialResult->Read( filename ),
                      "Failed to read partial result at: " +  filename );

    }
  }
  
  //------------------------------------------------------------------------
  //
  //  Save
  //
  //------------------------------------------------------------------------
  Bool ReconstructionSetup::Save  ( CSerializer   & oSerialBuf ) const
  {
    Bool bSuccess;
    bSuccess = oSerialBuf.InsertCompactObj                 ( FPeakFilter );
    bSuccess = bSuccess && oSerialBuf.InsertCompactObj     ( FPeakIntensityAccept );
    bSuccess = bSuccess && oSerialBuf.InsertCompactVector  ( oFZOrientationList );
    
    bSuccess = bSuccess && oExpSetup.Save                  ( oSerialBuf );
    bSuccess = bSuccess && oReconstructionLayer.Save       ( oSerialBuf );
    return bSuccess;
  }

  //------------------------------------------------------------------------
  //
  //  Restore
  //
  //------------------------------------------------------------------------
  Bool ReconstructionSetup::Restore( CDeserializer & oSerialBuf )
  {
    Bool bSuccess;
    oFZOrientationList.clear();
    bSuccess = oSerialBuf.GetCompactObj                     (  &FPeakFilter );
    bSuccess = bSuccess && oSerialBuf.GetCompactObj         ( & FPeakIntensityAccept );
    bSuccess = bSuccess && oSerialBuf.GetCompactVector      ( oFZOrientationList );

    bSuccess = bSuccess && oExpSetup.Restore                ( oSerialBuf );

    //--------- HACK  - need to separate Sample from MicIO
    oReconstructionLayer.InitializeMic ( InputParameters().nMicGridType );
    //--------- END HACK  - need to separate Sample from MicIO

    bSuccess = bSuccess && oReconstructionLayer.Restore     ( oSerialBuf );
 
    return bSuccess;
  }

}
