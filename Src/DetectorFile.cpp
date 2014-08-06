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
#include "DetectorFile.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//  Class  DetectorFile
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace InitFileIO
{
  //-------------------------------------------------------------------------------------
  //  GetDetector
  //
  //  TODO:  Reorganize.  Manually calculating coordinate origin seems dangerous, especially 
  //         when this has to happen at multiple places in the code.
  //-------------------------------------------------------------------------------------
  CDetector CDetectorInfo::GetDetector() const
  {
    CDetector oDetector =  CXDMDetectorFactory::MakeDetector( nNumJPixels, nNumKPixels, oLabFrameLocation,
                                                              fBeamCenterJ, fBeamCenterK,
                                                              vJUnitVector, vKUnitVector,
                                                              fPixelWidth, fPixelHeight,
                                                              oLabFrameOrientMatrix );
    return oDetector;
  }
  
  //-------------------------------------------------------------------------------------
  //
  //  ParseDetectorBlock
  //
  //  
  //
  //-------------------------------------------------------------------------------------
  Bool CDetectorInfo::ParseDetectorBlock( const string & sBuf )
  {

    vector< vector<string> > vsTokens;
    GeneralLib::Tokenize( vsTokens, sBuf, " \t\n");
    const char* sKeywordStrings[] = {
      // comment
      "#",

      "JUnitVector",
      "KUnitVector",

      "BeamCenterJ",
      "BeamCenterK",

      "LabFrameLocation",
      "LabFrameOrientation",

      "NumJPixels",
      "NumKPixels",

      "PixelJLength",
      "PixelKLength"       
    };

     enum EKeyword{
       // comment
       eComment,

       eJUnitVector,
       eKUnitVector,
       
       eBeamCenterJ,
       eBeamCenterK,
       
       eLabFrameLocation,
       eLabFrameOrientation,
       
       eNumJPixels,
       eNumKPixels,
       
       ePixelJLength,
       ePixelKLength,
       
       // error check
       eNumKeywords
     };

     // actual parsing
     for(Size_Type i = 0; i < vsTokens.size(); i ++)
     {
       Size_Type iFirstToken;


       if ( vsTokens[i].size() == 0)
       {
         iFirstToken = eComment;
       }
       else if ( vsTokens[i][0].find_first_of( sKeywordStrings[ eComment ] ) == 0) // if first token is comment
       {
         iFirstToken = eComment;
       }
       else
       {
         // Identify the keyword in the beginning of the line (i.e., vertex? texture?
         for(iFirstToken = 0; iFirstToken < eNumKeywords; iFirstToken ++)
         {
           if(strcmp(sKeywordStrings[iFirstToken], vsTokens[i][0].c_str()) == 0)
             break;
         }
         CONFIG_DEBUG( std::cout << i + 1 << " " << sKeywordStrings[iFirstToken] << " "
                       << vsTokens[i][0].c_str() << " Token " << iFirstToken << endl); 
       }

       //  Actual processing
       //
       switch ( iFirstToken )
       {
         case eComment:
           break;
         case eJUnitVector:
           vJUnitVector = ExtractVector( vsTokens[i], i );
           break;
         case eKUnitVector:
           vKUnitVector = ExtractVector( vsTokens[i], i );
           break;
         case eBeamCenterJ:
           fBeamCenterJ = ExtractReal( vsTokens[i], i );
           break;
         case eBeamCenterK:
           fBeamCenterK = ExtractReal( vsTokens[i], i );
           break;
         case eLabFrameLocation:
           oLabFrameLocation = ExtractVector( vsTokens[i], i );
           break;
         case eLabFrameOrientation:
           {
             vLabFrameOrientation = ExtractVector( vsTokens[i], i );
             vLabFrameOrientation = DegreeToRadian( vLabFrameOrientation );
             oLabFrameOrientMatrix.BuildActiveEulerMatrix( vLabFrameOrientation.m_fX,
                                                           vLabFrameOrientation.m_fY,
                                                           vLabFrameOrientation.m_fZ );
             break;
           }
         case eNumJPixels:
           nNumJPixels = ExtractInt( vsTokens[i], i );
           break;
         case eNumKPixels:
           nNumKPixels = ExtractInt( vsTokens[i], i );
           break;
         case ePixelJLength:    // note:  this follows the convention of the 
                                // program of having J as the "horizontal" direction
           fPixelWidth = ExtractReal( vsTokens[i], i );
           break;
         case ePixelKLength:
           fPixelHeight = ExtractReal( vsTokens[i], i );;
           break;
         default:
           {
             cerr << "[DetectorFile] Error: syntax not recognized:  Line " << i  << endl;
             cerr <<  vsTokens[i][0] << endl;
             exit(0);
             return false;
           }
           // insert error stuff here
           
       }
     }
     return true;
  }

  //-------------------------------------------------------------------------------------
  //
  //  Print
  //
  //-------------------------------------------------------------------------------------
  void CDetectorInfo::Print( ostream &os ) const
  {
    os << "JUnitVector " << vJUnitVector << endl;
    os << "KUnitVector " << vKUnitVector << endl; 
    
    os << "BeamCenterJ " << fBeamCenterJ << endl; 
    os << "BeamCenterK " << fBeamCenterK << endl;
    
    os << "LabFrameLocation  " << oLabFrameLocation << endl;
    os << "LabFrameOrientation  " << RadianToDegree( vLabFrameOrientation ) << endl;
    
    os << "NumJPixels  " << nNumJPixels << endl;
    os << "NumKPixels  " << nNumKPixels << endl;
    
    os << "PixelJLength " << fPixelWidth << endl;
    os << "PixelKLength  " << fPixelHeight << endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //  Class CDetectorFile
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////

  
  //-------------------------------------------------------------------------------------
  //
  //  FindDetectorBlock
  //
  //-------------------------------------------------------------------------------------
  Bool CDetectorFile::FindDetectorBlock( vector<string> &vBlocks, const string & sBuf )
  {
    return FindBlocks( vBlocks, sBuf );
  }

  
  //-------------------------------------------------------------------------------------
  //
  //  Parse
  //
  //  Ghetto parser for detector files
  //
  //-------------------------------------------------------------------------------------
  Bool CDetectorFile::Parse( vector<CDetector> &oDetectorList, const string & filename )
  {
  
    char *pBuffer = NULL;
    Size_Type nBufSize = 0;
    pBuffer = InitFileIO::ReadFileToBuf( nBufSize, filename );
    string sBufStr( pBuffer, nBufSize );

    vector<string> vsBlocks;
    if (! FindDetectorBlock( vsBlocks, sBufStr ) )
      return false;

    if ( vsBlocks.size() < 1 )
      return false;
    
    for ( Size_Type i = 0; i < vsBlocks.size(); i ++ )
    {
      CDetectorInfo oBlockInfo;
      std::cout << " In Detector Block " << i <<  " (Line numbers are specified within the block)" << std::endl;
      oBlockInfo.ParseDetectorBlock( vsBlocks[i] );

      CDetector oNewDetector = oBlockInfo.GetDetector();
      oDetectorList.push_back( oNewDetector );
    }
  
    delete [] pBuffer;
    return true;
  }

  
  //-------------------------------------------------------------------------------------
  //
  //  Parse
  //
  //  Ghetto parser for detector files
  //
  //-------------------------------------------------------------------------------------
  Bool CDetectorFile::ParseOptimizationInfo( vector<CDetectorInfo> &oDetStepSizeInfo, const string & filename )
  {
    char *pBuffer = NULL;
    Size_Type nBufSize = 0;
    pBuffer = InitFileIO::ReadFileToBuf(nBufSize, filename);
    string sBufStr( pBuffer, nBufSize );

    vector<string> vsBlocks;
    if ( ! FindDetectorBlock( vsBlocks, sBufStr ) )
      return false;
    
    if ( vsBlocks.size() < 1 )
    {
      RUNTIME_ASSERT(0, "Optimization file not in correct format\n");
      return false;
    }
    for ( Size_Type i = 0; i < vsBlocks.size(); i ++ )
    {
      CDetectorInfo oBlockInfo;
      std::cout << " In Detector Block " << i <<  " (Line numbers are specified within the block)" << std::endl;
      oBlockInfo.ParseDetectorBlock( vsBlocks[i] );
      oDetStepSizeInfo.push_back( oBlockInfo );
    }
  
    delete [] pBuffer;
    return true;
  }
  
}
