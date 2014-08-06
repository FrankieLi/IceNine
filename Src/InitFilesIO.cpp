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
#include "InitFilesIO.h"
#include "ConfigFile.h"
#include "Pixel.h"

using std::cerr;

namespace InitFileIO{
  
  
  //-------------------------------------------------------------------------------------
  // ExtractVector
  //-------------------------------------------------------------------------------------
  SVector3 ExtractVector( const vector< string > & vBuf, Size_Type nLineNumber )
  {
    SVector3 v;
    
    if ( vBuf.size() < 4 )
    {
      cerr << "[ExtractVector] Error: Incorrect number of parameters " << nLineNumber << endl;
      exit( 0 ); 
    }
    else
    {
      v.m_fX = atof( vBuf[1].c_str() );
      v.m_fY = atof( vBuf[2].c_str() );
      v.m_fZ = atof( vBuf[3].c_str() );
    }
    return v;
  }
  
  //-------------------------------------------------------------------------------------
  //  ExtractInt
  //-------------------------------------------------------------------------------------
  Int ExtractInt( const vector< string > & vBuf, Size_Type nLineNumber )
  {
    if ( vBuf.size() < 2 )
    {
      cerr << "[ExtractInt] Error: Incorrect number of parameters " << nLineNumber << endl;
      exit( 0 ); 
    }
    return atoi( vBuf[1].c_str() ) ;
  }

  //-------------------------------------------------------------------------------------
  //  ExtractReal
  //-------------------------------------------------------------------------------------
  Float ExtractReal( const vector< string > & vBuf, Size_Type nLineNumber )
  {
    if ( vBuf.size() < 2 )
    {
      cerr << "[ExtractReal] Error: Incorrect number of parameters " << nLineNumber << endl;
      exit( 0 ); 
    }
    return atof( vBuf[1].c_str() ) ;
  }
  
  //-------------------------------------------------------------------------------------
  //
  //  FindBlocks
  //
  //-------------------------------------------------------------------------------------
  Bool FindBlocks( vector<string> &vBlocks, const string & sBuf )
  {
 
    Size_Type nSearchPos = 0;
    string::const_iterator pStringIt = sBuf.begin();

    while ( nSearchPos < sBuf.size() )
    {

      Size_Type nStartPos = sBuf.find_first_of( '{', nSearchPos );
      Size_Type nEndPos = sBuf.find_first_of( '}', nSearchPos );

      RUNTIME_ASSERT( !( nStartPos == string::npos && nEndPos != string::npos ), 
                      "[DetectorFile] ERROR: Invalid Detector Block:  { expected\n");
      RUNTIME_ASSERT( !( nStartPos != string::npos && nEndPos == string::npos ), 
                      "[DetectorFile] ERROR: Invalid Detector Block:  } expected\n");

      // if neither are found, then job's finished here
      if ( nStartPos == string::npos && nEndPos == string::npos )
      {
        break;
      }
      else
      {
        nSearchPos = nEndPos + 1;  // advance search position
        string sToInsert( pStringIt + nStartPos + 1, pStringIt + nEndPos );
        vBlocks.push_back( sToInsert );
      }
    
    }

    return ( vBlocks.size() > 0 );
  }

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //
  //
  //
  //
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  bool ReadPeakFile(vector<CPeak> &vPeaks, const string &filename)
  {

    char *pBuffer = NULL;
    Size_Type nBufSize = 0;
    vector< vector<string> > vsTokens;
	
    pBuffer = ReadFileToBuf(nBufSize, filename);
    if(pBuffer == NULL){
      return false;
    }
	
	
    GeneralLib::Tokenize(vsTokens, pBuffer, "# \t\n");
	
    if(vsTokens[0].size() > 3){
      cerr << "ERROR: " << filename << " is not a peak file." << endl;
      return false;
    }
	

    // for each line
    CPeak currentPeak;
    Pixel p;
    Int peakID;
    Int currentID = 0;

    // get current peak ID
    if(vsTokens[1].size() < 4 ){
      cerr << "[CMicFile::ReadMicFile] Error: Incorrect number of columns, line: "
           << 1  << endl;
      return false;
    }
    peakID =  atoi(vsTokens[1][3].c_str());
	
    for(UInt i = 1; i < vsTokens.size(); i ++){
      if(vsTokens[i].size() < 4 ){
        cerr << "[CMicFile::ReadMicFile] Error: Incorrect number of columns, line: "
             << i  << endl;
        return false;
      }else{


        p.x = atoi(vsTokens[i][0].c_str());
        p.y = atoi(vsTokens[i][1].c_str());
        p.fIntensity = atof(vsTokens[i][2].c_str());
        currentID = atoi(vsTokens[i][3].c_str());
      }
		
      // add when either new ID or last pixel
      if( ( currentID != peakID )|| (i == (vsTokens.size() -1 ))){
        vPeaks.push_back(currentPeak);
        currentPeak.vPixelList.clear();
        currentPeak.vPixelList.push_back(p);
        peakID = currentID;			
      }else{
        currentPeak.vPixelList.push_back(p);
      }
    }

    cerr << "NumPeaks Read: " << vPeaks.size() << endl;

    delete[] pBuffer;
	
    return true;

  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //  Returns a Null terminated C-string
  //
  //
  //
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  char* ReadFileToBuf(Size_Type & nBufferSize, string filename)
  {

    ifstream infile;
    infile.open(filename.c_str());
    FILE *pFile;
    char *pBuffer;

    if ( (pFile = fopen(filename.c_str(), "r" )) == NULL ){
      cerr << "[ReadFileToBuf]Cannot Open File: " << filename.c_str()  << endl;
      return NULL;
    }

    // Get the size of the file
    fseek( pFile, 0L, SEEK_END ); // Position to end of file
    nBufferSize = ftell( pFile );           // Get file length
    rewind( pFile );                                // Back to start of file

    if( nBufferSize <= 0 )
    {
      return NULL;
    }
    
    // Read in the entire file and close the file handle
    pBuffer = new char[nBufferSize + 1];

    if ( fread( pBuffer, nBufferSize, 1, pFile ) <= 0 ){
      fclose( pFile );
      cerr << "[ReadFileToBuf]:File read error" << endl;
      return NULL;
    }

    fclose( pFile );

    pBuffer[ nBufferSize ] = '\0'; // NULL termination of string


    return pBuffer;

  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //
  //  ReadCrystalStructureFile
  // 
  //  Note:  Units used in this file is angstrom as opposed to mm everywhere else.  This is because
  //  the actual measurement of the scattering vector is rarely used in conjunction with the geometric aspect.
  // (Typically unit vectors are used)
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool ReadCrystalStructureFile(CUnitCell &oCell, const string &filename)
  {
    char *pBuffer = NULL;
    Size_Type nBufSize = 0;
    vector< vector<string> > vsTokens;
	
    pBuffer = ReadFileToBuf(nBufSize, filename);
    if(pBuffer == NULL){
      return false;
    }

    string sBuffStr( pBuffer, nBufSize );
    
    // RUNTIME_ASSERT(oCell.oTranslationVector.size() == 0, "Vector screwed up\n");
    GeneralLib::Tokenize( vsTokens, sBuffStr, ",# \t\n");
    RUNTIME_ASSERT( vsTokens.size() >= 4,
                   "[ReadCrystalStructureFile] ERROR: Invalid file type (too short)\n" );
	
    // parse a, b, c unit cell lengths
    RUNTIME_ASSERT( vsTokens[0].size() >= 3,
                   "[ReadCrystalStructureFile] ERROR: failed to parse a, b, c\n" );

    oCell.SetUnitCellLength(  atof( vsTokens[0][0].c_str() ),
                              atof( vsTokens[0][1].c_str() ),
                              atof( vsTokens[0][2].c_str() ) );
                              
    // parse alpha, beta, gamma unit cell angles
    RUNTIME_ASSERT(vsTokens[1].size() >= 3,
                   "[ReadCrystalStructureFile] ERROR: failed to parse alpha, beta, gamma\n");

    oCell.SetUnitCellBasisAngles( DEGREE_TO_RADIAN( atof(vsTokens[1][0].c_str()) ),
                                  DEGREE_TO_RADIAN( atof(vsTokens[1][1].c_str()) ),
                                  DEGREE_TO_RADIAN( atof(vsTokens[1][2].c_str()) )  );
    
    // parse number of basis atoms
    RUNTIME_ASSERT(vsTokens[2].size() >= 1,
                   "[ReadCrystalStructureFile] ERROR: failed to parse number of atoms \n");
    //    oCell.nNumAtoms = atoi(vsTokens[2][0].c_str());	

    oCell.SetNumAtoms( atoi( vsTokens[2][0].c_str() ) );
    
    // parse Z, number of protons in the atom, and the three coordinates
    //          of the Premitive Translation Vectors
    for(Size_Type i = 3; i < vsTokens.size(); i ++){
      
      RUNTIME_ASSERT( (vsTokens[i].size() >= 4), 
                      "[ReadCrystalStructureFile] ERROR: Invalid translation vector\n");
      CAtom pVector;		

      SVector3 oPos ( atof( vsTokens[i][1].c_str() ), 
                      atof( vsTokens[i][2].c_str() ), 
                      atof( vsTokens[i][3].c_str() ) );
      // oCell.oTranslationVector.push_back( pVector );
  
      oCell.AddAtom( atoi( vsTokens[i][0].c_str() ), oPos );
    }
	
    oCell.InitializeCoordinateSystem();

    delete [] pBuffer;
    return true;
  }



  ///////////////////////////////////////////////////////////
  //
  //
  //  ReadRotationIntervalFile (a.k.a "omegafiles")
  //
  //  Given a filename, read in the range of data collection interval and
  //  the experimental file range
  //
  //  TODO:  Create a more flexible format to replace the omega files
  //
  bool ReadRotationIntervalFiles( vector<SRange> &oRotationRange,
                                  vector<SIntRange> &oFileRange,
                                  Size_Type nNumFileRange,
                                  const string & filename )
    
  {
    char *pBuffer = NULL;
    Size_Type nBufSize = 0;
    vector< vector<string> > vsTokens;
	
    pBuffer = ReadFileToBuf(nBufSize, filename);
    if(pBuffer == NULL){
      return false;
    }

    string sBuffStr( pBuffer, nBufSize );
    
    GeneralLib::Tokenize( vsTokens, sBuffStr, ",# \t\n");
    RUNTIME_ASSERT( vsTokens.size() >= 4,
                   "[ReadRotationIntervalFile] ERROR: Invalid file type (too short)\n");

    // conform to existing omega file format, the first line is ignored
    for ( Size_Type i = 1; i <= nNumFileRange; i ++)
    {
      SIntRange s;
      s.nLow  = atoi( vsTokens[i][0].c_str() );
      s.nHigh = atoi( vsTokens[i][1].c_str() );
      oFileRange.push_back( s );
    }

    //
    //  Get the omega ranges
    //
    for(Size_Type i = oFileRange.size() + 1; i < vsTokens.size(); i ++){
      
      RUNTIME_ASSERT( (vsTokens[i].size() >= 2), 
                      "[ReadRotationIntervalFiles] ERROR: Invalid file range\n");
      
      SRange s;
      s.fLow  = DEGREE_TO_RADIAN( atof( vsTokens[i][0].c_str() ) );
      s.fHigh = DEGREE_TO_RADIAN( atof( vsTokens[i][1].c_str() ) );
      oRotationRange.push_back( s );
    }


    delete [] pBuffer;
    return true;
    
  }


  //------------------------------
  //
  //  ReadFundamentalZoneFile
  //
  //------------------------------
  bool ReadFundamentalZoneFile( vector<SVector3> & vFZEulerAngleList,
                                const string & sFilename )
  {
    char *pBuffer = NULL;
    Size_Type nBufSize = 0;
    vector< vector<string> > vsTokens;
	
    pBuffer = ReadFileToBuf(nBufSize, sFilename);
    if(pBuffer == NULL){
      return false;
    }

    string sBuffStr( pBuffer, nBufSize );
    
    GeneralLib::Tokenize( vsTokens, sBuffStr, ",# \t\n");
    
    for(Size_Type i = 0; i < vsTokens.size(); i ++){
      
      RUNTIME_ASSERT( (vsTokens[i].size() >= 3), 
                      "[ReadFundamentalZoneFile] ERROR: Invalid Euler Angles\n");
      
      
      Float fPhi   = DEGREE_TO_RADIAN( atof( vsTokens[i][0].c_str() ) );
      Float fTheta = DEGREE_TO_RADIAN( atof( vsTokens[i][1].c_str() ) );
      Float fPsi   = DEGREE_TO_RADIAN( atof( vsTokens[i][2].c_str() ) );

      SVector3 newAngles;

      newAngles.Set( fPhi, fTheta, fPsi );

      vFZEulerAngleList.push_back( newAngles );
    }


    delete [] pBuffer;
    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //  NumToSuffix
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  string NumToSuffix( UInt n, UInt length )
  {
    string s;
    stringstream tmpSS;
    Int originalSize;
    tmpSS << n;
    s = tmpSS.str();
    originalSize = (Int) s.size();
	
    for(Int i = 0; i < (  (Int)length - originalSize ); i++ ){
      s = '0' + s;
    }

    return s;
  }
  

}


