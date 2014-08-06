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
//
//   InitFilesIO.h
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:  This file exist only as an organizational helper.  Most functions used here are only
//             called once during initialization.  However, their existence lead to cluter in their caller
//             objects.  
//             All input, output, and related routines go into this file. 
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _INIT_FILES_IO_H
#define _INIT_FILES_IO_H

#include "CrystalStructure.h"
#include "PhysicalConstants.h"
#include "Peak.h"
#include "Parser.h"
#include <string>
#include <fstream>
#include "Types.h"
#include <vector>
#include "3dMath.h"
#include "Debug.h"
#include "Error.h"
#include <sstream>


///
//  TODO:  Possibly convert this into an object oriented design
//
//
//
using std::vector;
using std::string;
using std::stringstream;
using namespace GeneralLib;





namespace InitFileIO
{
  
  //
  //  pre:  vBuf has enough elements to be assigned to a vector
  //  post:  element 1, 2, 3 are assigned to v.m_fX, v.m_fY, v.m_fZ
  //
  SVector3 ExtractVector( const vector< string > & vBuf, Size_Type nLineNumber );
  Int      ExtractInt   ( const vector< string > & vBuf, Size_Type nLineNumber );
  Float    ExtractReal  ( const vector< string > & vBuf, Size_Type nLineNumber );
  Bool     FindBlocks   ( vector<string> &vBlocks, const string & sBuf );
  
  ////////////////////////////////////////////////////////////
  //
  //  ReadFileToBuf:
  //  Return:  Pointer to a buffer created by ReadFileToBuf
  //      !!!  It is the responsibility of the caller to free the
  //           memory location with delete []
  char* ReadFileToBuf(Size_Type & nBufferSize, string filename);
  
  bool ReadCrystalStructureFile(CUnitCell &oCell, const string &filename);

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
                                  const string & filename );


  //------------------------------
  //
  //  ReadFundamentalZoneFile
  //
  //  Note that it is a set of euler angles that's being returned
  //
  //-----------------------------
  bool ReadFundamentalZoneFile( vector<SVector3> & vFZEulerAngleList,
                                const string & sFilename );
  
  //------------------------------
  //   ReadDetectorDistanceFiles
  //
  //   Parse detector distance file and return a set of detectors
  //------------------------------
  bool ReadDetectorDistanceFiles( );

  ////////////////////////////////////////////
  //
  //  NumToSuffix
  //
  //
  ////////////////////////////////////////////
  string NumToSuffix( UInt n, UInt length );
}

#endif
