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
///////////////////////////////////////////////////////
//
//  OptimizationConfig.cpp
//   Implementation of optimization config.
//
///////////////////////////////////////////////////////

#include "OptimizationConfig.h"

namespace OptimizationConfig
{

  namespace StrainOpt
  {
    //--------------------------------------
    //  ParseStrainFile
    //--------------------------------------
    EOptimizationMethod ProcessOptMethod( const vector< std::string > & Line,
                                          int nLineNumber )
    {
      if( Line.size() < 2 )
      {
        std::cerr << "Line " << nLineNumber << std::endl
                  << "StrainOpt: Optimization Method expects an argument" << std::endl
                  <<  "[ MonteCarlo | ConjugateGradient ]" << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      
      EOptimizationMethod Method;
      for( Size_Type j = 0; j < eNumOptimizationMethods; j++ )
      {
        if ( OptMethodKeywords[j] == Line[1] )
          Method = static_cast<EOptimizationMethod>( j );
      }

      if( Method >= eNumOptimizationMethods )
      {
        std::cerr << "Line " << nLineNumber << std::endl
                  << "StrainOpt:  Unknown strain optimization method, choices are [ MonteCarlo | ConjugateGradient ]"
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }

      return Method;
    }

    //--------------------------------------
    //  ParseStrainFile
    //--------------------------------------
    EStrainType ProcessStrainType( const vector< std::string > & Line,
                                   int nLineNumber)
    {
      if( Line.size() < 2 )
      {
        std::cerr << "Line " << nLineNumber <<  std::endl
                  << "StrainOpt: Strain Type expects an arguments \n [ Diag | Sym | Asym ]"
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      
      EStrainType Method;
      for( Size_Type j = 0; j < eNumStrainTypes; j++ )
      {
        if ( StrainTypeKeywords[j] == Line[1] )
          Method = static_cast<EStrainType>( j );
      }

      if( Method >= eNumStrainTypes )
      {
        std::cerr << "Line " << nLineNumber << std::endl
                  << "StrainOpt:  Unknown strain type, choices are [ Diag | Sym | Asym ]"
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      return Method;
    }

    //--------------------------------------
    //  ParseStrainFile
    //--------------------------------------
    ESearchType ProcessSearchType( const vector< std::string > & Line,
                                   int nLineNumber)
    {
      if( Line.size() < 2 )
      {
        std::cerr << "Line " << nLineNumber <<  std::endl
                  << "StrainOpt: Strain Search Type expects an argument \n [ Global | Local | Mix ]"
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      
      ESearchType Method;
      for( Size_Type j = 0; j < eNumStrainSearches; j++ )
      {
        if ( SearchTypeKeywords[j] == Line[1] )
          Method = static_cast<ESearchType>( j );
      }
      if( Method >= eNumStrainSearches )
      {
        std::cerr << "Line " << nLineNumber << std::endl
                  << "StrainOpt:  Unknown strain search type, choices are [ Global | Local | Mix ]"
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      return Method;
    }

    //--------------------------------------
    //  ParseStrainFile
    //--------------------------------------
    std::pair<EConvergenceMethod, Float> ProcessConvergenceMethod( const vector< std::string > & Line,
                                                                   int nLineNumber)
    {
      if( Line.size() < 3 )
      {
        std::cerr << "Line " << nLineNumber <<  std::endl
                  << "StrainOpt: ConvergenceMethod expects two arguments " << std::endl
                  << "[ Quality fQ | HitRatio fH | Differential fTol ]" << std::endl
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      EConvergenceMethod Method;
      for( Size_Type j = 0; j < eNumConvergenceMethods; j++ )
      {
        if ( ConvergenceMethodKeywords[j] == Line[1] )
          Method = static_cast<EConvergenceMethod>( j );
      }
      if( Method >= eNumConvergenceMethods )
      {
        std::cerr << "Line " << nLineNumber << std::endl
                  << "StrainOpt:  Unknown ConvergenceMethod, choices are [ Quality fQ | HitRatio fH | Differential fTol ]" << std::endl
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      return std::make_pair( Method, InitFileIO::ExtractReal( Line, nLineNumber ) );
    }

    //--------------------------------------
    //  ParseMatrix
    //--------------------------------------
    GeneralLib::SMatrix3x3 ParseMatrix( const vector< vector< string > > & vsTokens,
                                        int nLineNumber )
    {
      GeneralLib::SMatrix3x3 oRes;
      
      if( ( static_cast<int>( vsTokens.size() ) - nLineNumber ) < 3  )  // need 3 lines
      {
        std::cerr << "Line " << nLineNumber <<  std::endl
                  << "StrainOpt: Expecting a 3 x 3 matrix here, but end of file is encountered prematurely " 
                  << std::endl;
        RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
      }
      
      for( int i = 0; i < 3; i ++ )
      {
        int nCurrentLine = nLineNumber + i;
        const vector< string > & Line = vsTokens[ nCurrentLine ];

        for( int j = 0; j < 3; j ++ )
        {
          int nShift = 0;
          if( i == 0 )
            nShift = 1;

          if( j + nShift > static_cast<int>( Line.size() ) )
          {
            std::cerr << "Line " << nCurrentLine <<  std::endl
                      << "StrainOpt: Expecting a 3 x 3 matrix here, but end of line is encountered prematurely " 
                    << std::endl;
            RUNTIME_ASSERT( 0, " Parsing Strain Optimization File Failed ");
          }          
          oRes.m[i][j] = atof( Line[ j + nShift ].c_str() );
        }
      }
      return oRes;
    }
  }

  //--------------------------------------
  //  ParseStrainFile
  //--------------------------------------
  bool StrainConfigFile::Parse( const std::string & sBuf )
  {
    bool *vInitializationCheck;
    vInitializationCheck = new bool[ StrainOpt::eNumKeywords ]; 
    
    for ( Int i = 0; i < StrainOpt::eNumKeywords; i ++ )
      vInitializationCheck[i] = false;

    std::vector< std::vector< std::string> > vsTokens;
    GeneralLib::Tokenize( vsTokens, sBuf, " \t\n");
    for(Size_Type i = 0; i < vsTokens.size(); i ++)
    {
      Size_Type iFirstToken;
      
      if ( vsTokens[i].size() == 0 )
      {
        iFirstToken = StrainOpt::eComment;
      }
      else if ( vsTokens[i][0].find_first_of( StrainOpt::Keywords[ StrainOpt::eComment ] )
                == 0 ) // if first token is comment
      {
        iFirstToken = StrainOpt::eComment;
      }
      else
      {
        // Identify the keyword in the beginning of the line (i.e., vertex? texture?
        for( iFirstToken = 0;
             iFirstToken < StrainOpt::eNumKeywords; iFirstToken ++)
        {
          if( strcmp( StrainOpt::Keywords[iFirstToken], vsTokens[i][0].c_str()) == 0 )
            break;
        }
        CONFIG_DEBUG( std::cout << i + 1 << " " << StrainOpt::Keywords[iFirstToken] << " "
                      << vsTokens[i][0].c_str() << " Token " << iFirstToken << endl ); 
      }
      
      switch( iFirstToken )   // look at first token of each line
      {
        case StrainOpt::eComment:  
          break;
        case StrainOpt::eStrainSearchType:
          SearchType = StrainOpt::ProcessSearchType( vsTokens[i], i );
          break;
        case StrainOpt::eStrainType:
          StrainType = StrainOpt::ProcessStrainType( vsTokens[i], i );
          break;
        case StrainOpt::eOptimizationMethod:
          OptimizationMethod = StrainOpt::ProcessOptMethod( vsTokens[i], i );
          break;
        case StrainOpt::eStrainStepSize:
          StrainSteps = StrainOpt::ParseMatrix( vsTokens, i );
          i += 2; // need to skip the matrix
          break;
        case StrainOpt::eMaxStrainStep:
          nMaxSteps =  InitFileIO::ExtractInt( vsTokens[i], i );
          break;
        case StrainOpt::eConvergenceMethod:
          std::tie( ConvergenceMethod, fConvergenceThresh )
            = StrainOpt::ProcessConvergenceMethod( vsTokens[i], i );
          break;
          
        default:
          {
            cerr << "[StrainConfig] Error: syntax not recognized:  Line " <<  i  << endl;
            exit(0);
            return false;
          }
      }
      vInitializationCheck[ iFirstToken ] = true;
    }

    bool bMissingTokens = false;
    for ( Int i = 0; i < StrainOpt::eNumKeywords; i ++ )
      if( ! vInitializationCheck[ i ] )
      {
        cerr << "[StrainConfig] Missing variable: \'" << StrainOpt::Keywords[i] << "\' not optional"  << endl;
        bMissingTokens = true;
      }
    
    delete [] vInitializationCheck;
    if ( bMissingTokens )
      RUNTIME_ASSERT(0, "[StrainConfig] Failed to parse Strain config file\n" );
    
    return true;
  }

  //--------------------------------------
  //  Print
  //--------------------------------------
  bool StrainConfigFile::Read( const std::string & sFilename )
  {
    Size_Type  nBufferSize = 0;
    char *pBuffer = InitFileIO::ReadFileToBuf( nBufferSize, sFilename );
    
    if( nBufferSize<= 0 || !pBuffer )
    {
      std::cerr << "ERROR: Unable to read strain config file" << endl; 
      exit(0);
      return false;
    }
	
    if( ! Parse( string( pBuffer, nBufferSize ) ) )
    {
      std::cerr << "ERROR: Unable to parse strain config file" << endl;
      exit(0);
      return false;
    }
    
    if( pBuffer )
      delete[] pBuffer;
    
    return true;
  }
  
  //--------------------------------------
  //  Print
  //--------------------------------------
  void StrainConfigFile::Print( ostream & os ) const
  {
    os << StrainOpt::Keywords[StrainOpt::eStrainStepSize]    << " " << StrainSteps << std::endl;
    os << StrainOpt::Keywords[StrainOpt::eOptimizationMethod] << " " << StrainOpt::OptMethodKeywords[ OptimizationMethod ] << std::endl;
    os << StrainOpt::Keywords[StrainOpt::eStrainType] << " " << StrainOpt::StrainTypeKeywords[ StrainType ] << std::endl;
    os << StrainOpt::Keywords[StrainOpt::eStrainSearchType] << " " << StrainOpt::SearchTypeKeywords[ SearchType ] << std::endl;
    os << StrainOpt::Keywords[StrainOpt::eConvergenceMethod] << " "
       << StrainOpt::ConvergenceMethodKeywords[ ConvergenceMethod ] << " " << fConvergenceThresh <<  std::endl;
    os << StrainOpt::Keywords[StrainOpt::eMaxStrainStep]    << " " << nMaxSteps << std::endl;
  }
}
