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
//  OptimizationConfig.h
//
//  Classes:  OptimizatonConfig
//
//  Purpose:  Implmentation of parsers and containers for
//            different optimization methods
//
///////////////////////////////////////////////////////

#ifndef ICENINE_OPT_CONFIG_H
#define ICENINE_OPT_CONFIG_H
#include "Types.h"
#include "3dMath.h"
#include "Debug.h"
#include "Parser.h"
#include <string>
#include "InitFilesIO.h"
#include <vector>
#include <tuple>
#include <sstream>


namespace OptimizationConfig
{
  namespace StrainOpt
  {
    enum EOptimizationMethod{
      eMonteCarlo, eConjugateGradient,
      eNumOptimizationMethods
    };
    
    static const char* OptMethodKeywords[] =
      { "MonteCarlo", "ConjugateGradient" };
    
    enum EStrainType{
      eDiag, eSym, eAsym, eNumStrainTypes
    };
    
    static const char* StrainTypeKeywords[] =
      { "Diag", "Sym", "Asym" };
    
    enum ESearchType{
      eGlobal, eLocal, eMix, eNumStrainSearches
    };
    
    static const char* SearchTypeKeywords[] =
      { "Global", "Local", "Mix" };
    
    enum EConvergenceMethod{
      eQuality, eHitRatio, eDifferential, eNumConvergenceMethods 
    };

    static const char* ConvergenceMethodKeywords[] =
      { "Quality", "HitRatio", "Differential" };
    
    enum EKeywords{
      eStrainSearchType, eStrainType, eOptimizationMethod,
      eStrainStepSize, eMaxStrainStep, eConvergenceMethod,
      eComment, eNumKeywords
    };

    static const char* Keywords[] =
      { "StrainSearchType", "StrainType", "OptimizationMethod",
        "StrainStepSize", "MaxStrainStep", "ConvergenceMethod", "#" };

    EOptimizationMethod ProcessOptMethod( const vector< std::string > & Line,
                                          int nLineNumber );
    EStrainType         ProcessStrainType( const vector< std::string > & Line,
                                           int nLineNumber );
    ESearchType         ProcessSearchType( const vector< std::string > & Line,
                                           int nLineNumber );
    std::pair<EConvergenceMethod, Float> ProcessConvergenceMethod( const vector< std::string > & Line,
                                                                   int nLineNumber );
    GeneralLib::SMatrix3x3 ParseMatrix( const vector< vector< string > > & vsTokens,
                                        int nLineNumber );
  };
  
  //--------------------------------------
  //  StrainConfigFile
  //--------------------------------------
  class StrainConfigFile
  {
  private:    
    bool Parse( const std::string & sBuf );

  public:
    
    GeneralLib::SMatrix3x3         StrainSteps;
    StrainOpt::EOptimizationMethod OptimizationMethod;
    StrainOpt::EStrainType         StrainType;
    StrainOpt::ESearchType         SearchType;
    StrainOpt::EConvergenceMethod  ConvergenceMethod;

    Float fConvergenceThresh;
    Int   nMaxSteps;
    
    bool Read( const std::string & sFilename );
    void Print( ostream & os ) const;
  };
  
  
}

#endif
