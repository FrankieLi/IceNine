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
//   PeakFilters.h 
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//   Purpose:  Definition of various peak filters used in the simulation
//
//
////////////////////////////////////////////////////////////
#ifndef _PEAKFILTERS_H_
#define _PEAKFILTERS_H_

namespace HEDM
{

  struct XDMPeakFilter
  {
    void SetSin2Theta( Float f) {  }
  };
  //--------------------------------------------------------------------------------------------------------
  //  Angular Accept    Filters to be removed to another location soon
  //
  //--------------------------------------------------------------------------------------------------------
  struct XDMEtaAcceptFn : public XDMPeakFilter
  {
    Float fMinEta;
    Float fMaxEta;
    Float fFormIntensity;   // form factor contribution of the intensity
    Float fSin2Theta;
  
    XDMEtaAcceptFn(  ): fMinEta( 0 ), fMaxEta( 0 ), fSin2Theta(1){};
    XDMEtaAcceptFn( Float _fMinEta, Float _fMaxEta ): fMinEta( _fMinEta ), fMaxEta( _fMaxEta ){};
  
    std::pair<Bool, Float> operator() ( const SVector3 & oScatteringDir )
    {
      Float fEta = atan2( fabs( oScatteringDir.m_fY ), fabs( oScatteringDir.m_fZ ) );
      bool bAccept = false;
      if( fEta < fMaxEta )
        bAccept = true;
    
      Float fIntensity = fFormIntensity / ( fabs( sin( fEta ) ) * fSin2Theta );
    
      return std::pair<Bool, Float> ( bAccept, fIntensity );
    }

    void SetSin2Theta( Float f )
    {
      fSin2Theta = f;
    }
  
  };

  //--------------------------------------------------------------------------------------------------------
  //  Angular Accept
  //
  //--------------------------------------------------------------------------------------------------------
  struct XDMSimpleEtaFilter : public XDMPeakFilter
  {
    Float fMinEta;
    Float fMaxEta;
    std::pair<bool, Float> operator() ( const SVector3 & oScatteringDir )
    {
      Float fEta = atan2( fabs( oScatteringDir.m_fY ), fabs( oScatteringDir.m_fZ ) );
      if( fEta < fMaxEta )
        return std::make_pair( true, 1 );
      return std::make_pair( false, 0 );
    }
  };

}


#endif
