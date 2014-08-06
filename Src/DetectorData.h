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
//
//   DetectorData.h
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:  This is a (abstract) base class which specifies the interface of different
//             representation of the data.  In reality, the data is stored in
//             TIFF formatted images.  Each image has a specific size.  However,
//             based on the specific data reduction technique and reconstruction
//             algorithm, different representation of the data may be needed.  For
//             example, bounding boxes represented using a quadtree is useful for
//             searching, while sparse matrix will be perfect for memory limited 
//             storage.  DetectorData therefore specifies a few limited functions
//             that's implemented by the base class.
//
//
//             Note that we're adopting JPixel and KPixel notation to prevent confusion
//             with 'rows' and 'columns'
//
////////////////////////////////////////////////////////////


#ifndef _DETECTOR_DATA_H_
#define _DETECTOR_DATA_H_

#include "Types.h"


class CDetectorData
{
public:
  // ------------------------
  //  IsDark
  //  Return true if the pixel at the location specified is not lit
  //  (Not lit at all)
  // ------------------------
  virtual Bool IsDark( Int nJPixel, Int nKPixel ) const = 0;

  // ------------------------
  //  IsBright
  //  Return true if the pixel at the location specified is lit
  //  (This is just for code beautification)
  // ------------------------
  virtual Bool IsBright( Int nJPixel, Int nKPixel  ) const = 0;

  virtual ~CDetectorData() {};
};


#endif
