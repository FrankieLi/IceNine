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
#include "DiscreteSearch.h"


namespace OrientationSearch
{
  namespace Utilities
  {
    //------------------------------------------------------------------------
    //
    //  Private:  GenerateLocalGrid
    //
    //  Purpose:  Construct a local, approximately Euclidean grid for high resolution,
    //            local search.  
    //
    //  Parameters:  fAngularCoverage is the "diameter" of the grid generated.
    //               A local grid of SO(3) with approximately uniform sampling is generated
    //               with this diameter.  Note that the diameter can be viewed as the
    //               diagonal of a cube on S^3 centered around the identity element.
    //               The maximum distance between two points in this approximately
    //               given by fAngularCoverage.
    //
    //               nLevel specifies the resolution of the grid.  Grid spacing is approximately
    //               fAngularCoverage / 2^nLevel for neighbors.
    //
    //  Postcondition:  vLocalGridList will be filled with orientation matrices covering
    //                  a local area
    //
    //------------------------------------------------------------------------
    void GenerateLocalGrid( vector<SMatrix3x3> &vLocalGridList,
                            Float fAngularCoverage, Int nLevel )
    {
      UniformGrid::CQuaternionGrid oGridGen;
      //-----------------------
      //  Get a uniformly structured grid
      //-----------------------
      vector<SQuaternion> oUniformStructureGrid = oGridGen.GetStructuredLocalGrid( fAngularCoverage, nLevel );  
      vLocalGridList.reserve( oUniformStructureGrid.size() );
      for ( Size_Type i = 0; i < oUniformStructureGrid.size(); i ++ )
      {
        SQuaternion q = oUniformStructureGrid[i];
        SMatrix3x3 oDelta = q.GetRotationMatrix3x3();
        vLocalGridList.push_back( oDelta );
      }
    }

    //------------------------------------------------------------------------
    //  GenerateLocalGrid
    //------------------------------------------------------------------------
    void GenerateLocalGrid( vector<SMatrix3x3> &vLocalGridList,
                            Float fAngularCoverage, Int nMinLevel,
                            Int nMaxLevel )
    {
      for( Int i = nMinLevel; i <= nMaxLevel; i ++ )
        GenerateLocalGrid( vLocalGridList, fAngularCoverage, i );
    }
  }  // namespace Utilities

  
}
