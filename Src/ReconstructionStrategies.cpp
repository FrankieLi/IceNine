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
#include "ReconstructionStrategies.h"


namespace ReconstructionStrategies
{
  namespace Utilities
  {
    //----------------------------------
    //  SubdivideTriangles
    //----------------------------------
    vector<SVoxel> SubdivideTriangles( const SVoxel & ParentVoxel, int GenMax,
                                       SVector3 PrimaryDiag,
                                       SVector3 SecondaryDiag,
                                       Float fNewSideLength )
    {
      SVector3 LeftVertex = XDMUtility::CXDMVoxel::LeftVertex( ParentVoxel );
      Int nCurrentGen = ParentVoxel.nGeneration;
      Int nGenDiff = GenMax - nCurrentGen;


      PrimaryDiag /= pow( Float(2), Float( GenMax ) );
      SecondaryDiag /= pow( Float(2), Float( GenMax ) );
  
      vector<SVoxel> ChildrenVoxels;
  
      SVector3 XStep( fNewSideLength, 0, 0 );
  
      Int nPoints = pow( 2, nGenDiff );
 
      for( Int i = 0; i < nPoints; i ++ )
      {
        for( Int j = 0; j < ( nPoints - i ) ; j ++ )
        {
          SVector3 CurrentVertex = LeftVertex + j * XStep;      
          SVoxel NewVoxel = ParentVoxel;
          NewVoxel.nGeneration = GenMax;
          NewVoxel.fSideLength = fNewSideLength;
          NewVoxel.pVertex[0] = CurrentVertex;
          NewVoxel.pVertex[1] = CurrentVertex + XStep;
          NewVoxel.pVertex[2] = CurrentVertex + PrimaryDiag;
          ChildrenVoxels.push_back( NewVoxel );
      
          if( i > 0  )
          {
            NewVoxel = ParentVoxel;
            NewVoxel.bVoxelPointsUp = ! ParentVoxel.bVoxelPointsUp;
            NewVoxel.nGeneration = GenMax;
            NewVoxel.fSideLength = fNewSideLength;
            NewVoxel.pVertex[0] = CurrentVertex;
            NewVoxel.pVertex[1] = CurrentVertex + XStep;
            NewVoxel.pVertex[2] = CurrentVertex + SecondaryDiag;
            ChildrenVoxels.push_back( NewVoxel );
          }
        }
        LeftVertex += PrimaryDiag;
      }
      return ChildrenVoxels;
    }

    //----------------------------------
    //  SetMinimumResolution
    //----------------------------------
    template<>
    void SetMinimumResolution( CMic & oSampleMic, Float fMaxVoxelSideLength )
    {
      vector<SVoxel> & vVoxelList = oSampleMic.GetVoxels();
      if( vVoxelList.size() <= 0 )
        return;
  
      Float fGen = ceil( log( oSampleMic.GetInitialSideLength() / fMaxVoxelSideLength ) / log( Float(2) ) );
      Int nMaxGen = static_cast<Int>( fGen );
  
  
      if( nMaxGen > 10 )  // sanity warning
      {
        std::cerr <<  "DANGER!!!!  You are about to regrid this sturcture by an exponent of 10!" << std::endl; 
        std::cerr <<  "This will increase your number of voxel by a factor of 4^10!  Are you SURE you want to do this!?" << std::endl;
        std::cerr <<  "If you are absolutely sure, comment out this set of error and recompile." ;
        exit( 0 );
      }
      if( nMaxGen <= 0 )
        return;

      vector<SVoxel> NewVoxels;

      Float fSideLength = oSampleMic.GetInitialSideLength() ;
      SVector3 UpDiag  ( fSideLength / Float(2), fSideLength * sqrt( Float( 3 )) / Float(2), 0 );
      SVector3 DownDiag( fSideLength / Float(2), -fSideLength * sqrt( Float( 3 )) / Float(2), 0 );
      Float fNewSideLength = fSideLength / pow( Float(2), Float( nMaxGen ) );
      for( Int i = 0; i < vVoxelList.size(); i ++ )   // a bit easier to keep track backward
      {
        if( vVoxelList[i].nGeneration < nMaxGen )
        {
          vector<SVoxel> vNewChildren;
          if( vVoxelList[i].bVoxelPointsUp )
            vNewChildren = SubdivideTriangles( vVoxelList[i], nMaxGen, UpDiag, DownDiag, fNewSideLength );
          else
            vNewChildren = SubdivideTriangles( vVoxelList[i], nMaxGen, DownDiag, UpDiag, fNewSideLength );
          NewVoxels.insert( NewVoxels.end(), vNewChildren.begin(), vNewChildren.end() );
        }
        else
        {
          NewVoxels.push_back( vVoxelList[i] );
        }
      }
      
      vVoxelList = NewVoxels;
      vector<SVoxel>( vVoxelList ).swap( vVoxelList );  // swap trick
    }
  } // ------------------------------- end of namespace Utilities
  
} // namespace Reconstructor
