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
//  Utilities.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//  Utilities used by IceNine
//
//
////////////////////////////////////////////////////////////


#ifndef _I9_UTILITIES_H_
#define _I9_UTILITIES_H_

#include "MicGrid.h"
#include "3dMath.h"
#include "Symmetry.h"
#include "Voxel.h"


namespace XDMUtility
{

  //----------------------------------
  // ReduceToFundamentalZone
  //----------------------------------
  template< class MicFileType, class Symmetry >
  void ReduceToFundamentalZone( MicFileType & MicFile, const Symmetry & oSym )
  {
    typedef typename MicFileType::VoxelType_iterator VoxelType_iterator;
    for( VoxelType_iterator it = MicFile.VoxelListBegin(); it != MicFile.VoxelListEnd(); ++ it )
    {
      GeneralLib::SQuaternion qOrient;
      qOrient.Set( it->oOrientMatrix );
      qOrient = LatticeSymmetry::ReduceToFundamentalZone( oSym, qOrient );
      it->oOrientMatrix = qOrient.GetRotationMatrix3x3();
    }
  }
  
  //----------------------------------
  //  GetBoundaryVoxels
  //
  //
  //  fMaxCost        - maximum cost of the voxel to be considered
  //  fMisorientation - minimum misorientation to be counted
  //  fRadius         - Distance from center of mass of the sample
  //  
  //
  //----------------------------------
  template< class Symmetry >
  vector<SVoxel> GetBoundaryVoxels( const Symmetry & oSym,
                                    Float fInitialSideLength,
                                    const vector<SVoxel> & oVoxelList,
                                    Float fMisorientation,
                                    Float fMaxCost, Float fRadius )
  {
    MicAnalysis::CMicGrid oMicGrid;
    oMicGrid.InsertReplace( fInitialSideLength, oVoxelList.begin(), oVoxelList.end() );
    

    SVector3 oCenterOfMass( 0, 0, 0);
    Int nVoxelsCounted = 0;
    for( Size_Type i = 0; i < oVoxelList.size(); i ++ )
    {
      if( oVoxelList[i].fConfidence >= ( Float(1) - fMaxCost ) )
      {
        oCenterOfMass += oVoxelList[i].GetCenter();
        nVoxelsCounted ++;
      }
    }
    
    oCenterOfMass /= Float( nVoxelsCounted );
    
    typedef MicAnalysis::CMicGrid::ShapePtr ShapePtr;
    vector< SVoxel > oBndVoxels;
    for( Size_Type i = 0; i < oVoxelList.size(); i ++ )
    {
      Float fR =  ( oVoxelList[i].GetCenter() - oCenterOfMass ).GetLength();
      if( fR <= fRadius )
      {
        vector<ShapePtr> oNgbList = oMicGrid.GetNeighbors( oVoxelList[i] );
        for( Size_Type j = 0; j < oNgbList.size(); j ++ )
        {
          if( LatticeSymmetry::GetMisorientation( oSym, oNgbList[j]->oOrientMatrix,
                                                  oVoxelList[i].oOrientMatrix ) > fMisorientation )
          {
            oBndVoxels.push_back( oVoxelList[i] );
            break;
          }
        }
      }
    }
    return oBndVoxels;
  }


  //----------------------------------
  //  FilteredOutStream
  //----------------------------------
  class FilteredOutStream
  {
  private:
    bool          bWritable;     // if this is a file to be outputted
  public:
    FilteredOutStream(): bWritable( false ) {}

    bool Writable() const
    {
      return bWritable;
    }
    
    void SetWritable()
    {
      bWritable = true;
    }
    std::ostream & Get()
    {
      return std::cout;
    }
    
  };
  
  //----------------------------------
  //  LogStream 
  //----------------------------------
  class LogStream
  {
  private:
    std::string   sFilename;
    std::ofstream oOutStream;
    bool          bWritable;     // if this is a file to be outputted
  public:
    LogStream(): bWritable( false ) {}
    LogStream( const std::string & s_ ):sFilename( s_ ), oOutStream(), bWritable( false ) {}

    bool Writable() const
    {
      return bWritable;
    }
    
    void SetWritable()
    {
      bWritable = true;
      oOutStream.open( sFilename.c_str() );
    }

    void open( const std::string & s_ )
    {
      sFilename = s_;
      if( bWritable )
        oOutStream.open( sFilename.c_str() );
    }
    
    void close( )
    {
      if(  bWritable )
        oOutStream.close();
    }

    std::ostream & Get()
    {
      return oOutStream;
    }

    
  };

#ifndef GET_LOG
#define GET_LOG( os )                        \
  if( ( !os.Writable() ) ) ;                 \
  else  os.Get() 
#endif

}  // end namespace


#endif
