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
//  CostFunctions.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//   Purpose:  Definition of Cost functions and/or functionals
//
//
////////////////////////////////////////////////////////////


#ifndef _COST_FUNCTIONS_H_
#define _COST_FUNCTIONS_H_

#include "OrientationSearch.h"
#include "3dMath.h"
#include "OverlapInfo.h"
#include "SearchDetails.h"
#include "DiffractionCore.h"
#include "SimulationData.h"
#include "Sample.h"

namespace CostFunctions
{
  namespace Details
  {
    class DefaultReflectionMap
    {
    public:
      DefaultReflectionMap(){}

      template< class SamplePoint >
      const vector<CRecpVector> & operator() ( const CSample & s,
                                               const SamplePoint & v ) const
      {
        DEBUG_ASSERT( v.nPhase < static_cast<int>( s.GetStructureList().size() ),
                      "BOUND ERROR!  nPhase >= size of Structure list\n" );
        return s.GetStructureList()[ v.nPhase ].GetReflectionVectorList();
      }
    };
    
    
  }
  
  //--------------------------------
  //  CostFunction Utilities
  //--------------------------------
  namespace Utilities
  {
    //--------------------------------------------------
    //  Classification of convergence
    //--------------------------------------------------
    template< typename SamplePoint >
    Int ClassifyConvergence( SamplePoint & Voxel,
                             const SearchDetails::SSearchParameter & oSearchParam )
    {
      Int ResultCode;
      if ( oSearchParam.fMaxConvergenceCost > Voxel.fCost )
      {
        ResultCode    = SearchDetails::CONVERGED;
        Voxel.nPhase = 1;   
      }
      else if( oSearchParam.fMaxAcceptedCost > Voxel.fCost )
      {
        ResultCode  = SearchDetails::PARTIAL;
        Voxel.nPhase = 1;   
      }
      else
      {
        ResultCode  = SearchDetails::NON_ACCEPTED;
        Voxel.nPhase = 0;   
      }
      return ResultCode;
    }
    
    //--------------------------------------------------
    //  Classification of convergence
    //--------------------------------------------------
    template< typename CostFunction, typename SamplePoint >
    Int ClassifyConvergence( SamplePoint & Voxel, CostFunction  & oOverlapFn,
                             const SearchDetails::SSearchParameter & oSearchParam )
    {
      return ClassifyConvergence( Voxel, oOverlapFn( Voxel ), oSearchParam );
    }

    //--------------------------------------------------
    //  Classification of convergence
    //--------------------------------------------------
    template< typename CostFunction, typename SamplePoint >
    vector<Int>
    ClassifyConvergence( SamplePoint & Voxels, CostFunction  & oOverlapFn,
                         const SearchDetails::SSearchParameter & oSearchParam )
    {
      vector<Int> ReturnCodes;
      for( Size_Type i = 0; i < Voxels.size(); i ++ )
        ReturnCodes.push_back( ClassifyConvergence( Voxels[i], oOverlapFn, oSearchParam ) );

      return ReturnCodes;
    }
    //--------------------------------------------------
    //  GetConfidence
    //--------------------------------------------------
    Float GetConfidence( const CostFunctions::SOverlapInfo & o );
    
    //--------------------------------------------------
    //  GetHitRatio
    //--------------------------------------------------
    Float GetHitRatio( const CostFunctions::SOverlapInfo & o );
     
  } //  end Utilities
  
  struct VoxelToCenter
  {
    template< class SamplePoint >
    inline vector<SVector3> operator()( const SamplePoint & oVoxel ) const
    {
      return vector<SVector3>( 1, oVoxel.GetCenter() );
    }
  };

  //----------------------------------------------
  //  Leave as explicit
  //----------------------------------------------
  struct VoxelToVertices
  {
    inline vector<SVector3> operator()( const SVoxel & oVoxel ) const
    {
      vector<SVector3> oRes(3);
      oRes[0] = oVoxel.pVertex[0];
      oRes[1] = oVoxel.pVertex[1];
      oRes[2] = oVoxel.pVertex[2];
      return oRes;
    }

    inline vector<SVector3> operator()( const SquareVoxel & oVoxel ) const
    {
      vector<SVector3> oRes(4);
      oRes[0] = oVoxel.pVertex[0];
      oRes[1] = oVoxel.pVertex[1];
      oRes[2] = oVoxel.pVertex[2];
      oRes[3] = oVoxel.pVertex[3];
      return oRes;
    }
  };
  
  //---------------------------------------------------------------------------------------
  //  InterploatedVoxelCostFunction
  //
  //
  //  Purpose:  Calculate the cost for the voxel by interpolating the calculation of the diffraction
  //            angles.
  //
  //
  //---------------------------------------------------------------------------------------
  template< class SampleVertexOverlapFn, class VertexExtractorT,
            class SamplePoint  >
  class InterpolatedVoxelCostFunction
  {
  public:
    typedef SamplePoint SamplePointT;
    
  private:
    SampleVertexOverlapFn OverlapFn;
    VertexExtractorT VoxelToVerticesFn;

    InterpolatedVoxelCostFunction();
  public:
    InterpolatedVoxelCostFunction( SampleVertexOverlapFn _f, VertexExtractorT _v ):
      OverlapFn( _f ), VoxelToVerticesFn( _v ) { }
    
    InterpolatedVoxelCostFunction( SampleVertexOverlapFn _f ):
      OverlapFn( _f ), VoxelToVerticesFn() { }
    
    SOverlapInfo operator()( const SamplePoint & oVoxel,
                             CSample & oSample,
                             const DetectorListT & oDetList,
                             const XDMSimulation::CSimulationRange & oRangeToIndexMap,
                             const SearchDetails::ScatteringVectorListT & oScatteringCandidates,
                             const SMatrix3x3 & oSmallRotation,
                             const CSimulationData  & oExpData )
    {
      vector<SPeakInfo> oPeakInfoList;
      HEDM::InterpolateObservablePeaks( oPeakInfoList, oScatteringCandidates, oSmallRotation );
      vector<SVector3> oVertexList = VoxelToVerticesFn( oVoxel );
      return OverlapFn.CalculateDiffractionOverlap( oVertexList, oPeakInfoList,
                                                    oSample, oDetList,
                                                    oRangeToIndexMap, oExpData );
    }
  };
  
  //---------------------------------------------------------------------------------------
  //  GeneralVoxelCostFunction
  //
  //
  //  Purpose:  Calculate the cost for the voxel given using *exact* calculation of the diffraction
  //            angles, with generalization of Refelection (reciprocal vector) map.
  //
  //---------------------------------------------------------------------------------------
  template< class SampleVertexOverlapFn, class VertexExtractorT,
            class ReflectionVectorMap, class SamplePoint >
  class GeneralVoxelCostFunction
  {
  public:
    typedef SamplePoint SamplePointT;
    
  protected:
    SampleVertexOverlapFn  OverlapFn;
    VertexExtractorT  VoxelToVerticesFn;
    ReflectionVectorMap RecpVectorMap;
    GeneralVoxelCostFunction();
    
  public:
    GeneralVoxelCostFunction(  SampleVertexOverlapFn _f ):
      OverlapFn( _f ), VoxelToVerticesFn() { }
    
    virtual SOverlapInfo operator()( const SamplePoint & oVoxel, CSample & oSample,
                                     const DetectorListT & oDetList,
                                     const XDMSimulation::CSimulationRange & oRangeToIndexMap,
                                     const CSimulationData & oExpData ) = 0;
    
  };
  
  //---------------------------------------------------------------------------------------
  //  VoxelCostFunction
  //
  //
  //  Purpose:  Calculate the cost for the voxel given using *exact* calculation of the diffraction
  //            angles.
  //
  //  Note:     By implementing a different VoxelCostFunction, say, one with adjustable MaxQ,
  //            we can arbitarily choose to calculate cost based on different peaks.
  //---------------------------------------------------------------------------------------
  template< class SampleVertexOverlapFn, class VertexExtractorT, class SamplePoint >
  class VoxelCostFunction :
    public GeneralVoxelCostFunction< SampleVertexOverlapFn,
                                     VertexExtractorT,
                                     Details::DefaultReflectionMap,
                                     SamplePoint >
  {
  public:
    typedef SamplePoint SamplePointT;
    typedef  GeneralVoxelCostFunction< SampleVertexOverlapFn,
                                       VertexExtractorT,
                                       Details::DefaultReflectionMap, SamplePointT >  Base;
    
    VoxelCostFunction(  SampleVertexOverlapFn _f ):
      Base(  _f ) { }
    
    //-----------------------------------------
    //  Computing Operator
    //-----------------------------------------
    SOverlapInfo operator()( const SamplePoint & oVoxel, CSample & oSample,
                             const DetectorListT & oDetList,
                             const XDMSimulation::CSimulationRange & oRangeToIndexMap,
                             const CSimulationData & oExpData )
    {
      const vector<CRecpVector> & oRecipVectors = Base::RecpVectorMap( oSample, oVoxel );
      vector<SPeakInfo> oPeakInfoList;      
      oPeakInfoList.reserve( oRecipVectors.size() * 2 );
      Base::OverlapFn.Simulator().GetObservablePeaks( oPeakInfoList, oVoxel.oOrientMatrix,
                                                oRecipVectors, oDetList );
      
      vector<SVector3> oVertexList = Base::VoxelToVerticesFn( oVoxel );
      return Base::OverlapFn.CalculateDiffractionOverlap( oVertexList, oPeakInfoList,
                                                          oSample, oDetList,
                                                          oRangeToIndexMap, oExpData );
      
    }
  };

  
}// CostFunctions namespace


#endif
