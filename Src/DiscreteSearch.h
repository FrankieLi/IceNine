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
//  DiscreteSearch.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:   Implementation of the discrete orientation search.
//             
//
////////////////////////////////////////////////////////////
#ifndef _DISCRETE_SEARCH_H_
#define _DISCRETE_SEARCH_H_

#include "Sampling.h"
#include "Quaternion.h"
#include "ExperimentSetup.h"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include "MicGrid.h"
#include "Utilities.h"
#include "SearchDetails.h"
#include "OverlapInfo.h"
#include "CostFunctions.h"
#include "Symmetry.h"

namespace OrientationSearch
{
  using std::vector;
  
  namespace Utilities
  {
    //------------------------------------------------------------------------
    //
    //  Private:  GenerateGrid  (This should only be ran once)
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
    //               fAngularCoverage * sqrt(3) / 2^nLevel for neighbors.
    //
    //  Postcondition:  vLocalGridList will be filled with orientation matrices covering
    //                  a local area
    //
    //------------------------------------------------------------------------
    void GenerateLocalGrid( vector<SMatrix3x3> &vLocalGridList,
                            Float fAngularCoverage, Int nLevel );
    
    //------------------------------------------------------------------------
    //  Same as above, but generate grid from min level to max level
    //------------------------------------------------------------------------
    void GenerateLocalGrid( vector<SMatrix3x3> &vLocalGridList,
                            Float fAngularCoverage, Int nMinLevel, Int nMaxLevel );
  }  // namespace Utilities


  struct TrivialPropertyMap
  {
    template< typename T >
    inline T Get( T o ) const { return o; }
  };
  
  //------------------------------------------------------------------------
  //  LinearSearch
  //
  //  Search through a list of orientations.  An acceptance function is used
  //  to determine the best orientation.  The result is saved by a type
  //  defined by the acceptance function.
  //------------------------------------------------------------------------
  template< class SearchPointProcessor >
  class LinearSearch
  {
  private:
    SearchPointProcessor Processor;
  public:
    typedef typename SearchPointProcessor::ResultType  ResultType;
    LinearSearch( SearchPointProcessor _Processor ): Processor( _Processor ) { }
    
    //------------------------------------------------------------------------
    //  FindOptimal
    //
    //  Find optimal given the search point and local interpolation around them.
    //  The specific nature of interpolation requires this to be implemented
    //  separatly from the typical FindOptimal
    //
    //  oVoxel - input voxel to optimize the orientation
    // 
    //------------------------------------------------------------------------
    template< class IterT, class PropertyMap, class VoxelType >
    ResultType FindOptimal( const VoxelType & oVoxel, IterT pStart, IterT pEnd,
                            PropertyMap Map )
    {
      ResultType oRes;
      for ( IterT pCur = pStart; pCur != pEnd; ++pCur  )
        Processor( oRes, oVoxel, Map.Get( *pCur ) );
      return oRes;
    }
  };
  
  //------------------------------------------------------------------------
  //  Standard exhausitve search.  
  //------------------------------------------------------------------------
  template< class CostFnT, class AcceptFnT >
  class LocalAngularInterpProcess
  {
  public:
    typedef typename CostFnT::SamplePointT SamplePointT;
  private:
    CSample                  & oCurrentLayer;  // this is used for sample geometry
    const DetectorListT      & oDetectorList;
    const CSimulationRange   & oRangeToIndexMap;
    const CSimulationData    & oExpData; 
    
    CostFnT   CostFunction;    // Cost function dictates how the overlap is calculated
    AcceptFnT AcceptanceFn;
    LocalAngularInterpProcess();  // this is not legal

    const CSimulation & Simulator;

    typedef vector<CRecpVector>::iterator RecpIterT;
    
  public:
    typedef SearchDetails::SCandidate   SCandidate;
    typedef CostFunctions::SOverlapInfo SOverlapInfo;
    typedef SearchDetails::ScatteringVectorListT ScatteringVectorListT;

    typedef vector<SCandidate> ResultType;

    //-------------------------------------
    //  ctor
    //  -- Full constructor
    //-------------------------------------
    LocalAngularInterpProcess( CSample & _Sample, const DetectorListT & _DetList,
                               const CSimulationRange & _RMap,
                               const CSimulationData & _Data,
                               CostFnT   _CostFn,
                               AcceptFnT _AcceptFn,
                               const CSimulation & _Sim ):
      oCurrentLayer( _Sample ), oDetectorList( _DetList ), oRangeToIndexMap( _RMap ),
      oExpData( _Data ), 
      CostFunction( _CostFn ), AcceptanceFn( _AcceptFn ), Simulator( _Sim )
    {    }

    //----------------------------------------------
    //  Process the search center for voxel
    //
    //  Parameters:
    //
    //
    //   pFirst and pEnd are the first and last iterator to the list of
    //   reciprocal lattice vectors.
    //----------------------------------------------
    void operator()( vector<SCandidate> & oRes, const SamplePointT & oVoxel,
                     const SMatrix3x3 & oSearchCenter,
                     const vector<SMatrix3x3> & oLocalGridList,
                     RecpIterT pFirst, RecpIterT pEnd )
    {
      ScatteringVectorListT oScatteringCandidates;
      HEDM::PrecomputeScatteringVector( oScatteringCandidates,
                                        oSearchCenter, oDetectorList,
                                        oRangeToIndexMap, pFirst,
                                        pEnd, Simulator );
      
      for ( Size_Type j = 0; j < oLocalGridList.size(); j ++ )
      {
        SOverlapInfo oOverlapInfo;
        oOverlapInfo.Initialize();
  
        oOverlapInfo = CostFunction ( oVoxel, oCurrentLayer,
                                      oDetectorList,
                                      oRangeToIndexMap,
                                      oScatteringCandidates,
                                      oLocalGridList[j],
                                      oExpData );

        if ( AcceptanceFn( oOverlapInfo ) )   // If accepted
          oRes.push_back( SCandidate( oLocalGridList[j] * oSearchCenter, 0 ) );
        
      }
    }
    
  };

  //------------------------------------------------------
  //
  //   DiscreteSearch
  //   --  Definition of the discrete search algorithm
  //------------------------------------------------------
  template< class VoxelCostFn, class CandidateAcceptFn >
  class DiscreteSearchFn
  {
  public:
    typedef OrientationSearch::LocalAngularInterpProcess< VoxelCostFn, CandidateAcceptFn >      SearchPointProcessor;

    typedef typename VoxelCostFn::SamplePointT SamplePointT;
    typedef SearchDetails::SCandidate          SCandidate;
    typedef vector<SCandidate>::iterator       CandidateIter;
  protected:
    vector<CRecpVector>       oRecipVectors;
    const vector<SMatrix3x3>  & oLocalGrid;
    SearchPointProcessor      oSearchProcess;
    
    template< typename SymmetryT >
    bool Acceptable( CandidateIter pCur, CandidateIter pLast,
                     const SCandidate & c, Float fMinAngle,
                     const SymmetryT & oSym )
    {
      for(; pCur != pLast; ++ pCur )
      {
        Float fMis = LatticeSymmetry::GetMisorientation( oSym,
                                                         pCur->oOrientation,
                                                         c.oOrientation );
        if( fMis < fMinAngle  && pCur->fCost < c.fCost )  // if there exists a better mattch
          return false;
      }
      return true;
    }
    
  public:
    
    DiscreteSearchFn( const CSimulation & oSimulator,
                      CSample & _Sample, const DetectorListT & _DetList,
                      const CSimulationRange & _RMap,
                      const CSimulationData & _Data,
                      const vector<SMatrix3x3> & _LocalGrid,
                      Int nPhase, VoxelCostFn oVoxelCostFunction,
                      CandidateAcceptFn oAcceptFn )
      : oLocalGrid( _LocalGrid ),
        oSearchProcess( _Sample, _DetList, _RMap, _Data,
                        oVoxelCostFunction, oAcceptFn, oSimulator )
    {
      using boost::tuple;
      using namespace boost::lambda;
      using boost::get;

      //----------------------------------------
      //   Need to change this to get to strain compatible - need a template parameter
      const vector<CUnitCell> & oCryStructList = _Sample.GetStructureList();
      oRecipVectors = oCryStructList[ nPhase ].GetReflectionVectorList();  
      //----------------------------------------
      std::sort( oRecipVectors.begin(), oRecipVectors.end(),
                 [](const CRecpVector& a, const CRecpVector& b) { return a.fMag < b.fMag; } );
    }
    
    //--------------------------------------------------
    //  GetCandidates
    //--------------------------------------------------
    template < typename  SearchPointIterT >
    vector<SCandidate> GetCandidates( const SamplePointT & oVoxel,
                                      SearchPointIterT pSearchPoint,
                                      SearchPointIterT pEnd,
                                      Float fQMax )
    {
      using boost::tuple;
      using namespace boost::lambda;
      using boost::get;
      
      typedef vector<CRecpVector>::iterator RecpIter;
      RecpIter pRecpEnd = std::find_if( oRecipVectors.begin(),
                                        oRecipVectors.end(),
                                        [fQMax](const CRecpVector& v) { return v.fMag > fQMax; } );
      
      typedef OrientationSearch::TrivialPropertyMap  PropMap;
      vector<SCandidate> oRes;
      for(; pSearchPoint != pEnd; ++ pSearchPoint )
        oSearchProcess( oRes, oVoxel, *pSearchPoint,
                        oLocalGrid, oRecipVectors.begin(), pRecpEnd );
      return oRes;
    }


    //--------------------------------------------------
    //  GetSpacedCandidates
    //--------------------------------------------------
    template < typename SearchPointIterT,
               typename CompleteCostFn,  typename SymmetryT >
    vector< vector<SCandidate> > GetSpacedCandidates( const SamplePointT & oVoxel,
                                                      SearchPointIterT pSearchPoint,
                                                      SearchPointIterT pEnd,
                                                      Float fQMax,
                                                      CompleteCostFn CompleteVoxelCostFn,
                                                      Float fAngularRadius,
                                                      SymmetryT & oSym )
    {
      using boost::tuple;
      using namespace boost::lambda;
      using boost::get;
      
      typedef vector<CRecpVector>::iterator RecpIter;
      RecpIter pRecpEnd = std::find_if( oRecipVectors.begin(),
                                        oRecipVectors.end(),
                                        [fQMax](const CRecpVector& v) { return v.fMag > fQMax; } );
      
      typedef OrientationSearch::TrivialPropertyMap  PropMap;
      vector< vector< SCandidate> > oRes;
      
      int nClique = 0;
      int nMaxRangeID = 0;
      int nMinCostID  = 0;
      int nCurrentID  = 0;
      Float fMaxRange = std::numeric_limits<Float>::min();
      Float fMinCost  = std::numeric_limits<Float>::max();

      GeneralLib::SVector3 BestOrientation;
      GeneralLib::SVector3 MaxRangeOrientation;

      
      for(; pSearchPoint != pEnd; ++ pSearchPoint )
      {
        vector<SCandidate> oLocalCandidates;
        oSearchProcess( oLocalCandidates, oVoxel, *pSearchPoint,
                        oLocalGrid, oRecipVectors.begin(), pRecpEnd );
        
        for( Size_Type i = 0; i < oLocalCandidates.size(); i ++ )
        {
          SamplePointT oTmp = oVoxel;
          oTmp.oOrientMatrix = oLocalCandidates[i].oOrientation;
          CostFunctions::SOverlapInfo oInfo = CompleteVoxelCostFn( oTmp );
          oLocalCandidates[i].fCost = 1 - CostFunctions::Utilities::GetConfidence( oInfo ); 
        }
	
        if( oLocalCandidates.size() > 0 )
        {
          nClique ++;
          //	  std::cout << "| nClique " << nClique << " size | "  << oLocalCandidates.size() << std::endl;
          CandidateIter pCur      = oLocalCandidates.begin();
          CandidateIter pFirstGood = pCur;
          pCur ++;
          Float fMin = std::numeric_limits<Float>::max();
          Float fMax = std::numeric_limits<Float>::min();
          while( pCur != oLocalCandidates.end() )
          {
            if( ! Acceptable( pFirstGood, oLocalCandidates.end(),
                              *pCur, fAngularRadius, oSym ) )
            {
              std::swap( *pFirstGood, *pCur );
              ++ pFirstGood;
            }
	    
            if( pCur->fCost < fMinCost )
            {
              fMinCost = pCur->fCost;
              nMinCostID = nCurrentID;
              //              BestOrientation = RADIAN_TO_DEGREE( pCur->oOrientation.GetEulerAngles() );
            }
            // if acceptable
            fMin = std::min( pCur->fCost, fMin );
            fMax = std::max( pCur->fCost, fMax );
            ++pCur;
          }
          if( fMaxRange < (fMax - fMin) )
          {
            //            MaxRangeOrientation = RADIAN_TO_DEGREE( pFirstGood->oOrientation.GetEulerAngles() );
            nMaxRangeID = nCurrentID;
          }	  
          nCurrentID ++;
          fMaxRange = std::max( fMaxRange, (fMax - fMin ) );
          vector<SCandidate> oAccepted;   // new list of acceted candidates for this region
          oRes.push_back( oAccepted );
          oRes.back().insert( oRes.back().begin(), pFirstGood, pCur );
        }
      }
      
        std::cout << "Num Cliques " << nClique << " MaxRange " << fMaxRange 
		<< " MaxRangeCliqueID: " << nMaxRangeID 
		<< " MinCostCliqueID " << nMinCostID 
		<< " " << fMinCost 
		<< " | " << BestOrientation << " | " << MaxRangeOrientation <<  std::endl;
        
      //      exit(0);
      return oRes;
    }
    


    //--------------------------------------------------
    //  GetSpacedCandidates2  - uses a rigid grid to resample
    //--------------------------------------------------
    template < typename SearchPointIterT,
               typename CompleteCostFn,  typename SymmetryT >
    vector< vector<SCandidate> > GetSpacedCandidates2( const SamplePointT & oVoxel,
                                                       SearchPointIterT pSearchPoint,
                                                       SearchPointIterT pEnd,
                                                       Float fQMax,
                                                       CompleteCostFn CompleteVoxelCostFn,
                                                       Float fAngularRadius,
                                                       SymmetryT & oSym )
    {
      using boost::tuple;
      using namespace boost::lambda;
      using boost::get;
      
      typedef vector<CRecpVector>::iterator RecpIter;
      RecpIter pRecpEnd = std::find_if( oRecipVectors.begin(),
                                        oRecipVectors.end(),
                                        [fQMax](const CRecpVector& v) { return v.fMag > fQMax; } );
      
      typedef OrientationSearch::TrivialPropertyMap  PropMap;
      vector< vector< SCandidate> > oRes;
      
      
      for(; pSearchPoint != pEnd; ++ pSearchPoint )
      {
        vector<SCandidate> oLocalCandidates;
        oSearchProcess( oLocalCandidates, oVoxel, *pSearchPoint,
                        oLocalGrid, oRecipVectors.begin(), pRecpEnd );
        
        for( Size_Type i = 0; i < oLocalCandidates.size(); i ++ )
        {
          SamplePointT oTmp = oVoxel;
          oTmp.oOrientMatrix = oLocalCandidates[i].oOrientation;
          CostFunctions::SOverlapInfo oInfo = CompleteVoxelCostFn( oTmp );
          oLocalCandidates[i].fCost = 1 - CostFunctions::Utilities::GetConfidence( oInfo ); 
        }
        
        if( oLocalCandidates.size() > 0 )
        {
          CandidateIter pCur      = oLocalCandidates.begin();
          CandidateIter pFirstGood = pCur;
          pCur ++;
          Float fMin = std::numeric_limits<Float>::max();
          Float fMax = std::numeric_limits<Float>::min();
          while( pCur != oLocalCandidates.end() )
          {
            if( ! Acceptable( pFirstGood, oLocalCandidates.end(),
                              *pCur, fAngularRadius, oSym ) )
            {
              std::swap( *pFirstGood, *pCur );
              ++ pFirstGood;
            }
            ++pCur;
          }
          vector<SCandidate> oAccepted;   // new list of acceted candidates for this region
          oRes.push_back( oAccepted );
          oRes.back().insert( oRes.back().begin(), pFirstGood, pCur );
        }
      }
      return oRes;
    }
  };
  
}  // namespace Orientation Search


#include "DiscreteSearch.tmpl.cpp"

#endif
