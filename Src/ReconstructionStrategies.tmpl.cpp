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



namespace ReconstructionStrategies
{
  //-----------------------------------------------------------------
  //  RestrictedStratifiedGrid
  //-----------------------------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void RestrictedStratifiedGrid< SamplePointT, SamplePointGrid >
  ::Initialize( const Mic & RecRegion_,
                Float fMinSideLength_,
                SVector3 oSampleCenter_,
                Float fSampleRadius_ )
  {
    ReconstructionRegion = RecRegion_;
    Utilities::SetMinimumResolution( ReconstructionRegion, fMinSideLength_ );
    SolutionGrid.InitializeWithReference( ReconstructionRegion );

    fSampleRadius = fSampleRadius_;
    oSampleCenter = oSampleCenter_;
    this->InitializeMicGridID( SolutionGrid, MultiStagedDetails::NOT_VISITED );
    
    typedef typename Mic::VoxelType_iterator Iter;
    for( Iter pCur = ReconstructionRegion.VoxelListBegin();
         pCur !=  ReconstructionRegion.VoxelListEnd(); ++ pCur )
    {
      SamplePointPtr pToInsert( & (* pCur ), Utilities::Null_Deleter() );
      SamplePointPtrList.push_back( pToInsert );
    }
    
    std::random_shuffle( SamplePointPtrList.begin(), SamplePointPtrList.end() );
    nElementsLeft       = SamplePointPtrList.size();
    pCurrentSamplePoint = SamplePointPtrList.begin();
    bFirstPass          = true;
  }


   //-----------------------------------------------------------------
  //  RestrictedStratifiedGrid::InitializeRestart
  //-----------------------------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void RestrictedStratifiedGrid<SamplePointT, SamplePointGrid>
  ::InitializeRestart( const Mic & RecRegion_,
                       Mic & RestartPoint,
                       Float fMinAcceptThresh,   // temporarily not using this - may need to add this back to improve efficiency
                       Float fMinSideLength_,
                       SVector3 oSampleCenter_,
                       Float fSampleRadius_ )
  {
    Initialize( RecRegion_, fMinSideLength_, oSampleCenter_, fSampleRadius_ );
    
    typedef typename Mic::VoxelType_const_iterator  const_iter;
    int nPartialVoxelsAdded = 0;
    if( RestartPoint.GetVoxels().size() > 0 )
    {
      Utilities::SetMinimumResolution( RestartPoint, fMinSideLength_ );  // make sure that side length is correct
      for( const_iter pCur = RestartPoint.VoxelListBegin();
           pCur !=  RestartPoint.VoxelListEnd(); ++ pCur )
      {
        SamplePointPtr pCenter = SolutionGrid( *pCur );
        RUNTIME_ASSERT( SamplePointGrid::IsValid( pCenter ),
                        "[InitializeRestart:] Unexpected result from outside the solution region\n" );        
        *pCenter      = *pCur;
        //       pCenter->nID = MultiStagedDetails::NOT_VISITED;   // should allow option to set this
        pCenter->nID = MultiStagedDetails::NOT_VISITED;
        nPartialVoxelsAdded ++;
      }
    }
    
    std::cout << "PartialVoxel Added: " << nPartialVoxelsAdded << std::endl;
  }
  
  //------------------------------------
  //  Reset
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void RestrictedStratifiedGrid<SamplePointT, SamplePointGrid>::Reset()
  {
    this->InitializeMicGridID( SolutionGrid, MultiStagedDetails::NOT_VISITED );
    pCurrentSamplePoint = SamplePointPtrList.begin();
    bFirstPass = true;
  }

  //------------------------------------
  //  Reset
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void RestrictedStratifiedGrid<SamplePointT, SamplePointGrid>
  ::Push( const SamplePointT & oResult )
  {
    SamplePointPtr pVoxel = SolutionGrid( oResult );
    RUNTIME_ASSERT( SamplePointGrid::IsValid( pVoxel ),
                    "Unexpected voxel outside of intended fit area\n");

    if( pVoxel->nID == MultiStagedDetails::VISITED ||
        pVoxel->nID == MultiStagedDetails::NOT_VISITED )   // 0 confidence case
    {
      nElementsVisited ++;
      *pVoxel = oResult;
      SolutionPointPtrList.push_back( pVoxel );
    }
    else if( pVoxel->fConfidence < oResult.fConfidence )
    {   
      *pVoxel = oResult;
      
      if( pVoxel->nID == MultiStagedDetails::FITTED )
        nElementsFitted ++;
    }
  }
  
  //------------------------------------
  //  Pop  (removal of first of the queue)
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void RestrictedStratifiedGrid<SamplePointT, SamplePointGrid>::Pop( )
  {
    if( pCurrentSamplePoint != SamplePointPtrList.end() )    // remove current voxelxs
    {
      (*pCurrentSamplePoint)->nID = MultiStagedDetails::VISITED;      
      ++ pCurrentSamplePoint;
      nElementsLeft --;
    }

    int nVisited = 0;
    int nOutside = 0;
    // skip over all visited voxels and unrestricted voxels
    while( pCurrentSamplePoint != SamplePointPtrList.end()
           && ( (* pCurrentSamplePoint )->nID != MultiStagedDetails::NOT_VISITED
                || ! this->IsInRestrictedRegion( *pCurrentSamplePoint,
                                           oSampleCenter, fSampleRadius ) ) )
    {
      nElementsLeft --;
      ++ pCurrentSamplePoint;
    }
    
    if( pCurrentSamplePoint == SamplePointPtrList.end() )
    {
      if( bFirstPass )
      {
        nElementsLeft = SamplePointPtrList.size();
        bFirstPass = false;
        pCurrentSamplePoint = SamplePointPtrList.begin();     // one time reset
        while( pCurrentSamplePoint != SamplePointPtrList.end()
               && ( (* pCurrentSamplePoint )->nID != MultiStagedDetails::NOT_VISITED
                    || ! this->IsInRestrictedRegion( *pCurrentSamplePoint,
                                               oSampleCenter, fSampleRadius ) ) )
        {
	  nElementsLeft --;
          ++ pCurrentSamplePoint;
        }
      }
      else
      {
        nElementsLeft = 0;
      }
    }
  }

  //------------------------------------
  //  RefreshSampleCenter()
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void RestrictedStratifiedGrid<SamplePointT, SamplePointGrid>::RefreshSampleCenter()
  {
    SVector3 oNewCenter( 0, 0, 0 );
    Float fWeight = 0;
    if( SolutionPointPtrList.size() > 0 )
    {
      Int nCounted = 0;
      for( Size_Type i = 0; i < SolutionPointPtrList.size(); ++ i)
      {
        Float fConf = SolutionPointPtrList[i]->fConfidence;
        if( fConf > 0 )
        {
          nCounted ++;
          fWeight += fConf;
          oNewCenter += ( SolutionPointPtrList[i]->GetCenter() * fConf );
        }
      }
      if( nCounted > 0)
        oNewCenter = oNewCenter / fWeight;

      Float fShiftDist = ( oSampleCenter - oNewCenter ).GetLength();
      if( fShiftDist > ( fSampleRadius / Float( 50 ) ) )
        bFirstPass = true;
      oSampleCenter = oNewCenter;

    }
  }
  
  //------------------------------------
  //  Solution
  //   Get current result
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  typename RestrictedStratifiedGrid<SamplePointT, SamplePointGrid>::Mic
  RestrictedStratifiedGrid<SamplePointT, SamplePointGrid>::Solution() const
  {
    Mic PartialSolution;
    PartialSolution.InitializeSampleLimits( & ReconstructionRegion );
    for( Int i = 0; i < SolutionGrid.Size1(); i ++ )
      for( Int j = 0; j < SolutionGrid.Size2(); j ++ )
        if( SamplePointGrid::IsValid( SolutionGrid( i, j ) ) )
          if( SolutionGrid(i, j)->nID != MultiStagedDetails::NOT_VISITED )
            PartialSolution.AddVoxel( *SolutionGrid(i, j) );
    
    return PartialSolution;
  }

  
  //---------------------------------------------------------------
  //
  //  BreadthFirstStrategy
  //
  //---------------------------------------------------------------
  
  //---------------------------------------------------------------
  //  Initialize
  //---------------------------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void BreadthFirstStrategy<SamplePointT, SamplePointGrid>
  ::Initialize( const MicFile<SamplePointT> & RecRegion_,
                Float fMinSideLength_ )
  {
    std::cout << "[LazyBFS Client]: Initialize begin " << std::endl;
    ReconstructionRegion = RecRegion_;
    std::cout << "[LazyBFS Client]: Set min resolution " << std::endl;
    Utilities::SetMinimumResolution( ReconstructionRegion, fMinSideLength_ );
    std::cout << "[LazyBFS Client]: initialize with reference " << std::endl;
    SolutionGrid.InitializeWithReference( ReconstructionRegion );
    std::cout << "[LazyBFS Client]: InitializeMicGridID " << std::endl;
    this->InitializeMicGridID( SolutionGrid, MultiStagedDetails::NOT_VISITED );
  }

  //---------------------------------------------------------------
  //  Push  
  //---------------------------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void BreadthFirstStrategy<SamplePointT, SamplePointGrid>
  ::Push( const SamplePointT & oResult )
  {
    SamplePointPtr pCenter = SolutionGrid( oResult );    
    RUNTIME_ASSERT( SamplePointGrid::IsValid( pCenter ),
                    "Unexpected result from outside the solution region[ invalid pointer ]\n" );
    
    if( pCenter->nID == MultiStagedDetails::VISITED
        || pCenter->nID == MultiStagedDetails::NOT_VISITED )  // when there's no fit status
    {
      *pCenter = oResult;
      oSolution.push_back( pCenter );
    }
    else if( pCenter->fConfidence <= oResult.fConfidence )
    {
      *pCenter = oResult;
    }
  }

  //---------------------------------------------------------------
  //  Pop
  //---------------------------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void BreadthFirstStrategy<SamplePointT, SamplePointGrid>::Pop( )
  {
   
    // remove first
    if( SamplePointPtrQueue.front()->nID == MultiStagedDetails::NOT_VISITED )
      SamplePointPtrQueue.front()->nID = MultiStagedDetails::VISITED;
    SamplePointPtrQueue.pop();

    //-----------------
    //  New - remove voxels that are visited and well fitted
    //  CHECK!
    //-----------------
    while( SamplePointPtrQueue.size() > 0 &&
           SamplePointPtrQueue.front()->nID == MultiStagedDetails::FITTED )
      SamplePointPtrQueue.pop();
  }

  //---------------------------------------------------------------
  //  Reset
  //
  //  Solution points are unique sample points that have been visited.
  //  They have to be reset to NON_VISITED.
  //---------------------------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void BreadthFirstStrategy<SamplePointT, SamplePointGrid>::Reset( )
  {
    while( ! SamplePointPtrQueue.empty() )
      SamplePointPtrQueue.pop();    
    oSolution.clear();
  }

  //---------------------------------------------------------------
  //  InsertSeed
  //
  //  This is the "breadth first" step
  //---------------------------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  void BreadthFirstStrategy<SamplePointT, SamplePointGrid>::InsertSeed( const SamplePointT & oSeedPoint )
  {    
    SamplePointPtr pCenter = SolutionGrid( oSeedPoint );
    
    if( pCenter->nID == MultiStagedDetails::NOT_VISITED )
    {
      pCenter->nID = MultiStagedDetails::VISITED;
      *pCenter = oSeedPoint;
      SamplePointPtrQueue.push( pCenter );
    }
    
    vector<SamplePointPtr> oNgbs = SolutionGrid.GetNeighbors( *pCenter );
    for( Size_Type i = 0; i < oNgbs.size(); ++ i )
    {
      if( oNgbs[i]->nID == MultiStagedDetails::NOT_VISITED )
      {
        oNgbs[i]->nID           = MultiStagedDetails::VISITED;
        oNgbs[i]->oOrientMatrix = oSeedPoint.oOrientMatrix;
        SamplePointPtrQueue.push( oNgbs[i] );
      }
    }
  }
}
