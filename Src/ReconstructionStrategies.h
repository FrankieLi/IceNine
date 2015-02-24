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
//  ReconstructonStrategies.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//  Purpose:   Implementation of strategies to subdivide and select a voxel
//             out of a layer or volume.  This includes selecting how
//             and where points in the sample space are to be
//             selected and reconstructed.  Some of these methods have a
//             feedback mechanism, where results are fed back in to improve
//             the selection of the next sample point.
//
////////////////////////////////////////////////////////////


#ifndef _RECONSTRUCTION_STRAGETIES_H_
#define _RECONSTRUCTION_STRAGETIES_H_


#include "MicIO.h"
#include "MicGrid.h"
#include "Voxel.h"
#include <queue>

namespace ReconstructionStrategies
{

  namespace Utilities
  {
    template< class GeneralMicFile >
    void SetMinimumResolution( GeneralMicFile & oSampleMic, Float fMaxVoxelSideLength )
    {
      std::cout << " |----Setting min resolution " << fMaxVoxelSideLength << " " << std::endl;
      oSampleMic.SetMinResolution( fMaxVoxelSideLength );
      std::cout << " |---- resulted in " << oSampleMic.GetNumVoxels() << std::endl;
    }

    template<>
    void SetMinimumResolution( CMic & oSampleMic, Float fMaxVoxelSideLength );
    
    
    vector<SVoxel> SubdivideTriangles( const SVoxel & ParentVoxel, int GenMax,
                                       SVector3 PrimaryDiag,
                                       SVector3 SecondaryDiag,
                                       Float fNewSideLength );
    struct Null_Deleter { void operator()( void const*) const{} };
  }

  //------------------------------------
  //  ReconstructionStrategy
  //------------------------------------
  template< class SamplePointT >
  class ReconstructionStrategy
  {

  public:

    //------------------------------------
    //  Reset
    //------------------------------------
    virtual void Reset() = 0;
    
    //------------------------------------
    //  Push - or save a result sample point
    //------------------------------------
    virtual void Push(const SamplePointT & oResult ) = 0;

    //------------------------------------
    // Accessors
    //------------------------------------
    
    //------------------------------------
    //  Pop - Discard the "top most" sample poin
    //------------------------------------
    virtual void Pop() = 0;
    
    //------------------------------------
    //  Current - return current sample point
    //------------------------------------
    virtual SamplePointT &       First() = 0;
    virtual const SamplePointT & First() const = 0;
    
    //------------------------------------
    //  Size - return the current number of sample points
    //------------------------------------
    virtual Size_Type Size() const = 0;
    //------------------------------------
    //  Empty - return true if empty
    //------------------------------------
    virtual bool Empty()     const = 0;
    
    //------------------------------------
    //  Solution
    //   Get current result
    //------------------------------------
    virtual MicFile<SamplePointT> Solution() const = 0;

  };
  //------------------------------------
  //  UniformSelection
  //  -- The simpliest way to select voxel
  //     by uniformly griding the sample space,
  //     then selecting voxels in sequential
  //     order.
  //
  //  The function is similar to a queue.
  //
  //------------------------------------
  class UniformSubidivisonSequentialSelection
  {
  public:
    typedef SVoxel SamplePointT;
    typedef vector<SVoxel>::iterator SamplePointIter;
  private:
    CMic ReconstructionRegion;
    CMic SolutionRegion;
    SamplePointIter pCurrentSamplePoint;
  public:
    
    //------------------------------------
    //  Initialize
    //------------------------------------
    void Initialize( const CMic & RecRegion_,
                     Float fMinSideLength )      
    {
      ReconstructionRegion = RecRegion_;
      Utilities::SetMinimumResolution( ReconstructionRegion, fMinSideLength );
      //SolutionRegion.SetInitialSideLength( ReconstructionRegion.GetInitialSideLength() );
      SolutionRegion.InitializeSampleLimits( & ReconstructionRegion );
      pCurrentSamplePoint = ReconstructionRegion.VoxelListBegin();
    }

    void Reset()
    {
      pCurrentSamplePoint = ReconstructionRegion.VoxelListBegin();
    }
    
    //------------------------------------
    //------------------------------------
    void Push(const SamplePointT & oResult )  // Insert Result
    {
      SolutionRegion.AddVoxel( oResult );
    }
    
    //------------------------------------
    // Accessors
    //------------------------------------

    //------------------------------------
    //  GetVoxelList - Deprecated
    //
    //  This is here to provide backward compatibility
    //  with some of the existing code.
    //------------------------------------
    vector<SVoxel> & GetVoxelList_DEP( )
    {
      return ReconstructionRegion.GetVoxels();
    }
    
    //------------------------------------
    //  Pop - Discard the "top most" sample point
    //------------------------------------
    virtual void Pop()
    {
      if( pCurrentSamplePoint != ReconstructionRegion.VoxelListEnd() )
        ++ pCurrentSamplePoint;
    }
    
    //------------------------------------
    //  Current - return current sample point
    //------------------------------------
    SamplePointT First()
    {
      return *pCurrentSamplePoint;
    }
    
    //------------------------------------
    //  Size - return the current number of sample points
    //------------------------------------
    Size_Type Size() const
    {
      return ReconstructionRegion.VoxelListEnd() - pCurrentSamplePoint;
    }

    //------------------------------------
    //  Empty - return true if empty
    //------------------------------------
    bool Empty() const
    {
      return ( ReconstructionRegion.VoxelListEnd() - pCurrentSamplePoint ) <= 0;
    }

    //------------------------------------
    //  ClearSolution
    //------------------------------------
    void ClearSolution()
    {
      SolutionRegion.GetVoxels().clear();
    }
    
    //------------------------------------
    //  Solution
    //   Get current result
    //------------------------------------
    const CMic & Solution() const
    {
      return SolutionRegion;
    }

    //------------------------------------
    //  Solution
    //   Get current result
    //------------------------------------
    CMic Solution()
    {
      return SolutionRegion;
    }
  };

  namespace MultiStagedDetails
  {
    enum VisitMarkerType
      {
        NOT_VISITED = -1,
        VISITED     = 0,
        FITTED      = 1,
        REFIT       = 2
      };
  }
  
  //------------------------------------
  //  GridStrategyTools -- should be template-template
  //
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  class GridStrategyTools
  {
  public:
    typedef boost::shared_ptr< SamplePointT > SamplePointPtr;
    
  protected:
    

    void InitializeMicGridID( SamplePointGrid & oMicGrid,
                              Int  nIDValue )
    {
      int n = 0;
      std::cout << "----InitializeMicGridID  " << std::endl;
      for( Int i = 0; i < oMicGrid.Size1(); i ++ )
        for( Int j = 0; j < oMicGrid.Size2(); j ++ )
          if( SamplePointGrid::IsValid( oMicGrid( i, j ) ) )
          {
            n ++;
            oMicGrid(i, j)->nID         = nIDValue;
            oMicGrid(i, j)->fConfidence = 0;
          }

      std::cout << "----initialized " << n << " Voxels " << std::endl;
    }

    bool IsInRestrictedRegion( const SamplePointPtr & pPoint,
                               const SVector3 oCenter, Float fRadius )
    {
      Float R = ( pPoint->GetCenter() - oCenter ).GetLength();
      return ( R < fRadius );
    }
  };
  
  //------------------------------------
  //  RestrictedStratifiedGrid
  //  
  //
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid  >
  class RestrictedStratifiedGrid
    : public ReconstructionStrategy< SamplePointT >,
      public GridStrategyTools< SamplePointT, SamplePointGrid >
  {
  public:
    typedef typename GridStrategyTools< SamplePointT, SamplePointGrid >::SamplePointPtr SamplePointPtr;
    typedef MicFile<SamplePointT> Mic;
    typedef boost::shared_ptr< Mic > MicPtr;
    
    Mic         ReconstructionRegion;
    
    vector< SamplePointPtr > SamplePointPtrList;
    vector< SamplePointPtr > SolutionPointPtrList;  
    typedef typename vector<SamplePointPtr>::iterator SamplePointPtrIter;
    SamplePointPtrIter pCurrentSamplePoint;
    
    SamplePointGrid SolutionGrid;
    
    bool     bFirstPass;    // use to figure out if we have already went around the
                            // SamplePointPtrList once already without a
                            // reset of sample radius or center
    Int      nElementsLeft;
    Int      nElementsFitted;
    Int      nElementsVisited;
    
    SVector3 oSampleCenter;
    Float    fSampleRadius;
        
  public:

    RestrictedStratifiedGrid()
      : bFirstPass( false ), nElementsLeft( 0 ),
        nElementsFitted( 0 ), nElementsVisited( 0 ),
        SamplePointPtrList(), SolutionPointPtrList(),
        pCurrentSamplePoint( SamplePointPtrList.begin() )
    {}
    
    virtual void Initialize( const Mic & ReconstructionRegion,
                             Float fMinSideLength,
                             SVector3 oSampleCenter,
                             Float fSampleRadius );

    //-----------------------------
    //  Initialize for Restart -- may need to change
    //  Values of RestartPoint will be assigned to Reconstruction region
    //
    //-----------------------------
    virtual void InitializeRestart( const Mic & ReconstructionRegion,
                                    Mic &  RestartPoint,
                                    Float fMinAcceptThresh,
                                    Float fMinSideLength,
                                    SVector3 oSampleCenter,
                                    Float fSampleRadius );
    
    //-----------------------------
    //  Restricted Region Mutators
    //-----------------------------
    void SetSampleCenter( const SVector3 & c )  { oSampleCenter = c; bFirstPass = true; }
    void SetSampleRadius( Float R )             { fSampleRadius = R; bFirstPass = true; }
    void RefreshSampleCenter();
    
    const SVector3 & SampleCenter() const { return oSampleCenter; }
    const Float & SampleRadius() const { return fSampleRadius; }

    //------------------------------------
    //  Reset
    //------------------------------------
    void Reset();
    
    //------------------------------------
    //  Push - or save a result sample point
    //------------------------------------
    virtual void Push( const SamplePointT & oResult );
    
    //------------------------------------
    // Accessors
    //------------------------------------
    
    //------------------------------------
    //  Pop
    //  - "remove" the top most element.
    //------------------------------------
    virtual void Pop();
    
    SamplePointT       & First()       { return *( * pCurrentSamplePoint ); }
    const SamplePointT & First() const { return *( * pCurrentSamplePoint ); }

    Size_Type Size() const { return nElementsLeft; }
    bool Empty()     const { return pCurrentSamplePoint == SamplePointPtrList.end(); }
    
    //------------------------------------
    //  Solution
    //   Get current result
    //------------------------------------
    Mic Solution() const;
    const vector< SamplePointPtr > & SolutionVector() const { return SolutionPointPtrList; }
  };

  //------------------------------------
  //  BreadthFirstStrategy
  //------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  class BreadthFirstStrategy
    : public ReconstructionStrategy< SamplePointT >,
      public GridStrategyTools< SamplePointT, SamplePointGrid >
  {
  public:
    typedef typename GridStrategyTools< SamplePointT, SamplePointGrid >::SamplePointPtr SamplePointPtr;
   
    
  private:
    MicFile<SamplePointT>            ReconstructionRegion;
    std::queue< SamplePointPtr >     SamplePointPtrQueue;
    vector< SamplePointPtr >         oSolution;
    
    typedef typename vector<SamplePointPtr>::iterator SamplePointPtrIter;
    SamplePointGrid SolutionGrid;
    
    //------------------------------------
    //  Solution
    //   Get current result
    //------------------------------------
    MicFile<SamplePointT> Solution() const
    {
      RUNTIME_ASSERT(0, "BreadthFirstStrategy::Solution() should not be called\n");
      return ReconstructionRegion;
    }   // This version should never be called
    
  public:
    
    void Initialize( const MicFile<SamplePointT> & ReconstructionRegion,
                     Float fMinSideLength );

    //------------------------------------
    //
    //------------------------------------
    void InsertSeed( const SamplePointT & oSeedPoint );
    
    //------------------------------------
    //  Reset
    //------------------------------------
    void Reset();
    
    //------------------------------------
    //  Push - or save a result sample point
    //------------------------------------
    void Push( const SamplePointT & oResult );

    //------------------------------------
    // Accessors
    //------------------------------------

    //------------------------------------
    //  Pop
    //  - "remove" the top most element.
    //------------------------------------
    virtual void Pop();
    
    SamplePointT       & First()       { return * SamplePointPtrQueue.front(); }
    const SamplePointT & First() const { return * SamplePointPtrQueue.front(); }
    
    Size_Type NumFitted() const { return oSolution.size(); }
    Size_Type Size()      const { return SamplePointPtrQueue.size(); }
    bool      Empty()     const { return SamplePointPtrQueue.empty() ; }
    
    const vector< SamplePointPtr > & SolutionVector() const { return oSolution; }
  };


  //--------------------------------------------
  //  RandomizeQueue
  //--------------------------------------------
  template< class SamplePointT, class SamplePointGrid >
  class RandomizedQueue :
    public RestrictedStratifiedGrid< SamplePointT, SamplePointGrid >
  {
    
  public:
    typedef typename GridStrategyTools< SamplePointT, SamplePointGrid >::SamplePointPtr SamplePointPtr;
    typedef MicFile<SamplePointT> Mic;
    typedef RestrictedStratifiedGrid< SamplePointT, SamplePointGrid > Base;

    void RandomizeReset()
    {
      std::random_shuffle( Base::SamplePointPtrList.begin(), Base::SamplePointPtrList.end() );
      Base::pCurrentSamplePoint = Base::SamplePointPtrList.begin();
    }

    template<class OutputIterator >
    void Get( int nVoxels, OutputIterator OutIter )
    {
      while( ! Base::Empty() && nVoxels > 0 )
      {
        OutIter = Base::First();
        nVoxels --;
        Base::Pop();
      }
    }

    void ClearSolution()
    {
      Base::SolutionPointPtrList.clear();
    }

    void RandomizedSetMaxElements( int n )
    {
      n = std::min( n, static_cast<int>( Base::SamplePointPtrList.size() ) );
      RandomizeReset();
      Base::SamplePointPtrList.erase( Base::SamplePointPtrList.begin(),
                                     Base::SamplePointPtrList.begin() + n );
      Base::pCurrentSamplePoint = Base::SamplePointPtrList.begin();
    }
  };
  
}


#include "ReconstructionStrategies.tmpl.cpp"

#endif
