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
//  SearchTraints.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//   Purpose:  Technically speaking, this is not a triat class.
//             The point of this file is to implement the composition
//             rule for typical orientation search algorithms.
//             This composition is true in most cases even
//             if the physics are drmatically different.
//
//   Note:     A word on terminology.  The term "vertex" is defined
//             in sample space, where as "pixel" is defined in the
//             detector space.  A "VertexCostFn" is therefore a function
//             that calculates the cost on a sample space vertex.  In
//             another words, "vertex" is an abstract representation of
//             a sample point on the sample space.  A "voxel" is then
//             an abstraction of a "region" in the sample space.  Taking
//             this to the other direction, a "peak" is a region in the
//             detector space.
////////////////////////////////////////////////////////////

#ifndef _SEARCH_TRAITS_
#define _SEARCH_TRAITS_

#include "CostFunctions.h"
#include "OverlapInfo.h"
#include "Simulation.h"

namespace OrientationSearch
{

  //---------------------------------------------------
  //
  //  GeneralSearchCostFn
  //   - PeakOverlapFn - Calculates the overlap of a "peak"
  //     in the detector space.  This peak could very well be
  //     be sampled by pixels.  In another words, this decides
  //     if a peak is overlapping, and the cost of this overlap.
  //
  //   - DetOverlapFn - Function that sums up all of the overlaps
  //                    within a detector.  This function decides
  //                    which peaks to send to PeakOverlapFn.  Since
  //                    this function calls PeakOverlapFn, the cost
  //                    for a detector could be weighted according to
  //                    whoever wrote DetOverlapFn.
  //
  //
  //   - VertexOverlapFn - As mentioned before, a vertex is a sample
  //                       space object, and this calculates the cost
  //                       of a single "vertex."  This accounts for all
  //                       possible diffractions specified by the writer
  //                       of VertexOverlapFn.
  //
  //   - VoxelCostFn     - A voxel is a region of sample space which
  //                       contains a collection of vertices.  This
  //                       function decides the cost of a voxel.  Similar
  //                       to all the functions above, VoxelCostFn calls
  //                       VertexOverlapFn.
  //
  //
  //  Note:  One can think of this as the definition rule for building an
  //         objective function for the search algorithm.
  //---------------------------------------------------
  template< class PeakFilterT,  class SimulatorT,
            class PeakOverlapFnT,
            template< class T > class DetOverlapFnT,
            template< class U1, class U2, class U3 > class VertexOverlapFnT,
            class VoxelToVertexFnT,
            template< class V1, class V2, class V3 > class VoxelCostFnT,
            class SamplePoint
  >
  class GeneralSearchCostFn
  {
  public:
    typedef PeakOverlapFnT                                                         PeakOverlapFn;
    typedef DetOverlapFnT< PeakOverlapFn >                                         DetOverlapFn;
    typedef VertexOverlapFnT< PeakFilterT, DetOverlapFn, SimulatorT >              VertexOverlapFn;
    typedef VoxelToVertexFnT                                                       VoxelToVertexFn;
    typedef VoxelCostFnT< VertexOverlapFn, VoxelToVertexFn, SamplePoint >          VoxelOverlapFn;
    
  protected:
    const SimulatorT    & _Simulator;   // perhaps check into const correctness

    PeakOverlapFn   PeakCounter;
    DetOverlapFn    DetCounter;
    VertexOverlapFn VertexOverlapCounter;
    VoxelToVertexFn VertexExtractor;
    VoxelOverlapFn  VoxelCostCounter;

    GeneralSearchCostFn();
  public:
    
    GeneralSearchCostFn( PeakFilterT PeakFilter, const SimulatorT & oSim )  
      : _Simulator( oSim ),
        PeakCounter(),
        DetCounter( PeakCounter ),
        VertexOverlapCounter( PeakFilter, DetCounter, _Simulator ),
        VertexExtractor(),
        VoxelCostCounter( VertexOverlapCounter  )
    {}

    //----------------------------------------
    //  Lazy man's way of initialization for PeakCounter - in case
    //  it has parameters, which happens sometimes.
    //----------------------------------------
    template< class ParamT  >
    GeneralSearchCostFn( PeakFilterT PeakFilter, const SimulatorT & oSim, ParamT Param )  
      : _Simulator( oSim ),
        PeakCounter( Param ),
        DetCounter( PeakCounter ),
        VertexOverlapCounter( PeakFilter, DetCounter, _Simulator ),
        VertexExtractor(),
        VoxelCostCounter( VertexOverlapCounter  )
    {}

    // Accessor
    PeakOverlapFn   & PeakCostFn()        { return PeakCounter; }
    DetOverlapFn    & DetectorCostFn()    { return DetCounter;  }
    VertexOverlapFn & VertexCostFn()      { return VertexOverlapCounter; }
    VoxelToVertexFn & VertexExtractorFn() { return VertexExtractor;   }
    VoxelOverlapFn  & VoxelCostFn()       { return VoxelCostCounter; }
    const SimulatorT& Simulator() const   { return _Simulator; }
  }; 

} // end orientation search

#endif
