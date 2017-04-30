/////////////////////////////////////////////////////////////////
//
//  File:    GrainAnalysis.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Implementation of grain analysis routines
//
/////////////////////////////////////////////////////////////////


#include "GrainAnalysis.h"


namespace GrainAnalysis
{
  //--------------------------------------------------------------------------------------------------------
  //  PointToPointMisorientation
  //
  //
  // find this point in map2
  // check for cases:  [triangle misorientation]
  //    Overlap - easy
  //    map1.triangle inside map2.tirangle - easy
  //    map1.triangle contains map2.triangle - hard    (regrid)
  //    map1.triangle intersects map2.triangle - hard  (regrid till contain)
  //
  // calculate misorientation
  // insert triangle into oMisorientationMap
  //
  //
  //  TODO:  Change Find -> Overlap, which is more robust.  Right now overlapping voxels will not be compared
  //         if one voxel doesn't lie inside another.
  //
  //--------------------------------------------------------------------------------------------------------
  CVoxelTree PointToPointMisorientation( const CVoxelTree & oMap1,  CVoxelTree oMap2,
                                         const CRefineCriterion &oStoppingCriterion )
  {
    CDiffOrientation<LatticeSymmetry::CCubicSymmetry> oDiffFun;
    return PointToPointCalculate( oMap1, oMap2, oStoppingCriterion, oDiffFun );
  }


  CVoxelTree PointToPointMisorientationHex( const CVoxelTree & oMap1, CVoxelTree oMap2, const CRefineCriterion &oStoppingCriterion )
  {
    CDiffOrientation<LatticeSymmetry::CHexagonalSymmetry> oDiffFun;
    return PointToPointCalculate( oMap1, oMap2, oStoppingCriterion, oDiffFun );
  }
  
  
  //--------------------------------------------------------------------------------------------------------
  //  PointToPointConfidence
  //
  //--------------------------------------------------------------------------------------------------------
  CVoxelTree PointToPointConfidence( const CVoxelTree & oMap1, CVoxelTree oMap2, const CRefineCriterion &oStoppingCriterion )
  {
    CDiffConfidence oDiffFun;
    return PointToPointCalculate( oMap1, oMap2, oStoppingCriterion, oDiffFun );
  }
 
  

  
}
