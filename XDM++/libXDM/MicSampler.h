/////////////////////////////////////////////////////////////////
//
//
//  MicSampler
//
//  Purpose:  Simple sampler for the MIC formats
//
//
//
/////////////////////////////////////////////////////////////////


#ifndef _MIC_SAMPLER_H
#define _MIC_SAMPLER_H


#include "Quadtree.h"
#include <boost/shared_ptr.hpp>
#include "BBox.h"
#include <vector>
#include "Sampling.h"

namespace Sampler2D
{
  using std::vector;
  template< typename VoxelT >
  struct FVoxelCollector
  {
    vector<VoxelT> oVoxelList;
    void operator() ( const PBRMath::Point & p, const VoxelT &oData )
    {
      oVoxelList.push_back( oData );
    }
  };
  
  //------------------------------------------------
  //
  //  Template Class
  //  CVoxelSampler
  //
  //  Purpose:  A simple sampler for generalized voxel
  //------------------------------------------------
  template< typename VoxelT, typename DataT >
  class CVoxelSampler
  {

  private:
    typedef typename RangeSearch::CQuadtree< VoxelT > VoxelTreeT;
    VoxelTreeT oVoxelTree;
    GeneralLib::CUniformRandomReal oRandGen;
  public:

    //------------------------------------------------
    //  Default Constructor
    //------------------------------------------------
    CVoxelSampler(): oVoxelTree(), oRandGen() {}

    //------------------------------------------------
    //  Initialize
    //------------------------------------------------
    void Initialize( Int nMaxTreeDepth, const PBRMath::BBox2D & oBBox )
    {
      oVoxelTree.SetMaxDepth( nMaxTreeDepth );
      oVoxelTree.TreeBound  ( oBBox         );
    }
    
    //------------------------------------------------
    //  AddVoxel
    //------------------------------------------------
    void AddVoxel( const VoxelT & oNewVoxel )
    {
      PBRMath::BBox2D oBox = oNewVoxel.GetBoundingBox();
      oVoxelTree.Add( oNewVoxel, oBox );
    }

    //------------------------------------------------
    //  AddVoxels
    //  --  Add a list of voxels to the structure
    //------------------------------------------------
    void AddVoxels( typename vector<VoxelT>::const_iterator & pFirst,
                    typename vector<VoxelT>::const_iterator & pEnd  )
    {
      for( typename vector<VoxelT>::iterator pIter = pFirst; pIter != pEnd; ++pIter )
        AddVoxel( *pIter );
    }

    //------------------------------------------------
    //  UniformRandomSample
    //
    //  Purpose:  Uniformly sample across the bounding box given
    //            in oBBox, which is defined to by the two points
    //
    //            pMin( x0, y0), pMax(x1, y1)
    //
    //            The box is therefore give inclusively as:
    //
    //             [x0, x1] x [y0, y1]
    //
    //  
    //             N sample points are taken randomly from this
    //             region.  The average is returned.  This function
    //             satisfies the requirement that as the area of
    //             oBBox -> 0, the returned value becomes the value
    //             of oMic at this point.
    //
    //  Return value:  A tuple of type (DataT, bool) is returned.
    //                 A point of type DataT, which represents the average
    //                 value is returned as the first value.  A flag is returned
    //                 to indicate if the requisite point is well defined.
    //                 False is returned if ill-defined points are found during
    //                 the sampling.
    //
    //------------------------------------------------
    template< typename AverageFn >
    std::pair< DataT, bool > UniformRandomSample( const PBRMath::BBox2D & oBBox,
                                                  Int nSamplePoints,
                                                  AverageFn oAverager ) const
    {
      bool bException = false;
      vector< DataT > oValueList;
      for( Int i = 0; i < nSamplePoints; i ++ )
      {
        Point oTestPoint;
        oTestPoint.x = oRandGen.GetRandomVariable( oBBox.pMin.x, oBBox.pMax.x );
        oTestPoint.y = oRandGen.GetRandomVariable( oBBox.pMin.y, oBBox.pMax.y );
        FVoxelCollector<VoxelT> oCollector;
        oVoxelTree.LookUp( oTestPoint, oCollector );
        if( oCollector.oVoxelList.size() == 0 )
          bException = true;

        for( Size_Type j = 0; j < oCollector.oVoxelList.size(); j ++ )
          if( oCollector.oVoxelList[j].Intersects( oTestPoint ) )
            oValueList.push_back( oCollector.oVoxelList[j].GetData() );
      }
      DataT oRes = oAverager( oValueList.begin(), oValueList.end() );
      return std::pair< DataT, bool > ( oRes, bException );
    }
  };
  
};

#include "MicSampler.tmpl.cpp"
#endif
