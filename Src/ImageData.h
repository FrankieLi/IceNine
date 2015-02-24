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
//
//   ImageData.h
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:  This is the simpliest representation of detector
//             data.  We are simply representing the CCD output
//             as a TIFF image, where location (j, k) may be
//             accessed directly.  Note that the choice to not use Raster
//             directly, but instead to write a wrapper around that.  This
//             is to separate the search function of ImageData away from
//             rasterization and other lower level access.  (Note that rasterization
//             in itself is non-trivial.
//
//   Note:    J is the horizontal direction, K is the vertical direction.
//            This convention came form the 3dXDM coordinate systems.
//
////////////////////////////////////////////////////////////


#ifndef _IMAGE_DATA_H_
#define _IMAGE_DATA_H_

#define BOOST_NO_HASH   // for boost graph

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp>

#include <boost/shared_ptr.hpp>

#include "DetectorData.h"
#include "XDMRaster.h"
#include "Parser.h"
#include "InitFilesIO.h"
#include "Peak.h"
#include "Pixel.h"
#include "Quadtree.h"
#include "Voxel.h"
#include "Serializer.h"
#include "PeakFile.h"

using namespace GeneralLib;
using std::string;
using std::ofstream;
using namespace PBRMath;
//-------------------------------------------------------------------------------------
//
//  Class CImageData
//
//-------------------------------------------------------------------------------------
class CImageData : public CDetectorData
{
  
protected:
  
  CXDMRaster oImage;

  
  struct SPixelCoordCmp
  {
    Int nWidth;
    Int nHeight;
    SPixelCoordCmp( Int w_, Int h_ )
      : nWidth( w_ ), nHeight( h_ ) {}
    
    bool operator() ( const Pixel & oLHS,
                      const Pixel & oRHS ) const
    {
      return ( oLHS.x * nWidth + oLHS.y )  < (oRHS.x * nWidth + oRHS.y );
    }
  };
  
  //-------------------------------------------------
  //  PixelToIndex
  //-------------------------------------------------
  Size_Type PixelToIndex( const Pixel & p ) const
  {
    return p.x * oImage.GetWidth() + p.y;
  }

  
public:

  //
  //  C O N S T R U C T O R S
  //
  CImageData();
  ~CImageData();
  CImageData( Int nWidth, Int nHeight );

  CImageData( const CImageData & oImage );

  CImageData & operator=( const CImageData & oRHS );

  //--------------------------
  //  Resize();
  //--------------------------
  void Resize( Int nWidth, Int nHeight );
  
  //--------------------------
  //  SetPixel();
  //--------------------------
  void SetPixel( Int nJPixel, Int nKPixel, Float fValue );
  
  //--------------------------
  //  AddToImage:  Add pixel to image
  //
  //  Note that this *adds* the value to an image location
  //  
  //--------------------------
  void AddToPixel( Int nJPixel, Int nKPixel, Float fValue );

  //--------------------------
  //   Add Triangle
  //--------------------------
  void AddTriangle( const Point &v0, const Point & v1, const Point & v2, Float fIntensity );

  //--------------------------
  //  Add Polygon
  //       Vertices must be in winding order. (i.e, right hand rule).
  //--------------------------
  void AddPolygon( const vector<Point> & oVertices, Float fFillValue );

  //--------------------------
  //   GetTriangleOverlap
  //
  //   Purpose:  Calculate the number of pixels overlapping
  //             between the triangle specified by the vertices
  //             and the image.
  //--------------------------
  Int GetTriangleOverlap( const Point &v0, const Point & v1, const Point & v2 );
  
  //--------------------------
  //  GetTriangleOverlapProperty
  //--------------------------
  std::pair<Int, Int>
  GetTriangleOverlapProperty( const Point &v0, const Point & v1, const Point & v2 );


  std::pair<Int, Int>
  GetSquareOverlapProperty( const Point &v0, const Point & v1,
                            const Point & v2, const Point & v3 );
  

  //--------------------------
  //  GetNumPixelsLit
  //
  //  Purpose:  count the number of pixels that the triangle specified by
  //            the three verticies lit on the detector.
  //--------------------------
  Int GetNumPixelsLit( const Point &v0, const Point & v1, const Point & v2 );
  
  //--------------------------
  //   ClearImage
  //--------------------------
  void ClearImage();


  //
  //  F I L E   I / O
  //

  //--------------------------
  //  ReadCXDMSimulationDataFile
  //--------------------------
  bool ReadCXDMSimulationDataFile  ( const string &filename );
    
  //--------------------------
  //   WritePixelsToASCII
  //--------------------------
  bool WritePixelsToASCII( const string & filename );

  //--------------------------
  //   WritePixelsToBin
  //--------------------------
  bool WritePixelsToBin( const string & filename );
  
  
  //
  //  A C C E S S O R S
  //
  
  Bool IsDark( Int nJPixel, Int nKPixel ) const;
  Bool IsBright( Int nJPixel, Int nKPixel ) const;
  Bool IsInBound( Int nJPixel, Int nKPixel) const;

  BBox2D GetImageBBox() const; 

  IntensityT At( Int nJPixel, Int nKPixel ) const
  { return oImage( nJPixel, nKPixel ); }
  
  //--------------------------
  // Save and Restore
  //--------------------------
  Bool Save   ( CSerializer & oSerialBuf ) const;
  Bool Restore( CDeserializer & oSerialBuf );
  
};


//-------------------------------------------------------------------------------------
//
//  Specialized 2D geometry namespace for simple things like triangle bounding box intersection,
//  separation axis test  (TODO:  move to another file in the future?) 
//
//-------------------------------------------------------------------------------------
namespace Geometry2D
{

  //-------------------------------------------------------------------------------------
  //
  //  CheckSeparatingAxis
  //
  //  Verticse of a rectangle is given as four vertices. 
  //
  //  Given separating axis specified by vertices v1 and v2 containing only the X and Y
  //  coordinates,  (Note that this is 2D)
  //
  //-------------------------------------------------------------------------------------
  bool CheckSeparatingAxis( const BBox2D & oBox, const SVoxel & oVoxel,
                            const SVector3 &v1, const SVector3 &v2 );
  
  //-------------------------------------------------------------------------------------
  //
  //  VoxelAABBOverlap()
  //
  //  Return:    True     if the AABB (axis aligned bounding box) overlaps with the voxel
  //             False    otherwise
  //
  //  NOTE:  This function only works in 2D.  Specifically, it only works with the x and y
  //         coordinates.
  //
  //-------------------------------------------------------------------------------------
  bool VoxelAABBOverlap( const BBox2D &oBox, const SVoxel &oVoxel );
  bool TriangleAABBOverlap( const BBox2D &oBox, const Point & p1, const Point & p2, const Point & p3);
  bool TriangleAABBOverlap( const BBox2D &oBox, const SVector3 &p1, const SVector3 &p2, const SVector3 &p3 );
  bool TriangleAABBOverlapEstimate( const BBox2D &oBox, const Point & p1, const Point & p2, const Point & p3);
}

namespace RangeSearch
{

  typedef boost::shared_ptr<CDetectorPeak> CPeakPtr;
  
  //-------------------------------------------------------------------------------------
  //  Class CPeakQuadtree -- a specialized class that has an early exit for
  //                         range search
  //-------------------------------------------------------------------------------------
  class CPeakQuadtree : public CQuadtree< CPeakPtr >
  {
  public:
    typedef CQuadNode<CPeakPtr> PeakTreeNodeT;
    
    //---------------------
    //  SRectOverlapProc
    //  Purpose:  function to return
    //  overlap with rectangular range
    //---------------------
    struct SAllRectOverlap
    {
      std::map< CPeakPtr, CPeak> oPeakMap;

      //  Default constructor
      SAllRectOverlap():oPeakMap(){}

      // Function operator
      inline bool operator() (  const PBRMath::BBox2D oRange, const CPeakPtr & pPeak) 
      {
        if( oRange.Overlaps( pPeak->GetBoundingBox() ) )
        {
          oPeakMap[ pPeak ] =  *pPeak;
          return true;
        }
        return false;
      }
    };

    //---------------------
    // SAnyRectOverlap
    // A short circuit is performed when
    // any overlap is found
    //---------------------
    struct SAnyRectOverlap
    {
      inline bool operator() (  const PBRMath::BBox2D oRange, const CPeakPtr & pPeak) const
      {
        if( oRange.Overlaps( pPeak->GetBoundingBox() ) )
          return true;
        return false;
      }
    };
    
  private:

    //---------------------
    //  Return true of there is any
    //  overlap between the the specified
    //  range and list of peaks on the detector.
    //  Note that the search will quit early
    //  if any is found.
    //
    //  NOTE:  Code taken directly from range search (maybe a merge in the
    //         not so distant future?)
    //---------------------
    template< class LookupProc >
    Bool FindOverlap( const PeakTreeNodeT   *node, 
                      const PBRMath::BBox2D &nodeBound,
                      const PBRMath::BBox2D &oSearchRange, 
                      LookupProc &process ) const
    {
      for( Size_Type i = 0; i < node->vData.size(); ++i )
      {
        Bool bFound =  process( oSearchRange, node->vData[i] );
        if( bFound )
	  {
	    std::cerr <<  i << std::endl;
          return true;
	  }
      }
      
      Bool overlap[4];
      PBRMath::Point pMid = 0.5 * nodeBound.pMin  +  nodeBound.pMax * 0.5;
      // node 0 and 2 are the ones with the lower x values => x in [x_min, x_mid]
      // node 1 and 3 are the ones with the upper x values =? x in [x_mid, x_max]
      overlap[0] = overlap[2] = oSearchRange.pMin.x < pMid.x;
      overlap[1] = overlap[3] = oSearchRange.pMax.x >= pMid.x;
	
      // y lower
      overlap[0] &= ( oSearchRange.pMin.y < pMid.y );  
      overlap[1] &= ( oSearchRange.pMin.y < pMid.y );

      // y upper
      overlap[2] &= ( oSearchRange.pMax.y >= pMid.y );
      overlap[3] &= ( oSearchRange.pMax.y >= pMid.y );
  
      for( Int nChild = 0; nChild < 4; ++ nChild ) 
      {
        if ( overlap[ nChild ] && node->opChildren[ nChild ] )
        {
          PBRMath::BBox2D childBound;
      
          // x direction  ture => upper, false => lower
          childBound.pMin.x = ( 1 & nChild ) ? pMid.x : nodeBound.pMin.x;
          childBound.pMax.x = ( 1 & nChild ) ? nodeBound.pMax.x : pMid.x;
      
          // y direction: true => lower, false => upper
          childBound.pMin.y = ( 2 & nChild ) ? pMid.y : nodeBound.pMin.y;
          childBound.pMax.y = ( 2 & nChild ) ? nodeBound.pMax.y : pMid.y;
          
          Bool bFound = FindOverlap( node->opChildren[ nChild ], childBound,
                                     oSearchRange, process );
          if( bFound )
            return true;
        }
      }
      return false;
    }

  public:

    //--------------------------
    //  NumPeakOverlap
    //--------------------------
    Size_Type NumPeakOverlaps( const PBRMath::BBox2D & oSearchRange ) const
    {
      SAllRectOverlap FProc;
      FindOverlap( &oRoot, oTreeBound, oSearchRange, FProc );
      return FProc.oPeakMap.size();
    }
    
    //--------------------------
    // Overlaps
    // Purpose:  Return true if there is any overlap
    //           at all.  False otherwise
    //--------------------------
    Bool Overlaps( const PBRMath::BBox2D & oSearchRange ) const
    {
      SAnyRectOverlap FProc;
      return FindOverlap( &oRoot, oTreeBound, oSearchRange, FProc );
    }

    //--------------------------
    //  ProcessOverlaps
    //--------------------------
    template< class ProcessFn >
    void ProcessOverlaps( ProcessFn & FProc, const PBRMath::BBox2D & oSearchRange ) const
    {
      FindOverlap( &oRoot, oTreeBound, oSearchRange, FProc );
    }
    
  };
  
}

//-------------------------------------------------------------------------------------
//
//  Class CSearchableImageData
//
//  Purpose:  Includes a bounding box that allows quick intersection searches
//
//-------------------------------------------------------------------------------------
class CSearchableImageData : public CImageData
{

private:
  

  
  // For testing purpose before installing kdtree   // moved to public to testing purposes
  // TODO:  Remove this -- but someone's has to write a
  //        serializer for the tree...
  vector<CDetectorPeak> vDetectorPeakList;

  //------------------------------
  //  need an int-> pointer map.
  //------------------------------
  typedef boost::shared_ptr<BBox2D> BBoxPtr;
  std::map< Int, BBoxPtr > IndexToPeakMap;

  RangeSearch::CPeakQuadtree oPeakTree;

  //--------------------------
  // Null_Deleter -- doesn't do anything
  //--------------------------
  struct Null_Deleter
  {
    void operator()(void const *) const { }
  };
  
  //--------------------------
  // BuildSearhTree
  //--------------------------
  void BuildSearchTree( );

  //--------------------------
  // BuildIndexToPeakPtrMap()
  //--------------------------
  void BuildIndexToPeakPtrMap( );
  
public:

  //--------------------------
  // Default constructor
  //--------------------------
  CSearchableImageData(): CImageData(), vDetectorPeakList(), oPeakTree() {};

  //--------------------------
  // Override the parent function to Extract Peaks
  //--------------------------
  bool ReadCXDMSimulationDataFile( const string &filename );

  //--------------------------
  //  ReadCXDMSimulationUFFFile
  //  Note:  When bServerMode is switched to true, the actual raster
  //         will NOT be built.  This is designed for memory saving in
  //         the case of server/client apps.  
  //--------------------------
  bool ReadCXDMSimulationUFFFile( const string &filename,
                                  bool bServerMode = false );

  //--------------------------
  //  WriteToUFF format
  //--------------------------
  bool WritePixelsToUFF( const string & filename );
  
  //--------------------------
  //   CheckOverlap  --
  //
  //   TODO: evolve this function into returning area of overlap?
  //
  //  Return true if there is any overlap
  //--------------------------
  Bool HasOverlap( const Point &v0, const Point & v1, const Point & v2 ) const;
  Bool HasOverlap( const vector<Point> & vPoints ) const;

  //--------------------------
  //
  //  Return the number of peaks overlapping
  //
  //--------------------------
  Int  GetNumOverlapPeaks( const Point &v0, const Point & v1, const Point & v2 ) const;
  Int  GetNumOverlapPeaks( const vector<Point> & vPoints ) const;

  //--------------------------
  //   GetNumPeakOverlap  --
  //
  //   Return the number of peaks (that is distinct intersections of bounding box)
  //
  //--------------------------
  Int GetNumPeaksOverlap( const Point &v0, const Point & v1, const Point & v2 ) const;

  //--------------------------
  //  HasOverlapPeak
  //
  //  Return true if there exist one peak that overlaps the bounding
  //  box oBox
  //--------------------------
  Bool HasOverlapPeak( const BBox2D & oBox ) const;

  //--------------------------
  // ProcessOverlappingBBox
  //--------------------------
  template< class FProcess >
  void ProcessOverlaps( FProcess & Processor, const BBox2D & oSearchBox ) const
  {
    oPeakTree.ProcessOverlaps( Processor, oSearchBox );
  }
  
  //--------------------------
  // Accessor
  //--------------------------
  void InitializePeakList();

  //--------------------------
  //
  //  GetNumPeaks
  //
  //--------------------------
  inline Size_Type GetNumPeaks() const
  {
    return vDetectorPeakList.size();
  }

  //--------------------------
  // GetPeakList
  //--------------------------
  inline const vector<CDetectorPeak> & GetPeakList() const
  {
    return vDetectorPeakList;
  }
  
  //--------------------------
  // Save and Restore
  //--------------------------
  Bool Save   ( CSerializer & oSerialBuf ) const ;
  Bool Restore( CDeserializer & oSerialBuf,
                bool bBuildSearchTree = true);
//   Bool Restore( CDeserializer & oSerialBuf )
//   { return Restore( oSerialBuf, true ); }

  template< typename CustomRasterType >
  Bool Overwrite( const CustomRasterType & Image )
  {
    return oImage.Overwrite( Image );
  }
  
};
#endif
