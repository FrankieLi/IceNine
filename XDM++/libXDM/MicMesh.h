/////////////////////////////////////////////////////////////////
//
//  File:    MicMesh.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Purpose: Special purpose triangular mesh class for MIC file
//           data analysis.  Note that we only support MIC file's
//           equilateral triangles.
//
//
//
//
//
/////////////////////////////////////////////////////////////////

#ifndef _MIC_MESH_H
#define _MIC_MESH_H

#include "BBox.h"
#include "3dMath.h"

#include "Voxel.h"
#include "MicIO.h"
#include "MathDebug.h"

#include "Types.h"
#include "Geometry.h"
#include "SimpleGraph.h"

#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "UniqueVertexMap.h"

#include "UltraLightWeightGraph.h"


#include <limits>
#include <vector>
#include <map>
#include <iostream>

#define BOOST_NO_HASH   // for boost graph
// Boost Graphics Library
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace XDMUtility
{
  
  using std::vector;
  using namespace GeneralLib;

  //---------------------------------------------------
  //  The purpose of this is to allow
  //---------------------------------------------------
  template < typename DataT >
  struct SDataPoint
  {
    DataT oData;
    Int nID;
    SDataPoint(): oData(){}
    
    SDataPoint( const DataT & oData_ ):
      oData( oData_ ){}

    SDataPoint( const SDataPoint<DataT> & oRHS ) :
      oData( oRHS.oData ){}
    
    
    //  Make this a template parameter 
    bool EquivilentCmp( const DataT & oRHS )
    {
      return oData.PerfectMatch( oRHS );
    }
  };
  
  //-----------------------------------------------------
  //  Null_Deleter   ( for Boost::shared_ptr )
  //-----------------------------------------------------
  struct Null_Deleter
  {
    void operator()(void const *) const { }
  };

  //-----------------------------------
  //  FNumericsError
  //  -- Definition of error provided by the user
  //
  //  Clearly, this could get a whole lot more complicated.
  //-----------------------------------
  struct FNumericsError
  {
    Float fError;

    FNumericsError( ): fError( std::numeric_limits<Float>::epsilon() ) {}
    FNumericsError( const FNumericsError & oFErr ):fError( oFErr.fError ){}
    FNumericsError( Float fError_ ):fError( fError_ ){}
    
    Float operator()() const
    {
      return fError;
    }
  };
 
  //-----------------------------------
  //  FNgbCollector -- collect all neighbors
  //  in the shape tree
  //-----------------------------------
  template< typename ObjT, typename BBoxT >
  struct FNgbCollector
  {
    vector<ObjT> oNgbList;
    void operator()( const BBoxT &b, const ObjT & o )
    {
      oNgbList.push_back( o );
    }
  };

  template< typename ObjT >
  struct FDuplicateCheck 
  {
    ObjT pKnown;
    FDuplicateCheck( const ObjT & p_ ): pKnown( p_ ){}
    bool operator() ( const ObjT & o ) const
    {
      return pKnown->nID == o->nID;
    }
  };

  template< typename ObjT >
  struct FShapePtrLessCmp
  {
    bool operator() ( const ObjT & oLHS, const ObjT & oRHS ) const
    {
      return oLHS->nID < oRHS->nID;
    }
  };
  
  template< typename ObjT >
  struct FShapePtrEqCmp
  {
    bool operator() ( const ObjT & oLHS, const ObjT & oRHS ) const
    {
      return oLHS->nID == oRHS->nID;
    }
  };

  template< typename ObjPtr >
  struct FDisconnectedCmp
  {
    const SVoxel & oCmp;
    FDisconnectedCmp( const SVoxel & o_ ): oCmp( o_ ){}
    bool operator() ( const ObjPtr & oRHS ) const
    {
      return ! oCmp.Overlaps( *oRHS );
    }
  };
  
  //-----------------------------------
  //  CMicMesh -- To be templated
  //
  //  Specialized mesh for MIC files
  //
  //  TODO:  change to MicMesh< ShapeT, ErrorDefinition >
  //
  //  Requirement on ShapeT
  //
  //  Must implement the following functions
  //  BBox2D    GetBoundingBox()
  //  bool      Overlap( ShapeT )   -- Definition of overlap
  //  VertexT   Vertex( Int i )     -- VertexT must be a comparable object
  //  bool      Connected( ShapeT ) -- Definition of connectedness
  //     bool      operator< ()
  //     OR   a comperator function
  //  vector< ShapeT >  Refine()
  //
  //-----------------------------------
  class CMicMesh
  {
    
  public:

    //-----------------------------------
    // The name "shape" comes from shape function from
    // FEM community.
    //-----------------------------------
 
    typedef SVoxel                              ShapeT; // to be removed to make MicMesh<ShapeT>
    typedef boost::shared_ptr< ShapeT >         ShapePtr;
    typedef FShapePtrLessCmp< ShapePtr >        ShapeWeakOrderCmp;
    
    typedef RangeSearch::CQuadtree< ShapePtr >          ShapeTreeT;
    typedef std::vector< ShapePtr >                     ShapePtrList;
    typedef ShapePtrList::iterator                      ShapePtrIter;
    typedef ShapePtrList::const_iterator                ConstShapePtrIter;
    
    typedef PBRMath::BBox2D                             BBox;
    typedef PBRMath::Point                              Point;
    typedef FNgbCollector<ShapePtr, BBox>               FNgbLookup;

  private:

    //  could be generalized to take equal operator
    struct SDeleteCmp
    {
      const ShapePtr & oToDelete;

      SDeleteCmp( const ShapePtr & oToDelete_ )
        :oToDelete( oToDelete_ ){}
      SDeleteCmp( const ShapeT & oToDelete_ )
        :oToDelete( ShapePtr( const_cast< ShapeT* >( &oToDelete_ ), Null_Deleter() ) ) {}

      bool operator() ( const ShapePtr & oRHS ) const
      {
        return ( oRHS->PerfectMatch( *oToDelete ) );
      }
    };
    
    ShapeTreeT      oShapeLocator;
    ShapePtrList    oShapePtrList;
    FNumericsError  FErrorEstimator;

    //-----------------------------------
    // UpdateErrorLimit
    //  Purpose:  Reset the limit of the error
    //            based on the new error.
    //            If the new error is less than the
    //            preset error estimate, a new error
    //            will be calculated.
    //-----------------------------------
    void UpdateErrorLimit( const BBox oBoundingBox );
    
  public:
    
    static const Int DefaultMaxDepth = 10;
    //-----------------------------------
    //  CMicMesh
    //-----------------------------------
    CMicMesh()
    {
      oShapeLocator.SetMaxDepth( DefaultMaxDepth );
      FErrorEstimator.fError = std::numeric_limits<Float>::epsilon();
    }

    CMicMesh( const BBox & oBox, Int nMaxDepth, Float fError)
      : oShapeLocator( oBox, nMaxDepth )
    {
      FErrorEstimator.fError = fError;
    }
    
    //-----------------------------------
    //  Insert
    //-----------------------------------
    ShapePtr Insert( const ShapeT & v );
    ShapePtr Insert( const ShapePtr & p );

    //-----------------------------------
    //  Delete
    //-----------------------------------
    void Delete( const ShapeT & v );
    
    //-----------------------------------
    //  Refine
    //-----------------------------------
    void Refine( const ShapePtr & pVoxel );
    
    //-----------------------------------
    //  Coarsen
    //-----------------------------------
    void Coarsen( const ShapePtr & pVoxel );
    
    ////////////////////////////////////////////////////
    //   A C C E S S O R S
    ////////////////////////////////////////////////////

    //-----------------------------------
    //  Access iterators
    //-----------------------------------
    ConstShapePtrIter Begin() const { return oShapePtrList.begin(); }
    ShapePtrIter      Begin()       { return oShapePtrList.begin(); }

    ConstShapePtrIter End() const   { return oShapePtrList.end(); }
    ShapePtrIter      End()         { return oShapePtrList.end(); }

    //-----------------------------------
    //  Find(&ShapeT)
    //
    //  Purpose:  Find a shape that's exactly the same as oSearchShape,
    //            within tolerence.
    //-----------------------------------
    
    //-----------------------------------
    //  FindOverlap(const &SVoxel)
    //
    //  Purpose:  Find all voxels within the voxel tree that overlaps the voxel, oCenter 
    //            All (Data) voxels in oRes must overlap oCenter.
    //
    //-----------------------------------
    ShapePtrList FindOverlap( const BBox & oCenter ) const;

    //-----------------------------------
    //  FindOverlap
    //  Purpose:  Find all overlapping voxel with oCenter, given a relative
    //            error limit.  (Fractional)   Note that this limits the resolution
    //            to fRelError * fSideLength
    //-----------------------------------
    ShapePtrList FindOverlap( const ShapeT & oCenter, Float fRelError = 0.001 ) const;
    
    //-----------------------------------
    //  GetNeighbors
    //  Precondition:  
    //  Purpose:  Return all voxels that shares an edge with the current voxel
    //
    //  Note:     Use this if you don't have the original voxel
    //
    //-----------------------------------
    ShapePtrList GetNeighbors( const ShapeT & oCenter ) const;
    
    void Find( vector<ShapePtr> & oRes, const BBox & oSearchBox ) const;
    bool Exists( const ShapeT & oCenter ) const;
  };
  

  //-----------------------------------------------------------------------------------
  //
  //   ShapeLocator   -- running range search and point location on shapes in 3D
  //
  //   TODO:  template-ize on shape and comparator
  //
  //-----------------------------------------------------------------------------------
  template< class ShapeContainer >
  class CGeneralShapeLocator
  {
  protected:
    
    struct SHeightCmp
    {
      Float fError;
      SHeightCmp( Float f_ ):fError(f_){}
      
      bool operator()( Float fLHS, Float fRHS ) const
      {
        if( fabs( fLHS - fRHS ) <= fError )
          return false;
        return fLHS < fRHS;
      }
    };
    
  public:

    typedef SVoxel                     UserDataT;
    typedef ShapeContainer             ShapeTreeT;
    typedef typename ShapeTreeT::ShapeT    ShapeT; // to be removed to make MicMesh<ShapeT>
    typedef typename ShapeTreeT::ShapePtr  ShapePtr;
    
    typedef boost::shared_ptr< ShapeTreeT >           ShapeTreePtr;
    typedef typename ShapeTreeT::BBox                 BBox2D;
    typedef GeneralLib::SBoundingBox                  BBox3D;
    typedef FNgbCollector<ShapePtr, BBox2D>           FNgbLookup2D;
    
    typedef std::map< Float, ShapeTreePtr, SHeightCmp > ShapeTree3D;
    typedef std::pair< Float, ShapeTreePtr >            HeightLayerPair;            
    typedef typename ShapeTree3D::iterator              LayerIter;
    typedef typename ShapeTree3D::const_iterator        ConstLayerIter;
    
  protected:

    BBox2D          oPlanarBBox;
    FNumericsError  FErrorEstimator;
    ShapeTree3D     oObjLocator;
    Int             nMaxDepth;
    
    //-----------------------------------
    //  UpdateErrorLimit  -- feels like ShapeLocator and MicMesh
    //                       should inherit from the same interface.
    //-----------------------------------
    void UpdateErrorLimit( const BBox2D oBoundingBox );
    
  public:
    

    static const Int DefaultMaxDepth = 10;
    
    //-----------------------------------
    //  CShapeLocator (default ctor)
    //-----------------------------------
    CGeneralShapeLocator() : oObjLocator(  SHeightCmp( Float( 100 ) * std::numeric_limits<Float>::epsilon() ) ),
                             nMaxDepth( DefaultMaxDepth )
    {
      FErrorEstimator.fError = Float( 100 ) * std::numeric_limits<Float>::epsilon();
    }

    //-----------------------------------
    //  CGeneralShapeLocator
    //-----------------------------------
    explicit CGeneralShapeLocator( const BBox2D & oPlanarBounds,
                                   Float fError, Int nMaxDepth_ ):
      oPlanarBBox( oPlanarBounds ), oObjLocator( fError ),
      nMaxDepth( nMaxDepth_ )
    {
      FErrorEstimator.fError = fError;
    }

    //-----------------------------------
    // Insert 
    //-----------------------------------
    ShapePtr Insert( const ShapeT   & oShape );
    ShapePtr Insert( const ShapePtr & pShape );

    //-----------------------------------
    //  Inserting iterators, [pFirst, pEnd )
    //-----------------------------------
    template< class IteratorT >
    void Insert( IteratorT pFirst, IteratorT pEnd )
    {
      for( IteratorT pCur = pFirst; pCur != pEnd; ++ pCur )
        Insert( *pCur );
    }
    
    vector<ShapePtr> Find( const BBox3D & oBox ) const;
    vector<ShapePtr> Find( const BBox2D & oBox, Float fLayerLoc ) const;

    //-----------------------------------
    //   GetNeighbors
    //    Purpose:  Get neighbors in 3D, i.e., from in layer and neighboring
    //              layers.
    //-----------------------------------
    vector<ShapePtr> GetNeighbors( const ShapePtr & pShape ) const;

    //-----------------------------------
    //     A C C E S S O R S
    //-----------------------------------
    LayerIter LayerBegin() { return oObjLocator.begin(); }
    LayerIter LayerEnd()   { return oObjLocator.end(); }

    ConstLayerIter LayerBegin() const { return oObjLocator.begin(); }
    ConstLayerIter LayerEnd()   const { return oObjLocator.end(); }


    // DEBUG function;
    Size_Type NumLayers() const { return oObjLocator.size(); }
    
  };


  typedef CGeneralShapeLocator<CMicMesh> CShapeLocator;
  
}  // end XDMUtilty

#include "MicMesh.tmpl.cpp"

#endif
