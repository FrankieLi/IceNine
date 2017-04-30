//-----------------------------------------------------------------------
//
//  Filename:  Polygon.h
//
//  Purpose:   General polygon related functions. (i.e., normal calculations
//             and so forth.)
//
//  Author:    Frankie Li
//  email:     sfli@cmu.edu
//
//-----------------------------------------------------------------------


#ifndef _POLYGON_H_
#define _POLYGON_H_

#include "Types.h"
#include "3dMath.h"
#include "BBox.h"

using namespace GeneralLib;

namespace Polygon
{
  static const Size_Type NUM_VERTICES = 3;
  //-------------------------------
  //
  //  Base class of Triangle
  //
  //-------------------------------
  template< typename VertexT >
  class CTriangle
  {
  public:
    static const int NumVertices = NUM_VERTICES;
    VertexT pVertex[NUM_VERTICES];
    
    //------------
    // ctor
    //------------
    CTriangle(){}
    CTriangle( const CTriangle & v_ )
    {
      for( Size_Type i = 0; i < NUM_VERTICES; i ++ )
        pVertex[i] = v_.pVertex[i];
    }
    
  };

  //-----------------------------------
  //  Generalized Shape abstract class, used as
  //  interface.
  //-----------------------------------
  template< class ShapeTrait >
  class CShape
  {
  public:

    typedef typename ShapeTrait::VertexT      VertexT;
    typedef typename ShapeTrait::BBoxT        BBoxT;
    typedef const VertexT                     Const_VertexT;
    typedef typename ShapeTrait::VertexIndexT VertexIndexT;

    template < typename ObjT, typename FNumericLimit >
    bool Connected( const ObjT & o, FNumericLimit FError ) const;

    template < typename ObjT >
    bool Overlaps( const ObjT & o ) const;

    template < typename ObjT >
    bool Subdivide( std::vector<ObjT> &vObj );

    
    virtual BBoxT           GetBoundingBox() const = 0;
    virtual Const_VertexT & Vertex( VertexIndexT i )  const = 0;
    virtual VertexT &       Vertex( VertexIndexT i )  = 0;
   
  };
};

#endif

