//
//  This file is modified from AABB_polygon_triangle_primitive.h
//  The current adaptation allows us to use facet handles from the Complex_3_in_3_dimensions
//
//**************************************************************************/
// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     : Stéphane Tayeb, Pierre Alliez
//
//******************************************************************************/
// File Description :
//
//******************************************************************************/
  
#ifndef AABB_Mesh_3_TRIANGLE_PRIMITIVE_H_
#define AABB_Mesh_3_TRIANGLE_PRIMITIVE_H_

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/iterator.h>

namespace CGAL {


  namespace Pandora
  {
    template< class C3T3, class Iterator >
    struct BoundarySelectionPredicate
    {
      typedef typename C3T3::Surface_index Surface_index;
      C3T3 * m_MeshPtr;
      BoundarySelectionPredicate( C3T3  *mesh )
        : m_MeshPtr( mesh ) { }
      
      bool operator()( const Iterator & it ) const
      {
        Surface_index SurfID = m_MeshPtr->surface_index( *it );
        return SurfID.first != SurfID.second;
      }
      
    };
  }
  
  

  /**
   * @class AABB_mesh_3_triangle_primitive
   *
   *
   */
  template< class C3T3, class AABB_Kernel, class Iterator >
  class AABB_mesh_3_triangle_primitive
  {
  public:
    typedef typename C3T3::Triangulation Triangulation;
    /// AABBPrimitive types
    typedef AABB_Kernel               GeomTraits;
    typedef typename GeomTraits::Point_3        Point;
    typedef typename GeomTraits::Triangle_3     Datum;
    typedef Iterator                   Id;
    typedef typename Triangulation::Finite_facets_iterator Finite_facet_iterator;

    typedef typename Triangulation::Geom_traits Tr_Kernel;
    typedef typename Tr_Kernel::Triangle_3      Tr_Triangle_3;
    
    typedef CGAL::Cartesian_converter< Tr_Kernel, GeomTraits>    TrKernel_To_AABB;
    /// Constructors
    AABB_mesh_3_triangle_primitive() {}
    AABB_mesh_3_triangle_primitive(const AABB_mesh_3_triangle_primitive& primitive)
      :  m_facet_handle( primitive.id() ), m_Datum( primitive.datum() )
    {
      
    }
    
    AABB_mesh_3_triangle_primitive( Id  IteratorID )
      : m_facet_handle( IteratorID )
    {
      TrKernel_To_AABB To_AABB;
      Finite_facet_iterator handle = *IteratorID;
      m_Datum = To_AABB( Tr_Triangle_3 ( handle->first->vertex( (handle->second+1) & 3)->point(),
                                         handle->first->vertex( (handle->second+2) & 3)->point(),
                                         handle->first->vertex( (handle->second+3) & 3)->point() )  );
    }
    

    // Default destructor, copy constructor and assignment operator are ok
    
    /// Returns by constructing on the fly the geometric datum wrapped by the primitive
    Datum datum() const
    {
      return m_Datum;
    }
    
    /// Returns a point on the primitive
    Point reference_point() const
    {
      return m_Datum.vertex(0);
    }
    
    
    /// Returns the identifier
    const Id& id() const { return m_facet_handle; }
    Id& id() { return m_facet_handle; }
    
  private:
    /// The id, here a polyhedron facet handle
    Id m_facet_handle;
    Datum m_Datum;
  };  // end class AABB_mesh_3_triangle_primitive
  

}  // end namespace CGAL


#endif // AABB_POLYHEDRON_TRIANGLE_PRIMITIVE_H_
