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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.5-branch/Mesh_3/include/CGAL/make_mesh_3.h $
// $Id: make_mesh_3.h 51555 2009-08-27 13:10:21Z stayeb $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description : make_mesh_3 function definition.
//******************************************************************************

#ifndef XDM_MAKE_MESH_3_H
#define XDM_MAKE_MESH_3_H

#include <iostream>
namespace CGAL {


/**
 * @brief This function meshes the domain defined by mesh_traits
 * (respecting criteria), and outputs the mesh to c3t3
 *
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if it is set to \c true, an exudation step will be done at
 *   the end of the Delaunay refinment process
 *
 * @return The mesh as a C3T3 object
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
C3T3 XDM_make_mesh_3(const MeshDomain&   domain,
                     const MeshCriteria& criteria,
                     bool exude = true)
{
  typedef typename MeshDomain::Point_3 Point_3;
  typedef typename MeshDomain::Weighted_point_3 Weighted_point_3;
  typedef typename MeshDomain::Index Index;
  typedef std::vector<std::pair<Point_3, Index> > Initial_points_vector;
  typedef std::vector<std::pair<Weighted_point_3, Index> > Constrained_points_vector;
  
  typedef typename Initial_points_vector::iterator Ipv_iterator;
  typedef typename Constrained_points_vector::iterator Cpv_iterator;
  typedef typename C3T3::Vertex_handle Vertex_handle;
  
  // Mesh initialization : get some points and add them to the mesh
  Initial_points_vector initial_points;
  Constrained_points_vector constrained_points;
  
  domain.construct_initial_constrained_points_object()( std::back_inserter( constrained_points ) );
  domain.construct_initial_points_object()( std::back_inserter( initial_points ) );
  domain.construct_initial_topology_preserving_points_object()( std::back_inserter( initial_points ) );;
  C3T3 c3t3;

  for ( Cpv_iterator it = constrained_points.begin() ;
        it != constrained_points.end() ;
        ++it )
  {
    Vertex_handle v = c3t3.triangulation().insert( it->first );

    v->SetProtectingRadius( it->first.weight() );
    c3t3.set_dimension(v,1); // by construction, points bnd points
    c3t3.set_index(v,it->second);
    
  }
  std::cout << "vertices facets edges " << std::endl;
  std::cout << c3t3.triangulation().number_of_vertices() << " "
            << c3t3.triangulation().number_of_facets() << " "  
            << c3t3.triangulation().number_of_edges() << std::endl;
  
  
  // Insert points and set their index and dimension
  for ( Ipv_iterator it = initial_points.begin() ;
        it != initial_points.end() ;
        ++it )
  {
    bool bCloseToConstraints = false;
    for ( Cpv_iterator it2 = constrained_points.begin() ;
          it2 != constrained_points.end() ;
          ++it2 )
    {
      float r = it2->first.weight();

      float dx = it->first.x() - it2->first.x();
      float dy = it->first.y() - it2->first.y();
      float dz = it->first.z() - it2->first.z(); 

      if( r*r > dx * dx + dy * dy + dz * dz )
        bCloseToConstraints = true;
    }
    
    if( !bCloseToConstraints )
    {
      Weighted_point_3 wp( it->first, 0 );
      //    Vertex_handle v = c3t3.triangulation().insert(it->first);
      Vertex_handle v = c3t3.triangulation().insert( wp );
      c3t3.set_dimension(v,2); // by construction, points are on surface
      c3t3.set_index(v,it->second);
      
    }
  }

   std::cout << "vertices facets edges " << std::endl;
   std::cout << c3t3.triangulation().number_of_vertices() << " "
             << c3t3.triangulation().number_of_facets() << " "  
             << c3t3.triangulation().number_of_edges() << std::endl;
   
  // Build mesher and launch refinement process
  // Don't reset c3t3 as we just created it
   //   refine_mesh_3(c3t3, domain, criteria, true, false);  // try to exude
   refine_mesh_3(c3t3, domain, criteria, false, false);  // turned of exudation - don't need it
  std::cout << "vertices facets edges " << std::endl;
  std::cout << c3t3.triangulation().number_of_vertices() << " "
            << c3t3.triangulation().number_of_facets() << " "  
            << c3t3.triangulation().number_of_edges() << std::endl;
  return c3t3;
};


}  // end namespace CGAL


#endif // CGAL_MAKE_MESH_3_H
