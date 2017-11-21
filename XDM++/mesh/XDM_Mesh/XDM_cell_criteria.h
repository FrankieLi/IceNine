// Copyright (c) 2004-2009  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.5-branch/Mesh_3/include/CGAL/Mesh_cell_criteria_3.h $
// $Id: Mesh_cell_criteria_3.h 51094 2009-08-06 13:11:07Z stayeb $
// 
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb


#ifndef XDM_MESH_CELL_CRITERIA_3_H
#define XDM_MESH_CELL_CRITERIA_3_H

#include <iostream>
#include <CGAL/Mesh_3/mesh_standard_cell_criteria.h>
#include "XDM_protecting_ball.h"


namespace CGAL
{
  template <typename Tr, typename Visitor_ = Mesh_3::Cell_criterion_visitor<Tr> >
  class Mesh_cell_criteria_3
  {
    typedef Visitor_ Visitor;
    typedef Mesh_3::Criteria<Tr,Visitor> Criteria;
    typedef typename Tr::Cell_handle Cell_handle;
    typedef typename Tr::Geom_traits GT;
    typedef typename Tr::Geom_traits::FT FT;
    typedef typename Tr::Geom_traits::Weighted_point Weighted_point;
    typedef typename Tr::Geom_traits::Bare_point     Point;
  
    typedef Mesh_cell_criteria_3<Tr> Self;
  
  public:
    typedef typename Visitor::Cell_quality Cell_quality;
    typedef typename Visitor::Cell_badness Cell_badness;
  
    /**
     * @brief Constructor
     * @param radius_edge_bound the radius-edge bound
     * @param radius_bound the radius bound (tet sizing)
     */
    Mesh_cell_criteria_3(const FT& radius_edge_bound,
                         const FT& radius_bound )
    {
      typedef Mesh_3::Cell_radius_criterion<Tr,Visitor> Radius_criterion;
      typedef Mesh_3::Cell_radius_edge_criterion<Tr,Visitor> Radius_edge_criterion;
    
      if ( 0 != radius_bound )
        criteria_no_constraint_.add(new Radius_criterion(radius_bound));
    
      if ( 0 != radius_edge_bound )
      {
        criteria_no_constraint_.add(new Radius_edge_criterion(radius_edge_bound));
        criteria_constraint_.add(new Radius_edge_criterion(radius_edge_bound));
      }
    }
  
    /// Destructor
    ~Mesh_cell_criteria_3() { }
  
    /**
     * @brief returns the badness of cell \c cell
     * @param cell the cell
     * @return the badness of \c cell
     */
    Cell_badness operator()(const Cell_handle& cell) const
    {
      int nOverlaps = Number_of_protecting_ball_overlap( cell );
    
      
      if( nOverlaps == 3 )
        if( Refinement_point_overlaps_protecting_balls( cell ) )
          return criteria_constraint_(cell);
        else
          return Cell_badness();
      else if( nOverlaps == 1 || nOverlaps == 2 )
        return criteria_constraint_(cell);
      else
        return criteria_no_constraint_(cell);
     
    }
  
  private:
    
    bool Refinement_point_overlaps_protecting_balls( const Cell_handle& cell ) const
    {
      typedef typename Tr::Vertex_handle Vertex_handle;

      //  clearly - this is a cheat so that I don't have to rewrite Refine_cells.h
      Point pToInsert =
        GT().construct_weighted_circumcenter_3_object()( cell->vertex(0)->point(),
                                                         cell->vertex(1)->point(),
                                                         cell->vertex(2)->point(),
                                                         cell->vertex(3)->point() );
      for( int i = 0; i < 4; i ++ )
      {
        const Vertex_handle & v = cell->vertex(i);
        const Weighted_point & p = v->point();
        float dx = p.x() - pToInsert.x();
        float dy = p.y() - pToInsert.y();
        float dz = p.z() - pToInsert.z();
        float r = v->ProtectingRadius();
        if( dx * dx + dy * dy + dz * dz < r * r )
          return true;
      }
      return false;
    }
    
    //----------------------------------------------
    //  Get_number_of_protecting_ball_overlap
    //  Given the facet f, return the number of protecting ball
    //  that'd lead to overlap
    //----------------------------------------------
    int Number_of_protecting_ball_overlap( const Cell_handle& f ) const
    {
      typedef typename Tr::Vertex_handle Vertex_handle;
      const Vertex_handle& v1 = f->vertex(0);
      const Vertex_handle& v2 = f->vertex(1);
      const Vertex_handle& v3 = f->vertex(2);
      const Vertex_handle& v4 = f->vertex(3);

      
      int pOverlapIndex[4] = { 0, 0, 0, 0 };
      if( XDM_test::IsOverlapping( v1, v2 ) )
      {
        pOverlapIndex[0] = 1; pOverlapIndex[1] = 1;
      }
      if( XDM_test::IsOverlapping( v1, v3 ) )
      {
        pOverlapIndex[0] = 1; pOverlapIndex[2] = 1;
      }
      if( XDM_test::IsOverlapping( v1, v4 ) )
      {
        pOverlapIndex[0] = 1; pOverlapIndex[3] = 1;
      }
      if( XDM_test::IsOverlapping( v2, v3 ) )
      {
        pOverlapIndex[1] = 1; pOverlapIndex[2] = 1;
      }
      if( XDM_test::IsOverlapping( v2, v4 ) )
      {
        pOverlapIndex[1] = 1; pOverlapIndex[3] = 1;
      }
      if( XDM_test::IsOverlapping( v3, v4 ) )
      {
        pOverlapIndex[2] = 1; pOverlapIndex[3] = 1;
      }
      int nOverlaps = 0 ;
      for( int i = 0; i < 4; i ++ )
        if( pOverlapIndex[i] > 0 )
          nOverlaps ++;
      
      return nOverlaps;
    }

    Criteria criteria_no_constraint_;
    Criteria criteria_constraint_;
  };  // end class Mesh_cell_criteria_3

}  // end namespace CGAL


#endif
