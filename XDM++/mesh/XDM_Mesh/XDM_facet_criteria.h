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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.5-branch/Mesh_3/include/CGAL/Mesh_facet_criteria_3.h $
// $Id: Mesh_facet_criteria_3.h 51555 2009-08-27 13:10:21Z stayeb $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Mesh_default_facet_criteria_3 class. See class description.
//******************************************************************************

#ifndef XDM_MESH_FACET_CRITERIA_3_H
#define XDM_MESH_FACET_CRITERIA_3_H


#include <CGAL/Surface_mesher/Standard_criteria.h>
//#include <CGAL/Mesh_3/Facet_on_surface_criterion.h>

//#include <boost/optional.hpp>
#include <CGAL/Mesh_3/mesh_standard_facet_criteria.h>
#include "XDM_protecting_ball.h"
namespace CGAL
{

    namespace BYB   // criterions from the paper Boltcheva et. al.'s paper
  {
    //--------------------------------------------------------------
    //
    //  A modified criterion that takes into account edges as well as surface
    //
    //--------------------------------------------------------------
    template <typename Tr, typename Visitor_>
    class Facet_on_surface_or_edge_criterion :
      public Mesh_3::Abstract_criterion<Tr, Visitor_>
    {
    private:
      typedef typename Tr::Facet Facet;

      typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
      typedef typename Base::Quality Quality;
      typedef typename Base::Badness Badness;

      typedef Facet_on_surface_or_edge_criterion<Tr,Visitor_> Self;

    public:
      /// Constructor
      Facet_on_surface_or_edge_criterion() {};
      /// Destructor
      ~Facet_on_surface_or_edge_criterion() {};

    protected:
      virtual void do_accept(Visitor_& v) const
      {
        v.visit(*this);
      }

      virtual Self* do_clone() const
      {
        // Call copy ctor on this
        return new Self(*this);
      }

      virtual Badness do_is_bad (const Facet& f) const
      {
        typedef typename Tr::Vertex_handle Vertex_handle;
        typedef typename Tr::Cell_handle Cell_handle;

        const Cell_handle& ch = f.first;
        const int i = f.second;
        const Vertex_handle& v1 = ch->vertex((i+1)&3);
        const Vertex_handle& v2 = ch->vertex((i+2)&3);
        const Vertex_handle& v3 = ch->vertex((i+3)&3);

        // Look if vertex are on surface
        if ( ( (v1->in_dimension() != 2) && (v1->in_dimension() != 1) ) ||
             ( (v2->in_dimension() != 2) && (v2->in_dimension() != 1) ) ||
             ( (v3->in_dimension() != 2) && (v3->in_dimension() != 1) ) 
             )
        {
          return Badness(Quality(1));
        }
        else
          return Badness();
      }
    }; // end class Facet_on_surface_or_edge_criterion
  } // end namespace BYB


  
  template<typename Tr, typename Visitor_ = Mesh_3::Facet_criterion_visitor<Tr> >
  class Mesh_facet_criteria_3
  {
  private:
    typedef Visitor_ Visitor;
    typedef Mesh_3::Criteria<Tr,Visitor> Criteria;
    
    typedef typename Tr::Facet Facet;
    typedef typename Tr::Geom_traits::FT FT;

    typedef Mesh_facet_criteria_3<Tr> Self;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Cell_handle Cell_handle;
  public:
    typedef typename Visitor::Facet_quality Facet_quality;
    typedef typename Visitor::Facet_badness Facet_badness;
    
    /**
     * @brief Constructor
     */
    Mesh_facet_criteria_3(const FT& angle_bound,
                          const FT& radius_bound,
                          const FT& distance_bound)
    {
      typedef Mesh_3::Aspect_ratio_criterion<Tr,Visitor> Aspect_criterion;
      typedef Mesh_3::Uniform_size_criterion<Tr,Visitor> Uniform_size_criterion;
      typedef Mesh_3::Curvature_size_criterion<Tr,Visitor> Curvature_criterion;
      //typedef Mesh_3::Facet_on_surface_criterion<Tr,Visitor> On_surface_criterion;
      typedef BYB::Facet_on_surface_or_edge_criterion<Tr,Visitor> On_surface_criterion;
      
      if ( 0 != angle_bound )
        criteria_no_constraint_.add(new Aspect_criterion(angle_bound));
      
      if ( 0 != radius_bound )
      {
        criteria_no_constraint_.add(new Uniform_size_criterion(radius_bound));
        criteria_constraint_.add(new Uniform_size_criterion(radius_bound));
      }
      if ( 0 != distance_bound )
      {
        criteria_no_constraint_.add(new Curvature_criterion(distance_bound));
        criteria_constraint_.add(new Curvature_criterion(distance_bound));
      }
      criteria_no_constraint_.add(new On_surface_criterion());
      criteria_constraint_.add(new On_surface_criterion());
    }

    /// Destructor
    ~Mesh_facet_criteria_3() { }

    /**
     * @brief returns the badness of facet \c facet
     * @param facet the facet
     * @return the badness of \c facet
     */
    Facet_badness operator()(const Facet& facet) const
    {
      int nOverlapping = Number_of_protecting_ball_overlap( facet );
      if( nOverlapping == 0 )
        return criteria_no_constraint_(facet);
      else if ( nOverlapping < 3 )
        return criteria_constraint_(facet);
      else
        return Facet_badness();
    }

  private:

    
    //----------------------------------------------
    //  Get_number_of_protecting_ball_overlap
    //  Given the facet f, return the number of protecting ball
    //  that'd lead to overlap
    //----------------------------------------------
    int Number_of_protecting_ball_overlap( const Facet & f ) const
    {
      const Cell_handle& ch = f.first;
      const int i = f.second;
      const Vertex_handle& v1 = ch->vertex((i+1)&3);
      const Vertex_handle& v2 = ch->vertex((i+2)&3);
      const Vertex_handle& v3 = ch->vertex((i+3)&3);

      int nOverlaps = 0;
      if( XDM_test::IsOverlapping( v1, v2 ) )
        nOverlaps ++;
      if( XDM_test::IsOverlapping( v2, v3 ) )
        nOverlaps ++;
      if( XDM_test::IsOverlapping( v1, v3 ) )
        nOverlaps ++;

      return nOverlaps;
    }
    
    Criteria criteria_constraint_;
    Criteria criteria_no_constraint_;
  };  // end class Mesh_facet_criteria_3

}  // end namespace CGAL

#endif // XDM_MESH_FACET_CRITERIA_3_H
