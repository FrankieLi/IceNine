//
//  Modified version of Labeld_mesh_domain_3.h to fit the requirement of HEDM.
//  Most of the code is copied right out of Labeled_mesh_domain_3.h.
//
//
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
//
//
// Author(s)     : Frankie Li
//******************************************************************************/
// File Description :
// class XDM_mesh_domain_3. See class description.
//******************************************************************************/

#ifndef XDM_MESH_3_DOMAIN_H
#define XDM_MESH_3_DOMAIN_H

#include <CGAL/Bbox_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Mesh_3/Creator_weighted_point_3.h>

#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include "XDM_data_function_wrapper.h"
#include <boost/variant.hpp>
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/multi_array.hpp>
#include "XDM_Data.h"
#include <CGAL/Image_3.h>
#include <iostream>
#include <vector>

#include <CGAL/Regular_triangulation_filtered_traits_3.h>
namespace CGAL {
  
  /**
   * @class XDM_mesh_domain
   *           is a Labled_mesh_domain_3
   *
   */
  template<class Image,
           class BGT,
           class Wrapper = Mesh_3::XDM_data_function_wrapper<Image, BGT> >
  class XDM_mesh_domain
    : public Mesh_3::Labeled_mesh_domain_3<Wrapper, BGT>
  {
    
  public:
    typedef typename CGAL::Regular_triangulation_filtered_traits_3<BGT>    RegTriTraits;
    
    //    typedef typename BGT::Weighted_point Weighted_point_3;
    typedef typename RegTriTraits::Weighted_point Weighted_point_3;
    typedef Mesh_3::Labeled_mesh_domain_3<Wrapper, BGT> Base;
    typedef typename BGT::Point_3 Point_3;
   
    typedef typename Base::Sphere_3 Sphere_3;
    typedef typename Base::FT FT;
    typedef BGT Geom_traits;
    typedef CGAL::Bbox_3 Bbox_3;

    /// Constructor
    XDM_mesh_domain( const Image& image,
                     const FT& error_bound = FT(1e-3) )
      : Base(Wrapper(image),
             compute_bounding_box(image),
             error_bound),
        function_( Wrapper( image ) ){ };

    /// Destructor
    virtual ~XDM_mesh_domain() { };

    
    //--------------------------------
    //  Initial Point Construction object
    //--------------------------------
    struct Construct_initial_constrained_points
    {
      Construct_initial_constrained_points(const XDM_mesh_domain & domain)
        : r_domain_(domain) {}
      
      template<class OutputIterator>
      OutputIterator operator()(OutputIterator pts, const int n = 12 )
      {
        int nIndex = -3;
        std::vector< XDM_test::TestPoint > oConstrainedPnts = r_domain_.GetBndPoints();
        
        for( std::vector< XDM_test::TestPoint >::const_iterator
               pIter = oConstrainedPnts.begin();
             pIter != oConstrainedPnts.end();
             ++ pIter )
          *pts = std::make_pair( Weighted_point_3( Point_3( pIter->x, pIter->y, pIter->z ),
                                                   pIter->w  ), nIndex );
        
        return pts;
      }
      
    private:
      const XDM_mesh_domain & r_domain_;
    };

    //--------------------------------
    //  Construct_initial_topology_preserving_points
    //
    //  Insert at least a point at the center-of-mass of each
    //  domain so that there's a gaurantee that each of the
    //  grains remains
    //--------------------------------
    struct Construct_initial_topology_preserving_points
    {
      Construct_initial_topology_preserving_points(const XDM_mesh_domain & domain)
        : r_domain_(domain) {}
      
      template<class OutputIterator>
      OutputIterator operator()(OutputIterator pts, const int n = 12 )
      {
        int nIndex = -3;
        std::vector< XDM_test::TestPoint > oPts = r_domain_.GetTopologyPoints();
        
        for( std::vector< XDM_test::TestPoint >::const_iterator
               pIter = oPts.begin(); pIter != oPts.end(); ++ pIter )
          *pts = std::make_pair( Weighted_point_3( Point_3( pIter->x, pIter->y, pIter->z ),
                                                   pIter->w  ), nIndex );
        
        return pts;
      }
    private:
      const XDM_mesh_domain & r_domain_;
    };
    
 
    /// Returns Construct_initial_points object
    Construct_initial_constrained_points
    construct_initial_constrained_points_object() const
    {
      return Construct_initial_constrained_points(*this);
    }

    /// Returns Construct_initial_points object
    Construct_initial_topology_preserving_points
    construct_initial_topology_preserving_points_object() const
    {
      return Construct_initial_topology_preserving_points(*this);
    }


    
    
  private:
    /// Returns a box enclosing image \c im
    Bbox_3 compute_bounding_box(const Image& im) const
    {
      return im.compute_bounding_box();
    }

    
    //--------------------------------
    //--------------------------------
    std::vector< XDM_test::TestPoint > GetBndPoints() const
    {
      return function_.GetBndPoints();
    }
      
      
    //--------------------------------
    //--------------------------------
    std::vector< XDM_test::TestPoint > GetTopologyPoints() const
    {
      return function_.GetTopologyPoints();
    }

  private:
    // Disabled copy constructor & assignment operator
    typedef XDM_mesh_domain<Image, BGT, Wrapper> Self;
    XDM_mesh_domain(const Self& src);
    Self& operator=(const Self& src);

    Wrapper function_;
  };  // end class XDM_mesh_domain


  
  
}  // end namespace CGAL

#endif // XDM_mesh_domain_3
