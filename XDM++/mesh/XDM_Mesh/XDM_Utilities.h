#ifndef XDM_UTILITIES_H
#define XDM_UTILITIES_H

#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/grid_simplify_point_set.h>

namespace CGAL
{
namespace XDM_test
{
  struct TestPoint
  {
    TestPoint() {}
    TestPoint( float x_, float y_,
               float z_, float w_ ):
        x( x_ ), y( y_ ), z( z_ ), w( w_ ){}
    
    float x;
    float y;
    float z;
    float w;
  };
}

}

namespace CGAL
{
  // modifications of CGAL's point processing method so that it
  // is aware of special points (quad points vs. triple points)
  namespace XDM_PointProcessing
  {
    //-------------------------------------------------------------------------------------
    //
    //  This is directly taken out of CGAL 3.6.1's grid_simplify_point_set and modified
    //
    //-------------------------------------------------------------------------------------
    /// Merges points which belong to the same cell of a grid of cell size = epsilon.
    ///
    /// This method modifies the order of input points so as to pack all remaining points first,
    /// and returns an iterator over the first point to remove (see erase-remove idiom).
    /// For this reason it should not be called on sorted containers.
    ///
    /// @commentheading Precondition: epsilon > 0.
    ///
    /// @commentheading Template Parameters:
    /// @param ForwardIterator iterator over input points.
    /// @param PointPMap is a model of boost::ReadablePropertyMap with a value_type = Point_3<Kernel>.
    ///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
    /// @param Kernel Geometric traits class.
    ///        It can be omitted and deduced automatically from PointPMap value_type.
    ///
    /// @return iterator over the first point to remove.

    // This variant requires all parameters.
    template <typename ForwardIterator,
              typename PointPMap,
              typename PointLabelPMap>
              //         typename Kernel>
    ForwardIterator grid_simplify_point_set(
                                            ForwardIterator first,  ///< iterator over the first input point.
                                            ForwardIterator beyond, ///< past-the-end iterator over the input points.
                                            PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
                                            PointLabelPMap label_pmap,  ///< -- property map to indicate label
                                            double epsilon) ///< tolerance value when merging 3D points.
                                            //const Kernel& kernel) ///< geometric traits.
    {
      // actual type of input points
      typedef typename std::iterator_traits<ForwardIterator>::value_type Enriched_point;
      
      CGAL_point_set_processing_precondition(epsilon > 0);
      
      // Merges points which belong to the same cell of a grid of cell size = epsilon.
      // points_to_keep[] will contain 1 point per cell; the others will be in points_to_remove[].
      Epsilon_point_set_3<Enriched_point, PointPMap> points_to_keep(epsilon, point_pmap);
      std::deque<Enriched_point> points_to_remove;
      std::deque<Enriched_point> special_points;
      
      for (ForwardIterator it=first ; it != beyond ; it++)
      {
        std::pair<typename Epsilon_point_set_3<Enriched_point, PointPMap>::iterator,bool> result;
        result = points_to_keep.insert(*it);
        if (!result.second) // if not inserted
        {
     //      int nLabel = get( label_pmap, it );
//           if( nLabel >= 4 )  // Quadruple point
//           {
//             if( get( label_pmap, result.first ) < 4 )   // not dual special point
//             {
//               points_to_remove.push_back( *( result.first ) );   // replacing with special point 
//               points_to_keep.erase( result.first );
//               points_to_keep.insert( *it );
//             }
//             else
//             {
//               special_points.push_back( *it );
//             }
//           }
          points_to_remove.push_back(*it);
        }
      }
      
      // Replaces [first, beyond) range by the content of points_to_keep, then points_to_remove.
      ForwardIterator first_point_to_remove =
        std::copy(points_to_keep.begin(), points_to_keep.end(), first);

      std::copy(special_points.begin(), special_points.end(), first_point_to_remove);
      std::copy(points_to_remove.begin(), points_to_remove.end(), first_point_to_remove);
      
      return first_point_to_remove;
    }
  }
}
#endif
