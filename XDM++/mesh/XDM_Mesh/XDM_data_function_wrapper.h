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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.5-branch/Mesh_3/include/CGAL/Mesh_3/Image_to_labeled_function_wrapper.h $
// $Id: Image_to_labeled_function_wrapper.h 51555 2009-08-27 13:10:21Z stayeb $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Image_to_labeled_function_wrapper declaration and implementation. See
// class description.
//******************************************************************************

#ifndef XDM_DATA_FUNCTION_WRAPPER_H
#define XDM_DATA_FUNCTION_WRAPPER_H

#include <CGAL/Image_3.h>
#include "XDM_Data.h"
namespace CGAL {

namespace Mesh_3 {

/**
 * @class Image_to_labeled_function_wrapper
 *
 * Wraps a labeled image into a labeled function which takes his values into
 * N. Uses trilinear interpolation.
 * Note: Image has to be labeled with unsigned char
 */
template<class XDM_Data,
         class BGT,
         typename Return_type = int,
         bool use_trilinear_interpolation=true>
         
class XDM_data_function_wrapper
{
public:
  typedef Return_type return_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  XDM_data_function_wrapper( const XDM_Data& image )
    : oData(image) { };

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~XDM_data_function_wrapper() { };

 
  std::vector< XDM_test::TestPoint> GetBndPoints( ) const
  {
    return oData.GetBndPoints( );
  }

   
  std::vector< XDM_test::TestPoint> GetTopologyPoints( ) const
  {
    return oData.GetTopologyPoints( );
  }

  /**
   * Returns an int corresponding to the label at point \c p
   * @param p the input point
   * @return the label at point \c p
   */
  return_type operator()( const Point_3& p, const bool = true) const
  {
//     if ( use_trilinear_interpolation )
//     {
    return
        static_cast<return_type>( oData.trilinear_interpolation( CGAL::to_double(p.x()),
                                                                 CGAL::to_double(p.y()),
                                                                 CGAL::to_double(p.z()) ) );
   //  }
//     else
//     {
//       return oData( CGAL::to_double(p.x()),
//                     CGAL::to_double(p.y()),
//                     CGAL::to_double(p.z()) );
//     }
  }



  //  need to overload construction of initial point estimate
private:
  
  const XDM_Data & oData;

};  // end class Image_to_labeled_function_wrapper


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // 
