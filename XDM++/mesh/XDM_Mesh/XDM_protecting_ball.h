//---------------------------------------------------------
//
//  XDM_protecting_ball
//---------------------------------------------------------


#ifndef XDM_PROTECTING_BALL_H_
#define XDM_PROTECTING_BALL_H_


namespace CGAL
{
namespace XDM_test
{
  template< class Vertex_handle>
  bool IsOverlapping( Vertex_handle v1, Vertex_handle v2 )
  {
    if( v1->in_dimension() != 1 || v2->in_dimension() != 1 )
      return false;
    
    float rMax = v1->ProtectingRadius() + v2->ProtectingRadius();
    
    float dx = v1->point().x() - v2->point().x();
    float dy = v1->point().y() - v2->point().y();
    float dz = v1->point().z() - v2->point().z();
    
    return ( rMax * rMax < dx * dx + dy * dy + dz * dz );
  }
}
}

#endif
