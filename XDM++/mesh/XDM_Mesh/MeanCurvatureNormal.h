//----------------------------------------------------
// MeanCurvatureNormal.h
//
//  Author:   Frankie Li
//  Purpose:  Implementation of MeanCurvatureNormal from Schroeder's paper
//            
//----------------------------------------------------

#ifndef MEAN_CURVATURE_NORMAL_H_
#define MEAN_CURVATURE_NORMAL_H_

#include <vector>

namespace Pandora
{
  namespace Details
  {
    template< class C3T3 >
    class MeanCurvatureNormal
    {
    private:
      typedef typename C3T3::Vertex_handle                           Vertex_handle;
      typedef typename C3T3::Triangulation::Geom_traits::Vector_3    Vector_3;
      typedef typename C3T3::Triangulation::Geom_traits::RT          RT;
      typedef typename C3T3::Facet                                   Facet;
      typedef typename C3T3::Triangulation::Geom_traits::Bare_point  Point3;
      typedef typename C3T3::Surface_index                           Surface_index;


    public:

      //--------------------------------
      //  operator() 
      //
      //  The way this is implemented is slightly out of order than the sum presented in Meyer.
      //  The reason is to simply the algorithm.  We are splitting the original:
      //
      //   A_mix = 1/8 \sum_{ j in ngb(i) } ( cot( alpha_ij ) + cot( beta_ij ) ) | x_i - x_j |^2
      //
      //   K = 1/( 2 A_mix )  \sum_{ j in ngb(i) } ( cot( alpha_ij ) + cot( beta_ij ) ) ( x_i - x_j )
      //
      //
      //   Note that the *mean curvature* is given by |K|/2,  
      //
      //
      //--------------------------------
      std::pair<Vector_3, double> operator() ( const Vertex_handle & vh, const C3T3 & Mesh,
                                               Surface_index BndSurfaceIndex, bool RunDebug = false )
      {
        std::vector<Facet> FacetList;
        Mesh.triangulation().finite_incident_facets( vh, std::back_inserter( FacetList ));      
        
        double A_obtuse_half  = 0;
        double A_obtuse_quart = 0;
        double A_acute        = 0;
        
        Vector_3  K( Point3(0, 0, 0), Point3(0, 0, 0) );

        if( RunDebug )
          std::cout << " Start: Num Facets " << FacetList.size() << std::endl; 
        for( int i = 0; i < FacetList.size(); i ++ )
        {
          if( Mesh.is_in_complex( FacetList[i] )
              &&  BndSurfaceIndex == Mesh.surface_index( FacetList[i] ) )
          {
            if( RunDebug )
              std::cout << " Facet IDs " << Mesh.surface_index( FacetList[i] ).first
                        << " " << Mesh.surface_index( FacetList[i] ).second << " ";
            //--------------------------
            int vID[2];
            int n = 0;
            for( int k = 0; k < 4; k ++ )
            {
              if ( k != FacetList[i].second && FacetList[i].first->vertex(k) != vh  )
              {
                vID[n] = k;
                ++ n;
              }
            }
            RUNTIME_ASSERT( n == 2, " Number of facet index traversed is not exactly 2" );
            Point3 p0 = vh->point();
            Point3 p1 = FacetList[i].first->vertex(vID[0])->point();
            Point3 p2 = FacetList[i].first->vertex(vID[1])->point();

            Vector_3 r12 = p2 - p1;
            Vector_3 r10 = p0 - p1;   // flipped on purpose
            Vector_3 r02 = p2 - p0;

            //-------------------------------------------  DEBUG
            typedef typename C3T3::Triangulation::Geom_traits::Plane_3    Plane_3;
            Plane_3 FacetPlane( p0, p1, p2 );
            Vector_3 PlaneNormal = FacetPlane.orthogonal_vector(); 
            if( RunDebug )
              std::cout << " [ " << PlaneNormal / CGAL::sqrt( PlaneNormal.squared_length() ) << " ] ";
            

            RT cos_theta1 =  (r10 * r12)  /  CGAL::sqrt( r10.squared_length() * r12.squared_length() );
            RT cos_theta2 =  (r02 * r12)  /  CGAL::sqrt( r02.squared_length() * r12.squared_length() );
            
            RT cot_theta1 = cos_theta1 / CGAL::sqrt(  1. - (cos_theta1 * cos_theta1) );
            RT cot_theta2 = cos_theta2 / CGAL::sqrt(  1. - (cos_theta2 * cos_theta2) );             

            //--------------------------
            // check obtuse, calculate mixed area
            if( CGAL::OBTUSE == CGAL::angle( p1, p0, p2 ) )
            {
              A_obtuse_half += CGAL::sqrt( Mesh.
                                           triangulation().
                                           triangle( FacetList[i] ).
                                           squared_area() );      // to be divided by 2
            }
            else if(  CGAL::OBTUSE == CGAL::angle( p0, p1, p2 )
                      || CGAL::OBTUSE == CGAL::angle( p0, p2, p1 ) )
            {
              A_obtuse_quart +=  CGAL::sqrt( Mesh.
                                             triangulation().
                                             triangle( FacetList[i] ).
                                             squared_area() );    // to be divided by 1/2
            }
            else  // otherwise, no obtuse angle, voronoi area
            {
              A_acute += ( r10.squared_length() * cot_theta2
                           + r02.squared_length() * cot_theta1 );
            }
            
            K = K +  ( -r10 * cot_theta2 + r02 * cot_theta1 );  // negative sign on purpose
            
            if( RunDebug )
            {
              std::cout << A_obtuse_half << " " << A_obtuse_quart << " " << A_acute << " " << std::endl
                        << " |---- " <<   ( r10 * cot_theta2 - r02 * cot_theta1 ) << std::endl
                        << " |---- " <<   r10 << std::endl  
                        << " |---- " <<   r02 << std::endl 
                        << "   |____ " << K <<  " || " << K / CGAL::sqrt( K.squared_length() ) << std::endl;
            }
          }// end if on complex

        } // end For
        RT MixedArea = A_acute        / RT( 8 )
          + A_obtuse_quart / RT( 4 )
          + A_obtuse_half  / RT( 2 ) ; 
        K = K / ( MixedArea * RT( 2 ) );
        if( RunDebug )
          std::cout << " End: Normal: "  << K / CGAL::sqrt( K.squared_length() ) <<  " || " << MixedArea << std::endl;
        
        return std::make_pair( K, MixedArea );
      }  // end operator()
      
    
    };
   
    
    
  } // end namespace Details

  
} // end namespace Pandora

#endif
