//----------------------------------------------------
//  SpecialTripleLineDistance
//
//  Author:   Frankie Li
//  Purpose:  Calculate the distance between special triple line
//            to a given point.
//            This is a class built as an Oracle.
//----------------------------------------------------


#ifndef TRIPLE_LINE_DISTANCE_H
#define TRIPLE_LINE_DISTANCE_H

//#include <CGAL/basic.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include "XDM_Mesh/AABB_mesh_3_triangle_primitive.h"

#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>

#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <limits>
#include <CGAL/bounding_box.h>
#include "MeshUtilities.h"
#include <CGAL/Cartesian_converter.h>
#include "Symmetry.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


namespace Pandora
{

  //-------------------------------------------------------------
  //  SpecialTripleLineDistance
  //  Purpose:  Calculate the distance between a spetial triple line and
  //            a given point.
  //-------------------------------------------------------------
  template< class C3T3 >
  class SpecialTripleLineDistance
  {
  public:
        
    struct Kernel: public CGAL::Exact_predicates_inexact_constructions_kernel {};
    typedef typename C3T3::Triangulation::Vertex_iterator            Vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_vertices_iterator   Finite_vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_facets_iterator     Finite_facet_iterator;
    typedef typename C3T3::Triangulation::Vertex                     Vertex;
    typedef typename C3T3::Triangulation::Vertex_handle              Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits::Bare_point    Point3;

    typedef typename C3T3::Cell_handle                               Cell_handle;
    typedef typename C3T3::Surface_index  Surface_index;
    typedef typename C3T3::Subdomain_index                           Subdomain_index;


    typedef typename C3T3::Triangulation::Locate_type     Locate_type;
    typedef typename C3T3::Triangulation::Facet           Facet;

    typedef typename C3T3::Triangulation::Geom_traits::Triangle_3  Triangle_3;
    typedef typename C3T3::Triangulation                           Triangulation;

    typedef typename C3T3::Triangulation::Geom_traits::Vector_3    Vector_3;
  
    typedef typename C3T3::Facet_iterator                          Facet_iterator;
    
    typedef GeneralLib::SQuaternion                                SQuaternion;
    
    //-------------------------------------------------------
    //  Random sampling
    //-------------------------------------------------------
    typedef GeneralLib::SMatrix3x3 SMatrix3x3;
    typedef GeneralLib::SVector3   SVector3;
    typedef GeneralLib::UniformGrid::CQuaternionGrid CQuaternionGrid;

    mutable GeneralLib::CUniformRandomReal RandGen;
    
    //-------------------------------------------------------
    //  kdtree
    //-------------------------------------------------------
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;
    typedef typename Neighbor_search::Tree Tree;


    
    typedef typename C3T3::Triangulation::Geom_traits Tr_Kernel;
    typedef CGAL::Simple_cartesian<double> AABB_Kernel;
    
    typedef typename AABB_Kernel::Ray_3     AABB_Ray;
    typedef typename AABB_Kernel::Segment_3 AABB_Segment;
    typedef typename AABB_Kernel::Point_3   AABB_Point;
    typedef typename AABB_Kernel::FT        AABB_FT;
    
    typedef typename C3T3::Triangulation::Geom_traits::FT FT;

    typedef typename std::list< Finite_facet_iterator >                                 Filtered_finite_facet_list;
    typedef typename std::list< Finite_facet_iterator >::const_iterator                 Filtered_facet_iterator;

    //--------------------------------------
    typedef std::list<AABB_Segment>::iterator                                         Segment_Iterator;
    typedef typename CGAL::AABB_segment_primitive<AABB_Kernel,  Segment_Iterator>     AABB_Primitive;
    typedef typename CGAL::AABB_traits<AABB_Kernel, AABB_Primitive>                   AABB_Mesh_Segment_Traits;
    typedef typename CGAL::AABB_tree<AABB_Mesh_Segment_Traits>                       Mesh_AABB_Tree;
    
    typedef typename Mesh_AABB_Tree::Object_and_primitive_id                          Object_and_primitive_id;
    typedef typename Mesh_AABB_Tree::Primitive_id                                     Primitive_id;
    typedef typename Mesh_AABB_Tree::Point_and_primitive_id                           Point_and_primitive_id;
    typedef CGAL::Cartesian_converter< Tr_Kernel, AABB_Kernel>                        TrKernel_To_AABB;

    //--------------------------------------
  

    
  public:
    
    //----------
    //  Also need a sigma map
    //----------
    SpecialTripleLineDistance( const C3T3 & Mesh,
                               std::vector<SMatrix3x3> & IDOrientMap_  ):
      Mesh_( Mesh ),
      IDOrientMap( IDOrientMap_ )
    {}
    
    
    //--------------------------------------
    //  DefineTripleJunctionAsSink
    //--------------------------------------
    void DefineTripleJunctionAsSink()
    {
      typedef typename C3T3::Triangulation::Edge_iterator Edge_iterator;
      typedef typename C3T3::Triangulation::Facet_circulator  Facet_circulator;
      typedef typename C3T3::Triangulation::Edge    Edge;

      SegmentList.clear();
      TotalNJunctionLength = 0;

      for( Edge_iterator pEdge = Mesh_.triangulation().edges_begin();
           pEdge != Mesh_.triangulation().edges_end(); ++ pEdge )
      {
        Facet_circulator pEnd = Mesh_.triangulation().incident_facets( *pEdge );
        Facet_circulator pCur = Mesh_.triangulation().incident_facets( *pEdge );
        if( ! Mesh_.triangulation().is_infinite( *pEdge ) )
        {
          bool bInterior = true;
          std::set< Subdomain_index > DomainSet;
          if( pCur != 0 ) // construct of using Cell_circulator
          {
            do
            {
              if( Mesh_.is_in_complex( *pCur ) )
              {
                Surface_index SurfID = Mesh_.surface_index( *pCur );
                DomainSet.insert( SurfID.first );
                DomainSet.insert( SurfID.second );
                if( SurfID.first < 0 || SurfID.second < 0 )
                  bInterior = false;
              }
              ++ pCur;
            } while( pCur != pEnd );
          }
          
          if( bInterior && DomainSet.size() > 2 )
          {
            Edge e = *pEdge;
            Cell_handle ch = e.get<0>();

            Point3 BarePoint1 = ( ch->vertex( e.get<1>() ) )->point().point();
            Point3 BarePoint2 = ( ch->vertex( e.get<2>() ) )->point().point();
            AABB_Point p1( BarePoint1.x(), BarePoint1.y(), BarePoint1.z() );
            AABB_Point p2( BarePoint2.x(), BarePoint2.y(), BarePoint2.z() );
            AABB_Segment Segment( p1, p2 );
            SegmentList.push_back( Segment );
            TotalNJunctionLength += CGAL::sqrt( Segment.squared_length() );

          }
        }// end if infinte edge
      }// end for

      NLineAABBTree.rebuild( SegmentList.begin(), SegmentList.end() );
      
    }

    //--------------------------------------
    //  DistanceToNearestSink_orig
    //
    //  Purpose:  Get the distance to the nearest "sink" - the types of
    //            said sink is dictated by the initialization.  For example,
    //            sinks can be special (of specific types) boundaries, free
    //            surfaces, or even special "grains."
    //
    //
    //            Using AABB Tree only - see http://www-compsci.swan.ac.uk/~csmark/PDFS/dist.pdf
    //--------------------------------------
    double DistanceToNearestSink( const Point3 &SearchCenter )
    {
      AABB_Point AABB_SearchCenter( SearchCenter.x(),
                                    SearchCenter.y(),
                                    SearchCenter.z() );
      
      AABB_Point Closest_Point = NLineAABBTree.closest_point( AABB_SearchCenter );
      double AABB_Dist  = CGAL::sqrt( CGAL::squared_distance( Closest_Point, AABB_SearchCenter ) );
      return AABB_Dist;
    }
    
    //--------------------------------------
    //  GetMeshIDSet
    //--------------------------------------
    std::set<int> GetMeshIDSet( const C3T3 & c3t3_ )
    {
      std::set<int> IDMap;
      typedef typename C3T3::Surface_index   Surface_index;
      typedef typename C3T3::Cell_handle     Cell_handle;
      typedef typename C3T3::Cell_iterator   Cell_iterator;

      for( Cell_iterator pCur = c3t3_.cells_begin();
           pCur != c3t3_.cells_end(); ++ pCur )
      {
        Cell_handle ch = pCur;
        int nID = c3t3_.subdomain_index( ch );
        IDMap.insert( nID );
      }
      return IDMap;
    }

    //--------------------------------------
    //  GetGrainID
    //--------------------------------------
    int GetGrainID( const Point3 & QueryLoc )
    {
      Cell_handle ch = Mesh_.triangulation().locate( QueryLoc );
      Subdomain_index GrainID =  Mesh_.subdomain_index( ch );

      if( ! Mesh_.triangulation().is_infinite( ch ) ) // not infinite => typical
        return GrainID;
      return -1;
    }

    //--------------------------------------
    //  GetTotalNJunctionLength()
    //  Precondition:  Sink boundary must be already defined.
    //--------------------------------------
    double GetTotalNJunctionLength() const
    {
      return TotalNJunctionLength;
    }
    
    //--------------------------------------
    //  Precondition:  Grain Exists
    //--------------------------------------
    GeneralLib::SQuaternion GetGrainOrientation( int GrainID )
    {
      GeneralLib::SQuaternion q;
      DEBUG_ASSERT( GrainID >= 0, "Precondition failed: Requirement: GrainID >= 0 \n");
      q.Set( IDOrientMap[ GrainID ] );
      return q;
    }
      
  private:
    
    //--------------------------------------
    //  IsGrainBoundary
    //--------------------------------------
    bool IsGrainBoundary( const Facet & f )
    {
      Surface_index ID = Mesh_.surface_index( f );
      return ( ID.first >= 0 && ID.second >= 0 && ID.first != ID.second );
    }
    
    //---------------
    //  Data
    //---------------
    
    std::vector<SMatrix3x3> & IDOrientMap;
    
    const C3T3 & Mesh_;
    
    Mesh_AABB_Tree NLineAABBTree;
    double  TotalNJunctionLength;
    std::list<AABB_Segment> SegmentList;
  };


}

#endif
