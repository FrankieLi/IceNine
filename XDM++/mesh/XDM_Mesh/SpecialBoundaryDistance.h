//----------------------------------------------------
//  SpecialBoundaryDistance.h
//
//  Author:   Frankie Li
//  Purpose:  Calculate the distance between special boundary 
//            to a given point.
//            This is a class built as an Oracle.
//----------------------------------------------------


#ifndef BOUNDARY_ANALYSIS_H
#define BOUNDARY_ANALYSIS_H

//#include <CGAL/basic.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include "AABB_mesh_3_triangle_primitive.h"

#include <CGAL/AABB_triangle_primitive.h>

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
  //  BoundaryAnalysis
  //  Purpose:  Provide standard toolkit on calculating distances and shift
  //            between two different meshes
  //-------------------------------------------------------------
  template< class C3T3 >
  class SpecialBoundaryDistance
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
    typedef typename std::multimap<Subdomain_index, Subdomain_index> CrossMeshIDMap;  // maps from grain in mesh 1 to mesh 2
    typedef typename std::pair<Subdomain_index, Subdomain_index>     IDPair;

    typedef typename C3T3::Triangulation::Locate_type     Locate_type;
    typedef typename C3T3::Triangulation::Facet           Facet;
    typedef typename C3T3::Triangulation::Geom_traits::Ray_3       Ray3;
    typedef typename C3T3::Triangulation::Geom_traits::Triangle_3  Triangle_3;
    typedef typename C3T3::Triangulation                           Triangulation;
    typedef typename Kernel::Iso_cuboid_3                          Iso_cuboid;
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
    
    typedef typename CGAL::AABB_mesh_3_triangle_primitive<C3T3, AABB_Kernel, Filtered_facet_iterator>   AABB_Primitive;
    typedef typename CGAL::AABB_traits<AABB_Kernel, AABB_Primitive>                   AABB_Mesh_Triangle_Traits;
    typedef typename CGAL::AABB_tree<AABB_Mesh_Triangle_Traits>                       Mesh_AABB_Tree;
    
    typedef typename Mesh_AABB_Tree::Object_and_primitive_id                          Object_and_primitive_id;
    typedef typename Mesh_AABB_Tree::Primitive_id                                     Primitive_id;
    typedef typename Mesh_AABB_Tree::Point_and_primitive_id                           Point_and_primitive_id;
    typedef CGAL::Cartesian_converter< Tr_Kernel, AABB_Kernel>                        TrKernel_To_AABB;
  

    
  public:

    //----------
    //  Also need a sigma map
    //----------
    SpecialBoundaryDistance( const C3T3 & Mesh,
                             std::vector<SMatrix3x3> & IDOrientMap_  ):
      Mesh_( Mesh ),
      IDOrientMap( IDOrientMap_ )
    {}
    
    //--------------------------------------
    //  DefineSpecialBoundaries
    //--------------------------------------
    void DefineSpecialBoundaries()
    {
      // probably require a function object here.
    }
    
    //--------------------------------------
    //  DefineBoundaryAsSink
    //--------------------------------------
    void DefineBoundaryAsSink()
    {
      TotalSinkBoundaryArea = 0;
      //--------------------------------------
      //  Select boundaries that we want to include, define as sink
      //--------------------------------------
      const int EmptySpaceID = -1;
      for( Finite_facet_iterator pIter = Mesh_.triangulation().finite_facets_begin();
           pIter != Mesh_.triangulation().finite_facets_end(); ++ pIter )
      {
        Surface_index SurfaceID = Mesh_.surface_index( *pIter );
        
        if( SurfaceID.first != EmptySpaceID
            && SurfaceID.second != EmptySpaceID
            && SurfaceID.first != SurfaceID.second )
        {
          BoundaryTriangleList.push_back( pIter ); 
          TotalSinkBoundaryArea += CGAL::sqrt( Mesh_.triangulation().triangle( *pIter ).squared_area() );
        }
      }

      std::cout << "Number of boundary facets " << BoundaryTriangleList.size() << std::endl;
      BoundaryAABBTree.rebuild( BoundaryTriangleList.begin(),
                                BoundaryTriangleList.end() );
      
    }
    
    //--------------------------------------
    //  DistanceToNearestSink
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
      
      AABB_Point Closest_Point = BoundaryAABBTree.closest_point( AABB_SearchCenter );
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
    //  GetTotalBoundaryArea()
    //  Precondition:  Sink boundary must be already defined.
    //--------------------------------------
    double GetTotalSinkBoundaryArea() const
    {
      return TotalSinkBoundaryArea;
    }
    
    //--------------------------------------
    //  Precondition:  Grain Exists
    //--------------------------------------
    GeneralLib::SQuaternion GetGrainOrientation( int GrainID )
    {
      GeneralLib::SQuaternion q;
      if ( GrainID >= 0 ) {
         q.Set( IDOrientMap[ GrainID ] );
      }
      else {
         q.Set( 1.0, 0.0, 0.0, 0.0 ) ;
      }
      // DEBUG_ASSERT( GrainID >= 0, "Precondition failed: Requirement: GrainID >= 0 \n");
      return q;
    }
      
  private:
    
    //--------------------------------------
    //  GetCenterOfMass
    //--------------------------------------
    template< class PtrIter >
    SVector3 GetCenterOfMass( PtrIter pCur, PtrIter pEnd )
    {
      SVector3 CenterOfMass(0, 0,0);
      int nElements = pEnd - pCur;
      for( ; pCur != pEnd; ++ pCur )
        CenterOfMass += SVector3( pCur->x(), pCur->y(), pCur->z() );
        
      CenterOfMass /= static_cast<double>( nElements );
      return CenterOfMass;
    }
    
    //-------------------------------------------------
    //  Center
    //   Purpose:  Return the center of the triangle.
    //-------------------------------------------------
    Point3 Center( const Triangle_3 & T )
    {
      Point3 p( ( T[0].x() + T[1].x() + T[2].x() ) / 3.,
                ( T[0].y() + T[1].y() + T[2].y() ) / 3.,
                ( T[0].z() + T[1].z() + T[2].z() ) / 3. );
      return p;
    }
    
    //-------------------------------------------------
    //  PointToFacetDistance
    //   Purpose:  Calculate perpendicular distance between
    //             point and facet
    //-------------------------------------------------
    double PointToFacetDistance( const C3T3 & c3t3,
                                 const Point3 & p,
                                 const Cell_handle  &ch,
                                 int OppVertex )
    {
      Triangle_3 T = c3t3.triangulation().triangle( ch, OppVertex );
      Point3 PlanePoint = T.supporting_plane().projection( p );
      return  CGAL::sqrt( CGAL::squared_distance(PlanePoint, p ) );
    }
    
    //--------------------------------------
    //  IsGrainBoundary
    //--------------------------------------
    bool IsGrainBoundary( const Facet & f )
    {
      Surface_index ID = Mesh_.surface_index( f );
      return ( ID.first >= 0 && ID.second >= 0 && ID.first != ID.second );
    }


    //-------------------------------------------------
    //  PointToFacetDistance
    //   Purpose:  Calculate perpendicular distance between
    //             point and facet
    //-------------------------------------------------
    double PointToFacetDistance( const Point3 & p,
                                 const Facet & f )
    {
      Triangle_3 T = Mesh_.triangulation().triangle( f.first, f.second );
      Point3 PlanePoint = T.supporting_plane().projection( p );
      return  CGAL::sqrt( CGAL::squared_distance( PlanePoint, p ) );
    }
    
    
    //---------------
    //  Data
    //---------------
    
    std::vector<SMatrix3x3> & IDOrientMap;
    
    const C3T3 & Mesh_;
    
    Iso_cuboid MeshBBox;

    Filtered_finite_facet_list BoundaryTriangleList;
    Mesh_AABB_Tree BoundaryAABBTree;

    double TotalSinkBoundaryArea;
  };


}

#endif
