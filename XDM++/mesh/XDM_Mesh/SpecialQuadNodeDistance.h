//----------------------------------------------------
//  SpecialTripleLineDistance
//
//  Author:   Frankie Li
//  Purpose:  Calculate the distance between special triple line
//            to a given point.
//            This is a class built as an Oracle.
//----------------------------------------------------


#ifndef QUAD_NODE_DISTANCE_H
#define QUAD_NODE_DISTANCE_H

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>


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
  //  SpecialQuadNodeDistance
  //  Purpose:  Calculate the distance between a spetial triple line and
  //            a given point.
  //-------------------------------------------------------------
  template< class C3T3 >
  class SpecialQuadNodeDistance
  {
  public:
        
    struct Kernel: public CGAL::Exact_predicates_inexact_constructions_kernel {};
    typedef typename C3T3::Triangulation::Vertex_iterator            Vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_vertices_iterator   Finite_vertices_iterator;
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

    
    typedef typename C3T3::Triangulation::Edge                     Edge;
    typedef typename C3T3::Triangulation::Edge_iterator            Edge_iterator;
    
    typedef typename C3T3::Triangulation::Facet_circulator   Facet_circulator;

    
    typedef GeneralLib::SQuaternion                                SQuaternion;

    
    //-------------------------------------------------------
    //  Random sampling
    //-------------------------------------------------------
    typedef GeneralLib::SMatrix3x3 SMatrix3x3;
    typedef GeneralLib::SVector3   SVector3;
    typedef GeneralLib::UniformGrid::CQuaternionGrid CQuaternionGrid;

    mutable GeneralLib::CUniformRandomReal RandGen;
    
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;
    typedef typename Neighbor_search::Tree Tree;

    typedef typename C3T3::Triangulation::Geom_traits Tr_Kernel;

    
  public:
    
    //----------
    //  Also need a sigma map
    //----------
    SpecialQuadNodeDistance( const C3T3 & Mesh,
                             std::vector<SMatrix3x3> & IDOrientMap_  ):
      Mesh_( Mesh ),
      IDOrientMap( IDOrientMap_ )
    {}
    
    
    //--------------------------------------
    //  DefineTripleJunctionAsSink
    //--------------------------------------
    void DefineQuadNodesAsSink()
    {
      QuadNodeList.clear();
      for( Finite_vertices_iterator pVertex = Mesh_.triangulation().finite_vertices_begin();
           pVertex != Mesh_.triangulation().finite_vertices_end(); ++ pVertex )
      {
        if( JunctionType( pVertex, IDOrientMap.size() ) >= 4 )
          QuadNodeList.push_back( pVertex->point().point() );
      }
      
      NumNJunctions = QuadNodeList.size();
      NJunctionTree.insert( QuadNodeList.begin(), QuadNodeList.end() );
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
      Neighbor_search search( NJunctionTree, SearchCenter, 1 );
      if( search.begin() == search.end() )
        return -4;
      
      return CGAL::sqrt( search.begin()->second );
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
    //  GetNumJunctions
    //--------------------------------------
    int GetNumJunctions() const
    {
      return NumNJunctions;
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


    //----------------------------------------------------------------
    //  IsValidEdge
    //
    //----------------------------------------------------------------
    int IsValidEdge( const Edge & e, int MaxID )
    {  
      if( ! Mesh_.triangulation().is_infinite( e ) )
      {
        std::set< Subdomain_index > DomainSet;
        Facet_circulator pEnd = Mesh_.triangulation().incident_facets( e );
        Facet_circulator pCur = Mesh_.triangulation().incident_facets( e );
        bool bException = false;
        if( pCur != 0 ) // construct of using Cell_circulator
        {
          do
          {
            if( Mesh_.is_in_complex( *pCur ) )
            {   
              Surface_index SurfID = Mesh_.surface_index( *pCur );
              DomainSet.insert( SurfID.first );
              DomainSet.insert( SurfID.second );
          
              if( SurfID.first  <= 0 || SurfID.first  >= MaxID )
                bException = true;
              if( SurfID.second <= 0 || SurfID.second >= MaxID )
                bException = true;
            }
            ++ pCur;
          } while( pCur != pEnd );
        }
        if( ! bException &&  DomainSet.size() > 2 )   // not surface
        {
          return DomainSet.size();
        }
      }
      return -1;
  
    }


    //----------------------------------------------------------------
    //  JunctionType
    //
    //----------------------------------------------------------------
    int JunctionType( Vertex_handle vh, int MaxID )
    {
      int NumValidEdges = 0;
      vector<Edge>              NgbEdgeList;
      std::set<Facet>           NgbFacetSet;
      std::set<Subdomain_index> DomainSet;
      Mesh_.triangulation().incident_edges( vh, std::back_inserter( NgbEdgeList ) );


      bool ValidJunction = true;
      for( int i = 0; i < NgbEdgeList.size(); i ++ )
      {
        std::set< Subdomain_index > LocalDomainSet;
        Facet_circulator pEnd = Mesh_.triangulation().incident_facets( NgbEdgeList[i] );
        Facet_circulator pCur = Mesh_.triangulation().incident_facets( NgbEdgeList[i] );
        bool bException = false;
        if( pCur != 0 ) // construct of using Cell_circulator
        {
          do
          {
            if( Mesh_.is_in_complex( *pCur ) )
            {
              Surface_index SurfID = Mesh_.surface_index( *pCur );
              LocalDomainSet.insert( SurfID.first );
              LocalDomainSet.insert( SurfID.second );
          
              if( SurfID.first  <= 0 || SurfID.first  >= MaxID )
                bException = true;
              if( SurfID.second <= 0 || SurfID.second >= MaxID )
                bException = true;
            }
            ++ pCur;
          } while( pCur != pEnd );
        }

        if( bException )
          ValidJunction = false;
        if( !bException && LocalDomainSet.size() > 2 )
        {
          DomainSet.insert( LocalDomainSet.begin(), LocalDomainSet.end() );
          NumValidEdges ++;
        }
      }
      if( !ValidJunction )
        return -1;
  
      if( DomainSet.size() > 3 )
        return NumValidEdges;
      else                       // exception
        return -1;
    }
    
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

    std::list< Point3 > QuadNodeList;
    Tree    NJunctionTree;
    int     NumNJunctions;
  };


}

#endif
