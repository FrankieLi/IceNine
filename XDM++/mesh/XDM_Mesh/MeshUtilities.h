//----------------------------------------------------
//  MeshUtilities.h
//
//  Author:   Frankie Li
//  Purpose:  Utilities for mesh calculations
//
//----------------------------------------------------

#ifndef MESH_UTILITIES_
#define MESH_UTILITIES_

#include <set>
#include <map>
#include <CGAL/IO/File_medit.h>
#include <CGAL/random_selection.h>
#include "3dMath.h"
#include "Sampling.h"
#include "Quaternion.h"

#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include "XDM_mesh_domain_3.h"
#include "XDM_Data.h"
#include "XDM_make_mesh.h"

namespace Pandora
{
  namespace Details
  {
    
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    template< class VectorT >
    VectorT UnitVector( const VectorT & v )
    {
      return VectorT( v.x(), v.y(), v.z(),
                      CGAL::sqrt( v.squared_length() ) );
    }
    
    //-----------------------------------------------------------------------------
    // IsTripleLineEdge
    //  - This is actually *EDGE LENGTH* instead of triple line.
    //
    //-----------------------------------------------------------------------------
    template < class C3T3, class Edge >
    bool IsTripleLineEdge(  const C3T3 & c3t3,
                            const Edge & e )
    {
      typedef typename C3T3::Triangulation::Cell_circulator  Cell_circulator;
      typedef typename C3T3::Subdomain_index Subdomain_index;
      typedef typename C3T3::Cell_handle Cell_handle;
      typedef typename C3T3::Vertex_handle Vertex_handle;
      
      std::set< Subdomain_index > DomainSet;
      Cell_circulator pEnd = c3t3.triangulation().incident_cells( e );
      Cell_circulator pCur = c3t3.triangulation().incident_cells( e );
      Cell_handle ch = e.get<0>();
      Vertex_handle v1 = ch->vertex( e.get<1>() );
      Vertex_handle v2 = ch->vertex( e.get<2>() );
      if( pCur != 0 )
      {
        do
        {
          int nMatch = 0;
          for( int i = 0; i < 4; i ++ )
            if( pCur->vertex(i) == v1 || pCur->vertex(i) == v2 )
              nMatch ++;
          if( nMatch == 2 )
            DomainSet.insert( c3t3.subdomain_index( pCur ) );
          ++ pCur;
        } while( pCur != pEnd );
      }
      return DomainSet.size() >= 3;
    }


    //--------------------------------
    //  GetCurvatureSign
    //   Convention is that all "surface normal" 
    //   must point from the RefID to
    //   the other ID
    //
    //   Return 1 if the sign of the surface normal is
    //   correct according to the convension.  -1 otherwise.
    //--------------------------------
    template< class Facet, class Point3, class Vector_3,
              class Subdomain_index, class C3T3 >
    int CorrectCurvatureSign( const Facet & f,
                              const Point3 & RefPoint,
                              const Vector_3 & K,    // Surface Normal
                              const Subdomain_index & RefID,
                              const C3T3 & Mesh )
    {
      typedef typename C3T3::Vertex_handle                           Vertex_handle;
      typedef typename C3T3::Triangulation::Geom_traits::Plane_3     Plane_3;
      typedef typename C3T3::Cell_handle                             Cell_handle;
      
      Plane_3 P = Mesh.
        triangulation().
        triangle( f ).
        supporting_plane();
            
      const Cell_handle & ch = f.first;
      Vertex_handle vh0 = ch->vertex( f.second );
      
      if( Mesh.subdomain_index( ch ) == RefID  )
      {
        if( P.has_on_positive_side( vh0->point() ) )  // vh0 is in interior of grain volu
          P = P.opposite();
      }
      else
      {
        if( P.has_on_negative_side( vh0->point() ) )  // vh0 is in interior of grain volu
          P = P.opposite();
      }
      
      Point3 TestPoint = RefPoint + K;

      if( P.has_on_negative_side( TestPoint ) )
        return -1;
      else
        return 1;
    }

    //-----------------------------------------------------------------------------
    //  IsSurfaceVertex
    //   Return true if the vertex handle vh points to a vertex on the surface of C3T3
    //-----------------------------------------------------------------------------
    template< class C3T3, class Vertex_handle >
    bool IsSurfaceVertex( const Vertex_handle & vh,
                          const C3T3 & Mesh )
    {
      typedef typename C3T3::Subdomain_index                         Subdomain_index;
      typedef typename C3T3::Cell_handle                             Cell_handle;
      std::set<Subdomain_index> IndexSet;
      
      std::vector<Cell_handle> CellList;
      Mesh.triangulation().finite_incident_cells( vh, std::back_inserter( CellList ));      
      for( int i = 0; i < CellList.size(); i ++  )
        IndexSet.insert( Mesh.subdomain_index( CellList[i] ) );
      return IndexSet.size() == 2;
    }
    
    //--------------------------------------
    //  ExtractBoundaryPoints
    //   Purpose:  Extract or vertices on grain boundaries
    //             and sample surface.
    //--------------------------------------
    template< class OutputIterator, class C3T3 >
    void ExtractBoundaryPoints( OutputIterator OutIter,
                                const C3T3 & c3t3 )
    {
      std::cout << "Extracting Boundary Points " << std::endl;
      typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertex_iterator;
      typedef typename C3T3::Triangulation::Cell_handle Cell_handle;
      //---------------------------------
      //  Need to change to select boundar vertices only
      //---------------------------------
      for( Finite_vertex_iterator pIter = c3t3.triangulation().finite_vertices_begin();
           pIter != c3t3.triangulation().finite_vertices_end(); ++ pIter )
      {
        if( c3t3.in_dimension( pIter) <= 2 )
        {
          std::vector<Cell_handle> IncidentCellList;
          c3t3.triangulation().incident_cells( pIter,
                                               std::back_inserter( IncidentCellList ) );
          bool bIsSpaceBoundary = false;
          for( int i = 0; i < IncidentCellList.size(); i ++ )
          {
            if( c3t3.subdomain_index( IncidentCellList[i] ) <= 0 )
              bIsSpaceBoundary = true;
          }
          
          if( !bIsSpaceBoundary )
          {
            *OutIter =  pIter->point();
            ++ OutIter;
          }
        }
      }
    }
    
    //--------------------------------------
    //  ExtractSurfacePoints
    //    Parameter:  c3t3    - mesh
    //                dZLimit - the offset from the top
    //                          and the bottom that delineates
    //                          the section of volume used for
    //                          registeration.
    //    Test:  Output and visualize
    //--------------------------------------
    template< class OutputIterator, class C3T3, class Subdomain_index >
    void ExtractSurfacePoints( OutputIterator OutIter,
                               const C3T3 & c3t3,
                               const Subdomain_index & EmptySpaceID )
    {
      std::cout << "Extracting Surface Points " << std::endl;
      typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertex_iterator;
      typedef typename C3T3::Triangulation::Cell_handle Cell_handle;
      //---------------------------------
      //  Need to change to select boundar vertices only
      //---------------------------------
      for( Finite_vertex_iterator pIter = c3t3.triangulation().finite_vertices_begin();
           pIter != c3t3.triangulation().finite_vertices_end(); ++ pIter )
      {
        if( c3t3.in_dimension( pIter) <= 2 )
        {
          std::vector<Cell_handle> IncidentCellList;
          c3t3.triangulation().incident_cells( pIter,
                                               std::back_inserter( IncidentCellList ) );
          bool bAdjacentSample = false;
          bool bAdjacentSpace  = false;
          for( int i = 0; i < IncidentCellList.size(); i ++ )
          {
            if( c3t3.subdomain_index( IncidentCellList[i] )
                   == EmptySpaceID
                || c3t3.subdomain_index( IncidentCellList[i] ) == 0 )  // this is a hack because
              bAdjacentSpace = true;                                   // of the inconsistency
                                                                       // in empty space declaration
            if( c3t3.subdomain_index( IncidentCellList[i] ) > 0  )
              bAdjacentSample = true;
            
          }
          if( bAdjacentSample && bAdjacentSpace )
          {
            *OutIter =  pIter->point();
            ++ OutIter;
          }
        }
      }
    }

    //--------------------------------------
    //  GetSurfaceGrainIDs
    //    Parameter:  c3t3    - mesh
    //                EmptySpaceID - the ID that indicates empty space
    //                OutIter - best to use a set
    //
    //
    //--------------------------------------
    template< class IDSet, class C3T3 >
    void GetSurfaceGrainIDs( IDSet & SurfaceIDSet,
                             const C3T3 & c3t3,
                             const typename C3T3::Subdomain_index & EmptySpaceID )
    {
      
      typedef typename C3T3::Facet_iterator Facet_iterator;
      typedef typename C3T3::Surface_index  Surface_index;
      for( Facet_iterator pCur = c3t3.facets_begin();
           pCur != c3t3.facets_end(); ++pCur )
      {
        Surface_index facet_index = c3t3.surface_index( *pCur );
        if( facet_index.first == EmptySpaceID || facet_index.first == 0  )
          SurfaceIDSet.insert( facet_index.second );
        else if( facet_index.second == EmptySpaceID || facet_index.second == 0 )
          SurfaceIDSet.insert( facet_index.first );
        
      }
    }
    
    //-----------------------------------------------------------------------------
    // IsMultidomainEdge
    //  - This is actually *EDGE LENGTH* instead of triple line.
    //
    //-----------------------------------------------------------------------------
    /*    template < class C3T3, class Edge, class Subdomain_index >
    bool IsMultidomainEdge( const C3T3 & c3t3,
                            const Edge & e,
                            const Subdomain_index & ID )
    {
      typedef typename C3T3::Triangulation::Cell_circulator  Cell_circulator;
      typedef typename C3T3::Cell_handle Cell_handle;
      typedef typename C3T3::Vertex_handle Vertex_handle;
      
      std::set< Subdomain_index > DomainSet;
      Cell_circulator pEnd = c3t3.triangulation().incident_cells( e );
      Cell_circulator pCur = c3t3.triangulation().incident_cells( e );
      Cell_handle ch = e.get<0>();
      Vertex_handle v1 = ch->vertex( e.get<1>() );
      Vertex_handle v2 = ch->vertex( e.get<2>() );
      int nIncidents = 0;   // debug
      int nDoubleMatch = 0; // debug
      if( pCur != 0 )
      {
        do
        {
          nIncidents ++;
          int nMatch = 0;
          for( int i = 0; i < 4; i ++ )
            if( pCur->vertex(i) == v1 || pCur->vertex(i) == v2 )
              nMatch ++;
          if( nMatch == 2 )
          {
            nDoubleMatch ++;
          }
          DomainSet.insert( c3t3.subdomain_index( pCur ) );
          ++ pCur;
        } while( pCur != pEnd );
      }
      if( nIncidents != nDoubleMatch )      
      {
        std::cout << "Incident is not what you think" << std::endl;
        exit(0);

      }
      if( DomainSet.find( ID ) == DomainSet.end() )
      {
        std::cout << "Multimap is not right" << std::endl;
        exit(0);
      }
      return ( DomainSet.size() >= 3 && DomainSet.find( ID ) != DomainSet.end() ) ;
    }
    */
    //-----------------------------------
    // BuildGrainTetMap
    //-----------------------------------
    template< class GrainIDToTetMapT, class C3T3 >
    void BuildGrainTetMap( GrainIDToTetMapT & GrainIDToTetMap,
                           const C3T3 & c3t3_ )
    {
      typedef typename C3T3::Subdomain_index Subdomain_index;
      typedef typename C3T3::Cell_iterator   Cell_iterator;
      typedef typename C3T3::Cell_handle     Cell_handle;

      for( Cell_iterator pCur = c3t3_.cells_begin();
           pCur != c3t3_.cells_end(); ++pCur )
      {
        Subdomain_index nID = c3t3_.subdomain_index( pCur );
        Cell_handle ch = pCur;
        GrainIDToTetMap.insert( std::make_pair( nID, ch ) );
      }
    }
      
    //-----------------------------------
    // BuildGrainFacetMap
    //-----------------------------------
    template< class GrainIDToFacetMapT, class C3T3 >
    void BuildGrainFacetMap( GrainIDToFacetMapT & GrainIDToFacetMap,
                             const C3T3 & c3t3_ )
    {
      typedef typename C3T3::Facet_iterator  Facet_iterator;
      typedef typename C3T3::Surface_index   Surface_index;
      for( Facet_iterator pCur = c3t3_.facets_begin();
           pCur != c3t3_.facets_end(); ++pCur )
      {
        Surface_index facet_index = c3t3_.surface_index( *pCur );
        GrainIDToFacetMap.insert( std::make_pair( facet_index.first,  *pCur ) );
        GrainIDToFacetMap.insert( std::make_pair( facet_index.second, *pCur ) );
      }
    }
   
    template < class InputT, class ReturnT >
    struct TrivialReturn
    {
      ReturnT operator()( InputT o )
      { return o; }
    };
    
    template < class InputT, class ReturnT >
    struct IterToSecond
    {
      ReturnT operator()( InputT o )
      { return o->second; }
    };

    //-----------------------------------
    // BuildFiniteGrainNgbMap
    //
    //  Given GrainTetMap
    //        For each GrainID
    //           For each Tet
    //               Find ngb
    //                   Add ngb ID -> set
    //  std::map< ID, std::set< ID > >
    //-----------------------------------
    template< class C3T3 >
    void BuildFiniteGrainNgbMap( std::map< int, std::set< int > > & ID2NgbMap,
                                 const C3T3 & c3t3_, int nMinID = 1 )
    {
      typedef std::map< int, std::set<int> > IDNgbMap;
      typedef IDNgbMap::iterator NgbIter;
      typedef typename C3T3::Cell_iterator   Cell_iterator;
      typedef typename C3T3::Cell_handle     Cell_handle;
      typedef typename C3T3::Surface_index   Surface_index;

      for( Cell_iterator pCur = c3t3_.cells_begin();
           pCur != c3t3_.cells_end(); ++ pCur )
      {
        Cell_handle ch = pCur;
        int nID = c3t3_.subdomain_index( ch );
        if( nID >= nMinID )
        {
          for ( int i = 0; i < 4; i ++ )
          {
            int nNgbID = c3t3_.subdomain_index( ch->neighbor( i ));
            if( nNgbID >= nMinID )
            {
              ID2NgbMap[ nID ].insert( nNgbID  );
              ID2NgbMap[ nNgbID ].insert( nID  );
            }
          }
        }
      }
    }
    
    //-----------------------------------
    // BuildGrainNgbMap
    //
    //  Given GrainTetMap
    //        For each GrainID
    //           For each Tet
    //               Find ngb
    //                   Add ngb ID -> set
    //  std::map< ID, std::set< ID > >
    //-----------------------------------
    template< class C3T3 >
    void BuildGrainNgbMap( std::map< int, std::set< int > > & ID2NgbMap,
                           const C3T3 & c3t3_ )
    {
      typedef std::map< int, std::set<int> > IDNgbMap;
      typedef IDNgbMap::iterator NgbIter;
      typedef typename C3T3::Cell_iterator   Cell_iterator;
      typedef typename C3T3::Cell_handle     Cell_handle;
      typedef typename C3T3::Surface_index   Surface_index;

      for( Cell_iterator pCur = c3t3_.cells_begin();
           pCur != c3t3_.cells_end(); ++ pCur )
      {
        Cell_handle ch = pCur;
        int nID = c3t3_.subdomain_index( ch );
        if( nID > 0 )
        {
          for ( int i = 0; i < 4; i ++ )
          {
            int nNgbID = c3t3_.subdomain_index( ch->neighbor( i ));
            ID2NgbMap[ nID ].insert( nNgbID  );
            ID2NgbMap[ nNgbID ].insert( nID  );
          }
        }
      }
    }
    
    
    //-------------------------------------------------
    //  BuildGrainEdgeMap
    //
    //  ID2EdgeMap needs to be a multi-map
    //-------------------------------------------------
    template< class C3T3, class ID2EdgeMap >
    void BuildGrainSurfaceEdgeMap( ID2EdgeMap & Map,
                                   const C3T3 & c3t3_ )
    {
      typedef std::map< int, std::set<int> > IDNgbMap;
      typedef IDNgbMap::iterator NgbIter;
      typedef typename C3T3::Cell_handle     Cell_handle;
      typedef typename C3T3::Triangulation::Edge_iterator Edge_iterator;
      typedef typename C3T3::Triangulation::Facet_circulator  Facet_circulator;
      typedef typename C3T3::Subdomain_index Subdomain_index;
      typedef typename C3T3::Surface_index Surface_index;

      for( Edge_iterator pEdge = c3t3_.triangulation().edges_begin();
           pEdge != c3t3_.triangulation().edges_end(); ++ pEdge )
      {
        Facet_circulator pEnd = c3t3_.triangulation().incident_facets( *pEdge );
        Facet_circulator pCur = c3t3_.triangulation().incident_facets( *pEdge );
        if( ! c3t3_.triangulation().is_infinite( *pEdge ) )
        {
          std::set< Subdomain_index > DomainSet;
          if( pCur != 0 ) // construct of using Cell_circulator
          {
            do
            {
              Surface_index SurfID = c3t3_.surface_index( *pCur );
              DomainSet.insert( SurfID.first );
              DomainSet.insert( SurfID.second );
              ++ pCur;
            } while( pCur != pEnd );
          }
          
          if( DomainSet.size() >= 2 )
          {
            for( typename std::set< Subdomain_index >::iterator pIter = DomainSet.begin();
                 pIter != DomainSet.end(); ++ pIter )
            {
              Map.insert( std::make_pair( *pIter, *pEdge ) );
            }
          }
        }// end infinte edge
      }// end for
    }
    
    //-------------------------------------------------
    //  TransformPoint
    //-------------------------------------------------
    template< class Point3, class SMatrix3x3, class SVector3 >
    Point3 TransformPoint(  const Point3  & p,
                            const SMatrix3x3 & R,
                            const SVector3 & T )
    {
      SVector3 Temp( static_cast<float>( p.x() ),
                     static_cast<float>( p.y() ),
                     static_cast<float>( p.z() ) );
      Temp = (R * Temp) + T; 
      Point3 Result( Temp.m_fX, Temp.m_fY, Temp.m_fZ ) ;
      return Result;
    }

    //-------------------------------------------------
    //  InverseTransformPoint
    //-------------------------------------------------
    template< class Point3, class SMatrix3x3, class SVector3 >
    Point3 InverseTransformPoint(  const Point3  & p,
                                   const SMatrix3x3 & R,
                                   const SVector3 & T )
    {
      SVector3 Temp( static_cast<float>( p.x() ),
                     static_cast<float>( p.y() ),
                     static_cast<float>( p.z() ) );
      SMatrix3x3 R_T = R;
      R_T.Transpose();
      Temp = R_T *  (Temp - T); 
      Point3 Result( Temp.m_fX, Temp.m_fY, Temp.m_fZ ) ;
      return Result;
    }

    
    //-----------------------------------------------------------------------------
    // CalculateCircularOrientationError
    //
    //-----------------------------------------------------------------------------
    template< class GrainPropMap, class EdgeErrorMap,
              class C3T3, class TripleLineMapper,
              class ID2OrientMapT   >
    void CalculateCircularOrientationError( GrainPropMap &        IDToPropMap,
                                            EdgeErrorMap &        EdgeToErrorMap,
                                            const C3T3 &          c3t3_,
                                            TripleLineMapper      Put,
                                            const ID2OrientMapT & IDToOrientMap )
    {
      typedef typename C3T3::Subdomain_index                 Subdomain_index;
      typedef typename std::set<Subdomain_index>::iterator   SetIter;
      typedef typename C3T3::Triangulation::Cell_circulator  Cell_circulator;
      typedef typename C3T3::Triangulation::Finite_edges_iterator    Finite_edge_iterator;
      typedef typename C3T3::Triangulation::Cell_circulator  Cell_circulator;
      typedef typename GeneralLib::SQuaternion SQuaternion;

      std::ofstream DebugStream( "DebugOrient.txt" );
      int nEdgeIndex = 0;
      for( Finite_edge_iterator pEdge = c3t3_.triangulation().finite_edges_begin();
           pEdge != c3t3_.triangulation().finite_edges_end(); ++ pEdge )
      {
        Cell_circulator pCur = c3t3_.triangulation().incident_cells( *pEdge );
        Cell_circulator pEnd = c3t3_.triangulation().incident_cells( *pEdge );
        
        bool bEdgeIncidentsSpace = false;
        typename std::set<Subdomain_index> DomainSet;
        SQuaternion qTot;
        qTot.Set( 1, 0, 0, 0 ); // set to identity
        
        if( pCur != 0 )
        {
          do
          {
            Cell_circulator pNext = pCur; ++ pNext;
            Subdomain_index IDa = c3t3_.subdomain_index( pCur );
            DomainSet.insert( IDa );
            Subdomain_index IDb = c3t3_.subdomain_index( pNext );
            DomainSet.insert( IDb );

            if( IDa > 0 && IDb > 0 )  // neither are space or empty
            {
              SQuaternion qA, qB;
              if( IDToOrientMap.find( IDa ) != IDToOrientMap.end() )
                qA = IDToOrientMap.find( IDa )->second;
              if( IDToOrientMap.find( IDb ) != IDToOrientMap.end() )
                qB = IDToOrientMap.find( IDb )->second;
              
              //              qTot = qTot * ( qA.Inverse() * qB );
              qTot = qA * qTot;
              DebugStream << nEdgeIndex << " "
                //                          << IDa << " " << IDb << " "
                          << IDa << " "
                          << qTot << " " << qA << std::endl;
            }
            else
            {
              bEdgeIncidentsSpace = true;   // qTot is no longer useful
            }
            ++ pCur;
          }while( pCur != pEnd );
        }  // if circulator is valid
        
        for( typename std::set<Subdomain_index>::iterator pID = DomainSet.begin();
             pID != DomainSet.end(); ++ pID )
          Put( IDToPropMap, *pID, qTot );
        Put( EdgeToErrorMap, pEdge, qTot, bEdgeIncidentsSpace, DomainSet );

        nEdgeIndex ++;
      }
     
    }
    
  }  // end Details
  
  //-------------------------------------------------
  //  BuildIDToGrainAverageMap
  //-------------------------------------------------
  template< class ID2OrientMapT, class CellOrientAverageMap,
            class C3T3 >
  void BuildIDToGrainAverageMap( ID2OrientMapT & ID2OrientMap,
                                 const CellOrientAverageMap & CellOrientMap,
                                 const C3T3 & c3t3_ )
  {
    typedef typename ID2OrientMapT::iterator ID2OrientMapIter;
    typedef typename CellOrientAverageMap::const_iterator Iter;
    typedef typename C3T3::Subdomain_index Subdomain_index;
    typedef typename C3T3::Cell_handle     Cell_handle;
    // CGAL::Default_cell_index_pmap<C3T3> oCellToID( c3t3_ );
    for( Iter pCur = CellOrientMap.begin();
         pCur != CellOrientMap.end(); ++ pCur )
    {
      Cell_handle      ch = pCur->first;
      Subdomain_index nID = c3t3_.subdomain_index( ch );
      ID2OrientMapIter pFound = ID2OrientMap.find( nID );
      if( pFound == ID2OrientMap.end() )
        ID2OrientMap[ nID ] = pCur->second;
    }
  }

  //-------------------------------------------------
  //  TransformMesh
  //-------------------------------------------------
  template <class C3T3>
  void TransformMesh( C3T3 & c3t3,
                      const GeneralLib::SMatrix3x3 & Rotation,
                      const GeneralLib::SVector3   & Translation )
  {
    typedef typename C3T3::Triangulation Tr;
    typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename Tr::Geom_traits::Bare_point Point3;
    for(  Finite_vertices_iterator vit
            = c3t3.triangulation().finite_vertices_begin();
          vit != c3t3.triangulation().finite_vertices_end();
          ++ vit )
    {
      Point3 p = vit->point();
      p = Details::TransformPoint( p, Rotation, Translation );
      vit->set_point( p );
    }
  }
  
} // namespace Pandora


#endif
