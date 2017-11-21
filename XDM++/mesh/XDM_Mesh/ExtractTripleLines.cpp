#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "XDM_mesh_triangulation_3.h"
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include "XDM_mesh_criteria_3.h"

#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include "XDM_mesh_domain_3.h"
#include "XDM_Data.h"
#include "XDM_make_mesh.h"
#include <fstream>
#include <string>
#include <vector>
#include <sstream>


#include "GrainCurvatureEst.h"
#include "XDM_Mesh_IO.h"
#include "3dMath.h"
#include "SimpleMeshVTKWriter.h"
#include "VTKPropertyMap.h"
#include "GrainCurvatureEst.h"
#include "XDM_Mesh_IO.h"
#include "ProcessMic.h"
#include "Symmetry.h"
#include <boost/shared_ptr.hpp>
#include "MeshUtilities.h"

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

// Mesh Criteria
typedef unsigned char XDMIndexT;
using std::vector;
typedef vector< vector< vector< XDMIndexT > > >  DataMatrixT;

typedef CGAL::XDM_test::XDM_Data<K> XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;

typedef CGAL::Mesh_criteria_3<XDM_Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

typedef std::map< int, int > ID2IDMap;
typedef std::vector< GeneralLib::SMatrix3x3 > IDOrientMap;


typedef XDM_C3t3::Subdomain_index                         Subdomain_index;
typedef XDM_C3t3::Surface_index                           Surface_index;
typedef XDM_C3t3::Cell_handle                             Cell_handle;
typedef XDM_C3t3::Triangulation::Vertex_handle            Vertex_handle;
typedef XDM_C3t3::Cell_iterator                           Cell_iterator;
typedef XDM_C3t3::Triangulation::Edge_iterator            Edge_iterator;
typedef XDM_C3t3::Triangulation::Facet_circulator         Facet_circulator;
typedef XDM_C3t3::Triangulation::Geom_traits::Triangle_3  Triangle_3;
typedef std::map< Vertex_handle, int >  VertexIDMapT;

typedef XDM_C3t3::Triangulation::Finite_vertices_iterator  Finite_vertices_iterator;
typedef XDM_C3t3::Triangulation::Edge                      Edge;
typedef XDM_C3t3::Triangulation::Cell                      Cell;
typedef XDM_C3t3::Triangulation::Facet                     Facet;




//----------------------------------------------------------------
//  IsValidEdge
//
//----------------------------------------------------------------
int IsValidEdge( const  XDM_C3t3 & c3t3, const Edge & e, int MaxID )
{  
  if( ! c3t3.triangulation().is_infinite( e ) )
  {
    std::set< Subdomain_index > DomainSet;
    Facet_circulator pEnd = c3t3.triangulation().incident_facets( e );
    Facet_circulator pCur = c3t3.triangulation().incident_facets( e );
    bool bException = false;
    if( pCur != 0 ) // construct of using Cell_circulator
    {
      do
      {
        if( c3t3.is_in_complex( *pCur ) )
        {
          Surface_index SurfID = c3t3.surface_index( *pCur );
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
int JunctionType( const XDM_C3t3 & c3t3, Vertex_handle vh, int MaxID )
{
  int NumValidEdges = 0;
  vector<Edge>              NgbEdgeList;
  std::set<Facet>           NgbFacetSet;
  std::set<Subdomain_index> DomainSet;
  c3t3.triangulation().incident_edges( vh, std::back_inserter( NgbEdgeList ) );
  for( int i = 0; i < NgbEdgeList.size(); i ++ )
  {
    std::set< Subdomain_index > LocalDomainSet;
    Facet_circulator pEnd = c3t3.triangulation().incident_facets( NgbEdgeList[i] );
    Facet_circulator pCur = c3t3.triangulation().incident_facets( NgbEdgeList[i] );
    bool bException = false;
    if( pCur != 0 ) // construct of using Cell_circulator
    {
      do
      {
        if( c3t3.is_in_complex( *pCur ) )
        {
          Surface_index SurfID = c3t3.surface_index( *pCur );
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
    
    if( !bException && LocalDomainSet.size() > 2 )
    {
      DomainSet.insert( LocalDomainSet.begin(), LocalDomainSet.end() );
      NumValidEdges ++;
    }
  }
  
  if( DomainSet.size() > 3 )
    return NumValidEdges;
  else                       // exception
    return -1;
}

struct EdgeInfo
{
  Vertex_handle  v1;
  Vertex_handle  v2;
  int nEdge;
  double fLength;
  int nComplexFacet;
};


//-----------------------------------------------------------
template< class Circulator >
int GetNumFacets( const XDM_C3t3 & c3t3, Circulator pCur, Circulator pStop )
{
  if( pCur == 0 )
    return 0;
  
  int nFacets = 0;
  do
  {
    if( c3t3.is_in_complex( *pCur ) )
    {
      nFacets ++;
    }
    ++ pCur;
  }while ( pCur != pStop );

  return nFacets;
}


//-----------------------------------------------------------
void WriteTripleLineVtk( const std::string & oFilename,
                         VertexIDMapT & IDMap,
                         const vector<EdgeInfo> & TripleLines )
{
  std::ofstream os( oFilename.c_str( ) );
  
  os << "# vtk DataFile Version 2.0 \n"
     << "Volume Test \n"
     << "ASCII \n"
     << "DATASET UNSTRUCTURED_GRID \n"
     << "POINTS " << IDMap.size() << " float \n";
  
  typedef std::map<Vertex_handle, int>::const_iterator MapIter;
  vector< Vertex_handle > VertexList( IDMap.size() );
  for( MapIter pCur = IDMap.begin(); pCur != IDMap.end(); ++ pCur )
    VertexList[ pCur->second ] = pCur->first;

  for( int i = 0; i < VertexList.size(); i ++ )
  {
    os << VertexList[i]->point().x() << " "
       << VertexList[i]->point().y() << " "
       << VertexList[i]->point().z() << " "
       << std::endl;
  }
  
  os << "CELLS " << TripleLines.size() << " " << TripleLines.size() * 3 << std::endl;

  for( int i = 0; i < TripleLines.size(); i ++ )
    os << "2 " <<  IDMap[ TripleLines[i].v1 ] << " "
       << IDMap[ TripleLines[i].v2 ] << std::endl;
  
  os << "CELL_TYPES " << TripleLines.size() << std::endl;
  for( int i = 0; i < TripleLines.size(); i ++ )
    os << 3 << std::endl;


  //--------------------------------------------------
  os << "CELL_DATA " << TripleLines.size() << std::endl;
  os << "SCALARS nConnect int 1" << std::endl;
  os << "LOOKUP_TABLE default" << std::endl;
  for( int i = 0; i < TripleLines.size(); i ++ )
    os << TripleLines[i].nEdge << std::endl;
  
  os << "SCALARS nFacet int 1" << std::endl;
  os << "LOOKUP_TABLE default" << std::endl;
  for( int i = 0; i < TripleLines.size(); i ++ )
    os << TripleLines[i].nComplexFacet << std::endl;
  os << "SCALARS Length float 1" << std::endl;
  os << "LOOKUP_TABLE default" << std::endl;
  for( int i = 0; i < TripleLines.size(); i ++ )
    os << TripleLines[i].fLength << std::endl;

  os.close();
}

//----------------------------------------
//  ReadIDOrientMapx
//----------------------------------------
void ReadIDOrientMap( IDOrientMap & Map, const std::string & filename )
{
  std::ifstream MapInput( filename.c_str() );
  while( MapInput.good() )
  {
    GeneralLib::SVector3 EulerDegree;
    
    MapInput >> std::skipws >> EulerDegree.m_fX
             >> std::skipws >> EulerDegree.m_fY
             >> std::skipws >> EulerDegree.m_fZ;
    EulerDegree = GeneralLib::DegreeToRadian( EulerDegree );
    GeneralLib::SMatrix3x3 R;
    
    R.BuildActiveEulerMatrix( EulerDegree.m_fX,
                              EulerDegree.m_fY,
                              EulerDegree.m_fZ );
    Map.push_back(R);
  }
  std::cout << Map.size() << " number of lines read from euler angle file " << std::endl;
  MapInput.close();
}




int main()
{
  std::string MeshFilename;
  std::cout << "Enter binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename;
  
  //-------------------------------------------------------------------
  //      Start computing and meshing
  //-------------------------------------------------------------------
  std::cout << "Start reading mesh  " << std::endl;
  std::ifstream MeshIs( MeshFilename.c_str(),
                        std::ios_base::in|std::ios_base::binary);
  if( ! MeshIs.is_open() )
  {
    std::cout << "Input mesh not found " << std::endl;
    exit( 0 );
  }
  CGAL::set_binary_mode( MeshIs );
  XDM_C3t3 c3t3;
  MeshIs >> c3t3;
  MeshIs.close();

  IDOrientMap OrientationMap;

  std::string OrientMapFilename;
  std::cout << "Enter first ID to orientation map filename: " << std::endl;
  std::cin  >> OrientMapFilename;
  ReadIDOrientMap( OrientationMap, OrientMapFilename );
  
  std::cout << " Done reading binary mesh " << std::endl;
  std::cout << " Extracting triple line information " << std::endl;

  std::cout << " Triple line output file:  " << std::endl;
  std::string sTripleFilename;
  std::cin >> std::skipws >> sTripleFilename;
  std::ofstream TripleLineFile( sTripleFilename.c_str() );

  
  std::cout << " N-Vertex output file:  " << std::endl;
  std::string sNVertexFilename;
  std::cin >> std::skipws >> sNVertexFilename;
  std::ofstream NVertexFile( sNVertexFilename.c_str() );
  
  
  std::cout << " Triple line vtk output file:  " << std::endl;
  std::string sVtkFilename;
  std::cin >> std::skipws >> sVtkFilename;
  
  
  // ------------------------------------------
  //  Write N-edge
  // ------------------------------------------
  vector<EdgeInfo> EdgeList;
  VertexIDMapT     VertexToIDMap;
  int nVertexID = 0;
  int nEdgeID = 0;
  for( Edge_iterator pEdge = c3t3.triangulation().edges_begin();
       pEdge != c3t3.triangulation().edges_end(); ++ pEdge )
  {
    if( ! c3t3.triangulation().is_infinite( *pEdge ) )
    {
      int NumGrainNgbs = IsValidEdge( c3t3, *pEdge, OrientationMap.size() );      
      if( NumGrainNgbs > 2 )   // not surface
      {
        Vertex_handle v1 = pEdge->get<0>()->vertex( pEdge->get<1>() ) ;
        Vertex_handle v2 = pEdge->get<0>()->vertex( pEdge->get<2>() ) ;

        if( VertexToIDMap.find( v1 ) == VertexToIDMap.end() )   // this is for printing
        {
          VertexToIDMap[v1] = nVertexID;
          nVertexID ++;
        }
        if( VertexToIDMap.find( v2 ) == VertexToIDMap.end() )
        {
          VertexToIDMap[v2] = nVertexID;
          nVertexID ++;
        }        
        Facet_circulator pEnd = c3t3.triangulation().incident_facets( *pEdge );
        Facet_circulator pCur = c3t3.triangulation().incident_facets( *pEdge );

        int nNumFacets = GetNumFacets( c3t3, pCur, pEnd );
        int nComplexFacet = 0;
        do
        {
          if( c3t3.is_in_complex( *pCur ) )
          {            
            Surface_index SurfID = c3t3.surface_index( *pCur );
            SVector3 e1 = GeneralLib::RadianToDegree( OrientationMap[SurfID.first ].GetEulerAngles() );
            SVector3 e2 = GeneralLib::RadianToDegree( OrientationMap[SurfID.second].GetEulerAngles() );
            
            Triangle_3 T = c3t3.triangulation().triangle( pCur->first, pCur->second ); 
            TripleLineFile << nEdgeID << " " << NumGrainNgbs << " "
                           << T.supporting_plane().orthogonal_direction() << " "
                           << e1 << " "
                           << e2 << " "
                           << c3t3.triangulation().segment( *pEdge ) << " " 
                           << CGAL::sqrt( T.squared_area() ) << " "
                           << T.is_degenerate() << " "
                           << SurfID.first << " " << SurfID.second << " " << nNumFacets << std::endl;
            nComplexFacet ++;
          }
          ++ pCur;
        } while( pCur != pEnd );
        EdgeInfo eInfo;
        eInfo.v1 =  pEdge->get<0>()->vertex( pEdge->get<1>() );
        eInfo.v2 =  pEdge->get<0>()->vertex( pEdge->get<2>() );
        eInfo.nComplexFacet = nComplexFacet;
        eInfo.nEdge  = NumGrainNgbs;
        eInfo.fLength = CGAL::sqrt( c3t3.triangulation().segment( *pEdge ).squared_length() );
        EdgeList.push_back( eInfo );
        nEdgeID ++;
      }
    }// end infinte edge
  }// end for

  TripleLineFile.close();

  std::cout << "Writing to vtk file " << std::endl;
  
  WriteTripleLineVtk( sVtkFilename, VertexToIDMap, EdgeList );



  // ------------------------------------------
  //  Write N-Vertex
  // ------------------------------------------
  int NJunctionID = 0;
  for( Finite_vertices_iterator pVertex = c3t3.triangulation().finite_vertices_begin();
       pVertex != c3t3.triangulation().finite_vertices_end(); ++ pVertex )
  {
    if( JunctionType( c3t3, pVertex, OrientationMap.size() ) > 3  )  // need to check at least 3 valid incoming edge of > triple line
    {
      vector<Edge>              NgbEdgeList;
      std::set<Facet>           NgbFacetSet;
      std::set<Subdomain_index> DomainSet;
      c3t3.triangulation().incident_edges( pVertex, std::back_inserter( NgbEdgeList ) );

      int NumIncidentNJunction = 0;
      for( int i = 0; i < NgbEdgeList.size(); i ++ )
      {
        int NumIncidentGrains = IsValidEdge( c3t3, NgbEdgeList[i], OrientationMap.size() );
        if( NumIncidentGrains > 2 )  // just changed from 0 -> 2
          NumIncidentNJunction ++;
      }

      for( int i = 0; i < NgbEdgeList.size(); i ++ )
      {
        int NumIncidentGrains = IsValidEdge( c3t3, NgbEdgeList[i], OrientationMap.size() );
        
        if( NumIncidentGrains > 2 )
        {
          Facet_circulator pEnd = c3t3.triangulation().incident_facets( NgbEdgeList[i] );
          Facet_circulator pCur = c3t3.triangulation().incident_facets( NgbEdgeList[i] );
          int nNumFacets = GetNumFacets( c3t3, pCur, pEnd );
          if( pCur != 0 ) // construct for using Edge_circulator
          {
            do
            {
              if( c3t3.is_in_complex( *pCur ) )
              {
                Surface_index SurfID = c3t3.surface_index( *pCur );
                DomainSet.insert( SurfID.first );
                DomainSet.insert( SurfID.second );
                NgbFacetSet.insert( *pCur );
                
                Triangle_3 T = c3t3.triangulation().triangle( pCur->first, pCur->second ); 
                
                // output to N-Vertex File
                NVertexFile << NJunctionID        << " "
                            << NumIncidentGrains  << " "     // num incident grains of the triple line
                            << NumIncidentNJunction << " "   // number of incident N-junction to this point
                            << i << " "
                            << T.supporting_plane().orthogonal_direction() << " "
                            << c3t3.triangulation().segment( NgbEdgeList[i] ) << " " 
                            << CGAL::sqrt( T.squared_area() ) << " "
                            << T.is_degenerate() << " "
                            << SurfID.first << " " << SurfID.second << " " << nNumFacets << " "
                            << pVertex->point().point() << std::endl;
              }
              ++ pCur;
            } while( pCur != pEnd );
          } // end ngb facet traversal
        }
      }
      NJunctionID ++;
    }// end if vertex is N-Junction

  }
  
  return 0;
}

