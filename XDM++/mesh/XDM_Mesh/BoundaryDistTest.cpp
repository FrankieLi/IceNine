#include <CGAL/basic.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include "AABB_mesh_3_triangle_primitive.h"

#include <CGAL/AABB_triangle_primitive.h>

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
#include "GrainGeometryAnalysis.h"
#include "DynamicsAnalysis.h"

#include <fstream>
#include <string>
#include <vector>
#include "Quaternion.h"
#include "GrainCurvatureEst.h"
#include "XDM_Mesh_IO.h"
#include "Symmetry.h"
#include "MeshUtilities.h"
#include "MeshAlignment.h"
#include "MeshBoundaryAligner.h"

#include "SimpleMeshVTKWriter.h"
#include "VTKPropertyMap.h"

#include <CGAL/iterator.h>

#include "SpecialBoundaryDistance.h"

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::XDM_test::XDM_Data<K> XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;

typedef std::vector< GeneralLib::SMatrix3x3 > IDOrientMap;


//--------------------------------------
//  BuildHistoryMap
//--------------------------------------

typedef std::map<int, int> HistoryMapT;
void BuildHistoryMap(  HistoryMapT & Map1, HistoryMapT & Map2, 
                       const std::map<int, int> & Mesh1_ToMesh2_IDMap,
                       const std::set<int> & IDSet1,
                       const std::set<int> & IDSet2 )
{
  typedef std::set<int>::const_iterator SetIter;

  for( SetIter it = IDSet1.begin(); it != IDSet1.end(); ++ it)
    Map1[ *it ] = -5;
  for( SetIter it = IDSet2.begin(); it != IDSet2.end(); ++ it)
    Map2[ *it ] = -5;
  
  

  typedef HistoryMapT::const_iterator CIter;
  for( CIter it = Mesh1_ToMesh2_IDMap.begin();
       it != Mesh1_ToMesh2_IDMap.end(); ++ it)
  {
    int nID = (rand() % Mesh1_ToMesh2_IDMap.size() + 1);
    Map1[ it->first  ] = nID;
    Map2[ it->second ] = nID;
  }
}

void BuildInteriorMap(  std::map<int, int>  & Map, 
                        const std::set<int> & MeshIDs,
                        const std::set<int> & SurfaceIDSet)
{
  int n = 0;
  typedef std::set<int>::const_iterator SetIter;

  Map[ -1 ] = 0;
  Map[  0 ] = 0;
  for( SetIter it = MeshIDs.begin(); it != MeshIDs.end(); ++ it)
  {
    if( SurfaceIDSet.find( *it ) != SurfaceIDSet.end() )
      Map[ *it ] = 0;
    else
      Map[ *it ] = 1;  // 1 is interior
  }
}


typedef XDM_C3t3::Cell_iterator                Cell_iterator;
typedef XDM_C3t3::Cell_handle Cell_handle;
typedef std::map<Cell_handle, int> CellIDMapT;
CellIDMapT BuildCellIDMap( const XDM_C3t3 & c3t3_,
                           const HistoryMapT IDToHistoryMap )
{
  CellIDMapT IDMap;
  for( Cell_iterator it = c3t3_.cells_begin();
       it != c3t3_.cells_end(); ++ it )
  {
    int nID = c3t3_.subdomain_index( it );
    if( IDToHistoryMap.find( nID ) != IDToHistoryMap.end() )
      IDMap[ it ] = IDToHistoryMap.find( nID )->second;
    else
      IDMap[it]   = -5;
  }
  return IDMap;
}




typedef std::map<Cell_handle, SVector3> CellToSVector3Map;

//---------------------------------------------------
//  Return orientation mapped to RGB color
//---------------------------------------------------
CellToSVector3Map GetRodriguzVectorMap( const XDM_C3t3& c3t3_,
                                        const IDOrientMap & OrientMap )
{
  CellToSVector3Map CellToOrientMap;
  for( Cell_iterator it = c3t3_.cells_begin();
       it != c3t3_.cells_end(); ++ it )
  {
    int nID = c3t3_.subdomain_index( it );
    if( nID < 0)
      CellToOrientMap[ it ] = SVector3( 0, 0, 0 );
    else if ( nID >= OrientMap.size() )
      std::cerr << "Average orientation file and mesh are mismatching.  ID exceeds size." << std::endl;
    else
    {
      SQuaternion q;
      q.Set( OrientMap[ nID ] );
      SVector3 rod( q.m_fX, q.m_fY, q.m_fZ );
      rod /= q.m_fW;
      // ToColor

      const float fColorEdge = (sqrt( 2.0 ) - 1) + 0.001;
      rod += SVector3(  fColorEdge, fColorEdge, fColorEdge );
      rod /= ( 2.0 * fColorEdge );
      CellToOrientMap[ it ] = rod;
    }
  }
  return CellToOrientMap;
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

//----------------------------------------
//  WriteIDOrientMapx
//----------------------------------------
void WriteIDOrientMap( const IDOrientMap & Map, const std::string & filename )
{
  std::ofstream MapOutput( filename.c_str() );
  for( int i = 0; i < Map.size(); i ++ )
    MapOutput << GeneralLib::RadianToDegree( Map[i].GetEulerAngles() ) << std::endl;

  MapOutput.close();
}


struct CubicFZThreshold
{
  float fThresh;
  bool operator( )( const GeneralLib::SQuaternion & q1,
                    const GeneralLib::SQuaternion & q2 )
  {
    float fMis =
      LatticeSymmetry
      ::GetMisorientation( LatticeSymmetry
                           ::CCubicSymmetry
                           ::Get(), q1, q2 );
    return fThresh >= fMis;
  }
  
};


int main()
{
  //----------------------------------------------------------------------------
  typedef CGAL::Simple_cartesian<double> AABB_Kernel;
  typedef AABB_Kernel  Geom_traits;
  typedef Geom_traits::Ray_3 Ray;
  typedef Geom_traits::Point_3 Point;
  typedef Geom_traits::Triangle_3 Triangle;
  typedef Geom_traits::Segment_3 Segment;
  
  typedef XDM_C3t3::Facet_iterator  Facet_iterator;
  typedef XDM_C3t3::Triangulation::Finite_facets_iterator  Finite_facet_iterator;
  typedef CGAL::Pandora::BoundarySelectionPredicate<XDM_C3t3, Finite_facet_iterator >  BoundarySelectFilter;
  
  typedef CGAL::Filter_iterator<Finite_facet_iterator, BoundarySelectFilter> Filtered_facet_iterator;
   
  typedef CGAL::AABB_mesh_3_triangle_primitive< XDM_C3t3, AABB_Kernel, Filtered_facet_iterator > AABB_Primitive;
  typedef CGAL::AABB_traits<AABB_Kernel, AABB_Primitive> AABB_Mesh_Triangle_Traits;
  typedef CGAL::AABB_tree<AABB_Mesh_Triangle_Traits> Mesh_AABB_Tree;

  typedef Mesh_AABB_Tree::Object_and_primitive_id Object_and_primitive_id;
  typedef Mesh_AABB_Tree::Primitive_id Primitive_id;
  typedef Mesh_AABB_Tree::Point_and_primitive_id Point_and_primitive_id;
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //                      I N P U T
  XDM_C3t3 c3t3;
  std::string MeshFilename;
  std::cout << "Enter binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename;

  std::string OrientMapFilename;
  std::cout << "Enter ID to orientation map filename: " << std::endl;
  std::cin  >> OrientMapFilename;

  
  std::ifstream MeshIS( MeshFilename.c_str(),
                        std::ios_base::in|std::ios_base::binary );
  CGAL::set_binary_mode( MeshIS );
  
  MeshIS >> c3t3;
  MeshIS.close();
  
  IDOrientMap OrientMap;
  ReadIDOrientMap( OrientMap, OrientMapFilename );  
  std::cout << "Done Reading " << std::endl;
  //
  //----------------------------------------------------------------------------
    
  Pandora::SpecialBoundaryDistance<XDM_C3t3> SinkDistance( c3t3, OrientMap );
  
  
  typedef Pandora::SpecialBoundaryDistance<XDM_C3t3>::Point3 BarePoint;

  double x, y, z;
  std::cout << "Enter test point " << std::endl;
  std::cin >> std::skipws >> x
            >> std::skipws >> y
            >> std::skipws >> z; 

  SinkDistance.DefineBoundaryAsSink();
  std::cout << SinkDistance.DistanceToNearestSink( BarePoint( x, y, z ) ) << std::endl;;
  
  
  return 0;
}

