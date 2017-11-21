#include <CGAL/basic.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include "XDM_Mesh/AABB_mesh_3_triangle_primitive.h"

#include <CGAL/AABB_triangle_primitive.h>

#include "BoundaryAnalysis.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "XDM_Mesh/XDM_mesh_triangulation_3.h"
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include "XDM_Mesh/XDM_mesh_criteria_3.h"

#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include "XDM_Mesh/XDM_mesh_domain_3.h"
#include "XDM_Mesh/XDM_Data.h"
#include "XDM_Mesh/XDM_make_mesh.h"
#include "XDM_Mesh/GrainGeometryAnalysis.h"
#include "XDM_Mesh/DynamicsAnalysis.h"

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
  XDM_C3t3 c3t3_a, c3t3_b;

  std::string MeshFilename1;
  std::cout << "Enter first binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename1;

  std::string MeshFilename2;
  std::cout << "Enter second binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename2;

  std::string OrientMapFilename_a;
  std::cout << "Enter first ID to orientation map filename: " << std::endl;
  std::cin  >> OrientMapFilename_a;

  std::string OrientMapFilename_b;
  std::cout << "Enter second ID to orientation map filename: " << std::endl;
  std::cin  >> OrientMapFilename_b;

  std::ifstream MeshIS( MeshFilename1.c_str(),
                        std::ios_base::in|std::ios_base::binary );
  CGAL::set_binary_mode( MeshIS );
  
  MeshIS >> c3t3_a;
  MeshIS.close();

  MeshIS.open( MeshFilename2.c_str(),
               std::ios_base::in|std::ios_base::binary );
  CGAL::set_binary_mode( MeshIS );
  
  MeshIS >> c3t3_b;
  MeshIS.close();

  IDOrientMap OrientMap_a, OrientMap_b;
  ReadIDOrientMap( OrientMap_a, OrientMapFilename_a );
  ReadIDOrientMap( OrientMap_b, OrientMapFilename_b );
  std::cout << "Done Reading " << std::endl;
  std::cout << "Enter number of neighbors to use " << std::endl;

  int nNgbs;
  std::cin >> std::skipws >> nNgbs;
  
  float fAngleThreshold;
  std::cout << "Enter Misorientation threshold (degrees) " << std::endl;
  std::cin >> std::skipws >> fAngleThreshold;
  
  //typedef XDM_Tr::Geom_traits Geom_traits;
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

  //----------------------------------------------------------------
  
  //  Mesh_AABB_Tree TestTree( FilteredBegin, FilteredEnd );
  

  
  //
  //
  //----------------------------------------------------------------------------------------

  
  Pandora::BoundaryAnalysis<XDM_C3t3> BndAnalysis( c3t3_a, c3t3_b, OrientMap_a, OrientMap_b );
  typedef Pandora::BoundaryAnalysis<XDM_C3t3>::BndMotionInfo BndMotionInfo;

  std::cout << "Building MeshFacetMap " << std::endl;
  BndAnalysis.InitializeByBoundary();
  vector<BndMotionInfo> BndMotionList;// =  BndAnalysis.MeshFacetMap( nNgbs, DEGREE_TO_RADIAN( fAngleThreshold ) );
  //vector<BndMotionInfo> BndMotionList;
  std::cout << "Begin projected distance " << std::endl;
  vector<BndMotionInfo> BndMotionList_ProjectedDistance =  BndAnalysis.ProjectedDistanceBndMotion( DEGREE_TO_RADIAN( fAngleThreshold ) );
  // vector<BndMotionInfo> BndMotionList_ProjectedDistance;
    
  std::string OutputFilename;
  std::cout << "Enter output filename" << std::endl;
  std::cin >> std::skipws >> OutputFilename;
  std::ofstream Outfile( OutputFilename.c_str() );
  for( Size_Type i = 0; i < BndMotionList.size(); i ++ )
  {
    Outfile << BndMotionList[i].Bnd1Vertex[0].point().point() << " "
            << BndMotionList[i].Bnd1Vertex[1].point().point() << " "
            << BndMotionList[i].Bnd1Vertex[2].point().point() << " "
      //        << BndMotionList[i].BndFacet1[0] << " "    using the above method to make sure that curvature is associated with vertex
      //      << BndMotionList[i].BndFacet1[1] << " "
      //      << BndMotionList[i].BndFacet1[2] << " "

            << CGAL::sqrt( BndMotionList[i].BndFacet1.squared_area() ) << " "
            << BndMotionList[i].FacetID1.first << " "
            << BndMotionList[i].FacetID1.second << " "
            << BndMotionList[i].FacetID2.first << " "
            << BndMotionList[i].FacetID2.second << " "
            << BndMotionList[i].Pair1.first << " "
            << BndMotionList[i].Pair1.second << " "
            << BndMotionList[i].Pair2.first << " "
            << BndMotionList[i].Pair2.second << " "
            << BndMotionList[i].fDistance << std::endl;
  }
  Outfile.close();
  typedef XDM_C3t3::Subdomain_index Subdomain_index;
  std::set<Subdomain_index> SurfaceGrainIDSet_a;
  std::set<Subdomain_index> SurfaceGrainIDSet_b;
  Pandora::Details::GetSurfaceGrainIDs( SurfaceGrainIDSet_a, c3t3_a, -1 );
  Pandora::Details::GetSurfaceGrainIDs( SurfaceGrainIDSet_b, c3t3_b, -1 );

  std::cout << " Ouputting Surface Grain Index " << std::endl;
  std::stringstream SurfaceGrainFile_a;
  SurfaceGrainFile_a << OutputFilename << ".SurfaceGrainID_a.txt";
  std::ofstream SurfGrainFile_a( SurfaceGrainFile_a.str().c_str() );
  for( std::set<int>::iterator it = SurfaceGrainIDSet_a.begin();
       it != SurfaceGrainIDSet_a.end(); ++ it )
    SurfGrainFile_a << *it << std::endl;

  std::stringstream SurfaceGrainFile_b;
  SurfaceGrainFile_b << OutputFilename << ".SurfaceGrainID_b.txt";
  std::ofstream SurfGrainFile_b( SurfaceGrainFile_b.str().c_str() );
  for( std::set<int>::iterator it = SurfaceGrainIDSet_b.begin();
       it != SurfaceGrainIDSet_b.end(); ++ it )
    SurfGrainFile_b << *it << std::endl;

  SurfGrainFile_a.close();
  SurfGrainFile_b.close();
  
  //------------------------------------------------------------------
  std::stringstream ss;
  ss << OutputFilename << ".projected";
  std::ofstream ProjectedOutfile( ss.str().c_str() );
  for( Size_Type i = 0; i < BndMotionList_ProjectedDistance.size(); i ++ )
  {
    bool IsSurface_a = false;
    bool IsSurface_b = false;
    
    if( SurfaceGrainIDSet_a.find( BndMotionList_ProjectedDistance [i].FacetID1.first ) != SurfaceGrainIDSet_a.end()
        || SurfaceGrainIDSet_a.find( BndMotionList_ProjectedDistance [i].FacetID2.first ) != SurfaceGrainIDSet_a.end() )
      IsSurface_a = true;
    
    if( SurfaceGrainIDSet_b.find( BndMotionList_ProjectedDistance [i].FacetID1.second ) != SurfaceGrainIDSet_b.end()
        || SurfaceGrainIDSet_b.find( BndMotionList_ProjectedDistance [i].FacetID2.second ) != SurfaceGrainIDSet_b.end() )
      IsSurface_b = true;
    
    ProjectedOutfile << BndMotionList_ProjectedDistance [i].Bnd1Vertex[0].point().point() << " "
                     << BndMotionList_ProjectedDistance [i].Bnd1Vertex[1].point().point() << " "
                     << BndMotionList_ProjectedDistance [i].Bnd1Vertex[2].point().point() << " "
                     << CGAL::sqrt( BndMotionList_ProjectedDistance[i].BndFacet1.squared_area() ) << " "
                     << BndMotionList_ProjectedDistance [i].FacetID1.first << " "     // ID1 of state a
                     << BndMotionList_ProjectedDistance [i].FacetID1.second << " "    // ID1 of state b
                     << BndMotionList_ProjectedDistance [i].FacetID2.first << " "     // ID2 of state a
                     << BndMotionList_ProjectedDistance [i].FacetID2.second  << " "   // ID2 of state b
                     << BndMotionList_ProjectedDistance [i].Pair1.first << " "        // Orientations
                     << BndMotionList_ProjectedDistance [i].Pair1.second << " "
                     << BndMotionList_ProjectedDistance [i].Pair2.first << " "
                     << BndMotionList_ProjectedDistance [i].Pair2.second << " "
                     << BndMotionList_ProjectedDistance [i].fDistance << " "
                     << BndMotionList_ProjectedDistance [i].fNormalError << " "
                     << BndMotionList_ProjectedDistance [i].bDisappeared << " "
                     << BndMotionList_ProjectedDistance [i].BndFacet2[0] << " "  // left here for compatibility
                     << BndMotionList_ProjectedDistance [i].BndFacet2[1] << " "
                     << BndMotionList_ProjectedDistance [i].BndFacet2[2] << " "
                     << CGAL::sqrt( BndMotionList_ProjectedDistance[i].BndFacet2.squared_area() ) << " "
                     << BndMotionList_ProjectedDistance [i].MinDisplacement << " "
                     << IsSurface_a << " "
                     << IsSurface_b << " "  // 47 (48 in matlab)
                     << BndMotionList_ProjectedDistance[i].VertexCurvature[0] << " "
                     << BndMotionList_ProjectedDistance[i].VertexCurvature[1] << " "
                     << BndMotionList_ProjectedDistance[i].VertexCurvature[2] << " "
                     << BndMotionList_ProjectedDistance[i].VertexMixedArea[0] << " "
                     << BndMotionList_ProjectedDistance[i].VertexMixedArea[1] << " "
                     << BndMotionList_ProjectedDistance[i].VertexMixedArea[2] << " "
                     << BndMotionList_ProjectedDistance[i].Bnd1VertexType[0] << " "
                     << BndMotionList_ProjectedDistance[i].Bnd1VertexType[1] << " "
                     << BndMotionList_ProjectedDistance[i].Bnd1VertexType[2] << " "

                     <<  std::endl;
    
  }
  ProjectedOutfile.close();
  //------------------------------------------------------------------
  
  BndAnalysis.WriteMeshToMeshIDMap("MeshToMeshID.txt");
  
  std::set<int> IDSet1, IDSet2;
  BndAnalysis.GetMeshIDSets( IDSet1, IDSet2 );

  std::map<int, int> Mesh1Map, Mesh2Map;
  std::map<int, int> Mesh1ToMesh2Map = BndAnalysis.GetUniqueMesh1_ToMesh2_IDMapByMajority();
  BuildHistoryMap( Mesh1Map, Mesh2Map, Mesh1ToMesh2Map,
                   IDSet1, IDSet2 );

  std::map<int, int> Mesh1_InteriorIDMap, Mesh2_InteriorIDMap;
  BuildInteriorMap( Mesh1_InteriorIDMap, IDSet1, SurfaceGrainIDSet_a );
  BuildInteriorMap( Mesh2_InteriorIDMap, IDSet2, SurfaceGrainIDSet_b );


  typedef  std::map<int, int>::iterator MapIter;

  std::ofstream Debug0("Debug0.txt" );
  for( MapIter it = Mesh1ToMesh2Map.begin(); it != Mesh1ToMesh2Map.end(); ++ it )
    Debug0 << it->first << " " << it->second << std::endl;
  
  std::ofstream Debug1("Debug1.txt" );
  for( MapIter it = Mesh1Map.begin(); it != Mesh1Map.end(); ++ it )
    Debug1 << it->first << " " << it->second << std::endl;

  std::ofstream Debug2("Debug2.txt" );
  for( MapIter it = Mesh2Map.begin(); it != Mesh2Map.end(); ++ it )
    Debug2 << it->first << " " << it->second << std::endl;

  
  //-----------------------------------------------------------------------
  // Printing
  
  //  typedef XDM_C3t3::Cell_handle                  Cell_handle;
  typedef XDM_C3t3::Cell_iterator                Cell_iterator;
  
  typedef Pandora::VtkPropMap VtkPropMap;
  typedef Pandora::VoxelFieldSelector::IdentitySelector   IdentitySelector;
  
  typedef Pandora::HandleFieldMap< Cell_handle, int, Cell_iterator, IdentitySelector > CellToNewIDMapT;
  typedef Pandora::HandleFieldMap< Cell_handle, SVector3, Cell_iterator, IdentitySelector > CellToVectorMap;
  std::vector< VtkPropMap > PointDataMaps;
  std::vector< VtkPropMap > CellDataMaps_1;
  std::vector< VtkPropMap > CellDataMaps_2;
  std::vector< VtkPropMap > FieldDataMaps;
  
  CellToSVector3Map OrientColorMap_a = GetRodriguzVectorMap( c3t3_a, OrientMap_a );
  CellToSVector3Map OrientColorMap_b = GetRodriguzVectorMap( c3t3_b, OrientMap_b );
  
  CellIDMapT CellToIDMap1 = BuildCellIDMap( c3t3_a, Mesh1Map );
  CellIDMapT CellToIDMap2 = BuildCellIDMap( c3t3_b, Mesh2Map );
  
  CellIDMapT CellToIDInteriorMap1 = BuildCellIDMap( c3t3_a, Mesh1_InteriorIDMap );
  CellIDMapT CellToIDInteriorMap2 = BuildCellIDMap( c3t3_b, Mesh2_InteriorIDMap );
  CellDataMaps_1.push_back(  VtkPropMap ( CellToNewIDMapT( "SCALARS", "HistoryIndex", "int 1", 1,
                                                           CellToIDMap1,
                                                           IdentitySelector(),
                                                           c3t3_a.cells_begin(), c3t3_a.cells_end() )   ) );

  CellDataMaps_1.push_back(  VtkPropMap ( CellToNewIDMapT( "SCALARS", "SurfaceIndex", "int 1", 1,
                                                           CellToIDInteriorMap1,
                                                           IdentitySelector(),
                                                           c3t3_a.cells_begin(), c3t3_a.cells_end() )   ) );
    
  CellDataMaps_1.push_back(  VtkPropMap ( CellToVectorMap( "COLOR_SCALARS", "OrientationRF", "3", 3,
                                                           OrientColorMap_a,
                                                           IdentitySelector(),
                                                           c3t3_a.cells_begin(), c3t3_a.cells_end(), false )   ) );

  // -------------
  CellDataMaps_2.push_back(  VtkPropMap ( CellToNewIDMapT( "SCALARS", "HistoryIndex", "int 1", 1,
                                                           CellToIDMap2,
                                                           IdentitySelector(),
                                                           c3t3_b.cells_begin(), c3t3_b.cells_end() )   ) );
  CellDataMaps_2.push_back(  VtkPropMap ( CellToNewIDMapT( "SCALARS", "SurfaceIndex", "int 1", 1,
                                                           CellToIDInteriorMap2,
                                                           IdentitySelector(),
                                                           c3t3_b.cells_begin(), c3t3_b.cells_end() )   ) );

  CellDataMaps_2.push_back(  VtkPropMap ( CellToVectorMap( "COLOR_SCALARS", "OrientationRF", "3", 3,
                                                           OrientColorMap_b,
                                                           IdentitySelector(),
                                                           c3t3_b.cells_begin(), c3t3_b.cells_end(), false )   ) );
  
  std::string OutputVtkName_a = "Mesh_a.vtk";
  std::string OutputVtkName_b = "Mesh_b.vtk";
  std::cout << "Writing to " << OutputVtkName_a << std::endl;
  std::ofstream vtkFile1( OutputVtkName_a.c_str() );
  Pandora::WriteMeshToUnstructuredVtk( vtkFile1, c3t3_a,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps_1.begin(), CellDataMaps_1.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  vtkFile1.close();
  
  std::cout << "Writing to " << OutputVtkName_b << std::endl;
  
  std::ofstream vtkFile2( OutputVtkName_b.c_str() );
  
  Pandora::WriteMeshToUnstructuredVtk( vtkFile2, c3t3_b,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps_2.begin(), CellDataMaps_2.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  vtkFile2.close();
  
  std::cout << "Done writing " << std::endl;
  //-----------------------------------------------------------------------
  
  
  return 0;
}

