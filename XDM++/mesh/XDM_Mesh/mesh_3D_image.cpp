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
#include <fstream>
#include <string>
#include <vector>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Lapack/Linear_algebra_lapack.h>
#include "3dMath.h"
#include "SimpleMeshVTKWriter.h"
#include "VTKPropertyMap.h"
#include "GrainCurvatureEst.h"
#include "XDM_Mesh_IO.h"
#include "ProcessMic.h"
#include "Symmetry.h"
#include <boost/shared_ptr.hpp>

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
// Mesh Criteria

typedef unsigned char XDMIndexT;
using std::vector;
typedef vector< vector< vector< XDMIndexT > > >  DataMatrixT;

typedef CGAL::XDM_test::XDM_Data<K> XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::Monge_via_jet_fitting<K>   Monge_via_jet_fitting;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;

typedef CGAL::Mesh_criteria_3<XDM_Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

typedef XDM_C3t3::Facet_iterator Facet_iterator;
typedef XDM_C3t3::Cell_iterator  Cell_iterator;

typedef XDM_C3t3::Surface_index Surface_index;
typedef XDM_C3t3::Triangulation::Vertex_handle Vertex_handle;
typedef XDM_C3t3::Facet Facet;
typedef XDM_C3t3::Triangulation::Cell Cell;
typedef XDM_C3t3::Triangulation::Vertex::Point Point;
typedef XDM_C3t3::Cell_handle Cell_handle;
typedef XDM_C3t3::Triangulation::All_vertices_iterator Vertex_iterator;
typedef XDM_C3t3::Triangulation::Finite_vertices_iterator Finite_vertex_iterator;
typedef XDM_C3t3::Triangulation::Finite_edges_iterator Finite_edges_iterator;

typedef XDM_Mesh_Domain::Subdomain_index Subdomain_index;


typedef std::set< Vertex_handle > Vertex_set;
typedef std::multimap< Subdomain_index, Cell_handle > GrainTetMap;
typedef std::multimap< Subdomain_index, Facet >       GrainFacetMap;

typedef GrainTetMap::iterator GrainTetIter;
typedef GrainFacetMap::iterator GrainFacetIter;
typedef GrainTetMap::const_iterator ConstGrainTetIter;

typedef std::pair< Vertex_handle, Vertex_handle > VertexPairs;

struct SCurvatureData
{
  SCurvatureData(){}
  SCurvatureData( float fK1_, float fK2_, float fCondition_ ):
    fK1( fK1_ ), fK2( fK2_ ), fConditionNum( fCondition_ ) {}
  
  float fK1;
  float fK2;
  float fConditionNum;
};

// each grain has a VertexCurveMap
typedef std::map< Vertex_handle, SCurvatureData > VertexCurveMap ;
typedef std::map< Facet, SCurvatureData > FacetCurveMap; 


//-----------------------------------------------------------------------------
//  WriteTripleLineVtk
//-----------------------------------------------------------------------------
void WriteTripleLineVtk( const std::string & oFilename,
                         std::map< Vertex_handle, int > & IDMap,
                         const std::vector< VertexPairs > & TripleLines )
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
    os << "2 " <<  IDMap[ TripleLines[i].first ] << " "
       << IDMap[ TripleLines[i].second ] << std::endl;
  
  os << "CELL_TYPES " << TripleLines.size() << std::endl;
  for( int i = 0; i < TripleLines.size(); i ++ )
    os << 3 << std::endl;
}

//-----------------------------------------------------------------------------
//  WriteCurveVTK
//
//-----------------------------------------------------------------------------
void WriteCurveVtk( const std::string & oFilename,
                    const XDM_C3t3 & c3t3,
                    GrainFacetMap & FacetMap,
                    Subdomain_index nGrainID,
                    const VertexCurveMap & oMap )
{
  
  std::ofstream os( oFilename.c_str( ) );
  typedef std::map< Vertex_handle, int > IDMap;
  IDMap VertexIDMap;
  int nID = 0;   // following vtk style

  int nFacets = FacetMap.count( nGrainID );
  os << "# vtk DataFile Version 2.0 \n"
     << "Volume Test \n"
     << "ASCII \n"
     << "DATASET UNSTRUCTURED_GRID \n"
     << "POINTS " << oMap.size() << " float \n";
  
  typedef VertexCurveMap::const_iterator MapIter;

  for( MapIter pCur = oMap.begin(); pCur != oMap.end(); ++ pCur )
  {
    VertexIDMap.insert( std::make_pair( pCur->first, nID ) );
    nID ++;
    os << pCur->first->point().x() << " "
       << pCur->first->point().y() << " "
       << pCur->first->point().z() << " "
       << std::endl;
  }

  os << "CELLS " << nFacets << " " << nFacets * 4 << std::endl;

  
  std::pair<GrainFacetIter, GrainFacetIter> pRange;
  pRange = FacetMap.equal_range( nGrainID );
  for( GrainFacetIter pCur = pRange.first; pCur != pRange.second; ++ pCur )
  {
    os << "3 ";
    for( int n = 0; n < 4; n ++ )
    {
      if( n != pCur->second.second )
        os << VertexIDMap[ pCur->second.first->vertex( n ) ] << " ";
    }
    os << std::endl;
  }
  
  os << "CELL_TYPES " << nFacets << std::endl;
  for( int i = 0; i < nFacets; i ++ )
    os << 5 << std::endl;
  
  os << "POINT_DATA " << VertexIDMap.size() << std::endl;
  os << "SCALARS Mean_Curvature float 1\n";
  os << "LOOKUP_TABLE default \n";

  for( MapIter pCur = oMap.begin(); pCur != oMap.end(); ++ pCur )
    os << ( pCur->second.fK1 + pCur->second.fK2 ) / float( 2 ) << std::endl;

  os << "SCALARS k1 float 1\n";
  os << "LOOKUP_TABLE default \n";
  for( MapIter pCur = oMap.begin(); pCur != oMap.end(); ++ pCur )
    os << pCur->second.fK1 << std::endl;

  
  os << "SCALARS k2 float 1\n";
  os << "LOOKUP_TABLE default \n";
  for( MapIter pCur = oMap.begin(); pCur != oMap.end(); ++ pCur )
    os << pCur->second.fK2 << std::endl;
  
  os << "SCALARS Condition float 1\n";
  os << "LOOKUP_TABLE default \n";
  for( MapIter pCur = oMap.begin(); pCur != oMap.end(); ++ pCur )
    os << pCur->second.fConditionNum << std::endl;

  
  os.close();
}



//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

int NumIncidentingDomains( const XDM_C3t3 & c3t3,
                           const Vertex_handle & v )
{
  vector< Cell_handle > oIncidentCells;
  c3t3.triangulation().incident_cells( v, std::back_inserter( oIncidentCells ) );
  
  typedef std::set< Subdomain_index > DomainSet;
  DomainSet UniqueDomains;
  for( vector< Cell_handle >::iterator pCur = oIncidentCells.begin();
       pCur != oIncidentCells.end(); ++ pCur )
    UniqueDomains.insert( c3t3.subdomain_index( *pCur ) );
  
  return UniqueDomains.size();
}

//-----------------------------------------------------------------------------
//
//
//-----------------------------------------------------------------------------
bool IsTripleLineEdge(  const XDM_C3t3 & c3t3,
                        const Finite_edges_iterator & pEdge )
{
  typedef XDM_C3t3::Triangulation::Cell_circulator  Cell_circulator;

  std::set< Subdomain_index > DomainSet;
  Cell_circulator pEnd = c3t3.triangulation().incident_cells( *pEdge );
  Cell_circulator pCur = c3t3.triangulation().incident_cells( *pEdge );

  if( pCur != 0 )
  {
    do
    {
      DomainSet.insert( c3t3.subdomain_index( pCur ) );
      ++ pCur;
    } while( pCur != pEnd );
  }
  return DomainSet.size() >= 3;
}

template< class PMap >
struct ToRF
{
  template < class T >
  GeneralLib::SVector3 operator()( const T & o  ) const
  {
    float fColorEdge = (sqrt( 2.0 ) - 1) + 0.001;
    
    PMap Map;
    GeneralLib::SQuaternion q = Map( o );
    GeneralLib::SVector3 rod( q.m_fX, q.m_fY, q.m_fZ );
    rod /= q.m_fW;

    rod += GeneralLib::SVector3(  fColorEdge, fColorEdge, fColorEdge );
    rod /= ( 2.0 * fColorEdge );
    return rod;
  }
};


int main()
{
  XDM_Data oTestData;

  std::string oFilename;
  int nX, nY, nZ;
  float fX, fY, fZ;
  std::cout << "Enter File name: " << std::endl;
  std::cin >> oFilename;
  std::cout << "Enter dimension (x, y, z): " << std::endl;
  std::cin >> nX >> std::skipws >> nY >> std::skipws >> nZ;
  std::cout << "Enter voxel dimension (fx, fy, fz): " << std::endl;
  std::cin >> fX >> std::skipws >> fY >> std::skipws >> fZ;
  
  XDMCGal::ReadDxFile( oTestData, oFilename, nX + 2, nY + 2 , nZ + 2, fX, fY, fZ );
  float fBndSmoothingLength, fPointSpacing, fMinEdgeLength;
  std::cout << "=== Boundary Point Parameters: " << std::endl;
  std::cout << "BndSmoothingLength, Point Spacing, Min Bnd Edge Length: " << std::endl;
  std::cin >> std::skipws >> fBndSmoothingLength
           >> std::skipws >> fPointSpacing
           >> std::skipws >> fMinEdgeLength;
  oTestData.MakeBndPoints( fBndSmoothingLength, fPointSpacing, fMinEdgeLength );
  XDM_Mesh_Domain oTestDomain( oTestData );

  float fAngle, fSize, fApproximation;

  std::cout << "Enter Angle (Degree)  Size  Approximation: " << std::endl;
  std::cin >> fAngle >> std::skipws >> fSize >> std::skipws >> fApproximation;
  std::cout << "Angle " << fAngle << " Size " << fSize << " Approximation " << fApproximation << std::endl;
  Facet_criteria facet_criteria( fAngle, fSize, fApproximation); // angle, size, approximation
  float fRatio;
  std::cout << "Enter Aspect Ratio and Tet size approximation: " << std::endl;
  std::cin >> fRatio >> std::skipws >> fSize;
  std::cout << "Ratio " << fRatio << " Size " << fSize << std::endl;
  
  Cell_criteria cell_criteria( fRatio, fSize ); // radius-edge ratio, size

  //Cell_criteria cell_criteria(1, 1); // radius-edge ratio, size
  
  Mesh_criteria criteria(facet_criteria, cell_criteria);
  std::cout << "started meshing " << std::endl;
  // Meshing
  XDM_C3t3 c3t3_orig = CGAL::XDM_make_mesh_3<XDM_C3t3>( oTestDomain, criteria );
  std::cout << "done " << std::endl;
  
  // Output
  std::ofstream medit_file("out.mesh");
  c3t3_orig.output_to_medit(medit_file);
   
  typedef MicAnalysis::CMicVolume::ShapePtr ShapePtr; 
  typedef Pandora::VoxelFieldSelector::ConfidenceSelector ConfSelector;
  typedef Pandora::VoxelFieldSelector::EulerAngleSelector EulerAngleSelector;
  typedef Pandora::VoxelFieldSelector::RodriguzSelector   RodriguzSelector;
  typedef Pandora::VoxelFieldSelector::RodColorSelector   RodColorSelector;
  typedef Pandora::VoxelFieldSelector::CubicFZRodColorSelector CubicFZRodColorSelector;
  typedef Pandora::VoxelFieldSelector::IdentitySelector   IdentitySelector;

    //----------------------------------------------------------------------------
  // To save:
  std::ofstream out_file("test-save-mesh.cgal.binary",
                         std::ios_base::out|std::ios_base::binary);
  CGAL::set_binary_mode(out_file);
  out_file << c3t3_orig;
  out_file.close();
  
  // To load:
  std::ifstream in_file("test-save-mesh.cgal.binary",
                        std::ios_base::in|std::ios_base::binary);
  CGAL::set_binary_mode(in_file);
  XDM_C3t3 c3t3;
  in_file >> c3t3;
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  //  Aligntment between mesh and mic files
  //
  GeneralLib::SVector3 MeshOrigin(0, 0, 0);
  GeneralLib::SVector3 MeshScaleVector( 1.  / 1000.0, 1. / 1000.0, 1.  / 1000.0 );
  MicAnalysis::CMicVolume oVolume = Pandora::ProcessMicVolume( MeshOrigin );
  std::cout << MeshOrigin << std::endl << MeshScaleVector << std::endl;
  //----------------------------------------------------------------------------

  typedef GeneralLib::SVector3                  SVector3;
  typedef GeneralLib::SQuaternion               SQuaternion;
  typedef std::map< Cell_handle, ShapePtr >     HandleShapePtrMap;
  typedef std::map< Cell_handle, SVector3 >     HandleVector3Map;
  typedef std::map< Cell_handle, float >        HandleFloatMap;
  
  typedef std::pair< SQuaternion, float >          MisorientationT;
  typedef std::map< Vertex_handle, ShapePtr >   VertexShapeMap;
  typedef std::map< Vertex_handle, MisorientationT >   VertexMisorientMap;
  typedef Pandora::VoxelFieldSelector::SecondSelector<SQuaternion, float>  SecondSelector;
  typedef Pandora::VoxelFieldSelector::FirstSelector<SQuaternion, float>   FirstSelector;
  
  HandleShapePtrMap CellShapePtrMap;
  HandleVector3Map  CellRFMap;
  HandleVector3Map  CellAveOrientMap;
  HandleFloatMap    CellMisorientMap;
  VertexShapeMap    VertexShapePtrMap;
  VertexMisorientMap VertexToMisorientPairMap;
  VertexMisorientMap VertexToAveMisorientPairMap;

  typedef Pandora
    ::VoxelFieldSelector
    ::ShapePtrListToOrientationGradient< LatticeSymmetry::CCubicSymmetry > ShapePtrListToOrientationGradient;

  typedef Pandora
    ::VoxelFieldSelector
    ::ShapePtrListToLocalAverageMisorientation< LatticeSymmetry::CCubicSymmetry > ShapePtrListToLocalAveMisor;


  std::cout << "Preparing output mesh --- " << std::endl;
  std::cout << "Enter Smoothing Lenghts (x, y, z) (mm) " << std::endl;
  float xSmooth, ySmooth, zSmooth;
  std::cin  >> xSmooth >> std::skipws >> ySmooth >> std::skipws >> zSmooth;
  
  std::cout << "Smoothing Lenghts (x, y, z) (mm) " << std::endl;
  std::cout  << xSmooth << " " << ySmooth << " " << zSmooth << std::endl;
  
  MicAnalysis::CMicVolume::BBox3D VertexSmoothingVolume;
  VertexSmoothingVolume.m_oBoxMin = SVector3( -xSmooth, -ySmooth, -zSmooth ) / float(2);
  VertexSmoothingVolume.m_oBoxMax = SVector3(  xSmooth,  ySmooth,  zSmooth ) / float(2);
  std::cout << "Constructing Cell Voxel Ptr Map" << std::endl;
  Pandora::Details::ConstructCellVoxelPtrMap( CellShapePtrMap, oVolume,
                                              c3t3, MeshOrigin, MeshScaleVector );
  std::cout << "Constructing Vertex Voxel Ptr Map" << std::endl;
  Pandora::Details::ConstructVertexVoxelPtrMap( VertexShapePtrMap, oVolume,
                                                c3t3, MeshOrigin, MeshScaleVector );
  std::cout << "Constructing Cell To Misorientation Map" << std::endl;
  Pandora::Details::ConstructCellToMisorientationMaps( CellMisorientMap, CellRFMap,
                                                       VertexShapePtrMap, c3t3,
                                                       MeshScaleVector,
                                                       LatticeSymmetry::CCubicSymmetry::Get() );
  std::cout << "Constructing Cell To Average Orientation Map" << std::endl;
  Pandora::Details::ConstructCellToGrainAverageMap( CellAveOrientMap, CellShapePtrMap,
                                                    c3t3, LatticeSymmetry::CCubicSymmetry::Get(),
                                                    Pandora::Details::QuatToCubicColor()  );

  std::cout << "Constructing Vertex to pair<Quaternion, float> map" << std::endl;
  
  Pandora::Details::ConstructVertexFieldMap( VertexToMisorientPairMap, oVolume,
                                             c3t3, MeshOrigin, MeshScaleVector,
                                             ShapePtrListToOrientationGradient(),
                                             VertexSmoothingVolume );

  Pandora::Details::ConstructVertexFieldMap( VertexToAveMisorientPairMap, oVolume,
                                             c3t3, MeshOrigin, MeshScaleVector,
                                             ShapePtrListToLocalAveMisor(),
                                             VertexSmoothingVolume );
  
  typedef Pandora::HandleFieldMap< Cell_handle, ShapePtr, Cell_iterator, ConfSelector > CellToConfMap;
  typedef Pandora::HandleFieldMap< Cell_handle, ShapePtr, Cell_iterator, EulerAngleSelector > CellToEulerMap;
  typedef Pandora::HandleFieldMap< Cell_handle, ShapePtr, Cell_iterator, RodriguzSelector >   CellToRodMap;

  typedef Pandora::HandleFieldMap< Cell_handle, SVector3, Cell_iterator, IdentitySelector >   CellToVectorMap;
  typedef Pandora::HandleFieldMap< Cell_handle, float, Cell_iterator, IdentitySelector >      CellToFloatMap;
  
  typedef Pandora::HandleFieldMap< Vertex_handle, MisorientationT, Finite_vertex_iterator,  SecondSelector >      VertexToFloatMap;
  typedef Pandora::HandleFieldMap< Vertex_handle, MisorientationT, Finite_vertex_iterator, ToRF<FirstSelector> > VertexToMisRFMap;
  typedef Pandora::HandleFieldMap< Vertex_handle, ShapePtr, Finite_vertex_iterator, RodColorSelector >   VertexToRodMap;
  typedef Pandora::HandleFieldMap< Vertex_handle, ShapePtr, Finite_vertex_iterator, CubicFZRodColorSelector>   VertexToCubRodMap;
   
  //----------------------
  //
  //  Instatiate maps
  //   Add property maps to map list
  //   
  //
  //----------------------
  typedef Pandora::VtkPropMap VtkPropMap;
  std::vector< VtkPropMap > PointDataMaps;
  std::vector< VtkPropMap > CellDataMaps;
  std::vector< VtkPropMap > FieldDataMaps;

  PointDataMaps.push_back(  VtkPropMap ( VertexToRodMap( "COLOR_SCALARS", "Vertex_Rodrigues", "3", 3,
                                                         VertexShapePtrMap, RodColorSelector(),
                                                         c3t3.triangulation().finite_vertices_begin(),
                                                         c3t3.triangulation().finite_vertices_end(), false ) ) );

  PointDataMaps.push_back(  VtkPropMap ( VertexToRodMap( "SCALARS", "Vertex_Rodrigues_test", "float 3", 3,
                                                         VertexShapePtrMap, RodColorSelector(),
                                                         c3t3.triangulation().finite_vertices_begin(),
                                                         c3t3.triangulation().finite_vertices_end() ) ) );
  
//   PointDataMaps.push_back(  VtkPropMap ( VertexToCubRodMap( "COLOR_SCALARS", "Cubic_Rodrigues", "3", 3,
//                                                             VertexShapePtrMap, CubicFZRodColorSelector(),
//                                                             c3t3.triangulation().finite_vertices_begin(),
//                                                             c3t3.triangulation().finite_vertices_end(), false ) ) );
  
//   PointDataMaps.push_back(  VtkPropMap ( VertexToFloatMap( "SCALARS", "Vertex_dTheta_dR", "float 1", 1,
//                                                            VertexToMisorientPairMap, SecondSelector(),
//                                                            c3t3.triangulation().finite_vertices_begin(),
//                                                            c3t3.triangulation().finite_vertices_end() ) ) );

  PointDataMaps.push_back(  VtkPropMap ( VertexToFloatMap( "SCALARS", "Ave_Misorient", "float 1", 1,
                                                           VertexToAveMisorientPairMap, SecondSelector(),
                                                           c3t3.triangulation().finite_vertices_begin(),
                                                           c3t3.triangulation().finite_vertices_end() ) ) );
  
  
  PointDataMaps.push_back(  VtkPropMap ( VertexToMisRFMap( "COLOR_SCALARS", "Vertex_MisorientOp", "3", 3,
                                                           VertexToMisorientPairMap, ToRF<FirstSelector>(),
                                                           c3t3.triangulation().finite_vertices_begin(),
                                                           c3t3.triangulation().finite_vertices_end(), false ) ) );
  
  CellDataMaps.push_back(  VtkPropMap ( CellToConfMap( "SCALARS", "Confidence", "float 1", 1,
                                                       CellShapePtrMap, ConfSelector(),
                                                       c3t3.cells_begin(), c3t3.cells_end() )   ) );

  CellDataMaps.push_back(  VtkPropMap ( CellToRodMap( "COLOR_SCALARS", "Rodrigues", "3", 3,
                                                      CellShapePtrMap, RodriguzSelector(),
                                                      c3t3.cells_begin(), c3t3.cells_end(), false ) ) );

  CellDataMaps.push_back(  VtkPropMap ( CellToVectorMap( "COLOR_SCALARS", "MisorientOp", "3", 3,
                                                         CellRFMap, IdentitySelector(),
                                                         c3t3.cells_begin(), c3t3.cells_end(), false ) ) );

  CellDataMaps.push_back(  VtkPropMap ( CellToVectorMap( "COLOR_SCALARS", "Average_Orientation", "3", 3,
                                                         CellAveOrientMap, IdentitySelector(),
                                                         c3t3.cells_begin(), c3t3.cells_end(), false ) ) );

//   FieldDataMaps.push_back(  VtkPropMap ( CellToEulerMap( "SCALARS", "Euler_Angles", "float", 3,
//                                                          CellShapePtrMap, EulerAngleSelector(),
//                                                          c3t3.cells_begin(), c3t3.cells_end() ) ) );
  
//   FieldDataMaps.push_back(  VtkPropMap ( CellToFloatMap( "SCALARS", "dTheta_dR", "float", 1,
//                                                          CellMisorientMap, IdentitySelector(),
//                                                          c3t3.cells_begin(), c3t3.cells_end() ) ) );


  std::cout << "Enter Output vtk name: " << std::endl;
  std::string OutputVtkName;
  std::cin >> OutputVtkName;
  std::cout << "Writing to " << OutputVtkName << std::endl;
  std::ofstream vtkFile( OutputVtkName.c_str() );
  Pandora::WriteMeshToUnstructuredVtk( vtkFile, c3t3,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps.begin(), CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  
  vtkFile.close();

  std::cout << "Done " << std::endl;
  std::ofstream BndIDToGrainIDFile( "BndID_2_GrainID_File.txt");
  Pandora::PrintBndIDToGrainID( BndIDToGrainIDFile, c3t3 );
  BndIDToGrainIDFile.close();


  return 0;
}
