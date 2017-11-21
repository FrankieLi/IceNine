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
#include "GrainGeometryAnalysis.h"

#include "M_EstimateWithIRLS.h"

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
typedef XDM_C3t3::Triangulation::Edge Edge;
typedef CGAL::Mesh_criteria_3<XDM_Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

// template < class CellIter, class C3T3, class Map >
// void BuildCellIDMap( Map & m, const C3T3 & c3t3_,
//                      CellIter pCur, CellIter pEnd )
// {
//   for(; pCur != pEnd; ++ pCur)
//     m[ pCur ] = c3t3_.subdomain_index( pCur );
// }


struct SetBBox
{
  template< class T1, class T2, class T3 >
  void operator()( T1 & map, const T2 & key, const T3 & value ) const
  {
    map[key].BBox = value;
  }
};
struct SetMeanwidth
{
  template< class T1, class T2, class T3 >
  void operator()( T1 & map, const T2 & key, const T3 & value ) const
  {
    map[key].MeanWidth = value;
  }
};
struct SetVolume
{
  template< class T1, class T2, class T3 >
  void operator()( T1 & map, const T2 & key, const T3 & value ) const
  {
    map[key].Volume = value;
  }
};

struct SetInterior
{
  template< class T1, class T2, class T3 >
  void operator()( T1 & map, const T2 & key, const T3 & value ) const
  {
    if( value )
      map[key].nIsInterior = 1;
    else
      map[key].nIsInterior = 0;
  }
};

int main()
{
  std::string MeshFilename;
  std::cout << "Enter binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename;

  std::string GeometryFilename;
  std::cout << "Enter filename for Geometry File " << std::endl;
  std::cin  >> GeometryFilename;

  std::cout << "Enter Output vtk name: " << std::endl;
  std::string OutputVtkName;
  std::cin >> OutputVtkName;
 
  
  //-------------------------------------------------------------------
  //      Start computing and meshing
  //-------------------------------------------------------------------
  std::cout << "Start reading mesh  " << std::endl;
  std::ifstream MeshIs( MeshFilename.c_str(),
                        std::ios_base::in|std::ios_base::binary);
  CGAL::set_binary_mode( MeshIs );
  XDM_C3t3 c3t3;
  
  MeshIs >> c3t3;
  MeshIs.close();

  typedef MicAnalysis::CMicVolume::ShapePtr ShapePtr; 
  
  typedef GeneralLib::SVector3                  SVector3;
  typedef GeneralLib::SQuaternion               SQuaternion;

  typedef XDM_C3t3::Subdomain_index              Subdomain_index;
  typedef XDM_C3t3::Surface_index                Surface_index;
  typedef XDM_C3t3::Cell_handle                  Cell_handle;
  typedef XDM_C3t3::Triangulation::Vertex_handle Vertex_handle;
  typedef XDM_C3t3::Cell_iterator                Cell_iterator;
  typedef XDM_C3t3::Facet_iterator               Facet_iterator;
  
  typedef std::map< Cell_handle, ShapePtr >     HandleShapePtrMap;
  typedef std::map< Cell_handle, SQuaternion >  HandleQuaternionMap;
  typedef std::map< Cell_handle, float >        HandleFloatMap;
  
  typedef std::pair< SQuaternion, float >          MisorientationT;
  typedef std::map< Vertex_handle, ShapePtr >   VertexShapeMap;
  typedef std::map< Vertex_handle, MisorientationT >   VertexMisorientMap;
  typedef Pandora::VoxelFieldSelector::SecondSelector<SQuaternion, float>  SecondSelector;
  typedef Pandora::VoxelFieldSelector::FirstSelector <SQuaternion, float>  FirstSelector;

  //-----------------------------------------------------------------------
  //  Printing
  typedef Pandora::VtkPropMap VtkPropMap;
  typedef Pandora::VoxelFieldSelector::IdentitySelector   IdentitySelector;
  
  std::vector< VtkPropMap > PointDataMaps;
  std::vector< VtkPropMap > CellDataMaps;
  std::vector< VtkPropMap > FieldDataMaps;

  std::cout << "Writing to " << OutputVtkName << std::endl;
  std::ofstream vtkFile( OutputVtkName.c_str() );
  Pandora::WriteMeshToUnstructuredVtk( vtkFile, c3t3,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps.begin(),  CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  vtkFile.close();
  
  //-----------------------------------------------------------------------
  typedef Pandora::Details::GeomProp        GP;
  typedef std::map< Subdomain_index, GP>    GeomPropMap;

  typedef Pandora::Details::GrainGeometry<XDM_C3t3> GrainGeometry;
  typedef GrainGeometry::IDNgbMapT    IDNgbMapT;
  typedef GrainGeometry::GrainTetMapT GrainTetMapT;
  typedef GrainGeometry::IDToEdgeMapT  IDToEdgeMapT;
  //-----------------------------------------------------------------------
  //  Build maps
  
  IDToEdgeMapT  IDToEdgeMap;
  IDNgbMapT    NgbMap;
  GrainTetMapT GrainTetMap;
  GeomPropMap  IDToGeomMap;
  std::cout << "Building grain surface edge map " << std::endl;
  Pandora::Details::BuildGrainSurfaceEdgeMap( IDToEdgeMap, c3t3 );
  std::cout << "Building grain tet map " << std::endl;
  Pandora::Details::BuildGrainTetMap( GrainTetMap, c3t3 );
  GrainGeometry GeometricExtractor( c3t3 );
  std::cout << "Calculating mean width  " << std::endl;
  GeometricExtractor.CalculateMeanwidth( IDToGeomMap, SetMeanwidth(), IDToEdgeMap );
  std::cout << "Calculating volume  " << std::endl;
  GeometricExtractor.CalculateVolume   ( IDToGeomMap, SetVolume(),    GrainTetMap );
  std::cout << "Building neighbor map  " << std::endl;
  Pandora::Details::BuildGrainNgbMap( NgbMap, c3t3 );
  std::cout << "Identify interior grains  " << std::endl;
  GeometricExtractor.SelectInterior    ( IDToGeomMap, SetInterior(),
                                         NgbMap, GrainTetMap );
  
  //-----------------------------------------------------------------------
  
  std::ofstream GeometryOS( GeometryFilename.c_str() );

  typedef GeomPropMap::iterator GeomIter;
  std::cout << "Number of grains in orientation map " << IDToGeomMap.size() << std::endl;

  double fTotalInteriorVolume = 0;
  double fTotalVolume         = 0;
  int nInteriorGrains = 0;
  for( GeomIter pIter = IDToGeomMap.begin(); pIter != IDToGeomMap.end(); ++ pIter )
  {
    GeometryOS << pIter->first << " "
               << pIter->second.MeanWidth << " "
               << pIter->second.Volume << " "
               << pIter->second.nIsInterior << " "
               << NgbMap[pIter->first].size() << " "
               << IDToEdgeMap.count( pIter->first ) << std::endl;
    if( pIter->first > 0 )
    {
      if( pIter->second.nIsInterior == 1 )
      {
        fTotalInteriorVolume += pIter->second.Volume;
        nInteriorGrains ++;
      }
      fTotalVolume += pIter->second.Volume;
    }
  }
  GeometryOS.close();

  double fInternalSurfaceArea = 0;
  double fTotalSurfaceArea    = 0;
  int nInternal = 0;
  int nTotal    = 0;

  //  Get boundary area
  for( Facet_iterator it = c3t3.facets_begin();
       it != c3t3.facets_end(); ++ it )
  {
    Surface_index sID = c3t3.surface_index( *it );

    if( sID.first > 0 && sID.second > 0 )
    {
      vector<int> V;
      
      for( int i = 0; i < 4; i ++ )
        if( i != it->second )
          V.push_back( i );
      K::Triangle_3 T( it->first->vertex( V[0] )->point(),
                       it->first->vertex( V[1] )->point(),
                       it->first->vertex( V[2] )->point() );
      double fPatchArea = static_cast<double>( std::sqrt( T.squared_area () ) );
      
      if( IDToGeomMap[sID.first].nIsInterior == 1
          || IDToGeomMap[sID.second].nIsInterior == 1 )
      {
        nInternal ++;
        fInternalSurfaceArea += fPatchArea;
      }
      nTotal ++;
      fTotalSurfaceArea += fPatchArea;
    }
  }
  std::cout << "Total Area, volume "    << fTotalSurfaceArea << ", "
            << fTotalVolume <<  " " << nTotal << " " << IDToGeomMap.size() << std::endl;
  std::cout << "Internal Area, volume " << fInternalSurfaceArea << ", "
            << fTotalInteriorVolume << " " << nInternal << " " << nInteriorGrains << std::endl;

  
  return 0;
}

