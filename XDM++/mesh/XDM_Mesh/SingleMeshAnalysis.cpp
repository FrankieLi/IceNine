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

// template < class CellIter, class C3T3, class Map >
// void BuildCellIDMap( Map & m, const C3T3 & c3t3_,
//                      CellIter pCur, CellIter pEnd )
// {
//   for(; pCur != pEnd; ++ pCur)
//     m[ pCur ] = c3t3_.subdomain_index( pCur );
// }

int main()
{
  std::string MeshFilename;
  std::cout << "Enter binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename;

  std::string GrainIDToOrientFilename;
  std::cout << "Enter filename for Grain ID to Orientation Map: " << std::endl;
  std::cin  >> GrainIDToOrientFilename;
  
  
  //-------------------------------------------------------------------
  //  Input for voxel correlations
  //-------------------------------------------------------------------
  std::cout << "Preparing output mesh --- " << std::endl;
  std::cout << "Enter Smoothing Lenghts (x, y, z) (mm) " << std::endl;
  float xSmooth, ySmooth, zSmooth;
  std::cin  >> xSmooth >> std::skipws >> ySmooth >> std::skipws >> zSmooth;
  std::cout << "Smoothing Lenghts (x, y, z) (mm) " << std::endl;
  std::cout  << xSmooth << " " << ySmooth << " " << zSmooth << std::endl;

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
  std::cout << " Done reading binary mesh " << std::endl;
  std::cout << " Start MicVolume calculation and alignment " << std::endl;

  //----------------------------------------------------------------------------
  //  Aligntment between mesh and mic files
  //
  GeneralLib::SVector3 MeshOrigin(0, 0, 0);
  GeneralLib::SVector3 MeshScaleVector( 1.  / 1000.0, 1. / 1000.0, 1.  / 1000.0 );
  MicAnalysis::CMicVolume oVolume = Pandora::ProcessMicVolume( MeshOrigin );
  std::cout << "MeshOrigin " <<  MeshOrigin << std::endl << "MeshScale " <<  MeshScaleVector << std::endl;
  //----------------------------------------------------------------------------
  
  typedef MicAnalysis::CMicVolume::ShapePtr ShapePtr; 
  
  typedef GeneralLib::SVector3                  SVector3;
  typedef GeneralLib::SQuaternion               SQuaternion;

  typedef XDM_C3t3::Subdomain_index              Subdomain_index;
  typedef XDM_C3t3::Cell_handle                  Cell_handle;
  typedef XDM_C3t3::Triangulation::Vertex_handle Vertex_handle;
  typedef XDM_C3t3::Cell_iterator                Cell_iterator;

  
  typedef std::map< Cell_handle, ShapePtr >     HandleShapePtrMap;
  typedef std::map< Cell_handle, SQuaternion >  HandleQuaternionMap;
  typedef std::map< Cell_handle, float >        HandleFloatMap;
  
  typedef std::pair< SQuaternion, float >            MisorientationT;
  typedef std::map< Vertex_handle, ShapePtr >        VertexShapeMap;
  typedef std::map< Vertex_handle, MisorientationT > VertexMisorientMap;
  typedef Pandora::VoxelFieldSelector::SecondSelector<SQuaternion, float>  SecondSelector;
  typedef Pandora::VoxelFieldSelector::FirstSelector<SQuaternion, float>   FirstSelector;
  
  HandleShapePtrMap CellShapePtrMap;
  HandleQuaternionMap  CellAveOrientMap;

  typedef Pandora
    ::VoxelFieldSelector
    ::ShapePtrListToOrientationGradient< LatticeSymmetry::CCubicSymmetry > ShapePtrListToOrientationGradient;

  typedef Pandora
    ::VoxelFieldSelector
    ::ShapePtrListToLocalAverageMisorientation< LatticeSymmetry::CCubicSymmetry > ShapePtrListToLocalAveMisor;

  MicAnalysis::CMicVolume::BBox3D VertexSmoothingVolume;
  VertexSmoothingVolume.m_oBoxMin = SVector3( -xSmooth, -ySmooth, -zSmooth ) / float(2);
  VertexSmoothingVolume.m_oBoxMax = SVector3(  xSmooth,  ySmooth,  zSmooth ) / float(2);
  std::cout << "Constructing Cell Voxel Ptr Map" << std::endl;
  Pandora::Details::ConstructCellVoxelPtrMap( CellShapePtrMap, oVolume,
                                              c3t3, MeshOrigin, MeshScaleVector );
  std::cout << "Constructing Cell To Average Orientation Map" << std::endl;
  Pandora::Details::ConstructCellToGrainAverageMap( CellAveOrientMap, CellShapePtrMap,
                                                    c3t3, LatticeSymmetry::CCubicSymmetry::Get(),
                                                    Pandora::Details::TrivialReturn<SQuaternion, SQuaternion>() );
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
                                       CellDataMaps.begin(), CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  vtkFile.close();
  //-----------------------------------------------------------------------
  
  typedef std::map< Subdomain_index, SQuaternion > ID2OrientMapT;
  ID2OrientMapT ID2OrientMap;
  Pandora::BuildIDToGrainAverageMap( ID2OrientMap, CellAveOrientMap, c3t3 );

  std::ofstream IDToOrientFile( GrainIDToOrientFilename.c_str() );
  typedef ID2OrientMapT::iterator ID2OrientIter;
  
  std::cout << "Number of cells  in CellAveOrientMap " << CellAveOrientMap.size() << std::endl; 
  std::cout << "Number of grains in orientation map " << ID2OrientMap.size() << std::endl;
  
//   for(  ID2OrientIter pCur = ID2OrientMap.begin();
//         pCur != ID2OrientMap.end(); ++ pCur )
//   {
//     IDToOrientFile << pCur->first << " " << pCur->second << std::endl;
//   }
  
//   IDToOrientFile.close();
  
  //----------------------------------------------
  //  Traditional debug output - to be phased out
  //----------------------------------------------
//   std::ofstream medit_file("out.mesh");
//   c3t3.output_to_medit(medit_file);
//   std::ofstream BndIDToGrainIDFile( "BndID_2_GrainID_File.txt");
//   Pandora::PrintBndIDToGrainID( BndIDToGrainIDFile, c3t3 );
//   BndIDToGrainIDFile.close();
  
  return 0;
}


