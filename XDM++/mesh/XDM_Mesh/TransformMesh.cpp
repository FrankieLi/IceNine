//-----------------------------------------------------------------
//
//  TransformMesh.cpp
//   Transform a binary mesh.  This is basically manual alignment.
//   It's dumb, but I don't have time to align the last few degrees.
//
//
//-----------------------------------------------------------------
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

#include "SimpleMeshVTKWriter.h"
#include "VTKPropertyMap.h"

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::XDM_test::XDM_Data<K> XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;

typedef std::map< int, GeneralLib::SQuaternion > IDOrientMap;

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
    return ( fThresh >= fMis );
  }
};


int main()
{

  std::string MeshFilenameIn;
  std::cout << "Enter inptut binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilenameIn;
  std::string MeshFilenameOut;
  std::cout << "Enter output binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilenameOut;

  XDM_C3t3 c3t3;
  std::ifstream MeshIS( MeshFilenameIn.c_str(),
                        std::ios_base::in|std::ios_base::binary );
  CGAL::set_binary_mode( MeshIS );
  MeshIS >> c3t3;
  MeshIS.close();
 

  SVector3 Trans;
  SMatrix3x3 Rotation;
  //-------------------------------------------------------------------------
  //  Read in transformations
  std::cout << "Define shift for mesh (x y z):" << std::endl;
  std::cin  >> std::skipws >> Trans.m_fX 
            >> std::skipws >> Trans.m_fY
            >> std::skipws >> Trans.m_fZ; 

  std::cout << "Translation of " << Trans << " in microns " <<  std::endl;
  SVector3 Eulers;
  std::cout << "Rotation for mesh in Euler angles (Degrees):" << std::endl;

  std::cin  >> std::skipws >> Eulers.m_fX 
            >> std::skipws >> Eulers.m_fY
            >> std::skipws >> Eulers.m_fZ; 
  std::cout << "Rotation of " << Eulers << " in degrees " <<  std::endl;
  Eulers = GeneralLib::DegreeToRadian( Eulers );
  Rotation.BuildActiveEulerMatrix( Eulers.m_fX, Eulers.m_fY, Eulers.m_fZ );
  //
  //-------------------------------------------------------------------------

  Pandora::TransformMesh( c3t3, Rotation, Trans );
  
  std::ofstream MeshOs( MeshFilenameOut.c_str(),
                        std::ios_base::out | std::ios_base::binary );
  CGAL::set_binary_mode( MeshOs );
  MeshOs << c3t3;
  MeshOs.close();
  std::cout << "Done outputting binary mesh " << std::endl;
  
  //----------------------------------------  DEBUG purpose

  std::string VtkName;
  std::cout << "VtkFilename for transformed mesh " << std::endl;
  std::cin >> VtkName;
  
                                                                   
  typedef Pandora::VtkPropMap VtkPropMap;
  typedef Pandora::VoxelFieldSelector::IdentitySelector   IdentitySelector;
  
  std::vector< VtkPropMap > PointDataMaps;
  std::vector< VtkPropMap > CellDataMaps;
  std::vector< VtkPropMap > FieldDataMaps;

  std::cout << "Writing to vtk " << std::endl;
  std::ofstream vtkFile( VtkName.c_str() );
  
  Pandora::WriteMeshToUnstructuredVtk( vtkFile, c3t3,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps.begin(),  CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  
  vtkFile.close();

  //----------------------------------------  DEBUG purpose
  return 0;
}
