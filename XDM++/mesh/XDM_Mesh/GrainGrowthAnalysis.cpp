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

//----------------------------------------
//   ReadIDOrientMap
//----------------------------------------
void ReadIDOrientMap( IDOrientMap & Map, const std::string & filename )
{
  std::ifstream MapInput( filename.c_str() );

  if( ! MapInput.is_open() )
  {
    std::cout << "IDOrientMap filename not found  " << std::endl;
    exit( 0 );
  }

  while( !MapInput.eof() )
  {
    GeneralLib::SQuaternion q;
    int n;
    MapInput >> std::skipws >> n
             >> std::skipws >> q.m_fW
             >> std::skipws >> q.m_fX
             >> std::skipws >> q.m_fY
             >> std::skipws >> q.m_fZ >> std::skipws;

    if( !MapInput.eof() )
      Map[n] = q;
  }
  MapInput.close();
}

//----------------------------------------
//   ReadEulerIDMap
//----------------------------------------
void ReadEulerIDMap( IDOrientMap & Map, const std::string & filename )
{
  std::ifstream MapInput( filename.c_str() );
  
  if( ! MapInput.is_open() )
  {
    std::cout << "IDOrientMap filename not found  " << std::endl;
    exit( 0 );
  }

  int n = 0;
  while( !MapInput.eof() )
  {
    GeneralLib::SVector3 Eulers;
  
    MapInput >> std::skipws >> Eulers.m_fX
             >> std::skipws >> Eulers.m_fY
             >> std::skipws >> Eulers.m_fZ;
    Eulers = GeneralLib::DegreeToRadian( Eulers );

    SQuaternion q;
    q.Set( Eulers.m_fX, Eulers.m_fY, Eulers.m_fZ );
    if( !MapInput.eof() )
      Map[n] = q;
    n++;
  }
  std::cout << "Read in " << n << " rows " << std::endl;
  MapInput.close();
}

//----------------------------------------
//   RotationOrientMap
//----------------------------------------
void RotateOrientMap( IDOrientMap & Map, const SMatrix3x3 & R )
{
  GeneralLib::SQuaternion qR;
  qR.Set( R );

  typedef IDOrientMap::iterator Iter;
  for( Iter pCur = Map.begin(); pCur != Map.end(); ++ pCur )
    pCur->second = qR * pCur->second;
  
}
 
struct CubicFZCost
{
  float operator( )( const GeneralLib::SQuaternion & q1,
                    const GeneralLib::SQuaternion & q2 )
  {
    float fMis =
      LatticeSymmetry
      ::GetMisorientation( LatticeSymmetry
                           ::CCubicSymmetry
                           ::Get(), q1, q2 );
    return fMis;
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

  XDM_C3t3 c3t3;
  std::ifstream MeshIS( MeshFilename1.c_str(),
                        std::ios_base::in|std::ios_base::binary );
  CGAL::set_binary_mode( MeshIS );
  if( ! MeshIS.is_open() )
  {
    std::cout << "Input mesh a file not found " << std::endl;
    exit( 0 );
  }
  
  MeshIS >> c3t3_a;
  MeshIS.close();

  MeshIS.open( MeshFilename2.c_str(),
               std::ios_base::in|std::ios_base::binary);
  if( ! MeshIS.is_open() )
  {
    std::cout << "Input mesh b file not found " << std::endl;
    exit( 0 );
  }
  CGAL::set_binary_mode( MeshIS );
  MeshIS >> c3t3_b;
  MeshIS.close();
  
  IDOrientMap OrientMap_a, OrientMap_b;

  ReadEulerIDMap( OrientMap_a, OrientMapFilename_a );
  ReadEulerIDMap( OrientMap_b, OrientMapFilename_b );
  std::cout << "Done reading input file " << std::endl;

  SVector3 Trans_a, Trans_b;
  SMatrix3x3 Rotation_a, Rotation_b;
  //-------------------------------------------------------------------------
  //  Read in transformations
  std::cout << "Define shift for mesh a (x y z):" << std::endl;
  std::cin  >> std::skipws >> Trans_a.m_fX 
            >> std::skipws >> Trans_a.m_fY
            >> std::skipws >> Trans_a.m_fZ; 

  SVector3 Eulers;
  std::cout << "Rotation for mesh a in Euler angles (Degrees):" << std::endl;
  std::cin  >> std::skipws >> Eulers.m_fX 
            >> std::skipws >> Eulers.m_fY
            >> std::skipws >> Eulers.m_fZ;
  Eulers = GeneralLib::DegreeToRadian( Eulers );
  Rotation_a.BuildActiveEulerMatrix( Eulers.m_fX, Eulers.m_fY, Eulers.m_fZ );

  std::cout << "Define shift for mesh b (x y z):" << std::endl;
  std::cin  >> std::skipws >> Trans_b.m_fX 
            >> std::skipws >> Trans_b.m_fY
            >> std::skipws >> Trans_b.m_fZ; 

  std::cout << "Rotation for mesh b in Euler angles (Degrees):" << std::endl;
  std::cin  >> std::skipws >> Eulers.m_fX 
            >> std::skipws >> Eulers.m_fY
            >> std::skipws >> Eulers.m_fZ; 
  Eulers = GeneralLib::DegreeToRadian( Eulers );
  Rotation_b.BuildActiveEulerMatrix( Eulers.m_fX, Eulers.m_fY, Eulers.m_fZ );
  //
  //-------------------------------------------------------------------------
  RotateOrientMap( OrientMap_a, Rotation_a );
  RotateOrientMap( OrientMap_b, Rotation_b );
  
  Pandora::TransformMesh( c3t3_a, Rotation_a, Trans_a );
  Pandora::TransformMesh( c3t3_b, Rotation_b, Trans_b );

  
  typedef std::map< int, Pandora::GrainDynamics::GrainHistoryProp > IDHistoryMap;
  IDHistoryMap GrainIDHistoryMap;
  float fGrainThresholdAngle;
  std::cout << "Enter minimum anglular (degree) deviation between grains" << std::endl; 
  std::cin >> std::skipws >> fGrainThresholdAngle;
  CubicFZCost CubicFZCost;
  typedef Pandora::GrainDynamics::GrainDynamicsAnalysis<XDM_C3t3, IDOrientMap> GrainGrowthAnalysis;
  GrainGrowthAnalysis GrainGrowthMapper(c3t3_a, c3t3_b, OrientMap_a, OrientMap_b );
  GrainGrowthMapper.GetGrainHistory( GrainIDHistoryMap, CubicFZCost,
                                     DEGREE_TO_RADIAN( fGrainThresholdAngle )  );

  std::string GrainIDFilename;
  std::cout << "Grain ID History filename " << std::endl;
  std::cin  >> GrainIDFilename;
  std::ofstream GrainHistoryOs( GrainIDFilename.c_str() );
  typedef IDHistoryMap::iterator IDMapIter;
  for( IDMapIter it = GrainIDHistoryMap.begin();
       it != GrainIDHistoryMap.end(); ++ it )
  {
    GrainHistoryOs << it->first << " "
                   << it->second.ID << " "
                   << it->second.fMisorient << " "
                   << it->second.IntersectingVolume << " "
                   << it->second.Volume_a << " "
                   << it->second.Volume_b << " "
                   << it->second.nNgb << " "
                   << it->second.MeanWidth << " "
                   << it->second.TripleLine << " "
                   << it->second.q_a << " "
                   << it->second.q_b << " "
                   << it->second.GrainCenter_a << " ";
    
    if ( it->second.bInternal )
      GrainHistoryOs << 1 << std::endl;
    else
      GrainHistoryOs << 0 << std::endl;
    
    //    GrainHistoryOs << std::endl;
  }
  GrainHistoryOs.close();
  
  //----------------------------------------  DEBUG purpose

  std::string VtkName_a, VtkName_b;
  std::cout << "VtkFilename for transformed mesh A " << std::endl;
  std::cin >> VtkName_a;
  
  std::cout << "VtkFilename for transformed mesh B " << std::endl;
  std::cin >> VtkName_b;
  
                                                                   
  typedef Pandora::VtkPropMap VtkPropMap;
  typedef Pandora::VoxelFieldSelector::IdentitySelector   IdentitySelector;
  
  std::vector< VtkPropMap > PointDataMaps;
  std::vector< VtkPropMap > CellDataMaps;
  std::vector< VtkPropMap > FieldDataMaps;

  std::cout << "Writing to vtk " << std::endl;
  std::ofstream vtkFile_a( VtkName_a.c_str() );
  std::ofstream vtkFile_b( VtkName_b.c_str() );
  Pandora::WriteMeshToUnstructuredVtk( vtkFile_a, c3t3_a,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps.begin(), CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  Pandora::WriteMeshToUnstructuredVtk( vtkFile_b, c3t3_b,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps.begin(), CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  
  vtkFile_a.close();
  vtkFile_b.close();
  //----------------------------------------  DEBUG purpose
  return 0;
}
