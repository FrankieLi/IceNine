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

#include "BoundaryAnalysis.h"
// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::XDM_test::XDM_Data<K> XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;

typedef std::vector< GeneralLib::SMatrix3x3 > IDOrientMap;

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

//----------------------------------------
//  TransformIDOrientMap
//----------------------------------------
void TransformIDOrientMap( IDOrientMap & Map, const GeneralLib::SMatrix3x3 & R  )
{
  for( int i = 0; i < Map.size(); i ++ )
    Map[i] = R * Map[i];
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
  
  std::cout << "Done Reading " << std::endl;
  using namespace GeneralLib;
  SMatrix3x3 InitialRotation;
  InitialRotation.SetIdentity();
  SVector3   InitialTranslation(0, 0, 0);
  SVector3   InitialStepSize(0, 0, 0);
  int nTestPoints  = 1000;
  int nRandomSteps = 1000;

  //-------------------------------------------------------------------------
  //  Read in transformations
  std::cout << "Define shift for mesh (x y z):" << std::endl;
  std::cin  >> std::skipws >> InitialTranslation.m_fX 
            >> std::skipws >> InitialTranslation.m_fY
            >> std::skipws >> InitialTranslation.m_fZ; 

  std::cout << "Translation of " << InitialTranslation << " in microns " <<  std::endl;


  std::cout << "Define shift for mesh (x y z):" << std::endl;
  std::cin  >> std::skipws >> InitialStepSize.m_fX 
            >> std::skipws >> InitialStepSize.m_fY 
            >> std::skipws >> InitialStepSize.m_fZ; 
  
  std::cout << "Step Size of " << InitialStepSize << " in microns " <<  std::endl;
  SVector3 Eulers;
  std::cout << "Rotation for mesh in Euler angles (Degrees):" << std::endl;

  std::cin  >> std::skipws >> Eulers.m_fX 
            >> std::skipws >> Eulers.m_fY
            >> std::skipws >> Eulers.m_fZ; 
  std::cout << "Rotation of " << Eulers << " in degrees " <<  std::endl;
  Eulers = GeneralLib::DegreeToRadian( Eulers );
  InitialRotation.BuildActiveEulerMatrix( Eulers.m_fX, Eulers.m_fY, Eulers.m_fZ );


  double SO3_Radius;
  std::cout << "Random rotation radius (degrees):  " << std::endl;
  std::cin  >> std::skipws >> SO3_Radius;

  std::cout << " SO3 Raidus (degree ): " << SO3_Radius << std::endl;
  
  std::cout << "Input number of test points " << std::endl;
  std::cin  >> std::skipws  >> nTestPoints;
  std::cout << "Test Points to be used " << nTestPoints << std::endl;
  
  std::cout << "Input number of iterations " << std::endl;
  std::cin  >> std::skipws  >> nRandomSteps;
  std::cout << "Number of iterations " << nRandomSteps << std::endl;

  int nRefinements;
  std::cout << "Number of refinement steps " << std::endl;
  std::cin  >> std::skipws  >> nRefinements;
  std::cout << "Number of refinements " << nRefinements << std::endl;

  float fAngleReduction, fTranslationReduction;
  std::cout << "[Anglular radius, Translation] reduction factors " << std::endl;
  std::cin  >> std::skipws  >> fAngleReduction >> std::skipws >> fTranslationReduction ;
  std::cout << "Reduction Factors [ " << fAngleReduction << " " << fTranslationReduction << " ] " << std::endl;

  
  std::string InitializeMethod; 
  std::cout << "Initialization type " << std::endl;
  std::cin  >> std::skipws >> InitializeMethod;
  std::cout << "Initialization method " << InitializeMethod << std::endl;
  
  //
  //-------------------------------------------------------------------------

  
  Pandora::BoundaryAligner<XDM_C3t3> SurfaceBndRegisterator( c3t3_a, c3t3_b );

  
  if( InitializeMethod.compare( "surface" ) == 0 )
    SurfaceBndRegisterator.InitializeBySurface();
  else if( InitializeMethod.compare( "boundary" ) == 0 )
    SurfaceBndRegisterator.InitializeByBoundary();
  else
  {
    std::cerr << "Error:  Unknown initialization type" << std::endl;
    exit(0);
  }

  string ResetOrigin;
  std::cout << " Reset origin for meshes 1 to (0, 0, 0) ? reset_com/[no_reset] " << std::endl;
  std::cin >> std::skipws >> ResetOrigin;
  
  SVector3 FinalTranslationShift(0, 0, 0);
  if( ResetOrigin.compare( "reset_com" ) == 0  )
  {
    std::cout << "Resetting origins" << std::endl;
    std::pair<SVector3, SVector3> CenterOfMassPair;
    CenterOfMassPair = SurfaceBndRegisterator.CenterOfMasses();
    SMatrix3x3 Identity;
    Identity.SetIdentity();
    Pandora::TransformMesh( c3t3_a, Identity, - CenterOfMassPair.first );

    InitialTranslation   -= CenterOfMassPair.first;
    std::cout << "   Reference Center of Mass " << CenterOfMassPair.first << std::endl;
    std::cout << "   To align mesh  Center of Mass " << CenterOfMassPair.second << std::endl;
    std::cout << "    New Initial Translation " << InitialTranslation << std::endl;
    if( InitializeMethod.compare( "surface" ) == 0 )
      SurfaceBndRegisterator.InitializeBySurface();
    else if( InitializeMethod.compare( "boundary" ) == 0 )
      SurfaceBndRegisterator.InitializeByBoundary();
    
    
    FinalTranslationShift = CenterOfMassPair.first;
  }

  

  
  SVector3 COMTranslation = SurfaceBndRegisterator.GetCenterOfMassTranslation();
  
  string UseTranslate;
  std::cout << "Use Center of mass translation " << COMTranslation << " use_com_trans/[no_com_trans] ? " << std::endl;
  std::cin >> std::skipws >> UseTranslate;
  
  if( UseTranslate.compare( "use_com_trans" ) == 0  )
  {
    InitialTranslation = COMTranslation;
    std::cout << "Using Center of mass" << std::endl;
  }

  double zShift;
  std::cout << " z shift for points to be used in alignment (microns) " << std::endl;
  std::cin >> std::skipws >> zShift;
  std::cout << " Shifting from top and bottom of mesh by " << zShift << std::endl;

  SurfaceBndRegisterator.RemoveEdgePointsByShift( zShift );
  string CostFilename;
  std::cin >> std::skipws >> CostFilename;
  std::cout << "Writing to " << CostFilename << std::endl;
  
  std::string VtkName;
  std::cout << "VtkFilename for transformed mesh " << std::endl;
  std::cin >> VtkName;

  std::ofstream CostFile( CostFilename.c_str() );
  double RefinementBestCost = std::numeric_limits<double>::max();
  for( int i = 0; i < nRefinements; ++ i )
  {
    double fCurrentBestCost;
    std::pair<SMatrix3x3, SVector3> oRes =   SurfaceBndRegisterator.AlignSample( fCurrentBestCost,
                                                                                 InitialRotation, InitialTranslation,
                                                                                 InitialStepSize,
                                                                                 nRandomSteps, nTestPoints,
                                                                                 DEGREE_TO_RADIAN( SO3_Radius), i );
    SO3_Radius      *= fAngleReduction;
    InitialStepSize *= fTranslationReduction;

    std::cout << i << " |  (COM shift included) "
              << FinalTranslationShift  + oRes.second << " |  "
              << GeneralLib::RadianToDegree( oRes.first.GetEulerAngles() )  << " | " 
              << SO3_Radius << " "
              << InitialStepSize << " "
              << fCurrentBestCost << " "
              << RefinementBestCost << std::endl;
    if( fCurrentBestCost < RefinementBestCost )
    {
      RefinementBestCost = fCurrentBestCost;
      InitialRotation    = oRes.first;
      InitialTranslation = oRes.second;

      CostFile << i << " |  (COM shift included) "
               << FinalTranslationShift  + InitialTranslation << " |  "
               << GeneralLib::RadianToDegree( InitialRotation.GetEulerAngles() )  << " | " 
               << SO3_Radius << " "
               << InitialStepSize << " "
               << fCurrentBestCost << std::endl;
      
      std::cout << "Step size "  << InitialStepSize << std::endl;
      std::cout << "SO3_Radius " << SO3_Radius << std::endl;

    }
    else
      std::cout << "not accepted " << std::endl;
  }
  CostFile.close();
  
  Pandora::TransformMesh( c3t3_b, InitialRotation, FinalTranslationShift + InitialTranslation );

  std::ofstream MeshOs( "RotatedMesh.bin",
                        std::ios_base::out | std::ios_base::binary );
  CGAL::set_binary_mode( MeshOs );
  MeshOs << c3t3_b;
  MeshOs.close();
  std::cout << "Done outputting binary mesh " << std::endl;


  IDOrientMap Map_b;
  ReadIDOrientMap( Map_b, OrientMapFilename_b );
  TransformIDOrientMap( Map_b, InitialRotation);
  WriteIDOrientMap(Map_b, "RotationGrainOrientations.txt");
  
  //----------------------------------------  DEBUG purpose

  
                                                                   
  typedef Pandora::VtkPropMap VtkPropMap;
  typedef Pandora::VoxelFieldSelector::IdentitySelector   IdentitySelector;
  
  std::vector< VtkPropMap > PointDataMaps;
  std::vector< VtkPropMap > CellDataMaps;
  std::vector< VtkPropMap > FieldDataMaps;

  std::cout << "Writing to vtk " << std::endl;
  std::ofstream vtkFile( "Rotated.vol.vtk" );
  
  Pandora::WriteMeshToUnstructuredVtk( vtkFile, c3t3_b,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps.begin(), CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  
  vtkFile.close();

  //----------------------------------------  DEBUG purpose


  
//   Pandora::BoundaryAnalysis<XDM_C3t3> BndAnalysis( c3t3_a, c3t3_b );



//   typedef Pandora::BoundaryAnalysis<XDM_C3t3>::BndMotionInfo BndMotionInfo;

//   vector<BndMotionInfo> BndMotionList =  BndAnalysis.MeshFacetMap( 20, DEGREE_TO_RADIAN( 5 ) );
    
  return 0;
}
