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

// #include "3dMath.h"
// #include "SimpleMeshVTKWriter.h"
// #include "VTKPropertyMap.h"
#include "GrainCurvatureEst.h"
#include "XDM_Mesh_IO.h"
// #include "ProcessMic.h"
// #include "Symmetry.h"
// #include <boost/shared_ptr.hpp>

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

// Mesh Criteria
typedef unsigned char XDMIndexT;
using std::vector;
typedef vector< vector< vector< XDMIndexT > > >  DataMatrixT;

typedef CGAL::XDM_test::XDM_Data XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;

typedef CGAL::Mesh_criteria_3<XDM_Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;



int main()
{

  //---------------------------------------------------------------
  //            Changing this to reading a mesh soon
  //
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
  oTestData.MakeBndPoints();
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
  Mesh_criteria criteria(facet_criteria, cell_criteria);
  
  std::cout << "started meshing " << std::endl;
  XDM_C3t3 c3t3 = CGAL::XDM_make_mesh_3<XDM_C3t3>( oTestDomain, criteria );

  //            
  //---------------------------------------------------------------
  
  typedef Pandora::CurveEstimator<XDM_C3t3> CurveEstimator;

  CurveEstimator Estimator( c3t3 );

  std::cout << "Enter Curve VTK file name: " << std::endl;
  std::string VtkFilename;
  std::cin >> VtkFilename;
  std::cout << "Writing to " << VtkFilename << std::endl;
  Estimator.WriteGrainCurveVtkFile( VtkFilename );
  return 0;
}
