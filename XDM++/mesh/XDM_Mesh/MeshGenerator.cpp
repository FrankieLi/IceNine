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
#include "XDM_Mesh_IO.h"
#include "XDM_Data.h"
#include "GrainCurvatureEst.h"

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





//---------------------------------------------------
// PrintCellIDToGrainIDMap
//---------------------------------------------------
void PrintCellIDToGrainIDMap( const string & OutFilename,  XDM_C3t3 & c3t3 )
{
  CGAL::Default_cell_index_pmap<XDM_C3t3>   oCellToID( c3t3 );
  typedef XDM_C3t3::Cell_iterator CellIterator;
  std::map<int, int> UniqueMap;
  for(XDM_C3t3::Cell_iterator cit = c3t3.cells_begin(); cit != c3t3.cells_end();
      ++ cit )
  {
    UniqueMap[ c3t3.subdomain_index( cit ) ] = get( oCellToID, cit );
  }


  std::ofstream os;
  os.open( OutFilename.c_str() );
  typedef std::map<int, int>::iterator MIter;
  for( MIter iter = UniqueMap.begin(); iter != UniqueMap.end(); ++ iter )
    os << iter->first << " " << iter->second << std::endl;
  os.close();
}

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
  oTestData.MakeTopologyPreservingPoints();
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
  Mesh_criteria criteria( facet_criteria, cell_criteria );

  std::string MeshFilename;
  std::cout << "Enter binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename;
  std::cout << MeshFilename << std::endl;
  //-------------------------------------------------------------------
  //      Start computing and meshing
  //-------------------------------------------------------------------
  std::cout << "started meshing " << std::endl;
  XDM_C3t3 c3t3 = CGAL::XDM_make_mesh_3<XDM_C3t3>( oTestDomain, criteria );
  std::cout << "Done Meshing " << std::endl;
  
  std::ofstream MeshOs( MeshFilename.c_str(),
                        std::ios_base::out | std::ios_base::binary );
  CGAL::set_binary_mode( MeshOs );
  MeshOs << c3t3;
  MeshOs.close();
  std::cout << "Done outputting binary mesh " << std::endl;

  //----------------------------------------------
  //  Traditional debug output - to be phased out
  //----------------------------------------------
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
  std::ofstream BndIDToGrainIDFile( "BndID_2_GrainID_File.txt");
  Pandora::PrintBndIDToGrainID( BndIDToGrainIDFile, c3t3 );
  PrintCellIDToGrainIDMap("GrainID2CellIDMap.txt", c3t3);
  BndIDToGrainIDFile.close();

  return 0;
}
