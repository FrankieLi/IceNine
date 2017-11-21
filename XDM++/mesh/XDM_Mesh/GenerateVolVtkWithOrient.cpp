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

typedef std::map< int, int > ID2IDMap;
typedef std::map< int, GeneralLib::SVector3 > ID2OrientMap;

//------------------------------------------------
//   ReadID2IDMap
//
//------------------------------------------------
void ReadID2IDMap( ID2IDMap & Map, const string & filename )
{
  std::ifstream MapInput( filename.c_str() );
  
  if( ! MapInput.is_open() )
  {
    std::cout << "ID2IDMap filename not found  " << std::endl;
    exit( 0 );
  }

  while( !MapInput.eof() )
  {
    int nOldID, nNewID;
    MapInput >> std::skipws >> nOldID
             >> std::skipws >> nNewID >> std::skipws;
    if( !MapInput.eof() )
      Map[ nOldID ] = nNewID;
  }
  MapInput.close();
}



//------------------------------------------------
//   ReadID2IDMap
//
//------------------------------------------------
void ReadID2OrientMap( ID2OrientMap & Map, const string & filename )
{
  std::ifstream MapInput( filename.c_str() );
  
  if( ! MapInput.is_open() )
  {
    std::cout << "ID2OrientMap filename not found  " << std::endl;
    exit( 0 );
  }

  int nID = 0;
  while( !MapInput.eof() )
  {
    GeneralLib::SVector3 v;
    MapInput >> std::skipws >> v.m_fX
             >> std::skipws >> v.m_fY
             >> std::skipws >> v.m_fZ >> std::skipws;
    if( !MapInput.eof() )
      Map[ nID ] = v;
  }
  MapInput.close();
}

template<class T, class U>
struct CompositSelector
{
  const Cell2OrientMap & _Map;
  CompositSelector(const CellID2Map & Map): _Map( Map )
  {
  }

  GeneralLib::SVector3 operator()( const std::pair<T, U> & o ) const
  {
    return _Map[ o.second ];
  }
};

//----------------------------------
//
//----------------------------------
template< class C3T3, class Cell2IDMap, class ID2IDMap >
void BuildCellToIDMap( Cell2IDMap & Map, const C3T3 & c3t3_,
                       const ID2IDMap & IDMap )
{
  typedef typename C3T3::Cell_iterator Iter;
  typedef typename C3T3::Cell_handle   CH;
  for( Iter ci = c3t3_.cells_begin(); ci != c3t3_.cells_end(); ++ ci )
  {
    CH ch = ci;

    if( IDMap.find( c3t3_.subdomain_index( ch ) ) != IDMap.end() )
      Map[ ch ] = IDMap.find( c3t3_.subdomain_index( ch ) )->second;
    else
      Map[ ch ] = -10;
  }
}

int main()
{
  std::string MeshFilename;
  std::cout << "Enter binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename;

  std::string ID2IDFilename;
  std::cout << "Enter filename for GrainID to NewID Map: " << std::endl;
  std::cin  >> ID2IDFilename;
  
  std::cout << "Enter Output vtk name: " << std::endl;
  std::string OutputVtkName;
  std::cin >> OutputVtkName;
 
  
  //-------------------------------------------------------------------
  //      Start computing and meshing
  //-------------------------------------------------------------------
  std::cout << "Start reading mesh  " << std::endl;
  std::ifstream MeshIs( MeshFilename.c_str(),
                        std::ios_base::in|std::ios_base::binary);
  if( ! MeshIs.is_open() )
  {
    std::cout << "Input mesh b file not found " << std::endl;
    exit( 0 );
  }
  CGAL::set_binary_mode( MeshIs );
  XDM_C3t3 c3t3;
  
  MeshIs >> c3t3;
  MeshIs.close();


  std::cout << " Done reading binary mesh " << std::endl;
  std::cout << " Start MicVolume calculation and alignment " << std::endl;

  typedef XDM_C3t3::Subdomain_index              Subdomain_index;
  typedef XDM_C3t3::Cell_handle                  Cell_handle;
  typedef XDM_C3t3::Triangulation::Vertex_handle Vertex_handle;
  typedef XDM_C3t3::Cell_iterator                Cell_iterator;

  typedef std::map< Cell_handle, Subdomain_index > CellToIndexMap;

  std::map<Cell_handle, int> CellToNewIDMap;
  std::map<int, int> OldIDToNewIDMap;
  ReadID2IDMap( OldIDToNewIDMap, ID2IDFilename );
  BuildCellToIDMap( CellToNewIDMap, c3t3, OldIDToNewIDMap );
  
  //-----------------------------------------------------------------------
  //  Printing
  typedef Pandora::VtkPropMap VtkPropMap;
  typedef Pandora::VoxelFieldSelector::IdentitySelector   IdentitySelector;
  typedef Pandora::HandleFieldMap< Cell_handle, int, Cell_iterator, IdentitySelector > CellToNewIDMapT;
  std::vector< VtkPropMap > PointDataMaps;
  std::vector< VtkPropMap > CellDataMaps;
  std::vector< VtkPropMap > FieldDataMaps;

//   CellDataMaps.push_back(  VtkPropMap ( CellToNewIDMapT( "SCALARS", "HistoryIndex", "int 1", 1,
//                                                          CellToNewIDMap, IdentitySelector(),
//                                                          c3t3.cells_begin(), c3t3.cells_end() )   ) );
  
  FieldDataMaps.push_back(  VtkPropMap ( CellToEulerMap( "SCALARS", "Euler_Angles", "float", 3,
                                                         CellToEulerMap, ID2EulerSelector(),
                                                         c3t3.cells_begin(), c3t3.cells_end() ) ) );
  

  std::cout << "Writing to " << OutputVtkName << std::endl;
  std::ofstream vtkFile( OutputVtkName.c_str() );
  Pandora::WriteMeshToUnstructuredVtk( vtkFile, c3t3,
                                       PointDataMaps.begin(), PointDataMaps.end(),
                                       CellDataMaps.begin(), CellDataMaps.end(),
                                       FieldDataMaps.begin(), FieldDataMaps.end() );
  vtkFile.close();
  //-----------------------------------------------------------------------

  return 0;
}

