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
#include "MeshUtilities.h"


// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::XDM_test::XDM_Data<K> XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;
typedef std::map< int, GeneralLib::SQuaternion > IDOrientMap;

// typedef GeneralLib::SVector3   SVector3;
// typedef GeneralLib::SQuaternion SQuaternion;

struct EdgeErrorProp
{
  EdgeErrorProp() {}
  EdgeErrorProp( float fMis, int nNgbs_, bool b ) :
    fMisorientError( fMis ),
    nNgbs( nNgbs_ ),
    bEdgeIncidentsSpace ( b ) {}
  
  float fMisorientError;
  int   nNgbs;
  bool  bEdgeIncidentsSpace;
};

struct EdgeOrientMapper
{
  template< class EdgeToErrorMapT, class EdgeIter, class DomainSetT >
  void operator()(  EdgeToErrorMapT   & ResultMap,
                    const EdgeIter    & pEdge,
                    const GeneralLib::SQuaternion & qTot,
                    bool bEdgeIncidentsSpace,
                    const DomainSetT  & DSet ) const
  {
    GeneralLib::SQuaternion qIdentity;
    qIdentity.Set( 1, 0, 0, 0 );
    float fMis = LatticeSymmetry::GetMisorientation( LatticeSymmetry::CCubicSymmetry::Get(),
                                                     qTot, qIdentity );
    
    
    EdgeErrorProp Res( fMis, DSet.size(), bEdgeIncidentsSpace );
    ResultMap.insert( std::make_pair( *pEdge, Res ) );
  }

  template< class IDToErrMap >
  void operator() ( IDToErrMap & Map, int ID,
                    const GeneralLib::SQuaternion qTot )
  {
    float fMis = float(2) * std::acos( std::min( float(1), qTot.m_fW ) ); 
    Map.insert( std::make_pair( ID, fMis )  );
  }
};


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
  std::string MeshFilename;
  std::cout << "Enter first binary filename for mesh saving: " << std::endl;
  std::cin  >> MeshFilename;

  std::string OrientMapFilename;
  std::cout << "Enter first ID to orientation map filename: " << std::endl;
  std::cin  >> OrientMapFilename;

  XDM_C3t3 c3t3;
  std::ifstream MeshIS( MeshFilename.c_str(),
                        std::ios_base::in|std::ios_base::binary );
  CGAL::set_binary_mode( MeshIS );
  if( ! MeshIS.is_open() )
  {
    std::cout << "Input mesh a file not found " << std::endl;
    exit( 0 );
  }

  MeshIS >> c3t3;
  MeshIS.close();
  
  IDOrientMap OrientMap;
  ReadEulerIDMap( OrientMap, OrientMapFilename );
  std::cout << "Done reading input file " << std::endl;

  EdgeOrientMapper Mapper;
  std::map< int,  float > IDToTotalMisorientMap;
  typedef XDM_C3t3::Triangulation::Finite_edges_iterator    Edge_iterator;
  typedef XDM_C3t3::Triangulation::Edge    Edge;
  std::map< Edge, EdgeErrorProp > EdgeToErrorMap;
  Pandora
    ::Details
    ::CalculateCircularOrientationError( IDToTotalMisorientMap,
                                         EdgeToErrorMap, c3t3,
                                         Mapper, OrientMap );

  std::string ErrorFilename;
  std::cout << "Enter output error filename: " << std::endl;
  std::cin  >> ErrorFilename;
  std::ofstream ErrorFile( ErrorFilename.c_str() );
  typedef std::map< Edge, EdgeErrorProp>::iterator EdgePropIter;

  typedef XDM_C3t3::Triangulation::Geom_traits::Bare_point Point3;

  for( EdgePropIter pCur = EdgeToErrorMap.begin();
       pCur != EdgeToErrorMap.end(); ++ pCur )
  {
    Point3 p1  = pCur->first.get<0>( )->vertex( pCur->first.get<1>() )->point();
    Point3 p2  = pCur->first.get<0>( )->vertex( pCur->first.get<2>() )->point();
    
    ErrorFile << pCur->second.nNgbs << " "
              << pCur->second.fMisorientError << " "
              << pCur->second.bEdgeIncidentsSpace << " "
              << p1 << " "
              << p2 << std::endl;
  }
  ErrorFile.close();
  return 0;
}
