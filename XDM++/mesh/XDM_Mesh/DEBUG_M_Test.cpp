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

#include "NJunctionPatchAnalysis.h"

#include "M_EstimateWithIRLS.h"
#include "ShapeGenerator.h"

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::XDM_test::XDM_Data<K> XDM_Data;
typedef CGAL::XDM_mesh_domain<XDM_Data, K> XDM_Mesh_Domain;
typedef CGAL::XDM_mesh_triangulation_3<XDM_Mesh_Domain>::type XDM_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<XDM_Tr> XDM_C3t3;

typedef std::vector< GeneralLib::SMatrix3x3 > IDOrientMap;


typedef CGAL::Mesh_criteria_3<XDM_Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;


//---------------------------------------------------------------
//  GenerateMesh
//---------------------------------------------------------------
XDM_C3t3 GenerateMesh( const XDM_Data & oData,
                       float fAngle, float fFacetSize, float fApproximation,
                       float fRatio, float fTetSize,
                       bool bWriteFile = false,
                       const std::string & MeditFilename = "")
{
  XDM_Mesh_Domain oTestDomain( oData );
  Facet_criteria facet_criteria( fAngle, fFacetSize,
                                 fApproximation);  // angle, size, approximation
  Cell_criteria cell_criteria( fRatio, fTetSize ); // radius-edge ratio, size
  Mesh_criteria criteria( facet_criteria, cell_criteria );
  std::cout << "Generating mesh " << std::endl;
  XDM_C3t3 c3t3 = CGAL::XDM_make_mesh_3<XDM_C3t3>( oTestDomain, criteria );
  std::cout << "Done " << std::endl;
  
  if( bWriteFile )
  {
    std::cout << "Writing to output mesh " << std::endl;
    std::ofstream medit_file( MeditFilename.c_str() );
    c3t3.output_to_medit(medit_file);
    std::cout << "Done outputting " << std::endl;
    medit_file.close();
    
  }
  return c3t3;
  
}


//----


//--------------------------------------
//  BuildHistoryMap
//--------------------------------------

typedef std::map<int, int> HistoryMapT;
void BuildHistoryMap(  HistoryMapT & Map1, HistoryMapT & Map2, 
                       const std::map<int, int> & Mesh1_ToMesh2_IDMap,
                       const std::set<int> & IDSet1,
                       const std::set<int> & IDSet2 )
{
  int n = 0;
  typedef std::set<int>::const_iterator SetIter;

  for( SetIter it = IDSet1.begin(); it != IDSet1.end(); ++ it)
    Map1[ *it ] = -5;
  for( SetIter it = IDSet2.begin(); it != IDSet2.end(); ++ it)
    Map2[ *it ] = -5;
  
  
  typedef HistoryMapT::const_iterator CIter;
  for( CIter it = Mesh1_ToMesh2_IDMap.begin();
       it != Mesh1_ToMesh2_IDMap.end(); ++ it)
  {
    Map1[ it->first  ] = n;
    Map2[ it->second ] = n;
    n++;
  }
}


typedef XDM_C3t3::Cell_iterator                Cell_iterator;
typedef XDM_C3t3::Cell_handle Cell_handle;
typedef std::map<Cell_handle, int> CellIDMapT;
CellIDMapT BuildCellIDMap( const XDM_C3t3 & c3t3_,
                           const HistoryMapT IDToHistoryMap )
{
  CellIDMapT IDMap;
  for( Cell_iterator it = c3t3_.cells_begin();
       it != c3t3_.cells_end(); ++ it )
  {
    int nID = c3t3_.subdomain_index( it );
    if( IDToHistoryMap.find( nID ) != IDToHistoryMap.end() )
      IDMap[ it ] = IDToHistoryMap.find( nID )->second;
    else
      IDMap[it]   = -5;
  }
  return IDMap;
}


typedef std::map<Cell_handle, SVector3> CellToSVector3Map;

//---------------------------------------------------
//  Return orientation mapped to RGB color
//---------------------------------------------------
CellToSVector3Map GetRodriguzVectorMap( const XDM_C3t3& c3t3_,
                                        const IDOrientMap & OrientMap )
{
  CellToSVector3Map CellToOrientMap;
  for( Cell_iterator it = c3t3_.cells_begin();
       it != c3t3_.cells_end(); ++ it )
  {
    int nID = c3t3_.subdomain_index( it );
    if( nID < 0)
      CellToOrientMap[ it ] = SVector3( 0, 0, 0 );
    else if ( nID >= OrientMap.size() )
      std::cerr << "Average orientation file and mesh are mismatching.  ID exceeds size." << std::endl;
    else
    {
      SQuaternion q;
      q.Set( OrientMap[ nID ] );
      SVector3 rod( q.m_fX, q.m_fY, q.m_fZ );
      rod /= q.m_fW;
      // ToColor

      const float fColorEdge = (sqrt( 2.0 ) - 1) + 0.001;
      rod += SVector3(  fColorEdge, fColorEdge, fColorEdge );
      rod /= ( 2.0 * fColorEdge );
      CellToOrientMap[ it ] = rod;
    }
  }
  return CellToOrientMap;
}

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
  std::cout << Map.size() << " number of lines read from euler angle file " << std::endl;
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
  XDM_C3t3 c3t3;

//   std::string MeshFilename;
//   std::cout << "Enter binary filename for mesh saving: " << std::endl;
//   std::cin  >> MeshFilename;

//   std::string OrientMapFilename;
//   std::cout << "Enter ID to orientation map filename: " << std::endl;
//   std::cin  >> OrientMapFilename;

//   std::ifstream MeshIS( MeshFilename.c_str(),
//                         std::ios_base::in|std::ios_base::binary );
//   CGAL::set_binary_mode( MeshIS );
  
//   MeshIS >> c3t3;
//   MeshIS.close();

  int ObjID = 5;
  int nBoxLength = 200;

  GeneralLib::SVector3 Center( 100, 100, 100 );
  GeneralLib::SVector3 PixelScale( 1, 1, 1 );
  float fRadius = 40;
  IDOrientMap OrientMap;
  GeneralLib::SMatrix3x3 Mat;

  Mat.SetIdentity();
  OrientMap.push_back( Mat );
  
//   ReadIDOrientMap( OrientMap, OrientMapFilename );

  std::cout << "Generate Sphere " << std::endl;
  XDM_Data DebugData;
  GeometricValidation::GenerateSphere( DebugData, nBoxLength, nBoxLength, nBoxLength,
                                       PixelScale, Center, fRadius, ObjID );

  
  float fMaxTetSize = 10;
  float fMaxTriSize = 6;
  float fDistanceApproximation = 1;
  
  c3t3 = GenerateMesh( DebugData,
                       15, fMaxTriSize, fDistanceApproximation,
                       10, fMaxTetSize,
                       true, "DEBUG_M_Test.medit" );
  
  //-----------------------------------------------------------------
  //  DEBUG M Estimates
  MEstimateIRLS::RobustCurvatureEstimate<XDM_C3t3> RCE;

  //   MeanCurvatureNormOp( it->first->vertex( k ),
//                        Mesh1_, Mesh1_.surface_index( *it ),
//                        RunDebug );

  typedef typename XDM_C3t3::Triangulation::Vertex                     Vertex;
  typedef typename XDM_C3t3::Vertex_handle                      Vertex_handle;
  typedef typename XDM_C3t3::Surface_index                      Surface_index;
  //  std::vector< Vertex_handle > VertexList;
 
  // Vertex_handle vh = c3t3.triangulation().vertices_begin();
//   Surface_index ID;
//   RCE.M_EstimateCurvature( vh, c3t3,
//                            ID,
//                            VertexList );
  
  //-----------------------------------------------------------------
  Pandora::NJunctionBoundaryExtractor<XDM_C3t3> SurfaceExtractor( c3t3, OrientMap );
  typedef Pandora::NJunctionBoundaryExtractor<XDM_C3t3>::VertexCurveMap   VertexCurveMapT;

   typedef Pandora::NJunctionBoundaryExtractor<XDM_C3t3>::SurfaceIndexBoundaryMap SurfaceIndexBoundaryMapT;
   typedef SurfaceIndexBoundaryMapT::iterator ID2BndIterator;
   SurfaceIndexBoundaryMapT TestMap;
   //   Pandora::AreaAveragedCurvature<XDM_C3t3, ID2BndIterator, VertexCurveMapT > ParallelCurveOp( TestMap.begin(), TestMap.end() );
//   ParallelCurveOp.GetPatchAveragedCurvature( 5, 0.02);
  
//   std::cout << "Extract boundary " << std::endl;
//   SurfaceExtractor.CalculateCurvature();
   std::cout << " Run M Estimate Curvature " << std::endl;
   //SurfaceExtractor.CalculateCurvature();
   SurfaceExtractor.MEstimateCurvature();

   std::cout << " Finished M Estimate  Curvature " << std::endl;
   SurfaceExtractor.PrintCurvatureMap("DEBUG_TestCurvatureMap.txt");
   SurfaceExtractor.PrintFacetMap("DEBUG_TestFacetMap.txt");
//   std::cout << "Finished extract boundary " << std::endl;
  
  return 0;
}
