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
#include "XDM_Mesh_IO.h"
#include "XDM_Mesh/XDM_Data.h"
#include "GrainCurvatureEst.h"

#include "MeshUtilities.h"
#include "GrainGeometryAnalysis.h"
#include "3dMath.h"
#include "ShapeGenerator.h"

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




//---------------------------------------------------------------
//---------------------------------------------------------------
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
struct SetEdgeLength
{
  template< class T1, class T2, class T3 >
  void operator()( T1 & map, const T2 & key, const T3 & value ) const
  {
    map[key].EdgeLength = value;
  }
};

//---------------------------------------------------------------
//  CalculateMeanWidth
//---------------------------------------------------------------
void CalculateMeanWidth( double & fMeanWidth, double & fVolume,
                         double & fEdgeLength,
                         XDM_C3t3 & c3t3, int ID )
{
  //-----------------------------------------------------------------------
  typedef XDM_C3t3::Subdomain_index              Subdomain_index;
  typedef Pandora::Details::GeomProp                GP;
  typedef std::map< Subdomain_index, GP>            GeomPropMap;
  typedef Pandora::Details::GrainGeometry<XDM_C3t3> GrainGeometry;
  typedef GrainGeometry::IDNgbMapT                  IDNgbMapT;
  typedef GrainGeometry::GrainTetMapT               GrainTetMapT;
  typedef GrainGeometry::IDToEdgeMapT               IDToEdgeMapT;
  //-----------------------------------------------------------------------
  //  Build maps
  
  IDToEdgeMapT  IDToEdgeMap;
  IDNgbMapT    NgbMap;
  GrainTetMapT GrainTetMap;
  GeomPropMap  IDToGeomMap;
  std::cout << " | Build surface - edge map | " << std::endl;
  Pandora::Details::BuildGrainSurfaceEdgeMap( IDToEdgeMap, c3t3 );
  std::cout << " | Build grain   - tet map | " << std::endl;
  Pandora::Details::BuildGrainTetMap        ( GrainTetMap, c3t3 );
  std::cout << " | Run geometric extractor | " << std::endl;
  GrainGeometry GeometricExtractor          ( c3t3 );
  std::cout << " | Calculating mean width | " << std::endl;
  GeometricExtractor.CalculateMeanwidth( IDToGeomMap, SetMeanwidth(), IDToEdgeMap );

  std::cout << " | Calculating volume | " << std::endl;
  GeometricExtractor.CalculateVolume   ( IDToGeomMap, SetVolume(),    GrainTetMap );
  std::cout << " | Getting Grain - Ngb map | " << std::endl;;
  Pandora::Details::BuildFiniteGrainNgbMap   ( NgbMap, c3t3 ); 

  std::cout << "Calculating triple line length  " << std::endl;
  GeometricExtractor.CalculateTripleLineLength( IDToGeomMap, SetEdgeLength(), 
                                                IDToEdgeMap );

 
  fMeanWidth  = IDToGeomMap[ ID ].MeanWidth;
  fVolume     = IDToGeomMap[ ID ].Volume;
  fEdgeLength = IDToGeomMap[ ID ].EdgeLength;
  std::cout << "Edge Length = " << fEdgeLength << std::endl;
  // exit(0);
}

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


//---------------------------------------------------------------
//  main
//---------------------------------------------------------------
int main()
{

  int ObjID = 5;
  int nBoxLength = 200;

  GeneralLib::SVector3 Center( 100, 100, 100 );
  GeneralLib::SVector3 PixelScale( 1, 1, 1 );
  
  vector<float> ScalingList;
  vector<float> ApproximationList;

  //  ScalingList.push_back( 10 );
  //  ScalingList.push_back( 5 );
  ScalingList.push_back( 1 );

  //  ApproximationList.push_back( 0.5 );
  ApproximationList.push_back( 1 );
  //  ApproximationList.push_back( 2 );
  // ApproximationList.push_back( 3 );

  //std::ofstream ScalingOutput("Cube_Scaling.txt");
  std::ofstream ScalingOutput;
  
  int nObjectType;
  std::cout << "Enter Object type: " << std::endl;
  std::cin >> std::skipws >> nObjectType;

  std::string OutputFilename;
  std::cout << "Enter Output filename " << std::endl;
  std::cin >> std::skipws >> OutputFilename;
  ScalingOutput.open( OutputFilename.c_str() );
  
  int nRotations = 1;
  int nSteps = 17;
    
  std::cout << "Enter Number of Rotations " << std::endl;
  std::cin >> std::skipws >> nRotations;
  std::cout << "Doing " << nRotations << " Rotation Per mesh parameter " << std::endl;
  std::cout << "Enter Number of Steps " << std::endl;
  std::cin >> std::skipws >> nSteps;
  std::cout << "Doing " << nSteps << " object sizes" << std::endl;
  
  for( int j = 0; j < ScalingList.size(); j ++ )
  {
    for( int k = 0; k < ApproximationList.size() ; k ++ )
    {

      float fMaxTetScaling = ScalingList[j];
      float fMaxTriScaling = ScalingList[j];
      float fApproximationScaling = ApproximationList[k];
      float fRadius = 40;
      
      for(int i = 0; i < nSteps; i ++ )
      {
        
        float fMaxTetSize = (fRadius * fRadius * fRadius / 6. ) / fMaxTetScaling;
        float fMaxTriSize = (fRadius * fRadius / 2.) / fMaxTriScaling;

        
        std::cout << "Step " << i << " "
                  << j << " " << k << " Max Tet " << fMaxTetSize
                  << " " << " Max Facet " << fMaxTriSize << " Radius " << fRadius <<  std::endl;

        GeneralLib::CRandomRotationGenerator RandQuatGen;
        double fAveMeanWidth  = 0;
        double fAveVolume     = 0;
        double fAveEdgeLength = 0;
        vector<double> MeanWidthList;
        vector<double> VolumeList;
        vector<double> EdgeLengthList;
        float fDistanceApproximation = PixelScale.m_fX * fApproximationScaling;


        for( int n = 0; n < nRotations; n ++ )
        {
          XDM_Data DebugData;
          SMatrix3x3 RotMat = RandQuatGen.GetRandomQuaternion().GetRotationMatrix3x3();

          if( nObjectType == 0 )
          {
            GeometricValidation::GenerateCube( DebugData, nBoxLength, nBoxLength, nBoxLength,
                                               PixelScale, Center, fRadius * 2,
                                               RotMat, ObjID );
          }
          else if ( nObjectType == 1)
          {
            GeometricValidation::GenerateSphere( DebugData, nBoxLength, nBoxLength, nBoxLength,
                                                 PixelScale, Center, fRadius, ObjID );
          }
          else if ( nObjectType == 2)
          {
            GeometricValidation::GenerateConstrainedCube( DebugData, nBoxLength, nBoxLength, nBoxLength,
                                                          PixelScale, Center, fRadius * 2,
                                                          RotMat, ObjID );
           
//             std::cin >> std::skipws >> fBndSmoothingLength
//                      >> std::skipws >> fPointSpacing
//                      >> std::skipws >> fMinEdgeLength;
            
//            DebugData.MakeBndPoints( 0, fRadius/4, fRadius/4 ); // no smoothing on purpose
            DebugData.MakeBndPoints( 2, 2, 1 ); // no smoothing on purpose
          }
          else
          {
            std::cout << "Unknown object type" << std::endl;
            exit(0);
          }
          
          std::stringstream MeditSS;
  
          MeditSS << "Cube." << i << "." << j << "." << k << ".medit";
          XDM_C3t3 c3t3 = GenerateMesh( DebugData,
                                        15, fMaxTriSize, fDistanceApproximation,
                                        10, fMaxTetSize,
                                        true, MeditSS.str() );
          double fMeanWidth;
          double fVolume;
          double fEdgeLength;
          std::cout << "Calculating Mean width " << std::endl;
          CalculateMeanWidth( fMeanWidth, fVolume, fEdgeLength, c3t3, ObjID );

          fAveEdgeLength += fEdgeLength;
          fAveVolume    += fVolume;
          fAveMeanWidth += fMeanWidth;
          MeanWidthList.push_back( fMeanWidth );
          VolumeList.push_back   ( fVolume );
          EdgeLengthList.push_back( fEdgeLength );
        }

        fAveVolume /= double( nRotations );
        fAveMeanWidth /= double( nRotations );
        fAveEdgeLength /= double( nRotations );
        double fVolumeStd    = 0;
        double fMeanWidthStd = 0;
        double fEdgeLengthStd = 0;
        for( int n = 0; n < VolumeList.size(); n ++ )
        {
          fVolumeStd     += ( VolumeList[n] - fAveVolume ) * ( VolumeList[n] - fAveVolume );
          fMeanWidthStd  += ( MeanWidthList[n] - fAveMeanWidth ) * ( MeanWidthList[n] - fAveMeanWidth );
          fEdgeLengthStd += ( EdgeLengthList[n] - fAveEdgeLength ) * ( EdgeLengthList[n] - fAveEdgeLength );
        }
        fVolumeStd     = std::sqrt( fVolumeStd )     / float( nRotations );
        fMeanWidthStd  = std::sqrt( fMeanWidthStd )  / float( nRotations );
        fEdgeLengthStd = std::sqrt( fEdgeLengthStd ) / float( nRotations );
        
        
        std::cout << "Done Calculating Mean width " << std::endl;
        ScalingOutput << i << " " << j << " " << k << " "
                      << PixelScale.m_fX << " "
                      << fDistanceApproximation << " "
                      << fMaxTetSize << " "
                      << fMaxTriSize << " "
                      << fRadius << " "
                      << fAveMeanWidth << " "
                      << fAveVolume <<  " "
                      << fMeanWidthStd << " "
                      << fVolumeStd << " "
                      << fAveEdgeLength << " "
                      << fEdgeLengthStd << std::endl;
        fRadius -=   2;
      }
    }
  }
  ScalingOutput.close();
  return 0;
}
