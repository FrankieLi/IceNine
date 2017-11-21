//----------------------------------------------------
//  NJunctionPatchAnalysis.h
//
//  Author:   Frankie Li
//  Purpose:  Extract boundary patches delinieated by n-junction 1-features
//            (n-junction lines, or in most cases triple lines).  This is
//            essentially a coarse grained analysis.
//            
//----------------------------------------------------


#ifndef N_JUNCTION_PATCH_ANALYSIS_H
#define N_JUNCTION_PATCH_ANALYSIS_H

#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <limits>
#include <CGAL/bounding_box.h>
#include "MeshUtilities.h"
#include <CGAL/Cartesian_converter.h>
#include "Symmetry.h"
#include "MeanCurvatureNormal.h"

#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Lapack/Linear_algebra_lapack.h>

#include "M_EstimateWithIRLS.h"
#include <Eigen/Dense>

//----------------------------------
//  Forward Declerations
//----------------------------------
namespace Pandora
{
  template< class C3T3, class Iterator, class VertexCurveMap >
  class AreaAveragedCurvature;
  
  template< class Vector_3, class RT >
  struct BoundaryVertexInfo;

  template< class VertexCurveMap, class FacetSet,  class RT, class Vector_3, class Point3>
  struct BoundaryPatchInfo;

  template< class C3T3 >
  class NJunctionBoundaryExtractor;
}

namespace Pandora
{

  //--------------------------------
  //  BndVertexInfo
  //--------------------------------
  template< class Vector_3, class RT >
  struct BoundaryVertexInfo
  {
    BoundaryVertexInfo( const Vector_3 & c, double a ):
      Curvature( c ), fMixedArea( a ) {}
    Vector_3       Curvature;
    RT             fCurvature;
    RT             fMixedArea;
    Vector_3       M_EstimateResidual;
    int            N_SampledPoints;
    
  };

  //--------------------------------
  //  BoundaryPatch
  //--------------------------------
  template< class VertexCurveMap, class FacetSet,  class RT, class Vector_3, class Point3>
  struct BoundaryPatchInfo
  {
    BoundaryPatchInfo()
      : CurvatureMap(), Facets(), TotalArea(0), TotalMixedArea(0) {}
    VertexCurveMap  CurvatureMap;
    FacetSet        Facets;
    RT              TotalArea;
    RT              TotalMixedArea;
    Vector_3        MeanCurvature;
    Vector_3        WeightedMeanCurvature;
    RT              LeastSquare_K1;
    RT              LeastSquare_K2;
    RT              LeastSquare_ConditionNum;
    Point3          FitOrigin;
    bool            Fitted;
    
    RT              fMeanCurvature;
    RT              fWeightedMeanCurvature;
    int             NumSurfacePoints;
  };

  
  template< class C3T3 >
  class NJunctionBoundaryExtractor
  {
  public:

    typedef typename Eigen::Matrix< Float, 3, 1 >       Vector3;
    typedef typename C3T3::Triangulation::Vertex_iterator          Vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_facets_iterator   Finite_facet_iterator;

    typedef typename C3T3::Triangulation::Vertex                   Vertex;
    typedef typename C3T3::Triangulation::Geom_traits::Bare_point  Point3;
    typedef typename C3T3::Triangulation::Geom_traits::RT          RT;
    typedef typename C3T3::Cell_handle                             Cell_handle;
    typedef typename C3T3::Facet                                   Facet;
    typedef typename C3T3::Vertex_handle                           Vertex_handle;

    typedef typename C3T3::Surface_index                           Surface_index;
    typedef typename C3T3::Subdomain_index                         Subdomain_index;
    typedef typename std::multimap<Subdomain_index, Subdomain_index>  CrossMeshIDMap;  // maps from grain in mesh 1 to mesh 2
    typedef typename std::pair<Subdomain_index, Subdomain_index>      IDPair;

    typedef typename C3T3::Triangulation::Locate_type              Locate_type;
    typedef typename C3T3::Triangulation::Geom_traits::Triangle_3  Triangle_3;
    typedef typename C3T3::Triangulation::Geom_traits::Plane_3     Plane_3;
    typedef typename C3T3::Triangulation::Geom_traits::Vector_3    Vector_3;
    typedef typename std::set<Vertex_handle>                       VertexSet;
    typedef typename std::set<Facet>                               FacetSet;
    
    // for fitting
    struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
    typedef typename CGAL::Monge_via_jet_fitting<K>   Monge_via_jet_fitting;
    
  
    
    typedef BoundaryVertexInfo< Vector_3, RT >                BndVertexInfo;  
    typedef typename std::map<Vertex_handle, BndVertexInfo>       VertexCurveMap;
    typedef BoundaryPatchInfo< VertexCurveMap,  FacetSet,  RT, Vector_3, Point3>  BoundaryPatch;

    typedef std::map< Surface_index, BoundaryPatch >    SurfaceIndexBoundaryMap;
    typedef typename SurfaceIndexBoundaryMap::iterator  SurfaceIndexBndMapIter;
    
    NJunctionBoundaryExtractor( const C3T3 & Mesh_,
                                std::vector<SMatrix3x3> & IDOrientMap_  ):
      Mesh( Mesh_ ), 
      IDOrientMap( IDOrientMap_ )
    {
      NJunctionSurface = ExtractBoundary();
    }

    //--------------------------------
    //  GetMajoritySign
    //
    //  - this needs to be checked
    //--------------------------------
    int GetMajoritySign( const Vertex_handle & ReferenceVertex,
                         const Vector_3 & K,
                         const Subdomain_index & RefID  )
    {
      std::vector<Facet> FacetList;
      Mesh.triangulation().finite_incident_facets( ReferenceVertex,
                                                   std::back_inserter( FacetList ));  // this might not be correct, may need filter      
      int nCount[2];
      
      nCount[0] = nCount[1] = 0;
      for( int i = 0; i < FacetList.size(); i ++ )
      {
        if ( Pandora::Details::CorrectCurvatureSign( FacetList[i],
                                                     ReferenceVertex->point(),
                                                     K,  RefID, Mesh ) > 0 )
          nCount[1] ++;
        else
          nCount[0] ++;
      }

      if( nCount[0] > nCount[1] )
        return -1;
      else
        return 1;
    }
    
    //--------------------------------
    //  PrintCurvatureMap
    //--------------------------------
    void PrintCurvatureMap( const std::string & filename )
    {
      std::ofstream os( filename.c_str() );
      for( SurfaceIndexBndMapIter PatchIter = NJunctionSurface.begin();
           PatchIter != NJunctionSurface.end(); ++ PatchIter  )
      {
        if( PatchIter->first.first != -1 && PatchIter->first.second != -1 )
        {
          os << PatchIter->first.first  << " "
             << PatchIter->first.second << " "
             << PatchIter->second.TotalArea << " "
             << PatchIter->second.TotalMixedArea << " "
             << PatchIter->second.MeanCurvature << " "
             << PatchIter->second.WeightedMeanCurvature << " "
             << PatchIter->second.LeastSquare_K1 << " "
             << PatchIter->second.LeastSquare_K2 << " "
             << PatchIter->second.LeastSquare_ConditionNum << " "
             << PatchIter->second.fMeanCurvature << " "
             << PatchIter->second.fWeightedMeanCurvature << " "
             << PatchIter->second.Facets.size() << " "
             << PatchIter->second.CurvatureMap.size() << " "
             << PatchIter->second.NumSurfacePoints << std::endl;
        }
      }
      os.close();
    }
    
    //--------------------------------
    //  PrintFacetMaps
    //--------------------------------
    void PrintFacetMap( const std::string & filename )
    {
      std::ofstream os( filename.c_str() );
      int nPatchID = 0;
      for( SurfaceIndexBndMapIter PatchIter = NJunctionSurface.begin();
           PatchIter != NJunctionSurface.end(); ++ PatchIter  )
      {
        if( PatchIter->first.first != -1 && PatchIter->first.second != -1 )
        {
          typedef typename FacetSet::iterator FIter;
          
          for( FIter it = PatchIter->second.Facets.begin();
               it != PatchIter->second.Facets.end(); ++ it )
          {
            os << nPatchID <<  " "  << PatchIter->first.first << " " << PatchIter->first.second << " ";
            for( int i = 0; i < 4; i ++ )
            {
              if( i != it->second )
              {
                os << it->first->vertex(i)->point().x() << " "
                   << it->first->vertex(i)->point().y() << " "
                   << it->first->vertex(i)->point().z() << " "
                   << PatchIter->second.CurvatureMap.find( it->first->vertex(i) )->second.fCurvature << " "
                   << PatchIter->second.CurvatureMap.find( it->first->vertex(i) )->second.M_EstimateResidual << " "
                  //<< PatchIter->second.CurvatureMap.find( it->first->vertex(i) )->second.Curvature << " "
                  //<< PatchIter->second.CurvatureMap.find( it->first->vertex(i) )->second.fMixedArea << " ";
                   << PatchIter->second.CurvatureMap.find( it->first->vertex(i) )->second.N_SampledPoints << " ";
                if( Pandora::Details::IsSurfaceVertex( it->first->vertex(i), Mesh ) )
                  os << " 1 ";
                else
                  os << " 0 ";

                if( Pandora::Details::IsSurfaceVertex( it->first->vertex(i), Mesh )
                    && PatchIter->second.CurvatureMap.find( it->first->vertex(i) )->second.N_SampledPoints == -1000 )
                {
                  std::cout << "Exception " << PatchIter->first.first << " " << PatchIter->first.second << std::endl;
                }
              }
            }

            os << std::endl;
          }
        }
        nPatchID ++;
      }
      os.close();
    }



    //--------------------------------
    //  MEstimateCurvature
    //--------------------------------
    void MEstimateCurvature()
    {
      std::cout << "Running MEstimate Curvature calculation - this will take a while " << std::endl;
      Pandora::Details::MeanCurvatureNormal<C3T3> MeanCurvatureNormalOp;
      MEstimateIRLS::RobustCurvatureEstimate<C3T3> RCE;
      
      typedef typename VertexCurveMap::iterator CurveMapIter;
      // iterator through all patches
      int nPatches = 0;
      for( SurfaceIndexBndMapIter PatchIter = NJunctionSurface.begin();
           PatchIter != NJunctionSurface.end(); ++ PatchIter  )
      {
        
        if( nPatches % 1000 == 0 )
          std::cout << nPatches << "/" << NJunctionSurface.size() << std::endl;
        nPatches ++;
        if( PatchIter->first.first != -1 && PatchIter->first.second != -1 )   // if not exterior surface
        {
          PatchIter->second.TotalMixedArea = 0;
          PatchIter->second.MeanCurvature  = Vector_3( Point3( 0, 0, 0 ), Point3( 0, 0, 0 ) );
          PatchIter->second.WeightedMeanCurvature   = Vector_3( Point3( 0, 0, 0 ), Point3( 0, 0, 0 ) );
          PatchIter->second.fMeanCurvature          = 0;
          PatchIter->second.fWeightedMeanCurvature  = 0;

          
          int nSurfacePoints = 0;
          for( CurveMapIter vit = PatchIter->second.CurvatureMap.begin();
               vit != PatchIter->second.CurvatureMap.end(); ++ vit )
          {
            if( Pandora::Details::IsSurfaceVertex( vit->first, Mesh ) )  // ignoring all n-junction lines
            {
              
              std::pair<Vector_3, double> oRes = MeanCurvatureNormalOp( vit->first, Mesh, PatchIter->first );
              
              int CurvatureSign =  GetMajoritySign( vit->first, oRes.first, PatchIter->first.first );
              oRes.first = CurvatureSign * oRes.first;
              
              PatchIter->second.TotalMixedArea += oRes.second;
              PatchIter->second.MeanCurvature   =  PatchIter->second.MeanCurvature + oRes.first;
              PatchIter->second.WeightedMeanCurvature
                =  PatchIter->second.WeightedMeanCurvature + oRes.first * oRes.second;
              
              PatchIter->second.fMeanCurvature         += CurvatureSign
                                                        * CGAL::sqrt( oRes.first.squared_length() );  
              PatchIter->second.fWeightedMeanCurvature += CurvatureSign * oRes.second * CGAL::sqrt( oRes.first.squared_length() );  
              
              vit->second.Curvature  = oRes.first;
              //    vit->second.fCurvature =  -1000;
              Vector3 ResidualErr;
              
              RCE.M_EstimateCurvature(vit->second.fCurvature,
                                      vit->second.N_SampledPoints,
                                      ResidualErr,
                                      vit->first, Mesh, PatchIter->first ); 

              vit->second.M_EstimateResidual = Vector_3( ResidualErr(0),
                                                         ResidualErr(1),
                                                         ResidualErr(2) );
              vit->second.fMixedArea = oRes.second;
              nSurfacePoints ++;
            }
            else
            {
              vit->second.fCurvature      =  -1000;
              vit->second.N_SampledPoints =  -1000;
              //              std::cout << "reject " << std::endl;
            }
            
          }
        }
      }
      
    }
    
    //--------------------------------
    //  CalculateCurvature
    //--------------------------------
    void CalculateCurvature()
    {
      Pandora::Details::MeanCurvatureNormal<C3T3> MeanCurvatureNormalOp;


      //      std::ofstream DebugCurvatureFile("DebugCurvature.csv");
      
      typedef typename Monge_via_jet_fitting::Monge_form   Monge_form;
      typedef typename VertexCurveMap::iterator CurveMapIter;
      // iterator through all patches
      for( SurfaceIndexBndMapIter PatchIter = NJunctionSurface.begin();
           PatchIter != NJunctionSurface.end(); ++ PatchIter  )
      {
        // std::cout << "Running Junfciton Surface  " << std::endl;
        if( PatchIter->first.first != -1 && PatchIter->first.second != -1 )
        {
          //          std::cout << "-------------------------" << std::endl;
          std::vector<Point3> PointList;
          PatchIter->second.TotalMixedArea = 0;
          PatchIter->second.MeanCurvature  = Vector_3( Point3( 0, 0, 0 ), Point3( 0, 0, 0 ) );
          PatchIter->second.WeightedMeanCurvature   = Vector_3( Point3( 0, 0, 0 ), Point3( 0, 0, 0 ) );
          PatchIter->second.fMeanCurvature          = 0;
          PatchIter->second.fWeightedMeanCurvature  = 0;

          int nSurfacePoints = 0;
          for( CurveMapIter vit = PatchIter->second.CurvatureMap.begin();
               vit != PatchIter->second.CurvatureMap.end(); ++ vit )
          {
            if( Pandora::Details::IsSurfaceVertex( vit->first, Mesh ) )  // ignoring all n-junction lines
            {              
              std::pair<Vector_3, double> oRes = MeanCurvatureNormalOp( vit->first, Mesh, PatchIter->first );

              //              DebugCurvatureFile << oRes.first << " " << vit->first->point() << std::endl;
              
              int CurvatureSign =  GetMajoritySign( vit->first, oRes.first, PatchIter->first.first );
              oRes.first = CurvatureSign * oRes.first;
              
              PatchIter->second.TotalMixedArea += oRes.second;
              PatchIter->second.MeanCurvature   =  PatchIter->second.MeanCurvature + oRes.first;
              PatchIter->second.WeightedMeanCurvature
                =  PatchIter->second.WeightedMeanCurvature + oRes.first * oRes.second;

              PatchIter->second.fMeanCurvature         += CurvatureSign
                                                        * CGAL::sqrt( oRes.first.squared_length() );  
              PatchIter->second.fWeightedMeanCurvature += CurvatureSign
                                                        * oRes.second
                                                        * CGAL::sqrt( oRes.first.squared_length() );  
              vit->second.Curvature  = oRes.first;
              vit->second.fCurvature = CurvatureSign * CGAL::sqrt( oRes.first.squared_length() );
              vit->second.fMixedArea = oRes.second;
              nSurfacePoints ++;
            }
            else
            {
              //              std::cout << "reject " << std::endl;
            }
            PointList.push_back( vit->first->point() );   // used for fitting regardless of surface v. n-junction
          }
          
          PatchIter->second.NumSurfacePoints = nSurfacePoints;
          PatchIter->second.MeanCurvature
            = PatchIter->second.MeanCurvature  / static_cast<RT>( nSurfacePoints ); 
          PatchIter->second.WeightedMeanCurvature 
            =    PatchIter->second.WeightedMeanCurvature / PatchIter->second.TotalMixedArea; 
          
          PatchIter->second.fMeanCurvature /= static_cast<RT>( nSurfacePoints ); 
          PatchIter->second.fWeightedMeanCurvature /= PatchIter->second.TotalMixedArea; 

          // Use Jet-fitting
          Monge_form FittedMorgeForm;
          Monge_via_jet_fitting MongeJetFit;
          if( PointList.size() > 15 )
          {
            FittedMorgeForm = MongeJetFit( PointList.begin(), PointList.end(), 2, 2 );
            
            PatchIter->second.LeastSquare_ConditionNum = MongeJetFit.condition_number();
            PatchIter->second.LeastSquare_K1           = FittedMorgeForm.principal_curvatures(0);
            PatchIter->second.LeastSquare_K2           = FittedMorgeForm.principal_curvatures(1);
            PatchIter->second.FitOrigin                = FittedMorgeForm.origin();
            PatchIter->second.Fitted                   = true;
          }
          else
          {
            PatchIter->second.LeastSquare_ConditionNum = -5;
            PatchIter->second.LeastSquare_K1           = -5;
            PatchIter->second.LeastSquare_K2           = -5;
            PatchIter->second.Fitted                   = false;
          }
        }
      }

      Pandora::AreaAveragedCurvature<C3T3, SurfaceIndexBndMapIter, VertexCurveMap >
        ParallelCurveOp( NJunctionSurface.begin(), NJunctionSurface.end(), Mesh );
      ParallelCurveOp.GetPatchAveragedCurvature( 20, 2 );
 
    }

  private:
    
    //--------------------------------
    //  ExtractBounary
    //--------------------------------
    SurfaceIndexBoundaryMap ExtractBoundary()
    {
      SurfaceIndexBoundaryMap         SurfaceMap;                                    

      int           InfiniteCount = 0;
      typedef typename C3T3::Facet_iterator FacetIter;
      for( FacetIter pIter = Mesh.facets_begin();
           pIter != Mesh.facets_end(); ++ pIter )
        
      {
        if( ! Mesh.triangulation().is_infinite( *pIter ) )
        {
          Surface_index ID = Mesh.surface_index( *pIter );
          
          if ( SurfaceMap.find( ID ) == SurfaceMap.end()  )
            SurfaceMap.insert( std::make_pair( ID, BoundaryPatch() ) );
          
          SurfaceMap[ID].TotalArea
            += CGAL::sqrt( Mesh.triangulation().triangle( *pIter ).squared_area() );
          SurfaceMap[ID].Facets.insert( *pIter );
          for (int i = 0; i < 4; i++)
          {
            if ( i != pIter->second )
            {
              const Vertex_handle& vh = pIter->first->vertex(i);
              BndVertexInfo Info( Vector_3( Point3(0,0,0), Point3(0,0,0) ), 0);
              Info.N_SampledPoints = 0;
              Info.fCurvature      = 0;
              SurfaceMap[ID].CurvatureMap.insert( std::make_pair( vh, Info ) );
            }
          }
        }
        else
        {
          InfiniteCount ++;
        }
      }
      std::cout << " Infinite Count " << InfiniteCount << std::endl;
      return SurfaceMap;
    }
    
    SurfaceIndexBoundaryMap         NJunctionSurface;
    const std::vector<SMatrix3x3> & IDOrientMap;
    const C3T3                    & Mesh;
 
  }; 
}



namespace Pandora
{
  template< class C3T3, class Iterator, class VertexCurveMap >
  class AreaAveragedCurvature
  {
  public:
    typedef typename C3T3::Triangulation::Vertex_iterator          Vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_facets_iterator   Finite_facet_iterator;
      
    typedef typename C3T3::Triangulation::Vertex                   Vertex;
    typedef typename C3T3::Triangulation::Geom_traits::Bare_point  Point3;
    typedef typename C3T3::Triangulation::Geom_traits::RT          RT;
    typedef typename C3T3::Cell_handle                             Cell_handle;
    typedef typename C3T3::Facet                                   Facet;
    typedef typename C3T3::Vertex_handle                           Vertex_handle;
      
    typedef typename C3T3::Surface_index                              Surface_index;
    typedef typename C3T3::Subdomain_index                            Subdomain_index;
    typedef typename std::multimap<Subdomain_index, Subdomain_index>  CrossMeshIDMap;  // maps from grain in mesh 1 to mesh 2
    typedef typename std::pair<Subdomain_index, Subdomain_index>      IDPair;
      
    typedef typename C3T3::Triangulation::Locate_type              Locate_type;
    typedef typename C3T3::Triangulation::Geom_traits::Triangle_3  Triangle_3;
    typedef typename C3T3::Triangulation::Geom_traits::Plane_3     Plane_3;
    typedef typename C3T3::Triangulation::Geom_traits::Vector_3    Vector_3;
    typedef typename std::set<Vertex_handle>                       VertexSet;
      
  public:
      
    AreaAveragedCurvature( Iterator Start, Iterator End,
                           const C3T3     & Mesh_ )
      : m_BoundaryStart( Start ), m_BoundaryEnd( End ), m_Mesh( Mesh_ )
    {
    }


    //---------------------------------------------------------------
    //  GetPatchAveragedCurvature
    //---------------------------------------------------------------
    void GetPatchAveragedCurvature( int NumPointsPerSide, RT Delta )
    {
      std::ofstream TestStream( "Average.txt" );
            
      typedef typename std::set<Facet> FacetSet;
      typedef typename FacetSet::const_iterator FacetSetIter;

      for( int i = -NumPointsPerSide; i <= NumPointsPerSide; ++ i)
      {
        TestStream << static_cast<RT>(i) * Delta << " ";
      }
      TestStream << std::endl;
      

      for( Iterator it = m_BoundaryStart; it != m_BoundaryEnd; ++ it  )
      {        
        TestStream << it->first.first << " " << it->first.second << " ";
        for( int i = -NumPointsPerSide; i <= NumPointsPerSide; ++ i)
        {
          RT TotalArea = 0;
          for( FacetSetIter FacetIter = it->second.Facets.begin();
               FacetIter != it->second.Facets.end(); ++ FacetIter )
          {
            TotalArea +=  GetParallelTranslatedArea( *FacetIter,
                                                     it->second.CurvatureMap,
                                                     it->first.first,
                                                     static_cast<RT>( i ) * Delta );
          }
          TestStream  << " " << TotalArea << " ";
        }
        TestStream << std::endl;
      }
    }
  private:
    
    //---------------------------------------------------------------
    //  GetParallelTranslatedArea
    //   Given Facet f represented by T = c3t3_.triangulation().triangle( f ),
    //   find the area of T after each vertex is translated 
    //
    //---------------------------------------------------------------
    RT GetParallelTranslatedArea( const Facet & f,
                                  const VertexCurveMap & CurvatureMap,
                                  Int InteriorID,
                                  RT delta ) const
    {
      Vertex_handle vh1 = f.first->vertex( ( f.second + 1 ) & 3 );
      Vertex_handle vh2 = f.first->vertex( ( f.second + 2 ) & 3 );
      Vertex_handle vh3 = f.first->vertex( ( f.second + 3 ) & 3 );
      
      if( CurvatureMap.find( vh1 ) != CurvatureMap.end()
          && CurvatureMap.find( vh2 ) != CurvatureMap.end()
          && CurvatureMap.find( vh3 ) != CurvatureMap.end() )
      {
        
        Vector_3 K1 = CurvatureMap.find( vh1 )->second.Curvature ;
        Vector_3 K2 = CurvatureMap.find( vh2 )->second.Curvature ;
        Vector_3 K3 = CurvatureMap.find( vh3 )->second.Curvature ;
        Point3 ReferencePoint = CGAL::centroid( m_Mesh.triangulation().triangle( f ) );
        
        //  just need a consistent direction
        K1 = K1 * Pandora::Details::CorrectCurvatureSign( f, ReferencePoint, K1, InteriorID, m_Mesh );
        K2 = K2 * Pandora::Details::CorrectCurvatureSign( f, ReferencePoint, K2, InteriorID, m_Mesh );
        K3 = K3 * Pandora::Details::CorrectCurvatureSign( f, ReferencePoint, K3, InteriorID, m_Mesh );
        
        
        Point3 v1 = vh1->point();
        Point3 v2 = vh2->point();
        Point3 v3 = vh3->point();
        
        
        
        //  translate in direction of K1, K2, and K3
        if( K1 != CGAL::NULL_VECTOR )
          K1 = Pandora::Details::UnitVector( K1 );
        if( K2 != CGAL::NULL_VECTOR )
          K2 = Pandora::Details::UnitVector( K2 );
        if( K3 != CGAL::NULL_VECTOR )
          K3 = Pandora::Details::UnitVector( K3 );

        v1 = v1 + K1 * delta;
        v2 = v2 + K2 * delta;
        v3 = v3 + K3 * delta;
        
        
        Triangle_3 T( v1, v2, v3 );
        return CGAL::sqrt( T.squared_area() );
        
      }
      else
      {
        bool bEdge = 
          ( Pandora::Details::IsSurfaceVertex( vh1 , m_Mesh ) ||
            Pandora::Details::IsSurfaceVertex( vh2 , m_Mesh ) ||
            Pandora::Details::IsSurfaceVertex( vh3 , m_Mesh ) );

        if( ! bEdge )
          std::cout << "Vertex not found, but facet is a surface facet " << std::endl;
      }
      return 0;
    }
    
    Iterator m_BoundaryStart, m_BoundaryEnd;
    const C3T3     & m_Mesh;
  };
}

#endif
