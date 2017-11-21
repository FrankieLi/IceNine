//-------------------------------------------------------
//
//  M_EstimateWithIRLS
//  M-estimate with iterative reweighting least square
//
//  Author:  Frankie Li
//  e-mail:  li31@llnl.gov
//
//  Implementation of robust curvature estimate:
//  "Robust statistical estimation of curvature on discretized
//   surfaces" - Evangelos Kalogerakis, Patricio Simari,
//               Derek Nowrouzezahrai, Karan Singh
//  Eurographics Symposium on Geometry Prossing (2007)
//
//
//  Requirement:  Eigen library
//
//-------------------------------------------------------

#ifndef M_EstimateWithIRLS_H
#define M_EstimateWithIRLS_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include "SimpleGraph.h"
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/multi_array.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include "MeanCurvatureNormal.h"

namespace MEstimateIRLS
{
  
  template< class C3T3 >
  class RobustCurvatureEstimate
  {
  private:

    typedef typename C3T3::Vertex_handle                           Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits::Vector_3    Vector_3;
    typedef typename C3T3::Triangulation::Geom_traits::RT          RT;
    typedef typename C3T3::Facet                                   Facet;
    typedef typename C3T3::Triangulation::Geom_traits::Bare_point  Point3;
    typedef typename C3T3::Triangulation::Geom_traits::Plane_3     Plane_3;
    typedef typename C3T3::Surface_index                           Surface_index;
    
    typedef CSimpleUGraph<Vertex_handle>                  SurfaceGraphT;
    typedef typename SurfaceGraphT::DataVertexMapT DataVertexMapT;
    

    typedef typename SurfaceGraphT::Edge                           EdgeT;
    typedef typename std::map<EdgeT, Float>                        WeightMapT;
    typedef typename std::vector<Vertex_handle>                    VertexListT;
    typedef typename std::set<Vertex_handle>                       VertexSetT;


    typedef typename Eigen::Matrix< Float, 3, 1 >       Vector3;
    typedef typename Eigen::Matrix< Float, 2, 2 >       Matrix22;
    typedef typename Eigen::Matrix< Float, 4, 1 >       Vector4;
    typedef typename Eigen::Matrix< Float, Eigen::Dynamic, 1 > VectorX;
    typedef typename Eigen::Matrix< Float, Eigen::Dynamic, 3 > Matrix3X;
    typedef typename Eigen::Matrix< Float, Eigen::Dynamic, Eigen::Dynamic > MatrixNM;
    

    //-------------------------------------------------
    //  Median  (dumbest possible way)
    //-------------------------------------------------
    inline Float GetMedian( vector<Float> & v )
    {
      std::sort( v.begin(), v.end() );

      return ( v[ v.size()/2 ] + v[ v.size()/2 + 1 ] ) / static_cast<Float>(2);
    }
    //-------------------------------------------------
    //  WeightFunction
    //-------------------------------------------------
    inline Float WeightFunction( Float r,  Float Sigma )
    {
      Float t = ( static_cast<Float>(1) + ( r * r / (Sigma * Sigma ) ));
      return static_cast<Float>( 2 ) / ( t * t );
    }
    
    //-------------------------------------------------
    //  ConstructDenseSampleSet
    //
    //  Input:   vertex set:  {v0, v1, v2...}, and normal set, {n0, n1, n2...}
    //  Output:  sample pairs (v_i, v_j, n_i, n_j) with distance cost
    //
    //  Return:  The A matrix, x, and b vector
    //
    //  Ordering:  (u, v)^T
    //
    //-------------------------------------------------
    void ConstructDenseSampleSet( Matrix3X & A,
                                  VectorX  & DataVector,
                                  const Vertex_handle & Center,
                                  const VertexListT & VertexList,
                                  const std::vector< Vector_3 > & NormalList,
                                  const Vector_3 & Tangent_u,
                                  const Vector_3 & Tangent_v )
    {
      int VectorSize = ( VertexList.size() * (VertexList.size() - 1) ); // Num sample points = n(n-1)/2
                                                                        // 2 piece of data per sample point
      A.resize( VectorSize, 3 );
      DataVector.resize( VectorSize ); 

      int nCurrentPos = 0;
      for( int i = 0; i < VertexList.size(); i ++ )
      {
        for( int j = i + 1; j < VertexList.size(); j ++ )
        {
          Vector_3 DeltaNormal = NormalList[j] - NormalList[i];
          Vector_3 DeltaP      = VertexList[j]->point() - VertexList[i]->point();
          
          DataVector( nCurrentPos ) = DeltaNormal * Tangent_u;
          A( nCurrentPos, 0 ) = DeltaP * Tangent_u;
          A( nCurrentPos, 1 ) = DeltaP * Tangent_v;
          A( nCurrentPos, 2 ) = 0;
          nCurrentPos ++;
          
          DataVector( nCurrentPos ) = DeltaNormal * Tangent_v;
          A( nCurrentPos, 1 ) = DeltaP * Tangent_u;
          A( nCurrentPos, 2 ) = DeltaP * Tangent_v;
          A( nCurrentPos, 0 ) = 0;
          
          nCurrentPos ++;
        }
      }
      //      std::cout << "VectorSize " << VectorSize << " nCurrentPos " << nCurrentPos << std::endl;
    }


    //-------------------------------------------------
    //  GetDiagonal_II_Form
    //-------------------------------------------------
    Matrix22 GetDiagonal_II_Form( const Vector3 & v )
    {
      Matrix22 M;
      M(0, 0) = v(0);
      M(0, 1) = v(1);
      M(1, 0) = v(1);
      M(1, 1) = v(2);

      Eigen::SelfAdjointEigenSolver<Matrix22> es(M);
      
      return es.eigenvalues().asDiagonal();
    }
    //-------------------------------------------------
    //
    //   BuildVertexList
    //   - Add up to approximate N-ring
    //
    //   - note:  Triple junction vertices are used for curvature
    //            estimates, even though curvature are not
    //            calculated on N-junctions
    //-------------------------------------------------
    void BuildVertexList( VertexSetT & VertexSet,
                          SurfaceGraphT & G, WeightMapT & EdgeWeightMap,
                          const Vertex_handle & CenterVertex,
                          const C3T3 & Mesh, Surface_index BndSurfaceIndex,
                          int NRingMax )
    {
      DEBUG_ASSERT( NRingMax >= 0, "Error:  M_EstimateWithIRLS::BuildVertex:  NRingMax < 0 -- must have at least 0-ring!" );
            
      VertexSet.insert( CenterVertex );
      std::vector<Facet> FacetList;
      Mesh.triangulation().finite_incident_facets( CenterVertex, std::back_inserter( FacetList ) );      
      
      if( NRingMax == 0 )  // base case
        return;
      
      for( int j = 0; j < FacetList.size(); j ++ )
      {
        int VertexIndex[3];
        int nIndex = 0;
        if( Mesh.is_in_complex( FacetList[j] )
            &&  BndSurfaceIndex == Mesh.surface_index( FacetList[j] ) )  // same surface ID
        {
          for( int n = 0; n < 4; n ++ )  // need all edges instead
            if( n != FacetList[j].second )
            {
              VertexIndex[nIndex] = n; 
              nIndex ++;
            }
          
          Vertex_handle vh0 = FacetList[j].first->vertex( VertexIndex[0] ); 
          Vertex_handle vh1 = FacetList[j].first->vertex( VertexIndex[1] ); 
          Vertex_handle vh2 = FacetList[j].first->vertex( VertexIndex[2] ); 
          
          //  Edges can be added to vertices that already exists
          bool  EdgeAdded;
          EdgeT e;
          G.AddEdge( &EdgeAdded, &e, vh0, vh1 );
          EdgeWeightMap[e] = CGAL::sqrt( (vh0->point() - vh1->point()).squared_length() );

          G.AddEdge( &EdgeAdded, &e, vh1, vh2 );
          EdgeWeightMap[e] = CGAL::sqrt( (vh1->point() - vh2->point()).squared_length() );
          
          G.AddEdge( &EdgeAdded, &e, vh2, vh0 );
          EdgeWeightMap[e] = CGAL::sqrt( (vh2->point() - vh0->point()).squared_length() );

          if( NRingMax > 0 )
          {
            BuildVertexList( VertexSet, G, EdgeWeightMap, vh0, Mesh, BndSurfaceIndex, NRingMax -1 );
            BuildVertexList( VertexSet, G, EdgeWeightMap, vh1, Mesh, BndSurfaceIndex, NRingMax -1 );
            BuildVertexList( VertexSet, G, EdgeWeightMap, vh2, Mesh, BndSurfaceIndex, NRingMax -1 );
          }
        }
      }
    }
    
 
    
  public:

    //-------------------------------------------------
    //  Given:  vertex point list {v0, v1, v2,...}
    //
    //
    //  0. for each vertex, build tangent plane
    //  1. for each vertex not on N-lines, c
    //       a.  Initial guess of curvature using Rusinkiewicz
    //       b.  Build residual list, find median (Sigma)
    //       [ Apply IRLS]
    //       d.  while n < N_Iter && Error > Epsilon
    //           -SolveWeightedLeastSequare
    //           -Update Residual List (weights)
    //           -Update Sigma
    //
    //       e. recompute normal
    //
    // a.  need to consider running all point distance
    // b.  need to consider number of components in surface patch
    // 
    //
    //-------------------------------------------------
    void M_EstimateCurvature( RT & MeanCurvature,
                              int   & NumSampledPoints,
                              Vector3 & TensorResidualErr,
                              const Vertex_handle & CenterVertex, const C3T3 & Mesh,
                              Surface_index         BndSurfaceIndex       )
    {
      SurfaceGraphT    BoundaryVertexGraph;
      WeightMapT       EdgeWeightMap;
      VertexListT      InitialVertexList;
      
      VertexSetT       VertexSet;
      BuildVertexList( VertexSet, BoundaryVertexGraph,
                       EdgeWeightMap,  CenterVertex,
                       Mesh, BndSurfaceIndex,  2);

      typedef boost::multi_array<Float, 2>                   DistanceMetric;
      Pandora::Details::MeanCurvatureNormal<C3T3>            MeanCurvatureNormalOp;
      typename boost::associative_property_map< WeightMapT > WeightPMap( EdgeWeightMap );
      
      InitialVertexList = BoundaryVertexGraph.GetVertices();
      int NumVertex = InitialVertexList.size();
      DistanceMetric D(boost::extents[NumVertex][NumVertex]);
      boost::johnson_all_pairs_shortest_paths( BoundaryVertexGraph.oBoostGraph, D,
                                               boost::get( boost::vertex_index, BoundaryVertexGraph.oBoostGraph ),
                                               WeightPMap, 0 );
      
      std::vector< Vector_3 >     NormalList;
      VertexListT                 VertexList;
      std::pair<Vector_3, double> CenterVertexCurvature = MeanCurvatureNormalOp( CenterVertex, Mesh, BndSurfaceIndex );
            
      for( int i = 0; i < InitialVertexList.size(); i ++ )
      {
        std::pair< Vector_3, double > Result = MeanCurvatureNormalOp( InitialVertexList[i], Mesh, BndSurfaceIndex );
        if( Result.second > 1e-8 )   // some epsilon
        {
          if( Result.first * CenterVertexCurvature.first > 0 )
            NormalList.push_back(   Result.first / CGAL::sqrt( Result.first.squared_length() ) );
          else
            NormalList.push_back( - Result.first / CGAL::sqrt( Result.first.squared_length() ) );   // point to the same dir
          VertexList.push_back( InitialVertexList[i] );
        }
      }

      Matrix3X A;
      VectorX  DataVector;
      Vector_3 Tangent_u, Tangent_v;            //  Construct Tangent Vectors...
      
      Plane_3 TangentPlane( Point3( 0, 0, 0 ), CenterVertexCurvature.first.direction() );
      Tangent_u = TangentPlane.base1();
      Tangent_v = TangentPlane.base2();
      Tangent_u = Tangent_u / CGAL::sqrt( Tangent_u.squared_length() );
      Tangent_v = Tangent_v / CGAL::sqrt( Tangent_v.squared_length() );
            
      ConstructDenseSampleSet( A, DataVector, CenterVertex, VertexList, NormalList, Tangent_u, Tangent_v );
      DataVertexMapT     DataMap = BoundaryVertexGraph.GetDataVertexMap();
      std::vector<Float> DistanceValues;
      for( int i = 0; i < VertexList.size(); i ++ )
      {
        int v1 = DataMap[ VertexList[i] ];
        for( int j = i + 1; j < VertexList.size(); j ++ )
        {
          int v2 = DataMap[ VertexList[j] ];
          DistanceValues.push_back( D[v1][v2] * D[v1][v2] );
          DistanceValues.push_back( D[v1][v2] * D[v1][v2] );   // deliberate duplicate
        }
      }
      NumSampledPoints = VertexList.size();
      
      Vector3 CurvatureSolution = A.fullPivLu().solve( DataVector );   

      //-----------------------------------
      // iterative step;
      //-----------------------------------
      int   MaxIterations = 50;
      Float ConvergenceEps = 1e-5;
      Float Epsilon = 1;
      int   iter = 0;
      int   NumRows = (VertexList.size() -1 ) * VertexList.size();

      std::vector<Float> ResidualValues;
      ResidualValues.resize( NumRows );
      MatrixNM WeightMatrix;
      WeightMatrix.resize( NumRows, NumRows );

      
      
      do
      {
        WeightMatrix = MatrixNM::Zero( NumRows, NumRows );
        VectorX Residual = A * CurvatureSolution - DataVector;
        for( int i = 0; i < NumRows; i ++ )
          ResidualValues[i] = std::abs( Residual(i) );
        Float Sigma = 1.4826 * GetMedian( ResidualValues );   // from Kalogerakis et al.

        for( int i = 0; i < NumRows; i += 2 )                 // recompute WeightMatrix, NumRows is by definition even
        {
          Float tmp = WeightFunction( Residual(i), Sigma ) + WeightFunction( Residual(i+1), Sigma );
          WeightMatrix( i,     i ) = tmp / DistanceValues[i];
          WeightMatrix( i+1, i+1 ) = tmp / DistanceValues[i];
        }

        
        MatrixNM LHS_Product = A.transpose() * WeightMatrix * A;
        Eigen::ColPivHouseholderQR<MatrixNM> DEC( LHS_Product);

        VectorX  A_DataVector_Prod = A.transpose() * WeightMatrix * DataVector;
        Vector3 NewSolution = DEC.solve( A_DataVector_Prod );
        TensorResidualErr = NewSolution - CurvatureSolution;
        Epsilon  = TensorResidualErr.norm()  / NewSolution.norm(); // norm() is L2 norm
        CurvatureSolution = NewSolution;
        iter ++;
        
      }while( iter < MaxIterations && Epsilon > ConvergenceEps ); 
      Matrix22 S = GetDiagonal_II_Form( CurvatureSolution );
      MeanCurvature = ( S(0, 0) + S(1, 1) ) / static_cast<Float>( 2 ) ;
    }
    
  };
}


#endif
