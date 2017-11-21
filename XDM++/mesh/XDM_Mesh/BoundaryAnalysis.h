//----------------------------------------------------
//  BoundaryAnalysis.h
//
//  Author:   Frankie Li
//  Purpose:  Boundary analysis for different meshes
//----------------------------------------------------


#ifndef BOUNDARY_ANALYSIS_H
#define BOUNDARY_ANALYSIS_H

//#include <CGAL/basic.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include "AABB_mesh_3_triangle_primitive.h"

#include <CGAL/AABB_triangle_primitive.h>

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
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "M_EstimateWithIRLS.h"
#include <Eigen/Dense>

namespace Pandora
{

  //-------------------------------------------------------------
  //  BoundaryAnalysis
  //  Purpose:  Provide standard toolkit on calculating distances and shift
  //            between two different meshes
  //-------------------------------------------------------------
  template< class C3T3 >
  class BoundaryAnalysis
  {
  public:
    typedef typename Eigen::Matrix< Float, 3, 1 >       Vector3;
    
    struct Kernel: public CGAL::Exact_predicates_inexact_constructions_kernel {};
    typedef typename C3T3::Triangulation::Vertex_iterator            Vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_vertices_iterator   Finite_vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_facets_iterator     Finite_facet_iterator;
    typedef typename C3T3::Triangulation::Vertex                     Vertex;
    typedef typename C3T3::Triangulation::Vertex_handle              Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits::Bare_point    Point3;

    typedef typename C3T3::Cell_handle                               Cell_handle;
    typedef typename C3T3::Surface_index  Surface_index;
    typedef typename C3T3::Subdomain_index                           Subdomain_index;
    typedef typename std::multimap<Subdomain_index, Subdomain_index> CrossMeshIDMap;  // maps from grain in mesh 1 to mesh 2
    typedef typename std::pair<Subdomain_index, Subdomain_index>     IDPair;

    typedef typename C3T3::Triangulation::Locate_type     Locate_type;
    typedef typename C3T3::Triangulation::Facet           Facet;
    typedef typename C3T3::Triangulation::Geom_traits::Ray_3       Ray3;
    typedef typename C3T3::Triangulation::Geom_traits::Triangle_3  Triangle_3;
    typedef typename C3T3::Triangulation                           Triangulation;
    typedef typename Kernel::Iso_cuboid_3                          Iso_cuboid;
    typedef typename C3T3::Triangulation::Geom_traits::Vector_3    Vector_3;
  
    typedef typename C3T3::Facet_iterator                          Facet_iterator;
    
    typedef GeneralLib::SQuaternion                                SQuaternion;
    
    //-------------------------------------------------------
    //  Random sampling
    //-------------------------------------------------------
    typedef GeneralLib::SMatrix3x3 SMatrix3x3;
    typedef GeneralLib::SVector3   SVector3;
    typedef GeneralLib::UniformGrid::CQuaternionGrid CQuaternionGrid;

    mutable GeneralLib::CUniformRandomReal RandGen;
    
    //-------------------------------------------------------
    //  kdtree
    //-------------------------------------------------------
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;
    typedef typename Neighbor_search::Tree Tree;


    
    typedef typename C3T3::Triangulation::Geom_traits Tr_Kernel;
    typedef CGAL::Simple_cartesian<double> AABB_Kernel;
    
    typedef typename AABB_Kernel::Ray_3     AABB_Ray;
    typedef typename AABB_Kernel::Segment_3 AABB_Segment;
    typedef typename AABB_Kernel::Point_3   AABB_Point;
    typedef typename AABB_Kernel::FT        AABB_FT;
    
    typedef typename C3T3::Triangulation::Geom_traits::FT FT;

//     typedef CGAL::Pandora::BoundarySelectionPredicate<C3T3, Finite_facet_iterator >  BoundarySelectFilter;
//     typedef CGAL::Filter_iterator<Finite_facet_iterator, BoundarySelectFilter>       Filtered_facet_iterator;

    typedef typename std::list< Finite_facet_iterator >                                 Filtered_finite_facet_list;
    typedef typename std::list< Finite_facet_iterator >::const_iterator                 Filtered_facet_iterator;
    
    typedef typename CGAL::AABB_mesh_3_triangle_primitive<C3T3, AABB_Kernel, Filtered_facet_iterator>   AABB_Primitive;
    typedef typename CGAL::AABB_traits<AABB_Kernel, AABB_Primitive>                   AABB_Mesh_Triangle_Traits;
    typedef typename CGAL::AABB_tree<AABB_Mesh_Triangle_Traits>                       Mesh_AABB_Tree;
    
    typedef typename Mesh_AABB_Tree::Object_and_primitive_id                          Object_and_primitive_id;
    typedef typename Mesh_AABB_Tree::Primitive_id                                     Primitive_id;
    typedef typename Mesh_AABB_Tree::Point_and_primitive_id                           Point_and_primitive_id;
    typedef CGAL::Cartesian_converter< Tr_Kernel, AABB_Kernel>                        TrKernel_To_AABB;
    
    //-------------------------------------
    //  BndMotionInfo
    //-------------------------------------
    struct BndMotionInfo
    {
      Finite_facet_iterator FacetIter1;
      Finite_facet_iterator FacetIter2;

      vector<Vertex>   Bnd1Vertex;
      vector<SVector3> VertexCurvature;
      vector<Float>    VertexMixedArea;
      vector<int>      Bnd1VertexType;

      Triangle_3 BndFacet1;
      Triangle_3 BndFacet2;
      Surface_index FacetID1;
      Surface_index FacetID2;
      std::pair< SQuaternion, SQuaternion > Pair1;
      std::pair< SQuaternion, SQuaternion > Pair2;
      double   fDistance;
      SVector3 MinDisplacement;
      //  these two are used for projected distance method only
      bool   bDisappeared;
      double fNormalError;
    };

    //  Get the Connected-ness of this vertex (in terms of IDs)
    int GetJunctionType( const C3T3& Mesh, const Vertex_handle & vh )
    {
      typename std::vector< Cell_handle > Cells;
      typename std::set< int > UniqueIDs;
      Mesh.triangulation().finite_incident_cells( vh, std::back_inserter(Cells) );
      for( int i = 0; i < Cells.size(); i ++ )
      { 
        int nID =  Mesh.subdomain_index( Cells[i] );
        if( Mesh.is_in_complex( Cells[i] ) &&  nID > 0 )
        {
          UniqueIDs.insert( nID );
        }
      }
      return UniqueIDs.size();
    }

    
  public:

    BoundaryAnalysis( const C3T3 & Mesh1,
                      const C3T3 & Mesh2,
                      std::vector<SMatrix3x3> & IDOrientMap1_,
                      std::vector<SMatrix3x3> & IDOrientMap2_  ):
      Mesh1_( Mesh1 ), Mesh2_( Mesh2 ),
      IDOrientMap1( IDOrientMap1_ ),
      IDOrientMap2( IDOrientMap2_ )
    {}
    
    //--------------------------------------
    //  InitializeByBoundary
    //--------------------------------------
    void InitializeByBoundary()
    {
      MeshSurfacePointList1.clear();
      MeshSurfacePointList2.clear();

      ExtractBoundaryPoints( std::back_inserter( MeshSurfacePointList1 ),
                             Mesh1_, -1 );
      ExtractBoundaryPoints( std::back_inserter( MeshSurfacePointList2 ),
                             Mesh2_, -1 );
    }

    //--------------------------------------
    //  CheckFullCompatibility
    //
    //  Purpose:  Checks to see if the surfaces are
    //            a)  Compatible and
    //            b)  Belong to the "same side," i.e., exterior or interior
    //                directions are consistent
    //--------------------------------------
    bool CheckFullCompatibility( const C3T3 & Mesh1_, const C3T3 & Mesh2_,
                                 const Facet & f1, const Facet & f2,
                                 double fThreshold, 
                                 double & Angle1, double & Angle2,
                                 IDPair & Pair1, IDPair & Pair2,
                                 FT & NormalAngleError )
    {
      if( IsCompatible( Mesh1_, Mesh2_, f1, f2, fThreshold,
                        Angle1, Angle2, Pair1, Pair2 ) )
      {
        Triangle_3 T1          = Mesh1_.triangulation().triangle( f1 );
        Vector_3   PatchNormal = T1.supporting_plane().orthogonal_vector();
        int nSign_1 =  Pandora::Details::CorrectCurvatureSign( f1, CGAL::centroid( T1 ),  PatchNormal,
                                                               Pair1.first, Mesh1_ );
        
        //----------------------------------------------------------------
        //  We are deliberately tossing away a lot of boundaries.  For
        //  example, when a grain can be corresponded via misorientation, and
        //  a neighboring grain has disappeared, the resulting boundary will
        //  not be selecting in this particular method.  
        //----------------------------------------------------------------
        Triangle_3 T2                     = Mesh2_.triangulation().triangle( f2 );
        Vector_3 Intersected_Patch_Normal = T2.supporting_plane().orthogonal_vector();
        int nSign_2                       = Pandora::Details::CorrectCurvatureSign( f2,
                                                                                    CGAL::centroid( T2 ),
                                                                                    Intersected_Patch_Normal,
                                                                                    Pair1.second, Mesh2_ );
        
        Vector_3 u1 = nSign_1 * Pandora::Details::UnitVector( PatchNormal              );
        Vector_3 u2 = nSign_2 * Pandora::Details::UnitVector( Intersected_Patch_Normal );
        
        FT u1Dotu2 = u1 * u2;
        u1Dotu2 = std::min( u1Dotu2, static_cast<FT>(  1 ) );  // clip at [0, 1]
        u1Dotu2 = std::max( u1Dotu2, static_cast<FT>( -1 ) );  // clip at [-1, 0]
        NormalAngleError  = acos( u1Dotu2  );

        if ( u1Dotu2 > 0 )   // non-negative indicates same side
          return true;
      }

      return false;
    }
    
    //--------------------------------------
    //  CenterOfMasses
    //  Purpose:  Return center of mass for both meshes
    //--------------------------------------
    std::pair<SVector3, SVector3> CenterOfMasses()
    {
      
      SVector3 CenterOfMass1 = GetCenterOfMass( MeshSurfacePointList1.begin(),
                                                MeshSurfacePointList1.end() );
      SVector3 CenterOfMass2 = GetCenterOfMass( MeshSurfacePointList2.begin(),
                                                MeshSurfacePointList2.end() );
      return std::make_pair( CenterOfMass1, CenterOfMass2 );
    }

    struct NotInBox
    {
      Iso_cuboid BBox;
      NotInBox( Iso_cuboid b ): BBox( b ){}
      bool operator()( const Point3 & p ) const
      {
        return BBox.has_on_unbounded_side( p );
      }
    };

    //--------------------------------------
    //  IterativeGetFirstIntersection
    //--------------------------------------
    bool IterativeGetFirstIntersection( int NumTries, AABB_FT delta,
                                        const Mesh_AABB_Tree & ReferenceTree,
                                        const AABB_Ray & QueryRay,
                                        std::list<Object_and_primitive_id > & ResultList )
    {
      while( NumTries > 0 )
      {
        NumTries --;

        AABB_Segment SegmentQuery( QueryRay.source(),
                                   QueryRay.source() + delta * Pandora::Details::UnitVector( QueryRay.to_vector() ) );
        if ( ReferenceTree.do_intersect( SegmentQuery ) )
        {
          ReferenceTree.all_intersections( SegmentQuery,
                                           std::back_inserter( ResultList ) );
          return true;
        }
        delta *= 2;
      }
      return false;
    }
    
    //--------------------------------------
    //  GetClosestIntersection
    //
    //--------------------------------------
    boost::tuple<AABB_FT, Finite_facet_iterator, GeneralLib::SVector3>
    GetClosestIntersection( const std::list<Object_and_primitive_id> & IntersectionList,
                            const AABB_Point & Origin, const C3T3 & Mesh, int EmptyID )
    {
      typedef typename std::list<Object_and_primitive_id>::const_iterator ResultIterator;

      bool    Success       = false;
      AABB_FT MinDistance   = std::numeric_limits< AABB_FT >::max();
      Finite_facet_iterator Nearest_Intersection;
      GeneralLib::SVector3  MinDisplacement;
      
      for( ResultIterator it = IntersectionList.begin(); it != IntersectionList.end(); ++ it  )
      {
        AABB_Point IntersectionPoint;
        AABB_FT    Distance;
        if( CGAL::assign( IntersectionPoint, it->first ) )
        {
          Surface_index ID = Mesh.surface_index( * (* (it->second) )  );
          if( ID.second != EmptyID && ID.first != EmptyID )
          {
            Success = true;
            Distance = CGAL::sqrt( CGAL::squared_distance( IntersectionPoint, Origin ) );
            if( Distance < MinDistance )
            {
              GeneralLib::SVector3 p( IntersectionPoint.x(), IntersectionPoint.y(), IntersectionPoint.z() );
              GeneralLib::SVector3 O( Origin.x(), Origin.y(), Origin.z() );
              MinDisplacement      = p - O;
              MinDistance          = Distance;
              Nearest_Intersection = *( it->second );
            }
          }
        }
      }
      if ( ! Success )
        MinDistance = static_cast< AABB_FT >( -1 );
      return boost::make_tuple( MinDistance, Nearest_Intersection, MinDisplacement );
    }
    
    //--------------------------------------
    //  ProjectedDistanceBndMotion
    //   Purpose:  Calculate the shift and create a facet-to-facet map between
    //             the two meshes using the local normal on the facet.
    //   Algorithm:  A ray is drawn from the local normal to the two directions (negative and
    //               positive normal direction) to find the closest surface patch with compatible
    //               facet.
    //
    //--------------------------------------
    vector<BndMotionInfo> ProjectedDistanceBndMotion( double fAngleThreshold, int EmptyID = -1  )
    {
      vector<BndMotionInfo> BndMotionList;
      
      Pandora::Details::MeanCurvatureNormal<C3T3> MeanCurvatureNormOp;
      MEstimateIRLS::RobustCurvatureEstimate<C3T3> RCE;
      std::cout << " Before building " << std::endl;

      Filtered_finite_facet_list FacetList;
      for( Finite_facet_iterator it = Mesh2_.triangulation().finite_facets_begin();
           it != Mesh2_.triangulation().finite_facets_end(); ++ it )
      {
        Surface_index SurfID = Mesh2_.surface_index( * it );
        if( SurfID.first != SurfID.second )
          FacetList.push_back( it );
      }
      
      Mesh_AABB_Tree ReferenceTree( FacetList.begin(), FacetList.end() );
      
      std::cout << " After building " << std::endl;
      
      TrKernel_To_AABB ToAABB;

      int nFacetsProcessed = 0;
      for( Facet_iterator it = Mesh1_.facets_begin();
           it != Mesh1_.facets_end(); ++ it )
      {
        nFacetsProcessed ++;
        if ( nFacetsProcessed % 5000 == 1 )
        {
          std::cout << " --- " << nFacetsProcessed << " / " << Mesh1_.triangulation().number_of_finite_facets() << std::endl;
        }
        Surface_index ID = Mesh1_.surface_index( *it );
        if( ID.second != EmptyID && ID.first != EmptyID )
        {
          //-----------------------------------------------------------
          //  NOTE: Since we are using ray intersection, by definition, the closest
          //  hit has to be the corresponding projected spot.  However, if the 
          //  projected spot is not a compatible facet, then the the original
          //  facet must have disappeared.
          //-----------------------------------------------------------
          Triangle_3 T           = Mesh1_.triangulation().triangle( *it ) ;
        
          Vector_3   PatchNormal = T.supporting_plane().orthogonal_vector();
          Point3     FacetCenter = CGAL::centroid( T );
        
          // construct ray from normal at vertex origin
          bool bPos = false;
          bool bNeg = false;
        
          AABB_Ray Positive_Ray_Query = ToAABB( Ray3( FacetCenter,  PatchNormal.direction() ) );  // may want to change to segment instead
          double Pos_Angle1, Pos_Angle2, Neg_Angle1, Neg_Angle2;
          IDPair Pos_Pair1,  Pos_Pair2,  Neg_Pair1,  Neg_Pair2;

          FT Pos_Distance, Neg_Distance;
          FT Pos_Normal_Angle, Neg_Normal_Angle;

          boost::tuple<AABB_FT, Finite_facet_iterator, GeneralLib::SVector3> Best_Pos_Intersection, Best_Neg_Intersection;
          
          std::list< Object_and_primitive_id > ResultList;
          if( IterativeGetFirstIntersection( 3, AABB_FT( 20 ),  ReferenceTree,
                                             Positive_Ray_Query,
                                             ResultList ) )
          {
            if( ResultList.size() > 0 )
            {
              Best_Pos_Intersection = GetClosestIntersection( ResultList, ToAABB( FacetCenter ),
                                                              Mesh2_, EmptyID );
              Pos_Distance = boost::get<0>( Best_Pos_Intersection );
              if( Pos_Distance > static_cast<AABB_FT>(-1) )
              {
                bPos = CheckFullCompatibility( Mesh1_, Mesh2_, *it,
                                               * boost::get<1>( Best_Pos_Intersection ),
                                               fAngleThreshold, Pos_Angle1, Pos_Angle2,
                                               Pos_Pair1, Pos_Pair2, Pos_Normal_Angle );
              }
            }
          
          }
          ResultList.clear();
          
          AABB_Ray Negative_Ray_Query = ToAABB( Ray3( FacetCenter, - PatchNormal.direction()  ));
          if(  IterativeGetFirstIntersection( 3, AABB_FT( 20 ),  ReferenceTree,
                                              Negative_Ray_Query,
                                              ResultList ) )
          {
            if( ResultList.size() > 0 )
            {
              Best_Neg_Intersection = GetClosestIntersection( ResultList, ToAABB( FacetCenter ),
                                                              Mesh2_, EmptyID );
              Neg_Distance = boost::get<0>( Best_Neg_Intersection );
              if( Neg_Distance > static_cast<AABB_FT>(-1) )
              {
                bNeg = CheckFullCompatibility( Mesh1_, Mesh2_, *it,
                                               *boost::get<1>(Best_Neg_Intersection ),
                                               fAngleThreshold, Neg_Angle1, Neg_Angle2,
                                               Neg_Pair1, Neg_Pair2, Neg_Normal_Angle );
              }
            }
          }
        

          //------------------------------------------
          //  Need to handle multiple cases:
          //
          //  1.  Check that the boundary normals are pointing in the same direction
          //      1a.  Record boundary normal (PatchNormal) angle as "error"
          //  
          //
          //------------------------------------------

          BndMotionInfo MotionInfo;
          if( bNeg && bPos )  // this really shouldn't happen...  not both sides could point to the same dir
          {
            std::cerr << "Unexpected error:  Multiple compatible boundaries" << std::endl;
            std::cerr << " Negative Normal Angle " << RADIAN_TO_DEGREE( Neg_Normal_Angle )
                      << " Positive Normal Angle " << RADIAN_TO_DEGREE( Pos_Normal_Angle ) << std::endl;
          }
          else if ( bNeg )
          {
            MeshToMeshIDMap.insert( Neg_Pair1 );
            MeshToMeshIDMap.insert( Neg_Pair2 );
            MotionInfo.BndFacet1 = T;
            MotionInfo.BndFacet2 = Mesh2_.triangulation().triangle( *boost::get<1>(Best_Neg_Intersection) );
            MotionInfo.FacetID1 = Neg_Pair1;
            MotionInfo.FacetID2 = Neg_Pair2;
          
            SQuaternion q1, q2;
            q1.Set( IDOrientMap1[ Neg_Pair1.first ] );
            q2.Set( IDOrientMap2[ Neg_Pair1.second ] );
            MotionInfo.Pair1 = std::make_pair( q1, q2 );
            MotionInfo.fDistance = Neg_Distance;
            MotionInfo.MinDisplacement = boost::get<2>( Best_Neg_Intersection );
            q1.Set( IDOrientMap1[ Neg_Pair2.first ] );
            q2.Set( IDOrientMap2[ Neg_Pair2.second ] );
            MotionInfo.Pair2 = std::make_pair( q1, q2 );
            MotionInfo.fNormalError = Neg_Normal_Angle;
            MotionInfo.bDisappeared = false;
          }
          else if ( bPos )
          {
            MeshToMeshIDMap.insert( Pos_Pair1 );
            MeshToMeshIDMap.insert( Pos_Pair2 );
            MotionInfo.BndFacet1 = T;
            MotionInfo.BndFacet2 = Mesh2_.triangulation().triangle( *boost::get<1>( Best_Pos_Intersection ) );
            MotionInfo.FacetID1 = Pos_Pair1;
            MotionInfo.FacetID2 = Pos_Pair2;
          
            SQuaternion q1, q2;
            q1.Set( IDOrientMap1[ Pos_Pair1.first ] );
            q2.Set( IDOrientMap2[ Pos_Pair1.second ] );
            MotionInfo.Pair1 = std::make_pair( q1, q2 );
            MotionInfo.fDistance = Pos_Distance;
            MotionInfo.MinDisplacement = boost::get<2>( Best_Pos_Intersection );
            q1.Set( IDOrientMap1[ Pos_Pair2.first ] );
            q2.Set( IDOrientMap2[ Pos_Pair2.second ] );
            MotionInfo.Pair2 = std::make_pair( q1, q2 );
            MotionInfo.fNormalError = Pos_Normal_Angle;
            MotionInfo.bDisappeared = false;
          }
          else   // Nothing is found - disappeared boundary - mark it?
          {
            MotionInfo.BndFacet1 = T;
            MotionInfo.BndFacet2 = T;
            MotionInfo.fDistance = -10;
            MotionInfo.bDisappeared = true;
          }

          // calculate vertex normal curvature
          for( int k = 0; k < 4; k ++ )
          {
            if ( k != it->second )
            {
              Surface_index ID = Mesh1_.surface_index( *it  );
              bool RunDebug = false;
              std::pair<Vector_3, double> Result =   MeanCurvatureNormOp( it->first->vertex( k ),
                                                                          Mesh1_, Mesh1_.surface_index( *it ),
                                                                          RunDebug );

              int JunctionType = GetJunctionType( Mesh1_, it->first->vertex(k) );
              
              if( Pandora::Details::IsSurfaceVertex( it->first->vertex( k ),
                                                     Mesh1_ ) )  // ignoring all n-junction lines
              {
                
                typedef typename C3T3::Triangulation::Geom_traits::RT            RT;
                RT RCE_Curvature;
                int NumSamplePoints;
                Vector3 ResidualErr;
                RCE.M_EstimateCurvature( RCE_Curvature,
                                         NumSamplePoints,
                                         ResidualErr,
                                         it->first->vertex( k ), 
                                         Mesh1_, Mesh1_.surface_index( *it ) ); 
                
                
                
                SVector3 Curvature( Result.first.x(),
                                    Result.first.y(),
                                    Result.first.z() );
                //MotionInfo.VertexCurvature.push_back( Curvature );   // really need to reorient normal
                // DEBUG
                Curvature.Normalize();
                MotionInfo.VertexCurvature.push_back( Curvature * RCE_Curvature );
                MotionInfo.VertexMixedArea.push_back( Result.second );
                MotionInfo.Bnd1Vertex.push_back( *( it->first->vertex(k) ) );
                MotionInfo.Bnd1VertexType.push_back( JunctionType );
              }
              else
              {
                MotionInfo.VertexCurvature.push_back( SVector3( 0, 0, 0 ) );
                MotionInfo.VertexMixedArea.push_back( Result.second );
                MotionInfo.Bnd1Vertex.push_back( *( it->first->vertex(k) ) );
                MotionInfo.Bnd1VertexType.push_back( JunctionType );
              }
              
              
            }
          }
          
          BndMotionList.push_back( MotionInfo );
        }
      }
      std::cout << "Done Projected distance " << std::endl;
      
      return BndMotionList;
    }
    
    //--------------------------------------
    //  MeshFacetMap
    //   Purpose:  Calculate the shift and create a facet-to-facet map between
    //             the two meshes.
    //   Parameter:  nMaxNearestNeighbors - determines the number of nearest neighbors
    //               to look at in the kdtree to find the best facet match in
    //               calculating the nearest shift
    //
    //               fAngleThreshold - determines the maximum misorientation between
    //                                 two grains to be counted as identical
    //
    //  Return:     Forest - indicating the correspondance between grains in mesh A to mesh B
    //                       represented by a Forest.  The forest is really a graph that where
    //                       the nodes are 
    //--------------------------------------
    vector<BndMotionInfo> MeshFacetMap( int nMaxNearestNeighbors, double fAngleThreshold, int EmptyID = -1  )
    {
      vector<BndMotionInfo> BndMotionList;
      
      ReferencePointTree = Tree( MeshSurfacePointList1.begin(),
                                 MeshSurfacePointList1.end() );

      std::cout << "Running MeshFacetMap with fAngleThreshold = " << fAngleThreshold << std::endl;
      
      for( Facet_iterator pIter = Mesh2_.facets_begin();
           pIter != Mesh2_.facets_end(); ++ pIter )
        
      {
        Surface_index ID = Mesh2_.surface_index( *pIter );
        if( ID.second != EmptyID && ID.first != EmptyID )
        {
          Triangle_3 T = Mesh2_.triangulation().triangle( pIter->first, pIter->second );
          
          Point3 Mesh2FacetCenter = Center( T );
          
          // locate this facet on the reference tree
          //
          Neighbor_search search( ReferencePointTree, Mesh2FacetCenter, nMaxNearestNeighbors );
          
          IDPair Pair1, Pair2;
          Facet BestFacet;
          bool bFound = false;
          double fBestDistance = std::numeric_limits<double>::max();
          for( Search_iterator pSearchPtr = search.begin(); pSearchPtr != search.end() && !bFound;
               ++ pSearchPtr )
          {
            int i, j;
            Locate_type LocT;
            Cell_handle ch = Mesh1_.triangulation().locate( pSearchPtr->first, LocT, i, j );
            double fDistance = FindClosestCompatibleFacet( *pIter, Mesh2FacetCenter,
                                                           ch, LocT, i, j,
                                                           fAngleThreshold,
                                                           Pair1, Pair2 );

            if( fDistance > 0 )
              bFound = true;
            BestFacet = Facet(ch, i) ;
            fBestDistance = fDistance;
          }
          
          if( bFound )
          {
            MeshToMeshIDMap.insert( Pair1 );
            MeshToMeshIDMap.insert( Pair2 );
            BndMotionInfo MotionInfo;
            MotionInfo.BndFacet1 = T;
            // MotionInfo.FacetID1 = Mesh1_.surface_index( BestFacet );
            //            MotionInfo.FacetID2 = Mesh2_.surface_index( *pIter );
            
            MotionInfo.FacetID1 = Pair1;
            MotionInfo.FacetID2 = Pair2;

            SQuaternion q1, q2;
            q1.Set( IDOrientMap1[ Pair1.first ] );
            q2.Set( IDOrientMap2[ Pair1.second ] );
            MotionInfo.Pair1 = std::make_pair( q1, q2 );
            MotionInfo.fDistance = fBestDistance;
            //             double test
            //               = LatticeSymmetry
            //               ::GetMisorientation( LatticeSymmetry
            //                                    ::CCubicSymmetry::Get(), q1, q2 );
            //             if( test > DEGREE_TO_RADIAN( 5 ) )
            //               std::cerr << "Orientation exception 1 " << std::endl;
            
            q1.Set( IDOrientMap1[ Pair2.first ] );
            q2.Set( IDOrientMap2[ Pair2.second ] );
            MotionInfo.Pair2 = std::make_pair( q1, q2 );
            
            //             if( test > DEGREE_TO_RADIAN( 5 ) )
            //               std::cerr << "Orientation exception 2 " << std::endl;

            BndMotionList.push_back( MotionInfo );
          }
        }
      }
      return BndMotionList;
    }
 
    //--------------------------------------
    //  WriteMeshToMeshIDMap
    //--------------------------------------
    void WriteMeshToMeshIDMap( const string & filename )
    {
      std::ofstream os( filename.c_str() );
      typedef typename CrossMeshIDMap::iterator Iter;

      for( Iter it = MeshToMeshIDMap.begin();
           it != MeshToMeshIDMap.end(); ++ it )
        os << it->first << " " << it->second << std::endl;

      os.close();
      
    }

    //--------------------------------------
    //  GetMeshIDSet
    //--------------------------------------
    std::set<int> GetMeshIDSet( const C3T3 & c3t3_ )
    {
      std::set<int> IDMap;
      typedef typename C3T3::Surface_index   Surface_index;
      typedef typename C3T3::Cell_handle     Cell_handle;
      typedef typename C3T3::Cell_iterator   Cell_iterator;

      for( Cell_iterator pCur = c3t3_.cells_begin();
           pCur != c3t3_.cells_end(); ++ pCur )
      {
        Cell_handle ch = pCur;
        int nID = c3t3_.subdomain_index( ch );
        IDMap.insert( nID );
      }
      return IDMap;
    }


    void GetMeshIDSets( std::set<int> & IDSet_1,  std::set<int> & IDSet_2 )
    {
      IDSet_1 = GetMeshIDSet( Mesh1_ );
      IDSet_2 = GetMeshIDSet( Mesh2_ );
    }
    
    //--------------------------------------
    //  GetUniqueMeshToMeshIDMapByMajority
    //--------------------------------------
    std::map<int, int> GetUniqueMesh1_ToMesh2_IDMapByMajority(  )
    {
      std::set<int> MeshID_1 = GetMeshIDSet( Mesh1_ );
      
      typedef typename std::map<int, int>  IDMap;
      
      IDMap CrossMeshIDMap;
      typedef typename std::set<int>::iterator  SetIter;
      typedef typename CrossMeshIDMap::iterator Iter;

      for( SetIter pID = MeshID_1.begin(); pID != MeshID_1.end(); ++ pID  )
      {
        if( *pID != -1 )
        {
          std::pair< Iter, Iter > Range = MeshToMeshIDMap.equal_range( *pID );
          IDMap RangeIDMap;
          for( Iter it = Range.first; it != Range.second; ++ it )
          {
            if( RangeIDMap.find( it->second ) != RangeIDMap.end() )
              RangeIDMap[ it->second ] ++;
            else
              RangeIDMap[ it->second ] = 0;
          }
          
          if( RangeIDMap.size() == 0 )
          {
            CrossMeshIDMap[ *pID ] = -5;
          }
          else if( RangeIDMap.size() > 1 )
          {
            int nBestID   = 0;
            int nMaxCount = 0;
            for( IDMap::iterator p = RangeIDMap.begin(); p != RangeIDMap.end(); ++ p )
            {
              if( p->second > nMaxCount )
              {
                nBestID   = p->first;
                nMaxCount = p->second;
              }
            }
            CrossMeshIDMap[ *pID ] = nBestID;
          }
          else
          {
            CrossMeshIDMap[ *pID ] = RangeIDMap.begin()->first;
          }
        }
      }
      return CrossMeshIDMap;
    }

    //--------------------------------------
    //  GetUniqueMeshToMeshIDMapByOrientation
    //--------------------------------------
    std::map<int, int> GetUniqueMesh1_ToMesh2_IDMapByOrientation(  )
    {
      std::set<int> MeshID_1 = GetMeshIDSet( Mesh1_ );
      
      typedef typename std::map<int, int>  IDMap;
      
      IDMap CrossMeshIDMap;
      typedef typename std::set<int>::iterator  SetIter;
      typedef typename CrossMeshIDMap::iterator Iter;

      for( SetIter pID = MeshID_1.begin(); pID != MeshID_1.end(); ++ pID  )
      {
        std::pair< Iter, Iter > Range = MeshToMeshIDMap.equal_range( *pID );
        std::set<int> RangeIDSet;
        for( Iter it = Range.first; it != Range.second; ++ it )
          RangeIDSet.insert( it->second );
        
        if( RangeIDSet.size() == 0 )
        {
          CrossMeshIDMap[ *pID ] = -5;
        }
        else if( RangeIDSet.size() > 1 )
        {
          int nBestID   = 0;
          double fMinMisorientation = std::numeric_limits<double>::max();
          for( typename std::set<int>::iterator p = RangeIDSet.begin();
               p != RangeIDSet.end(); ++ p )
          {
            double fMis = LatticeSymmetry
              ::GetMisorientation( LatticeSymmetry
                                   ::CCubicSymmetry::Get(),
                                   IDOrientMap1[ *pID ],
                                   IDOrientMap2[ *p ] );
            if( fMis < fMinMisorientation )
            {
              nBestID = *p;
              fMinMisorientation = fMis;
            }
          }
          CrossMeshIDMap[ *pID ] = nBestID;
        }
        else
        {
          CrossMeshIDMap[ *pID ] = *RangeIDSet.begin();
        }
      }
      return CrossMeshIDMap;
    }
    
  private:
    
    //--------------------------------------
    //  IsCompatiable
    //  Purpose:  Return true if f1 and f2 are
    //            compatible, i.e., the two sides
    //            have grains with the same orientations
    //
    //  Pair1, Pair2:  Pairs mapping IDs between Mesh1 -> Mesh2
    //
    //--------------------------------------
    bool IsCompatible( const C3T3 & c3t3_1, const C3T3 & c3t3_2,
                       const Facet & f1, const Facet & f2, double fThreshold,
                       double & Angle1, double & Angle2,
                       IDPair & Pair1, IDPair & Pair2  )
    {
      Surface_index ID1 = c3t3_1.surface_index( f1 );
      Surface_index ID2 = c3t3_2.surface_index( f2 );


      if( ID1.first == -1 || ID2.first == -1
          || ID1.second == -1 || ID2.second == -1 )
        return false;
      
      SMatrix3x3 f1_a, f1_b, f2_a, f2_b;

      f1_a = IDOrientMap1[ ID1.first ];
      f1_b = IDOrientMap1[ ID1.second ];

      f2_a = IDOrientMap2[ ID2.first ];
      f2_b = IDOrientMap2[ ID2.second ];

      Angle1 = LatticeSymmetry::GetMisorientation( LatticeSymmetry::CCubicSymmetry::Get(),
                                                   f1_a, f2_a );
      Angle2 = LatticeSymmetry::GetMisorientation( LatticeSymmetry::CCubicSymmetry::Get(),
                                                   f1_b, f2_b );
          
      if ( Angle1 < fThreshold  && Angle2 < fThreshold  )
      {
        Pair1 = std::make_pair( ID1.first,  ID2.first );
        Pair2 = std::make_pair( ID1.second, ID2.second );
        return true;
      }
      Angle1 = LatticeSymmetry::GetMisorientation( LatticeSymmetry::CCubicSymmetry::Get(),
                                                   f1_a, f2_b );
      Angle2 = LatticeSymmetry::GetMisorientation( LatticeSymmetry::CCubicSymmetry::Get(),
                                                   f1_b, f2_a );
      
      if(  Angle1 < fThreshold  && Angle2 < fThreshold  )
      {
        Pair1 = std::make_pair( ID1.first,  ID2.second );
        Pair2 = std::make_pair( ID1.second, ID2.first );
        return true;
      }
      
      return false;
    }
    
    //--------------------------------------
    //   FindClosestCompatibleFacet
    //   Note:  Search center may not be at the same location as the facet f
    //
    //--------------------------------------
    double FindClosestCompatibleFacet( const Facet & f,
                                       const Point3 & SearchCenter,
                                       const Cell_handle & ch,
                                       Locate_type LocT, int i, int j,
                                       double fThreshold,
                                       IDPair & PairToAdd1,
                                       IDPair & PairToAdd2 )
    {
      double Angle1, Angle2;
      switch( LocT )
      {
        case C3T3::Triangulation::FACET:  // handle at the end
          {
            if( IsCompatible( Mesh1_, Mesh2_, Facet(ch, i),
                              f, fThreshold, Angle1, Angle2,
                              PairToAdd1, PairToAdd2 ) )
            {
              double fDistance = PointToFacetDistance( Mesh1_, SearchCenter, ch, i );
              return fDistance;
            }
            return -1;
          }
        case C3T3::Triangulation::CELL :
          {
            // std::cout << "CELL " <<std::endl;
            double fDistance = std::numeric_limits<double>::max();
            bool bFound = false;
            double fBestAngle = std::numeric_limits<double>::max();
            
            for( int n = 0; n < 4; n ++ )
            {
              if( Mesh1_.is_in_complex( ch, n ) )
              {
                IDPair Pair1, Pair2;
                if( IsCompatible( Mesh1_, Mesh2_, Facet(ch, i), f,
                                  fThreshold, Angle1, Angle2,
                                  Pair1, Pair2 ) )
                {
                  double d = PointToFacetDistance( Mesh1_, SearchCenter, ch, n );
                  if( (Angle1 + Angle2) < fBestAngle )  // want best min misorientation
                  {
                    fBestAngle = Angle1 + Angle2;                    
                    fDistance = d;
                    PairToAdd1 = Pair1;
                    PairToAdd2 = Pair2;
                    bFound = true;
                  }
                }
              }
            }
            
            if( bFound )
            {
              return fDistance;
            }
            else
              return -1;
          }

        case C3T3::Triangulation::VERTEX:
        case C3T3::Triangulation::EDGE  :
        default:
          {
            std::cerr << "Exception: Degenerate triangle " << std::endl;
            return -10000;
          }
      }// end switch
    }
    
    //--------------------------------------
    //  GetCenterOfMass
    //--------------------------------------
    template< class PtrIter >
    SVector3 GetCenterOfMass( PtrIter pCur, PtrIter pEnd )
    {
      SVector3 CenterOfMass(0, 0,0);
      int nElements = pEnd - pCur;
      for( ; pCur != pEnd; ++ pCur )
        CenterOfMass += SVector3( pCur->x(), pCur->y(), pCur->z() );
        
      CenterOfMass /= static_cast<double>( nElements );
      return CenterOfMass;
    }
    

    //--------------------------------------
    //  ExtractBoundaryPoints
    //   Purpose:  Extract or vertices on grain boundaries
    //             and sample surface.
    //--------------------------------------
    template< class OutputIterator >
    void ExtractBoundaryPoints( OutputIterator OutIter,
                                const C3T3 & c3t3,
                                const Subdomain_index & EmptySpaceID )
    {
      std::cout << "Extracting Boundary Points " << std::endl;

      //---------------------------------
      //  Need to change to select boundar vertices only
      //---------------------------------
      // for( Finite_facet_iterator pIter = c3t3.triangulation().finite_facets_begin();
//            pIter != c3t3.triangulation().finite_facets_end(); ++ pIter )
      for( Facet_iterator pIter = c3t3.facets_begin();
           pIter != c3t3.facets_end(); ++ pIter )

      {
        Surface_index SurfaceID = c3t3.surface_index( *pIter );

        if( SurfaceID.first != EmptySpaceID && SurfaceID.second != EmptySpaceID )
        {
          Triangle_3 T = c3t3.triangulation().triangle( *pIter );
          *OutIter =  Center( T );
          ++ OutIter;
        }
      }
    }

    //-------------------------------------------------
    //  Center
    //   Purpose:  Return the center of the triangle.
    //-------------------------------------------------
    Point3 Center( const Triangle_3 & T )
    {
      Point3 p( ( T[0].x() + T[1].x() + T[2].x() ) / 3.,
                ( T[0].y() + T[1].y() + T[2].y() ) / 3.,
                ( T[0].z() + T[1].z() + T[2].z() ) / 3.     );
      return p;
    }
    
    //-------------------------------------------------
    //  PointToFacetDistance
    //   Purpose:  Calculate perpendicular distance between
    //             point and facet
    //-------------------------------------------------
    double PointToFacetDistance( const C3T3 & c3t3,
                                 const Point3 & p,
                                 const Cell_handle  &ch,
                                 int OppVertex )
    {
      Triangle_3 T = c3t3.triangulation().triangle( ch, OppVertex );
      Point3 PlanePoint = T.supporting_plane().projection( p );
      return  CGAL::sqrt( CGAL::squared_distance(PlanePoint, p ) );
    }
    
    std::vector<Point3> MeshSurfacePointList1, MeshSurfacePointList2;
    
    std::vector<SMatrix3x3> & IDOrientMap1;
    std::vector<SMatrix3x3> & IDOrientMap2;

    
    const C3T3 & Mesh1_;
    const C3T3 & Mesh2_;

    Iso_cuboid      MeshBBox1;
    Iso_cuboid      MeshBBox2;
    CrossMeshIDMap  MeshToMeshIDMap;
    Tree ReferencePointTree;
  };


}

#endif
