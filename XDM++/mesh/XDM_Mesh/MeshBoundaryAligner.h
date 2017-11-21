//----------------------------------------------------
//  MeshBoundaryAligner.h
//
//  Author:   Frankie Li
//  Purpose:  A generalized method that aligns two meshes
//            to minimize the total boundary mismatch.
//----------------------------------------------------


#ifndef MESH_BOUNDARY_ALIGNER_H
#define MESH_BOUNDARY_ALIGNER_H

#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <limits>
#include <CGAL/bounding_box.h>
#include "MeshUtilities.h"
#include <CGAL/Cartesian_converter.h>
namespace Pandora
{

  //-------------------------------------------------------------
  //  BoundaryAligner
  //
  //  (Should use BoundaryAnalysis class need refactor)
  //-------------------------------------------------------------
  template< class C3T3 >
  class BoundaryAligner
  {
  public:
 
    struct Kernel: public CGAL::Exact_predicates_inexact_constructions_kernel {};
    
    typedef typename C3T3::Triangulation::Vertex_iterator Vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertex_iterator;
    typedef typename C3T3::Triangulation::Vertex          Vertex;
    typedef typename C3T3::Triangulation::Geom_traits::Bare_point Point3;
    typedef typename C3T3::Triangulation::Geom_traits::RT RT;
    typedef typename C3T3::Cell_handle                    Cell_handle;
    typedef typename C3T3::Subdomain_index                Subdomain_index;
    typedef typename C3T3::Triangulation::Locate_type     Locate_type;
    
    typedef typename C3T3::Triangulation::Geom_traits::Triangle_3  Triangle_3;
    
    typedef typename Kernel::Iso_cuboid_3                 Iso_cuboid;
    

    typedef GeneralLib::SQuaternion                       SQuaternion;
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

    
  public:

    BoundaryAligner( const C3T3 & Mesh1, const C3T3 & Mesh2 ):
      Mesh1_( Mesh1 ), Mesh2_( Mesh2 )
    {}



    //--------------------------------------
    //  GetOVerlapZRegion
    //   Purpose:  Specialized overlap region
    //             selector for the typical HEDM
    //             geometry.  The region selected
    //             is equivilent to the "z-layers"
    //             overlapping between the two meshes.
    //             Note that this is mostly an approximation
    //             since only the z direction is considered.
    //--------------------------------------
    Iso_cuboid GetOverlapZRegion( const Iso_cuboid & b1,
                                  const Iso_cuboid & b2,
                                  const SMatrix3x3 & R,
                                  const SVector3 & T )
    {

      double fMinZ = std::numeric_limits<double>::max();
      double fMaxZ = std::numeric_limits<double>::min();
      
      for( int i = 0; i < 8; i ++ )
      {
        Point3 pTmp = Details::InverseTransformPoint( b1.vertex(i), R, T );
        fMinZ = std::min( pTmp.z(), fMinZ );
        fMaxZ = std::max( pTmp.z(), fMaxZ );
      }

      Iso_cuboid NewB2( Point3( b2.min().x(), b2.min().y(), fMinZ),
                        Point3( b2.max().x(), b2.max().y(), fMaxZ) );
      return NewB2;
    }
    //--------------------------------------
    //  GetOverlappingRandomizedPoints
    //
    //  parameters: b2 - bounding box that surrounds all points from pFirst -> pEnd
    //
    //--------------------------------------
    template< class PtrIterator, class BBox >
    PtrIterator RandomizeOverlapSamplePoints( PtrIterator pFirst,
                                              PtrIterator pEnd,
                                              const BBox & b1,
                                              const BBox & b2,
                                              const SMatrix3x3 & R,
                                              const SVector3   & T,
                                              int MaxPoints )
    {
      Iso_cuboid OverlapRegion = GetOverlapZRegion( b1, b2, R, T );
      PtrIterator pLastOverlapPoint = std::remove_if( pFirst, pEnd,
                                                      NotInBox( OverlapRegion )  );  // not actually deleted
      int nOverlapPoints = pLastOverlapPoint - pFirst;
      
      std::random_shuffle( pFirst, pLastOverlapPoint );
      int nPointsUsed = std::min( MaxPoints, nOverlapPoints );
      
      //      std::cout << "NPointsUsed " << nOverlapPoints << std::endl;
      return pFirst + nPointsUsed;
    }
    
    //--------------------------------------
    //  InitializeBySurface
    //--------------------------------------
    void InitializeBySurface()
    {
      MeshSurfacePointList1.clear();
      MeshSurfacePointList2.clear();
      Pandora::Details::ExtractSurfacePoints( std::back_inserter( MeshSurfacePointList1 ),
                                              Mesh1_, -1 );
      Pandora::Details::ExtractSurfacePoints( std::back_inserter( MeshSurfacePointList2 ),
                                              Mesh2_, -1 );
    }
    
    //--------------------------------------
    //  InitializeByBoundary
    //--------------------------------------
    void InitializeByBoundary()
    {
      MeshSurfacePointList1.clear();
      MeshSurfacePointList2.clear();

      Pandora::Details::ExtractBoundaryPoints( std::back_inserter( MeshSurfacePointList1 ),
                                               Mesh1_ );
      Pandora::Details::ExtractBoundaryPoints( std::back_inserter( MeshSurfacePointList2 ),
                                               Mesh2_ );
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
    
    //--------------------------------------
    //  GetCenterOfMassTranslation
    //   Purpose:  Return COM translation vector from Mesh 1 -> Mesh 2
    //--------------------------------------
    SVector3 GetCenterOfMassTranslation()
    {
      assert( MeshSurfacePointList1.size() > 0 && MeshSurfacePointList2.size() > 0 );
      SVector3 CenterOfMass1(0, 0, 0);
      SVector3 CenterOfMass2(0, 0, 0);

      CenterOfMass1 = GetCenterOfMass( MeshSurfacePointList1.begin(),
                                       MeshSurfacePointList1.end() );
      CenterOfMass2 = GetCenterOfMass( MeshSurfacePointList2.begin(),
                                       MeshSurfacePointList2.end() );
      return CenterOfMass2 - CenterOfMass1;
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
    //  RemoveEdgePointByShift
    //    Remove all points that 
    //
    //--------------------------------------
    void RemoveEdgePointsByShift( double zShift )
    {
      MeshBBox1 = CGAL::bounding_box( MeshSurfacePointList1.begin(),
                                      MeshSurfacePointList1.end()  );
      MeshBBox2 = CGAL::bounding_box( MeshSurfacePointList2.begin(),
                                      MeshSurfacePointList2.end()  );

      std::cout << "-------------------------------" << std::endl;
      std::cout << " Box 1 " << MeshBBox1.min() << " " << MeshBBox1.max() << std::endl;
      std::cout << " Box 2 " << MeshBBox2.min() << " " << MeshBBox2.max() << std::endl;
      MeshBBox1 = Iso_cuboid(
                             Point3( MeshBBox1.xmin(),
                                     MeshBBox1.ymin(),
                                     MeshBBox1.zmin() + zShift ),
                             Point3( MeshBBox1.xmax(),
                                     MeshBBox1.ymax(),
                                     MeshBBox1.zmax() - zShift )
                             );

      MeshBBox2 = Iso_cuboid(
                             Point3( MeshBBox2.xmin(),
                                     MeshBBox2.ymin(),
                                     MeshBBox2.zmin() + zShift ),
                             Point3( MeshBBox2.xmax(),
                                     MeshBBox2.ymax(),
                                     MeshBBox2.zmax() - zShift )
                             );

      std::cout << " Box 1 " << MeshBBox1.min() << " " << MeshBBox1.max() << std::endl;
      std::cout << " Box 2 " << MeshBBox2.min() << " " << MeshBBox2.max() << std::endl;
      std::cout << "-------------------------------" << std::endl;

      std::cout << "Before clipping " << MeshSurfacePointList1.size() << " " << MeshSurfacePointList2.size() << std::endl;
      MeshSurfacePointList1.erase( std::remove_if( MeshSurfacePointList1.begin(),
                                                   MeshSurfacePointList1.end(),
                                                   NotInBox( MeshBBox1 ) ),
                                   MeshSurfacePointList1.end() );
      MeshSurfacePointList2.erase( std::remove_if( MeshSurfacePointList2.begin(),
                                                   MeshSurfacePointList2.end(),
                                                   NotInBox( MeshBBox2 ) ),
                                   MeshSurfacePointList2.end() );

      std::cout << "After clipping " << MeshSurfacePointList1.size() << " " << MeshSurfacePointList2.size() << std::endl;
      //       std::remove_if( MeshSurfacePointList2.begin(), MeshSurfacePointList2.end(),
      //                       NotInBox( MeshBBox2 ) );
      
      MeshBBox1 = CGAL::bounding_box( MeshSurfacePointList1.begin(),
                                      MeshSurfacePointList1.end()  );
      MeshBBox2 = CGAL::bounding_box( MeshSurfacePointList2.begin(),
                                      MeshSurfacePointList2.end()  );
    }
    //--------------------------------------
    //  AlignBySampleSurface
    //   This is designed mostly cylindrical
    //   sample alignment
    //
    //   This is clearly a local optimization step with
    //   T = 0 Monte Carlo
    //--------------------------------------
    std::pair<SMatrix3x3, SVector3>
    AlignSample( double & fBestCost, const SMatrix3x3 & StartRotation,
                 const SVector3 &   StartTranslation,
                 const SVector3 &   TransStepSize,
                 int nMaxIterations,
                 int nTestPoints, float SO3_Radius, int nFileID = 0 )
    {
      assert( MeshSurfacePointList1.size() > 0 && MeshSurfacePointList2.size() > 0 );
      
      ReferencePointTree = Tree( MeshSurfacePointList1.begin(),
                                 MeshSurfacePointList1.end() );

      MeshBBox1 = CGAL::bounding_box( MeshSurfacePointList1.begin(),
                                      MeshSurfacePointList1.end()  );
      MeshBBox2 = CGAL::bounding_box( MeshSurfacePointList2.begin(),
                                      MeshSurfacePointList2.end()  );

      assert( (MeshSurfacePointList2.end() - MeshSurfacePointList2.begin() ) > nTestPoints );
      
      typedef typename std::vector<Point3>::iterator PointIter;
      PointIter pLastOverlappingPoint = RandomizeOverlapSamplePoints( MeshSurfacePointList2.begin(),
                                                                      MeshSurfacePointList2.end(),
                                                                      MeshBBox1, MeshBBox2,
                                                                      StartRotation,
                                                                      StartTranslation, nTestPoints );;
      
      std::cout << " STEP SIZE " << TransStepSize << std::endl;
      std::cout << " STEP SIZE " << SO3_Radius << std::endl;

      UniformGrid::CQuaternionGrid RandSO3Gen;
      std::vector<SQuaternion> RandomLocalPoints
        = RandSO3Gen.GetRandomLocalGrid( SO3_Radius, nMaxIterations );
      
      fBestCost = GetMisalignCostByTree( StartRotation, StartTranslation,
                                         MeshSurfacePointList2.begin(),
                                         pLastOverlappingPoint );
      fBestCost /= double( pLastOverlappingPoint - MeshSurfacePointList2.begin() );

      SMatrix3x3 Rotation    = StartRotation;
      SVector3   Translation = StartTranslation;
      std::stringstream ss;
      ss << "Cost." << nFileID << ".txt";
      std::ofstream outfile(ss.str().c_str() );
      for( int i = 0; i < nMaxIterations; i ++ )
      {
        float fDx = RandGen.GetRandomVariable( -TransStepSize.m_fX, TransStepSize.m_fX );
        float fDy = RandGen.GetRandomVariable( -TransStepSize.m_fY, TransStepSize.m_fY );
        float fDz = RandGen.GetRandomVariable( -TransStepSize.m_fZ, TransStepSize.m_fZ );
        
        SMatrix3x3 CurrentRot = RandomLocalPoints[i].GetRotationMatrix3x3() * Rotation ;
        SVector3   CurrentPos = Translation + SVector3( fDx, fDy, fDz );
        pLastOverlappingPoint = RandomizeOverlapSamplePoints( MeshSurfacePointList2.begin(),
                                                              MeshSurfacePointList2.end(),
                                                              MeshBBox1, MeshBBox2, CurrentRot,
                                                              CurrentPos, nTestPoints );
        double cost = GetMisalignCostByTree( CurrentRot, CurrentPos,
                                             MeshSurfacePointList2.begin(),
                                             pLastOverlappingPoint );
        int nOverlapPoints = pLastOverlappingPoint - MeshSurfacePointList2.begin();
        cost /= double( nOverlapPoints );
        
        outfile << " " << nOverlapPoints
                << " " << nTestPoints
                << " " << cost << " "
                << fBestCost << " "<< std::endl;
        if( fBestCost > cost )
        {
          fBestCost   = cost;
          Rotation    = CurrentRot;
          Translation = CurrentPos;
        }
      }
      outfile.close();
      std::cout << "Translation [internal] " << Translation << std::endl;
      std::cout << "Rotation [internal  " << GeneralLib::RadianToDegree( Rotation.GetEulerAngles() ) << std::endl;

      std::stringstream ss2;
      ss2 << "SinglePoint." << nFileID << ".csv";
      WriteMisalignCostByTree( Rotation, Translation,MeshSurfacePointList2.begin(),
                               MeshSurfacePointList2.end(),  ss2.str() );
      return std::make_pair( Rotation, Translation );
    }

    
    
  private:

    
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
    //  GetMisalignCostByTree
    //    Purpose:  Calculate the misalignment cost by tree
    //
    //--------------------------------------
    template< class PointIter >
    double GetMisalignCostByTree( const SMatrix3x3 & Rotation,
                                  const SVector3   & Translation,
                                  PointIter pCur, PointIter pEnd )
    {
      double cost = 0;
      //      std::vector<double> CostVector( pEnd - pCur );
      
      for(int n = 0; pCur != pEnd; ++ pCur, ++n )
      {
        Point3 p = Details::TransformPoint( *pCur, Rotation, Translation );
        Neighbor_search search( ReferencePointTree, p, 1 );
        double CurrentDistance = 0;
        if( search.begin() != search.end() )
        {
          CurrentDistance = search.begin()->second;
          //CostFile << CurrentDistance << " " << p <<  std::endl;
          cost += CurrentDistance;
        }
        else
        {
          std::cerr << "Exception " << std::endl;
          std::cout << "WTF ---------------- " << std::endl;
          
          cost += 500000000;   // else, extra penalty, gotta think about this
        }

      }

//       double StdSq = 0;
//       double MeanCost = cost / double( CostVector.size() );
//       for( int i = 0; i < CostVector.size(); i ++ )
//         StdSq += (CostVector[i] - cost) * ( CostVector[i] - cost );

//       StdSq /=  double( CostVector.size() );
      
//       exit(0);
      return cost;
    }

    //--------------------------------------
    //  GetMisalignCostByTree
    //    Purpose:  Calculate the misalignment cost by tree
    //
    //--------------------------------------
    template< class PointIter >
    void WriteMisalignCostByTree( const SMatrix3x3 & Rotation,
                                  const SVector3   & Translation,
                                  PointIter pCur, PointIter pEnd,
                                  const string & filename )
    {
      std::ofstream CostFile( filename.c_str() );

      for(int n = 0; pCur != pEnd; ++ pCur, ++n )
      {
        Point3 p = Details::TransformPoint( *pCur, Rotation, Translation );
        Neighbor_search search( ReferencePointTree, p, 1 );
        if( search.begin() != search.end() )
        {
          double CurrentDistance = search.begin()->second;
          CostFile << CurrentDistance << " " << p <<  std::endl;
        }
      }
      CostFile.close();
    }


    //--------------------------------------
    //  GetMisalignCostByMesh
    //    Purpose:  Calculate the misalginment cost by
    //              directly using the mesh.  Reference
    //              mesh is Mesh1
    //--------------------------------------
    template< class PointIter >
    double GetMisalignCostByMesh( const SMatrix3x3 & Rotation,
                                  const SVector3   & Translation,
                                  PointIter pCur, PointIter pEnd )
    {
      double cost = 0;
      for(; pCur != pEnd; ++ pCur )
      {
        Point3 p = Details::TransformPoint( *pCur, Rotation, Translation );
        Cell_handle ch;
        int i, j;
        Locate_type LocT;
        ch = Mesh1_.triangulation().locate( p, LocT, i, j );
        cost += PointToCellDistance( Mesh1_, ch, p, LocT, i, j );
      }
      return cost;
    }
    
    //--------------------------------------
    //  ExtractBoundaryPoints
    //   Purpose:  Extract or vertices on grain boundaries
    //             and sample surface.
    //--------------------------------------
//     template< class OutputIterator >
//     void ExtractBoundaryPoints( OutputIterator OutIter,
//                                 const C3T3 & c3t3,
//                                 const Subdomain_index & EmptySpaceID )
//     {
//       std::cout << "Extracting Boundary Points " << std::endl;

//       //---------------------------------
//       //  Need to change to select boundar vertices only
//       //---------------------------------
//       for( Finite_vertex_iterator pIter = c3t3.triangulation().finite_vertices_begin();
//            pIter != c3t3.triangulation().finite_vertices_end(); ++ pIter )
//       {
//         if( c3t3.in_dimension( pIter) <= 2 )
//         {
//           std::vector<Cell_handle> IncidentCellList;
//           c3t3.triangulation().incident_cells( pIter,
//                                                std::back_inserter( IncidentCellList ) );
//           bool bIsSpaceBoundary = false;
//           for( int i = 0; i < IncidentCellList.size(); i ++ )
//           {
//             if( c3t3.subdomain_index( IncidentCellList[i] ) <= 0 )
//               bIsSpaceBoundary = true;
//           }
          
//           if( !bIsSpaceBoundary )
//           {
//             *OutIter =  pIter->point();
//             ++ OutIter;
//           }
//         }
//       }
//     }
    
 //    //--------------------------------------
//     //  ExtractSurfacePoints
//     //    Parameter:  c3t3    - mesh
//     //                dZLimit - the offset from the top
//     //                          and the bottom that delineates
//     //                          the section of volume used for
//     //                          registeration.
//     //    Test:  Output and visualize
//     //--------------------------------------
//     template< class OutputIterator >
//     void ExtractSurfacePoints( OutputIterator OutIter,
//                                const C3T3 & c3t3,
//                                const Subdomain_index & EmptySpaceID )
//     {
//       std::cout << "Extracting Surface Points " << std::endl;
//       //---------------------------------
//       //  Need to change to select boundar vertices only
//       //---------------------------------
//       for( Finite_vertex_iterator pIter = c3t3.triangulation().finite_vertices_begin();
//            pIter != c3t3.triangulation().finite_vertices_end(); ++ pIter )
//       {
//         if( c3t3.in_dimension( pIter) <= 2 )
//         {
//           std::vector<Cell_handle> IncidentCellList;
//           c3t3.triangulation().incident_cells( pIter,
//                                                std::back_inserter( IncidentCellList ) );
//           bool bAdjacentSample = false;
//           bool bAdjacentSpace  = false;
//           for( int i = 0; i < IncidentCellList.size(); i ++ )
//           {
//             if( c3t3.subdomain_index( IncidentCellList[i] )
//                    == EmptySpaceID
//                 || c3t3.subdomain_index( IncidentCellList[i] ) == 0 )  // this is a hack because
//               bAdjacentSpace = true;                                   // of the inconsistency
//                                                                        // in empty space declaration
//             if( c3t3.subdomain_index( IncidentCellList[i] ) > 0  )
//               bAdjacentSample = true;
            
//           }
//           if( bAdjacentSample && bAdjacentSpace )
//           {
//             *OutIter =  pIter->point();
//             ++ OutIter;
//           }
//         }
//       }
//     }

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
      return  CGAL::squared_distance(PlanePoint, p );
    }
    
    //-------------------------------------------------
    //  PointToCellDistance
    //   Need to handle the returned cell based on
    //     CELL    - Find nearest facet that is a boundary
    //     FACET   - Find projected distance
    //     EDGE    - By definition (projected distance) = 0
    //     VERTEX  - By definition (vertex-on-vertex)   = 0
    //-------------------------------------------------
    double PointToCellDistance( const C3T3 & c3t3,
                                const Cell_handle & ch,
                                const Point3 & p,
                                Locate_type LocT, int i, int j)
    {
      switch(LocT)
      {
        case C3T3::Triangulation::VERTEX:
        case C3T3::Triangulation::EDGE  :
          return 0;
        case C3T3::Triangulation::FACET:  // handle at the end
          {
            double fDistance = PointToFacetDistance( c3t3, p, ch, i );
            std::cout << "  PointToFacetDistance " << fDistance << std::endl;
            //          return PointToFacetDistance( c3t3, p, ch, i );
            return fDistance;
          }
        case C3T3::Triangulation::CELL :
          {
            double fDistance = std::numeric_limits<double>::max();
            bool bFound = false;
            for( int n = 0; n < 4; n ++ )
            {
              if( c3t3.is_in_complex( ch, n ) )
              {
                fDistance = std::min( fDistance,
                                      PointToFacetDistance( c3t3, p, ch, n ) );
                bFound = true;
              }
            }
            
            if( bFound )
            {
              //              std::cout << "  PointToFacetDistance " << fDistance << std::endl;
              return fDistance;
            }
            // else go to default
          }
        default:
          {
            Neighbor_search search( ReferencePointTree, p, 1 );
            if( search.begin() != search.end() )
            {
              double fDistance =  CGAL::squared_distance(search.begin()->first,  p ); 
              //  std::cout << "  PointToFacetDistance " << fDistance << std::endl;
              return fDistance;
            }
            else
              std::cerr << "Search exception after PointToCellDistance " << std::endl;
            //std::cout << "  PointToFacetDistance " << 10000 << std::endl;
            return 10000;
          }
      }
      
    };
    
    
    std::vector<Point3> MeshSurfacePointList1, MeshSurfacePointList2;
    
    const C3T3 & Mesh1_;
    const C3T3 & Mesh2_;

    Iso_cuboid MeshBBox1;
    Iso_cuboid MeshBBox2;
    
    Tree ReferencePointTree;
  };


}


#endif
