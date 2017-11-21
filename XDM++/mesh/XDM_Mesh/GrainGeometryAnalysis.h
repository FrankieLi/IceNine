//----------------------------------------------------
//  TripleLineAnalysis.h
//
//  Author:   Frankie Li
//  Purpose:  Estimate triple line prop
//
//----------------------------------------------------

#ifndef GRAIN_GEOMETRY_H
#define GRAIN_GEOMETRY_H
#include "MeshUtilities.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Lapack/Linear_algebra_lapack.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/IO/File_medit.h>
#include "SimpleMeshVTKWriter.h"
#include "3dMath.h"
#include <map>
#include <vector>
#include <CGAL/Bbox_3.h>
#include <CGAL/Plane_3.h>
#include <limits>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>

namespace Pandora
{
  namespace Details
  {
    struct GeomProp
    {
      double        Volume;
      double        MeanWidth;
      double        EdgeLength;
      int           nIsInterior;
      CGAL::Bbox_3  BBox;
    };

    template<class IterT, class T>
    struct DereferenceReturn
    {
      T operator() ( IterT p ) const { return * p; }
    };
    //-----------------------------------------------------------------------------
    //  Calculate Geometry of Grain
    //
    //-----------------------------------------------------------------------------
    template< class C3T3>
    class GrainGeometry
    {
    public:
      mutable std::ofstream Normal_os;
      mutable std::ofstream Edge_os;
//       mutable std::ofstream Tet_os;
//       mutable bool DEBUG_FLAG;
      
      struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};  // move to template later
      struct K2: public CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt {};  // move to template later
      
//       typedef CGAL::Cartesian<K2> GeomKernel;
//       typedef typename GeomKernel::Plane_3        Plane3;
//       typedef typename GeomKernel::Direction_3    Direction3;

      typedef CGAL::Cartesian<double> GeomKernel;
      typedef typename CGAL::Plane_3< GeomKernel >        Plane3;
      typedef typename CGAL::Direction_3< GeomKernel >    Direction3;
      typedef typename CGAL::Vector_3< GeomKernel >        Vector3;

      typedef typename C3T3::Triangulation::Geom_traits::Bare_point Point3;
      typedef typename C3T3::Triangulation::Geom_traits::Tetrahedron_3 Tetrahedron3;
      typedef typename C3T3::Triangulation::Geom_traits::FT FT;
      typedef K  Kernel;

      typedef typename C3T3::Triangulation::Facet_circulator Facet_circulator;
      typedef typename C3T3::Triangulation::Cell_circulator  Cell_circulator;
      typedef typename C3T3::Cell_handle Cell_handle;

      typedef typename C3T3::Subdomain_index Subdomain_index;
      typedef typename C3T3::Surface_index   Surface_index;
      typedef typename C3T3::Triangulation::Vertex_handle Vertex_handle;
      typedef typename C3T3::Facet                  Facet;
      typedef typename C3T3::Triangulation::Edge    Edge;
      
      typedef typename C3T3::Facet_iterator  Facet_iterator;
      typedef typename C3T3::Cell_iterator   Cell_iterator;
      
      typedef std::multimap< Subdomain_index, Cell_handle > GrainTetMap;
      typedef std::multimap< Subdomain_index, Facet >       GrainFacetMap;

      typedef CGAL::Bbox_3 BBox;

      typedef std::map< int, std::set<int> >         IDNgbMapT;
      typedef std::multimap< Subdomain_index, Cell_handle >  GrainTetMapT;
      typedef std::multimap< Subdomain_index, Edge >         IDToEdgeMapT;
      
    private:
      const C3T3 & c3t3_;
      
    public:

      GrainGeometry( const C3T3 & c3t3 ): c3t3_( c3t3 )
      {
        //        Edge_os.open("Edge.txt");
        //       Normal_os.open("NormalDebug.txt");
//         Tet_os.open("TetDebug.txt");
      }
      
      //-----------------------------------------------------------------------------
      //  Calculate the center of mass of the list of tetrahedrons
      //-----------------------------------------------------------------------------
      template< class CellHandleIter, class HandleFn >
      Point3 CenterOfMass( CellHandleIter pFirst, CellHandleIter pEnd, HandleFn GetHandle ) const
      {
        FT x_c = 0;
        FT y_c = 0;
        FT z_c = 0;
        FT vTotal = 0;

        for( ; pFirst != pEnd; ++ pFirst )
        {
          Cell_handle ch = GetHandle( pFirst );
          Point3 p1 = ch->vertex(0)->point();
          Point3 p2 = ch->vertex(1)->point();
          Point3 p3 = ch->vertex(2)->point();
          Point3 p4 = ch->vertex(3)->point();
          Tetrahedron3 t ( p1, p2, p3, p4 );
          FT v = t.volume();
          FT x = p1.x() + p2.x() + p3.x() + p4.x();
          FT y = p1.y() + p2.y() + p3.y() + p4.y();
          FT z = p1.z() + p2.z() + p3.z() + p4.z();

          x_c += x * v / 4. ;
          y_c += y * v / 4. ;
          z_c += z * v / 4. ;
          vTotal += v;
        }
        x_c /= vTotal;
        y_c /= vTotal;
        z_c /= vTotal;

        return Point3( x_c, y_c, z_c );
      }

      //---------------------------------
      //  IsConvex
      //
      //  Assumed that P1 is oriented outward away from ID
      //---------------------------------
      bool IsConvex( const Plane3 & P1,
                     const Edge & E2, const Facet & F2 ) const
      {
        typedef typename CGAL::Point_3< GeomKernel >    SPoint3;

        int FarVertex = -1;
        Vertex_handle v1 = E2.get<0>()->vertex( E2.get<1>() );
        Vertex_handle v2 = E2.get<0>()->vertex( E2.get<2>() );
        Cell_handle ch = F2.first;
        for( int i = 0; i < 4; i ++ )
        {
          if( i != F2.second )
            if( ch->vertex( i ) != v1
                && ch->vertex( i ) != v2 )
              FarVertex = i;
        }

        if( FarVertex == -1 )
        {
          std::cout << "Can't find far vertex, don't know what to do" << std::endl;
          exit(0);
        }
        SPoint3 pTest( ch->vertex( FarVertex )->point().x(),
                       ch->vertex( FarVertex )->point().y(),
                       ch->vertex( FarVertex )->point().z() );
        return P1.has_on_negative_side( pTest ); 
      }
      
      //-----------------------------------------------------------------------------
      //  Calculate the center of mass of the list of tetrahedrons
      //-----------------------------------------------------------------------------
      template< class CellHandleIter, class HandleFn >
      BBox BoundingVolume( CellHandleIter pFirst, CellHandleIter pEnd, HandleFn GetHandle ) const
      {
        if( pFirst == pEnd )
          return BBox( 0, 0, 0, 0, 0, 0 );

        FT xMax = std::numeric_limits<FT>::min();
        FT yMax = std::numeric_limits<FT>::min();
        FT zMax = std::numeric_limits<FT>::min();

        FT xMin = std::numeric_limits<FT>::max();
        FT yMin = std::numeric_limits<FT>::max();
        FT zMin = std::numeric_limits<FT>::max();
        
        for( ; pFirst != pEnd; ++ pFirst )
        {
          Cell_handle ch = GetHandle( pFirst );
          for( int i = 0; i < 4; i ++ )
          {
            if( ! c3t3_.triangulation().is_infinite(ch->vertex(i)) )
            {
              Point3 p = ch->vertex(i)->point();
              xMax = std::max( xMax, p.x() );
              yMax = std::max( yMax, p.y() );
              zMax = std::max( zMax, p.z() );
              
              xMin = std::min( xMin, p.x() );
              yMin = std::min( yMin, p.y() );
              zMin = std::min( zMin, p.z() );
            }
          }
        }
        return BBox( xMin, yMin, zMin,
                     xMax, yMax, zMax );
      }

      //---------------------------------
      // GetOrientedPlane
      //---------------------------------
      Plane3 GetOrientedPlane( Cell_handle ch, int nOppIndex,
                               const Subdomain_index & ID ) const
      {
        typedef typename CGAL::Point_3< GeomKernel >    SPoint3;
        typedef typename C3T3::Triangulation::Geom_traits::Triangle_3 Triangle3;
        Triangle3 T = c3t3_.triangulation().triangle( ch, nOppIndex );   // by definition, points in-ward
        SPoint3 sp1( T[0].x(), T[0].y(), T[0].z() );
        SPoint3 sp2( T[1].x(), T[1].y(), T[1].z() );
        SPoint3 sp3( T[2].x(), T[2].y(), T[2].z() );
        Plane3 P( sp1, sp2, sp3 );

        if( c3t3_.subdomain_index( ch ) == ID )
          P = P.opposite();
        return P;
      }


      //-----------------------------------------------------------------------------
      // TurningAngle
      //-----------------------------------------------------------------------------
      FT TurningAngle( const Edge edge, const Facet & f1, const Facet & f2,
                       const Subdomain_index & ID ) const
      {
        std::vector< int > Vertices1;
        for( int i = 0; i < 4; i ++ )
          if( i != f1.second )
            Vertices1.push_back( i );

        Cell_handle c1 = f1.first;
        Plane3 P1 = GetOrientedPlane( c1, f1.second, ID );
        Point3 pTest1 = c1->vertex( Vertices1[0] )->point();

        std::vector< int > Vertices2;
        for( int i = 0; i < 4; i ++ )
          if( i != f2.second )
            Vertices2.push_back( i );
        Cell_handle c2 = f2.first;

        Plane3 P2 = GetOrientedPlane( c2, f2.second, ID );
        Point3 pTest2 = c2->vertex( Vertices2[0] )->point();

        Vector3 n1 = P1.orthogonal_vector();
        Vector3 n2 = P2.orthogonal_vector();
        Vector3 u1 = Pandora::Details::UnitVector( n1 );
        Vector3 u2 = Pandora::Details::UnitVector( n2 );

        FT u1Dotu2 = u1 * u2;
        u1Dotu2 = std::min( u1Dotu2, static_cast<FT>( 1 ) );  // clip at [0, 1]
        u1Dotu2 = std::max( u1Dotu2, static_cast<FT>( -1) );  // clip at [-1, 0]
        FT Angle  = acos( u1Dotu2  );
        if( !IsConvex( P1, edge, f2 ) )
          Angle = - Angle;
        
        return Angle;   // normalization by CGAL_PI not done yet
      }

      //-----------------------------------------------------------------------------
      // TripleLineLength
      //-----------------------------------------------------------------------------
      template< class EdgeIter, class HandleFn >
      FT TripleLineLength( EdgeIter pFirst, EdgeIter pEnd,
                           const Subdomain_index & ID,
                           HandleFn GetHandle ) const
      {
        if( pFirst == pEnd )
          return FT( 0 );
        
        FT Length = 0;
        for( ; pFirst != pEnd; ++ pFirst )
        {
          Edge e = GetHandle( pFirst );
          Cell_handle ch = e.get<0>();
          if( Details::IsTripleLineEdge( c3t3_, e) )
          {
            if( !c3t3_.triangulation().is_infinite( e ) )
            {
              Point3 p1 = ( ch->vertex( e.get<1>() ) )->point();
              Point3 p2 = ( ch->vertex( e.get<2>() ) )->point();
        
              FT EdgeLength = CGAL::sqrt( CGAL::squared_distance( p1, p2 ) );
              Length += EdgeLength;
            }
          }
        }
        return Length;
      }
      
      //-----------------------------------------------------------------------------
      //  Calculate Meanwidth of this grain
      //
      //  L( D ) = 1/ (2 * pi ) * sum_{ i = 1 } ^\nu ( \epsilon_i \beta_i )
      //-----------------------------------------------------------------------------
      template< class EdgeIter, class HandleFn >
      FT MeanWidth( EdgeIter pFirst, EdgeIter pEnd,
                    const Subdomain_index & ID,
                    HandleFn GetHandle ) const
      {
        if( pFirst == pEnd )
          return FT( 0 );

        FT fMeanWidth = 0;

        for( ; pFirst != pEnd; ++ pFirst )
        {
          Edge e = GetHandle( pFirst );
          
          if( ! c3t3_.triangulation().is_infinite( e ) )
          {
            Cell_handle ch = e.get<0>();
            Facet_circulator pEnd = c3t3_.triangulation().incident_facets( e );
            Facet_circulator pCur = c3t3_.triangulation().incident_facets( e );

            std::set< Facet >  FacetSet;
            if( pCur != 0 ) // construct of using Cell_circulator
            {
              do
              {
                if( ! c3t3_.triangulation().is_infinite( *pCur ) )
                {
                  Surface_index SurfID = c3t3_.surface_index( *pCur );
                  if( SurfID.first == ID || SurfID.second == ID )
                    FacetSet.insert( *pCur );
                }
                ++ pCur;
              } while( pCur != pEnd );
            }   

            if( FacetSet.size() == 2 )
            {
              typedef typename std::set<Facet>::iterator FIter;
              FIter p = FacetSet.begin();

              Facet f1 = *p; ++p;
              Facet f2 = *p;

              Vertex_handle vh1 =  ch->vertex( e.get<1>() );
              Vertex_handle vh2 =  ch->vertex( e.get<2>() ); 
            
              Point3 p1 = vh1->point();
              Point3 p2 = vh2->point();
              FT EdgeLength = CGAL::sqrt( CGAL::squared_distance( p1, p2 ) );
            
              FT TAngle = TurningAngle( e, f1, f2, ID );
            
              fMeanWidth += EdgeLength * TAngle; 
            
            }
            else if ( FacetSet.size() != 0  )  // non-manifold behavior
            {
              Point3 p1 = ( ch->vertex( e.get<1>() ) )->point();
              Point3 p2 = ( ch->vertex( e.get<2>() ) )->point();
              std::cout << ID << " " << FacetSet.size() << " " << p1 << " " << p2 << std::endl;
            }
          }
        }
        return fMeanWidth / ( 2 * CGAL_PI );
      }
  
      //-----------------------------------------------------------------------------
      //  Calculate volume
      //-----------------------------------------------------------------------------
      template< class CellHandleIter, class HandleFn >
      FT Volume ( CellHandleIter pFirst, CellHandleIter pEnd, HandleFn GetHandle  ) const
      {
        if( pFirst == pEnd )
          return 0;
        FT vTotal = 0;
        for( ; pFirst != pEnd; ++ pFirst )
        {
          Cell_handle ch = GetHandle( pFirst );
          Point3 p1 = ( ch->vertex(0) )->point();
          Point3 p2 = ( ch->vertex(1) )->point();
          Point3 p3 = ( ch->vertex(2) )->point();
          Point3 p4 = ( ch->vertex(3) )->point();
          Tetrahedron3 t ( p1, p2, p3, p4 );
          vTotal += t.volume();
        }
        return vTotal;
      }

      //-----------------------------------------------------------------------------
      //  CalculateGrainCenter
      //  Parameters:
      //    IDToPropMap: Maps from GrainID to Property.  This function does
      //                 not care about the exact implementation, but CenterMapper will
      //                 take care of how to assign to GrainPropMap
      //    AssignCenterFn:  Assign grain center to GrainPropMap (this is poor man's property map)
      //    IDToTetMap:  Map from GrainID to tetrahedron.  Tetrahedron is assumed be a model
      //                 of the CGAL Tets. (i.e., contains volume)  GrainTetMap is a model
      //                 of std::map
      //-----------------------------------------------------------------------------
      //template< class GrainPropMap, class CenterMapper, class GrainTetMap >
      template< class GrainPropMap, class CenterMapper >
      void CalculateGrainCenter( GrainPropMap &      IDToPropMap,
                                 CenterMapper        Put,
                                 const GrainTetMapT & IDToTetMap )
      {
        typedef typename GrainTetMapT::const_iterator MapIter;
        typedef typename std::set<Subdomain_index>::iterator SetIter;
        std::set<Subdomain_index> UniqueIDs;
        
        for( MapIter pCur = IDToTetMap.begin();
             pCur != IDToTetMap.end(); ++ pCur )
          UniqueIDs.insert( pCur->first );

        for( SetIter pCur = UniqueIDs.begin();
             pCur != UniqueIDs.end(); ++ pCur )
        {
          Subdomain_index ID = *pCur;
          std::pair< MapIter, MapIter > TetIters = IDToTetMap.equal_range( ID );      
          Point3 Center = CenterOfMass( TetIters.first, TetIters.second,
                                        IterToSecond< MapIter, Cell_handle >() );
          Put( IDToPropMap, ID, Center );
        }
      }
      
      //-----------------------------------------------------------------------------
      //  CalculateBoundingVolume
      //  Parameters:
      //    IDToPropMap: Maps from GrainID to Property.  This function does
      //                 not care about the exact implementation, but CenterMapper will
      //                 take care of how to assign to GrainPropMap
      //    AssignCenterFn:  Assign grain center to GrainPropMap (this is poor man's property map)
      //    IDToTetMap:  Map from GrainID to tetrahedron.  Tetrahedron is assumed be a model
      //                 of the CGAL Tets. (i.e., contains volume)  GrainTetMap is a model
      //                 of std::map
      //-----------------------------------------------------------------------------
      //  template< class GrainPropMap, class VolumeMapper, class GrainTetMap >
      template< class GrainPropMap, class VolumeMapper >
      void CalculateBoundingVolume( GrainPropMap &      IDToPropMap,
                                    VolumeMapper        Put,
                                    const GrainTetMapT & IDToTetMap )
      {
        typedef typename GrainTetMapT::const_iterator MapIter;
        typedef typename std::set<Subdomain_index>::iterator SetIter;
        std::set<Subdomain_index> UniqueIDs;

        for( MapIter pCur = IDToTetMap.begin();
             pCur != IDToTetMap.end(); ++ pCur )
          UniqueIDs.insert( pCur->first );
        
        for( SetIter pCur = UniqueIDs.begin();
             pCur != UniqueIDs.end(); ++ pCur )
        {
          Subdomain_index ID = *pCur;
          std::pair< MapIter, MapIter > TetIters = IDToTetMap.equal_range( ID );   
          BBox Box = BoundingVolume( TetIters.first, TetIters.second,
                                     IterToSecond< MapIter, Cell_handle >() );
          Put( IDToPropMap, ID, Box );
        }
      }
      
      //-----------------------------------------------------------------------------
      //  CalculateMeanwidth
      //  Parameters:
      //    IDToPropMap: Maps from GrainID to Property.  This function does
      //                 not care about the exact implementation, but CenterMapper will
      //                 take care of how to assign to GrainPropMap
      //    AssignCenterFn:  Assign grain center to GrainPropMap (this is poor man's property map)
      //    IDToTetMap:  Map from GrainID to tetrahedron.  Tetrahedron is assumed be a model
      //                 of the CGAL Tets. (i.e., contains volume)  GrainTetMap is a model
      //                 of std::map
      //-----------------------------------------------------------------------------
      template< class GrainPropMap, class VolumeMapper >
      void CalculateMeanwidth( GrainPropMap &      IDToPropMap,
                               VolumeMapper        Put,
                               const IDToEdgeMapT & IDToEdgeMap )
      {
        typedef typename IDToEdgeMapT::const_iterator MapIter;        
        typedef typename std::set<Subdomain_index>::iterator SetIter;
        std::set<Subdomain_index> UniqueIDs;

        for( MapIter pCur = IDToEdgeMap.begin();
             pCur != IDToEdgeMap.end(); ++ pCur )
          UniqueIDs.insert( pCur->first );

        for( SetIter pCur = UniqueIDs.begin();
             pCur != UniqueIDs.end(); ++ pCur )
        {

          Subdomain_index ID = *pCur;
          if( ID > 0 )
          {
            std::pair< MapIter, MapIter > EdgeIters = IDToEdgeMap.equal_range( ID );
          
            FT L = MeanWidth( EdgeIters.first, EdgeIters.second,
                              ID, IterToSecond< MapIter, Edge >() );
            Put( IDToPropMap, ID, L );
          }
        }
      }

      //-----------------------------------------------------------------------------
      //  CalculateTripleLineLength
      //
      //-----------------------------------------------------------------------------
      template< class GrainPropMap, class TripleLineMapper >
      void CalculateTripleLineLength( GrainPropMap &       IDToPropMap,
                                      TripleLineMapper     Put,
                                      const IDToEdgeMapT & IDToEdgeMap )
      {
        typedef typename IDToEdgeMapT::const_iterator MapIter;        
        typedef typename std::set<Subdomain_index>::iterator SetIter;
        std::set<Subdomain_index> UniqueIDs;

        for( MapIter pCur = IDToEdgeMap.begin();
             pCur != IDToEdgeMap.end(); ++ pCur )
          UniqueIDs.insert( pCur->first );

        for( SetIter pCur = UniqueIDs.begin();
             pCur != UniqueIDs.end(); ++ pCur )
        {

          Subdomain_index ID = *pCur;
          if( ID > 0 )
          {
            std::pair< MapIter, MapIter > EdgeIters = IDToEdgeMap.equal_range( ID );
            
            FT L = TripleLineLength( EdgeIters.first, EdgeIters.second,
                                     ID, IterToSecond< MapIter, Edge >() );
            Put( IDToPropMap, ID, L );
          }
        }
      }

      //-----------------------------------------------------------------------------
      //  CalculateVolume
      //-----------------------------------------------------------------------------
      //      template< class GrainPropMap, class VolumeMapper, class GrainTetMap >
      template< class GrainPropMap, class VolumeMapper >
      void CalculateVolume( GrainPropMap &      IDToPropMap,
                            VolumeMapper        Put,
                            const GrainTetMapT & IDToTetMap )
      {
        typedef typename GrainTetMapT::const_iterator MapIter;
        
        std::set<Subdomain_index> UniqueIDs;
        for( MapIter pCur = IDToTetMap.begin();
             pCur != IDToTetMap.end(); ++ pCur )
          UniqueIDs.insert( pCur->first );
        typedef typename std::set<Subdomain_index>::iterator SetIter;
        for( SetIter pCur = UniqueIDs.begin();
             pCur != UniqueIDs.end(); ++ pCur )
        {
          Subdomain_index ID = *pCur;
          std::pair< MapIter, MapIter > TetIters = IDToTetMap.equal_range( ID );      
          FT V = Volume( TetIters.first, TetIters.second,
                         IterToSecond< MapIter, Cell_handle >() );
          Put( IDToPropMap, ID, V );
        }
      }

      //-----------------------------------------------------------------------------
      //  SelectInterior
      //-----------------------------------------------------------------------------
      template< class GrainPropMap, class InteriorMapper>
      void SelectInterior( GrainPropMap       & IDToPropMap,
                           InteriorMapper       Put,
                           const IDNgbMapT    & IDNgbMap,
                           const GrainTetMapT & IDToTetMap )
      {
        typedef typename GrainTetMapT::const_iterator MapIter;
        
        std::set<Subdomain_index> UniqueIDs;
        for( MapIter pCur = IDToTetMap.begin();
             pCur != IDToTetMap.end(); ++ pCur )
          UniqueIDs.insert( pCur->first );
        typedef typename std::set<Subdomain_index>::iterator SetIter;
        for( SetIter pCur = UniqueIDs.begin();
             pCur != UniqueIDs.end(); ++ pCur )
        {
          Subdomain_index ID = *pCur;

          // find exterior
          bool bInterior = true;
          IDNgbMapT::const_iterator pIDSet = IDNgbMap.find( ID );
          if( pIDSet != IDNgbMap.end() )
          {
            for( std::set<int>::const_iterator pIter = pIDSet->second.begin();
                 pIter != pIDSet->second.end(); ++ pIter  )
            {
              if( *pIter <= 0 )
                bInterior = false;
            }
            Put( IDToPropMap, ID, bInterior );
          }
          else
          {
            Put( IDToPropMap, ID, false );
          }
        }
      }

      
    };
  }
}

#endif




     //    if( ID == 2 )
//         {
//           if( ( u1 * u2 ) < 0 )
//           {
//             std::cout << " First " << " ";
//             for( int i = 0; i < 4; i ++ )
//             {
//               std::cout << c1->vertex( i )->point().x() << " ";
//               std::cout << c1->vertex( i )->point().y() << " ";
//               std::cout << c1->vertex( i )->point().z() << " ";
//             }
//             std::cout << " [ " << f1.second << " ] " << u1 << std::endl;
            
//             std::cout << " Second " << " ";
//             for( int i = 0; i < 4; i ++ )
//             {
//               std::cout << c2->vertex( i )->point().x() << " ";
//               std::cout << c2->vertex( i )->point().y() << " ";
//               std::cout << c2->vertex( i )->point().z() << " ";
//             }
//             std::cout << " [ " << f2.second  << " ] " << u2 << std::endl;         

//             DEBUG_FLAG = true;
//             GetPlane( c1->vertex( Vertices1[0] )->point(),
//                       c1->vertex( Vertices1[1] )->point(),
//                       c1->vertex( Vertices1[2] )->point(),
//                       c1->vertex( f1.second ),
//                       ID );
//             GetPlane( c2->vertex( Vertices2[0] )->point(),
//                       c2->vertex( Vertices2[1] )->point(),
//                       c2->vertex( Vertices2[2] )->point(),
//                       c2->vertex( f2.second ),
//                       ID );
//             std::cout << "Done outputting " << std::endl;
//             std::cin.get();
//           }
//         }

     //    if( ID == 2 )
//           Normal_os << pTest1 << " "
//                     << pTest2 << " "
//                     << u1 << " " << u2 << " "<<  Angle << std::endl;




   //    //---------------------------------
//       // GetPlane
//       //
//       //  By definition, vh is going to be a handle that
//       //  points to an interior vertex.  (if it isn't interior,
//       //  then all vertices are neighboring multiple ID, which
//       //  means this is a single cell grain.)
//       //
//       //  Return a plane that's oriented outward away from the
//       //  ID.  
//       //---------------------------------
//       Plane3 GetPlane( const Point3 & p1,
//                        const Point3 & p2,
//                        const Point3 & p3,
//                        const Vertex_handle vh,
//                        const Subdomain_index & ID ) const
//       {
//         typedef typename CGAL::Point_3< GeomKernel >    SPoint3;

//         SPoint3 sp1( p1.x(), p1.y(), p1.z() );
//         SPoint3 sp2( p2.x(), p2.y(), p2.z() );
//         SPoint3 sp3( p3.x(), p3.y(), p3.z() );
//         SPoint3 spOpp( vh->point().x(), vh->point().y(), vh->point().z() );
//         Plane3 P( sp1, sp2, sp3 );

//         int Dimension = c3t3_.in_dimension( vh );

//         if( P.oriented_side( spOpp ) == CGAL::ON_ORIENTED_BOUNDARY )
//         {
//           std::cout << "Limits of numerical type " << std::endl;
//           exit(0);
//         }
        
//         typename std::vector< Cell_handle > Cells;
//         c3t3_.triangulation().finite_incident_cells(vh, std::back_inserter(Cells));
//         if( Cells.begin() != Cells.end() )    // since this vertex is interior, any ID will work
//         {            
//           if( c3t3_.subdomain_index( *Cells.begin() )== ID )          // if the vertex is in the ID region
//           {
//             if( P.oriented_side( spOpp ) != CGAL::ON_NEGATIVE_SIDE ) // interior should be at the back of facet
//             {
//               P = P.opposite();
//               if( DEBUG_FLAG ) std::cout << " In ID, spOpp on negative ";
//             }
//             else
//               if( DEBUG_FLAG ) std::cout << " In ID, No Flip ";
//           }
//           else
//           {
//             if( P.oriented_side( spOpp ) == CGAL::ON_NEGATIVE_SIDE ) // Exterior should be on positive side 
//             {
//               P = P.opposite();
//               if( DEBUG_FLAG ) std::cout << " In Ext, spOpp on negative ";
//               else
//                 if( DEBUG_FLAG ) std::cout << " In Ext, No Flip ";

//             }
//           }
//         }
//         else
//         {
//           if( DEBUG_FLAG ) std::cout << " No finite incident ";
//         }
        
//         if( DEBUG_FLAG ) std::cout << std::endl;
          
//         return P;
//       }
