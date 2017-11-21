//----------------------------------------------------
//  DynamicsAnalysis.h
//
//  Author:   Frankie Li
//  Purpose:  Estimate triple line prop
//
//----------------------------------------------------

#ifndef DYNAMICS_ANALYSIS_H
#define DYNAMICS_ANALYSIS_H

#include <CGAL/IO/File_medit.h>
#include "SimpleMeshVTKWriter.h"
#include <map>
#include <vector>
#include "GrainGeometryAnalysis.h"
#include "Quaternion.h"
#include "3dMath.h"
#include "boost/multi_array.hpp"
//----------------------------------------------------
//  For grain tracking
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/algorithm.h>
//----------------------------------------------------


namespace Pandora
{
  namespace Details
  {
    
  }

  namespace Stuff
  {
    int i;
  }
  namespace GrainDynamics
  {
    struct GrainHistoryProp
    {
      int ID;
      float fMisorient;
      GeneralLib::SQuaternion qMis;
      float IntersectingVolume;
      float Volume_a;
      float Volume_b;
      int   nNgb;
      float MeanWidth;
      float TripleLine;
      bool  bInternal;
      GeneralLib::SQuaternion q_a;
      GeneralLib::SQuaternion q_b;
      GeneralLib::SVector3    GrainCenter_a;
    };
    
    struct TrivialMapInsert
    {
      template< class T1, class T2, class T3 >
      void operator()( T1 & map, const T2 & key, const T3 & value ) const
      {  map[key] = value; }
    };
    

    //--------------------------------------------------------------
    // GetBestMatch
    //--------------------------------------------------------------
    template< class Map, class IDIter, class CostFnT  >
    std::pair<int, float>
    GetBestMatch( const Map & ID2OrientMap, IDIter pFirst, IDIter pEnd,
                  const SQuaternion qRef, CostFnT CostFn )
    {
      if( pFirst == pEnd )
        return std::make_pair( -1, 1000 );
      int nBestID = -1;
      float fBestCost = 1000;
      for( ; pFirst != pEnd; ++ pFirst )
      {
        if( ID2OrientMap.find( *pFirst ) == ID2OrientMap.end() )
        {
          std::cout << "ID2OrientMap is incomplete" << std::endl;
          exit(0);
        }
        float fCost = CostFn( ID2OrientMap.find( *pFirst )->second, qRef );
        if( fBestCost > fCost )
        {
          fBestCost = fCost;
          nBestID   = *pFirst;
        }
      }
      return std::make_pair( nBestID, fBestCost );
    }
    
    //--------------------------------------------------------------
    // GetEnclosingGrain
    //--------------------------------------------------------------
    //------------------------------------------------------------
    template< class C3T3, class ID2OrientMapT >
    class GrainDynamicsAnalysis
    {
    public:
      typedef typename C3T3::Cell_handle              Cell_handle;
      typedef typename std::map< int, std::set<int> > IDNgbMap;
      
      typedef typename Details::GrainGeometry<C3T3> GrainGeometry;
      typedef typename GrainGeometry::BBox          BBox;
      typedef typename GrainGeometry::FT            FT;
      typedef typename GrainGeometry::Point3        Point3;
      typedef typename GrainGeometry::Tetrahedron3  Tetrahedron3;
      typedef typename std::map< int, Point3 >      ID2CenterMapT;
      typedef typename std::map< int, BBox>         ID2BBoxMapT;
      typedef typename std::map< int, FT>           ID2VolumeMapT;
      typedef typename GrainGeometry::GrainTetMapT  GrainTetMap;
      
      typedef typename C3T3::Subdomain_index        Subdomain_index;
      typedef typename ID2CenterMapT::iterator      GrainCenterIter;
      typedef typename IDNgbMap::iterator           IDNgbMapIter;
      typedef typename ID2OrientMapT::iterator      OrientIter;
      
    private:

      bool GrainDynamicsDebugFlag;

      //--------------------------------------------------------------
      //--------------------------------------------------------------
      template< class IDIterator >
      std::vector<SQuaternion> GetNgbMisorinetation( IDIterator pCur,
                                                     IDIterator pEnd,
                                                     const SQuaternion & qRef,
                                                     const ID2OrientMapT Map )
      {
        std::vector<SQuaternion> MisorientList;
        for( ; pCur != pEnd; ++ pCur )
        {
          OrientIter qIter = Map.find( *pCur );
          if ( qIter == Map.end() )
            std::cerr << "Ngb ID does not lead to any orientation " << std::endl;
          MisorientList.push_back( 
                                  LatticeSymmetry
                                  ::GetMisorientationOp( LatticeSymmetry
                                                         ::CCubicSymmetry
                                                         ::Get(),
                                                         qRef,
                                                         qIter->second )
                                  );
        }
        return MisorientList;
      }
                                                     
      
      //--------------------------------------------------------------
      //  IsInternalGrain  -- could use some improvement 
      //--------------------------------------------------------------
      bool IsInternalGrain( Subdomain_index ID_a,
                            const GrainTetMap & TetMap_a,
                            const C3T3 & c3t3_ ) const
      {
        typedef typename GrainTetMap::const_iterator TetIter;
        std::pair< TetIter, TetIter > TetRange = TetMap_a.equal_range( ID_a );

        if( TetRange.first == TetRange.second )
          return false;
        
        typedef typename std::map< Subdomain_index, float > VolumeMap;
        VolumeMap IDVolumeMap;
        
        for( TetIter pCur = TetRange.first;
             pCur != TetRange.second; ++ pCur )
        {
          Cell_handle ch = pCur->second;
          if( c3t3_.subdomain_index( ch ) <= 0 )
            return false;
          for ( int i = 0; i < 4; i ++ )
          {
            int nNgbID = c3t3_.subdomain_index( ch->neighbor( i ));
            if( nNgbID <= 0 )
              return false;
          }
        }
        return true;
      }

      //--------------------------------------------------------------
      //  TopologicalDifference
      //    Purpose:  Given grain IDa from c3t3_a and IDb from c3t3_b,
      //              c3t3_a, c3t3_b are aligned, and IDa, and IDb persumed
      //              to refer to two grains that are supposedly the cloest
      //              match between c3t3_a and c3t3_b, return the topological
      //              difference in terms of number of neighbors and the
      //              difference in boundaries
      //  (Is this the same as graph isomotery?  I hope not.)
      //  
      //  Parameters:  EdgeAcceptThreshold - Angle in radian where an edge will 
      //               be accepted, i.e., d( q_a, q_b ) < EdgeAcceptThreshold )
      //               This parameter should be dictated by the grain ID threshold.
      //
      //--------------------------------------------------------------
      template< class Subdomain_index >
      int TopologicalDifference( Subdomain_index ID_a,
                                 Subdomain_index ID_b,
                                 double EdgeAcceptThreshold )
      {
        
        IDNgbMapIter pNgb_a = NgbMap_a.find( ID_a );
        IDNgbMapIter pNgb_b = NgbMap_b.find( ID_b );

        if( pNgb_a == NgbMap_a.end() || pNgb_b == NgbMap_b.end() )
          return -1;
        OrientIter qPtr_a = ID2OrientMap_a.find( ID_a );
        OrientIter qPtr_b = ID2OrientMap_b.find( ID_b );
        
        assert( qPtr_a != ID2OrientMap_a.end()
                && qPtr_b != ID2OrientMap_b.end() );
        
        typedef typename std::set<int>::iterator IDIter;
        std::vector<SQuaternion> MisorientList_a, MisorientList_b;

        MisorientList_a = GetNgbMisorientation( pNgb_a->second.begin(),
                                                pNgb_a->second.end(),
                                                qPtr_a->second,
                                                ID2OrientMap_a );

        MisorientList_b = GetNgbMisorientation( pNgb_b->second.begin(),
                                                pNgb_b->second.end(),
                                                qPtr_b->second,
                                                ID2OrientMap_b );

        // build matrix, but save it?  Look for non-uniqueness?
        boost::multi_array<double, 2> CostMatrix;
        CostMatrix.resize( boost::extents[ MisorientList_a.size() ][ MisorientList_b.size()  ] );
        for( int i = 0; i < MisorientList_a.size(); i ++ )
        {
          for( int j = 0; j < MisorientList_b.size(); j ++ )
          {
            CostMatrix[i][j] =
              LatticeSymmetry
              ::GetMisorientation( LatticeSymmetry
                                   ::CCubicSymmetry
                                   ::Get(),
                                   MisorientList_a[i],
                                   MisorientList_b[j] );
          }
        }
        
        return std::abs( int(pNgb_a->second.size()) - int(pNgb_b->second.size()) );
      }
      
      //--------------------------------------------------------------
      //  FindBestByTopologyMatch
      //--------------------------------------------------------------
      
      //--------------------------------------------------------------
      //  FindBestBoundaryMatch  -- need kdtree
      //--------------------------------------------------------------

      
      //--------------------------------------------------------------
      //  Calculate the intersecting volume of the two grains
      //--------------------------------------------------------------
      std::pair< Subdomain_index, double >
      FindBestVolumeOverlap( Subdomain_index ID_a,
                             const GrainTetMap & TetMap_a,
                             const C3T3 & c3t3_ ) const
      {
        typedef typename GrainTetMap::const_iterator TetIter;
        std::pair< TetIter, TetIter > TetRange = TetMap_a.equal_range( ID_a );

        if( TetRange.first == TetRange.second )
          return std::make_pair( -1, 0 );
        
        typedef typename std::map< Subdomain_index, float > VolumeMap;
        VolumeMap IDVolumeMap;
        
        for( TetIter pCur = TetRange.first;
             pCur != TetRange.second; ++ pCur )
        {
          Cell_handle ch = pCur->second;
          if( ! c3t3_.triangulation().is_infinite( ch ) ) // this is typically edge
          {
            Point3 p1 = ch->vertex(0)->point();
            Point3 p2 = ch->vertex(1)->point();
            Point3 p3 = ch->vertex(2)->point();
            Point3 p4 = ch->vertex(3)->point();
            //Tetrahedron3 t ( p1, p2, p3, p4 );
            Tetrahedron3 t = c3t3_.triangulation().tetrahedron( ch );
            FT v = t.volume();
            FT x = p1.x() + p2.x() + p3.x() + p4.x();
            FT y = p1.y() + p2.y() + p3.y() + p4.y();
            FT z = p1.z() + p2.z() + p3.z() + p4.z();
            x /= 4.;
            y /= 4.;
            z /= 4.;
            Cell_handle ch_b = c3t3_.triangulation().locate(  Point3( x, y, z ) );
            if( ! c3t3_.triangulation().is_infinite( ch_b ) ) // this is typically edge
            {
              Subdomain_index GrainID_b =  c3t3_.subdomain_index( ch_b );
              if( IDVolumeMap.find( GrainID_b ) == IDVolumeMap.end() )
                IDVolumeMap[ GrainID_b ]  = v; 
              else
                IDVolumeMap[ GrainID_b ] += v; 
            }
          }  // ignore infinite cells
        }
        typedef typename VolumeMap::iterator VIter;

        double fBestVolume = 0;
        Subdomain_index BestGrainID = -1;
        for( VIter pCur = IDVolumeMap.begin();
             pCur != IDVolumeMap.end();
             ++ pCur )
        {
          if( pCur->second > fBestVolume )
          {
            fBestVolume = pCur->second; 
            BestGrainID = pCur->first;
          }
        }
        return std::make_pair( BestGrainID, fBestVolume );
      }
      
      //--------------------------------------------------------------
      //--------------------------------------------------------------
      template <class BBox3 >
      double IntersectingVolume( const BBox3 & Box_a,
                                 const BBox3 & Box_b ) const
      {
        BBox3 Intersection( std::max( Box_a.xmin(), Box_b.xmin() ),
                            std::max( Box_a.ymin(), Box_b.ymin() ),
                            std::max( Box_a.zmin(), Box_b.zmin() ),
                            std::min( Box_a.xmax(), Box_b.xmax() ),
                            std::min( Box_a.ymax(), Box_b.ymax() ),
                            std::min( Box_a.zmax(), Box_b.zmax() )     );

        double dx = Intersection.xmax() - Intersection.xmin();
        double dy = Intersection.ymax() - Intersection.ymin();
        double dz = Intersection.zmax() - Intersection.zmin();

        return dx * dy * dz;
      }

      
      //--------------------------------------------------------------
      // GetBestMatch
      //--------------------------------------------------------------
      template< class IDIter, class CostFnT  >
      std::pair<int, float>
      GetBestMatch( const ID2OrientMapT & ID2OrientMap,
                    const ID2BBoxMapT & ID2BBoxMap,
                    IDIter pFirst, IDIter pEnd,
                    const SQuaternion qRef,
                    const BBox & BBox,
                    CostFnT CostFn ) const
      {
        if( pFirst == pEnd )
          return std::make_pair( -1, 1000 );
        
        int nBestID = -1;
        float fBestCost = 1000;
        float fMaxVolume = 0;
        for( ; pFirst != pEnd; ++ pFirst )
        {
          if( ID2BBoxMap.find( *pFirst ) != ID2BBoxMap.end() )
          {
            if( CGAL::do_overlap( ID2BBoxMap.find( *pFirst )->second, BBox ) )
            {
              double fVolume = IntersectingVolume( BBox, ID2BBoxMap.find( *pFirst )->second );
              if( fVolume > fMaxVolume )
              {
                fBestCost = CostFn( ID2OrientMap.find( *pFirst )->second, qRef );
                fMaxVolume = fVolume;
                nBestID = *pFirst;
              }

              float fCost = CostFn( ID2OrientMap.find( *pFirst )->second, qRef );
            
              if( GrainDynamicsDebugFlag )
              {
                std::cout << *pFirst << " " << fCost << " " << fVolume << std::endl;
              }

            }
          }
        }
        return std::make_pair( nBestID, fBestCost );
      }
      
      //--------------------------------------------------------------
      // GetEnclosingGrain
      //--------------------------------------------------------------
      bool GetEnclosingGrain( Subdomain_index & GrainID,
                              const C3T3 & c3t3_, Cell_handle ch ) const
      {
        if( c3t3_.triangulation().is_infinite( ch ) ) // this is typically edge
        {
          // should deal with neighbors here, ignoring it now.
          return false;
        }
        GrainID =  c3t3_.subdomain_index( ch );
        return true;
      }

      const C3T3 & c3t3_a;
      const C3T3 & c3t3_b;
      const ID2OrientMapT &ID2OrientMap_a;
      const ID2OrientMapT &ID2OrientMap_b; 

      ID2CenterMapT ID2CenterMap_a, ID2CenterMap_b;
      ID2BBoxMapT   ID2BBoxMap_a, ID2BBoxMap_b;
      ID2VolumeMapT ID2VolumeMap_a, ID2VolumeMap_b;
      
      GrainGeometry GeomExtractor_a;
      GrainGeometry GeomExtractor_b;
      IDNgbMap NgbMap_a, NgbMap_b;
      GrainTetMap GrainTetMap_a, GrainTetMap_b;
      
    public:

      GrainDynamicsAnalysis( const C3T3 & c3t3_a_,
                             const C3T3 & c3t3_b_,
                             const ID2OrientMapT &ID2OrientMap_a_,
                             const ID2OrientMapT &ID2OrientMap_b_ ) :
        c3t3_a( c3t3_a_ ), c3t3_b( c3t3_b_ ),
        ID2OrientMap_a( ID2OrientMap_a_ ),
        ID2OrientMap_b( ID2OrientMap_b_ ),
        GeomExtractor_a( c3t3_a_ ),
        GeomExtractor_b( c3t3_b_ ) {}
      //--------------------------------------------------------------
      //  MatchGrains
      //
      //  Purpose:  - Associate grains that are notcritical.  In
      //              in another words, look for grains that exists
      //              in both volumes.
      //  Direction of time:  a -> b
      // for each grain in a
      //     assert( GrainIDHistoryList.find( ID_a ) == GrainIDHistoryList.end() );
      //     locate center in c3t3_b
      //     get grain ID using CellToID_b
      //     get q_b, orientation
      //     if d( q_b, q_a ) < threshold )
      //        GrainIDHistoryList[ ID_a ].push_back( ID_b )
      //     else
      //        for each ID_ngb_b <- NgbMap_b[ ID_b ]
      //           if( d( q_ID_ngb_b, q_a ) < threshold )
      //               // perhaps check overlap?
      //             GrainIDHistoryList[ ID_a ].push_back( ID_ngb_b )
      //    //---- may add this in the future -----------------
      //         if ( NotFound )
      //         for each cell in Grain_a
      //             locate cell in c3t3_b
      //             if d(q_b, q_a) < threshold )
      //               ID_b -> Candidate_set
      //--------------------------------------------------------------
      template< class CostFnT  > 
      void GetGrainHistory( std::map< int, GrainHistoryProp >  & GrainIDHistoryList,
                            CostFnT CostFn, Float fThresh )
      {
        Details::BuildFiniteGrainNgbMap( NgbMap_a, c3t3_a );
        Details::BuildFiniteGrainNgbMap( NgbMap_b, c3t3_b );
      
        Details::BuildGrainTetMap( GrainTetMap_a, c3t3_a );
        Details::BuildGrainTetMap( GrainTetMap_b, c3t3_b );

        GeomExtractor_a.CalculateGrainCenter( ID2CenterMap_a, TrivialMapInsert(), GrainTetMap_a );
        GeomExtractor_b.CalculateGrainCenter( ID2CenterMap_b, TrivialMapInsert(), GrainTetMap_b );

        GeomExtractor_a.CalculateBoundingVolume( ID2BBoxMap_a, TrivialMapInsert(), GrainTetMap_a );
        GeomExtractor_b.CalculateBoundingVolume( ID2BBoxMap_b, TrivialMapInsert(), GrainTetMap_b );
        
        GeomExtractor_a.CalculateVolume( ID2VolumeMap_a, TrivialMapInsert(), GrainTetMap_a );
        GeomExtractor_b.CalculateVolume( ID2VolumeMap_b, TrivialMapInsert(), GrainTetMap_b );

        //-------------------------------------------------------------
        // Sneak peak debug
        typedef typename GrainGeometry::IDToEdgeMapT  IDToEdgeMapT;
        typedef typename std::map< Subdomain_index, double > IDToFloatMap;
        IDToFloatMap ID2MeanWidth_a;
        IDToFloatMap ID2TripLineSum_a;
        IDToEdgeMapT IDToEdgeMap;
        Pandora::Details::BuildGrainSurfaceEdgeMap( IDToEdgeMap, c3t3_a );
        std::cout << "Calculating mean width  " << std::endl;
        GeomExtractor_a.CalculateMeanwidth( ID2MeanWidth_a, TrivialMapInsert(), IDToEdgeMap );

        GeomExtractor_a.CalculateTripleLineLength( ID2TripLineSum_a, TrivialMapInsert(), 
                                                   IDToEdgeMap );
        std::cout << "Building neighbor map  " << std::endl;
        
        //-------------------------------------------------------------
        
        for( typename ID2CenterMapT::iterator pCur = ID2CenterMap_a.begin();
             pCur != ID2CenterMap_a.end(); ++ pCur )
        {
          Subdomain_index ID_a = pCur->first;
          if( ID_a > 0 )   // if not empty space
          {
            assert( GrainIDHistoryList.find( ID_a ) == GrainIDHistoryList.end() );
            Point3 GrainCenter_a = pCur->second;
            Cell_handle ch = c3t3_b.triangulation().locate( GrainCenter_a );
            Subdomain_index ID_b;
            bool bFinite = GetEnclosingGrain( ID_b, c3t3_b, ch );
            if( bFinite )
            {
              GrainDynamicsDebugFlag = false;

              IDNgbMapIter ID_b_Ngbs = NgbMap_b.find( ID_b );

              if( ID_b_Ngbs != NgbMap_b.end() )
              {
                SQuaternion q_a = ID2OrientMap_a.find( ID_a )->second;
                //               std::pair<int, float> Result
                //                 = GetBestMatch( ID2OrientMap_b, ID2BBoxMap_b,
                //                                 ID_b_Ngbs->second.begin(),
                //                                 ID_b_Ngbs->second.end(),
                //                                 q_a,
                //                                 ID2BBoxMap_a.find( ID_a )->second,
                //                                 CostFn );
                std::pair< Subdomain_index, double >  Result
                  = FindBestVolumeOverlap( ID_a, GrainTetMap_a, c3t3_b );

                Subdomain_index BestID = Result.first;
                Subdomain_index BestVolume = Result.second;
                
                if( BestID > 0 )
                {
                  float fCost = CostFn( ID2OrientMap_b.find( BestID )->second,
                                        ID2OrientMap_a.find( ID_a )->second   );
                  if( fCost< fThresh )
                  {
                    GrainIDHistoryList[ ID_a ].ID = BestID;
                    GrainIDHistoryList[ ID_a ].fMisorient = fCost;
                    GrainIDHistoryList[ ID_a ].IntersectingVolume = BestVolume;
                    GrainIDHistoryList[ ID_a ].Volume_a = ID2VolumeMap_a[ ID_a ];
                    GrainIDHistoryList[ ID_a ].Volume_b = ID2VolumeMap_b[ BestID ];
                    GrainIDHistoryList[ ID_a ].nNgb     = NgbMap_a[ ID_a ].size();
                    GrainIDHistoryList[ ID_a ].MeanWidth  = ID2MeanWidth_a[ ID_a ];
                    GrainIDHistoryList[ ID_a ].TripleLine = ID2TripLineSum_a[ ID_a ];
                    GrainIDHistoryList[ ID_a ].q_a        = ID2OrientMap_a.find( ID_a   )->second;
                    GrainIDHistoryList[ ID_a ].q_b        = ID2OrientMap_b.find( BestID )->second;
                    GrainIDHistoryList[ ID_a ].GrainCenter_a.m_fX = GrainCenter_a.x();
                    GrainIDHistoryList[ ID_a ].GrainCenter_a.m_fY = GrainCenter_a.y();
                    GrainIDHistoryList[ ID_a ].GrainCenter_a.m_fZ = GrainCenter_a.z();

                    GrainIDHistoryList[ ID_a ].bInternal  = IsInternalGrain( ID_a, GrainTetMap_a, c3t3_a );
                  }
                }
              }
            }  // else we'll just ignore for now - need to fix this
          }
        }
      }
      
    };  // end class  GrainDynamicsAnalysis
  }  // end GrainDynamics
}


#endif


//   if( ID_a == 2140 )
//           {
//             typedef typename std::set<int>::iterator NgbSetIter;
//             NgbSetIter pIter = ID_b_Ngbs->second.begin();
//             for( ; pIter != ID_b_Ngbs->second.end(); ++ pIter )
//             {
//               std::cout << ID2OrientMap_a.find( ID_a )->second << std::endl;
//               std::cout << ID2OrientMap_b.find( *pIter )->second << std::endl;
//               std::cout << *pIter << " " 
//                         << RADIAN_TO_DEGREE( CostFn( ID2OrientMap_a.find( ID_a )->second,
//                                                      ID2OrientMap_b.find(*pIter)->second ) ) << std::endl; 
//             }
//           }
     //  std::ofstream DebugIDStream( "NgbDebugIDa.txt" );
//       typedef IDNgbMap::iterator IDIter;
//       for( IDIter pCur = NgbMap_a.begin(); pCur != NgbMap_a.end();
//            ++ pCur )
//       {
//         DebugIDStream << pCur->first << " ";
//         typedef std::set<int>::iterator SetIter;
//         for( SetIter it = pCur->second.begin(); it != pCur->second.end(); ++ it )
//         {
//           DebugIDStream << *it << " ";
//         }
//         DebugIDStream << std::endl;
//       }
//       //  DEBUG
//       exit(0);
//       //  DEBUG 
