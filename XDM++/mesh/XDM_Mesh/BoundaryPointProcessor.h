//-------------------------------------------------------------------
//
//  BoundaryPointProcessor.h
//
//
//
//-------------------------------------------------------------------

#ifndef XDM_BOUNDARY_POINT_PROCESSOR_H
#define XDM_BOUNDARY_POINT_PROCESSOR_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/Memory_sizer.h>

#include <fstream>
#include <string>
#include <vector>
#include <ostream>
#include "XDM_Utilities.h"

namespace XDM
{
  template< class Kernel > 
  class BoundaryPointProcessor
  {
    
  public:
    // types for K nearest neighbors search structure
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;
    
    typedef typename Neighbor_search::Tree Tree;
    typedef typename Kernel::Point_3 Point;
    typedef typename std::vector<Point> PointList;
    
    typedef typename PointList::iterator PointIter;

    typedef CGAL::XDM_test::TestPoint WeightedPoint;
    typedef std::vector<WeightedPoint> WeightedPointList;
    typedef typename std::vector< std::pair< WeightedPoint, int > > WPointDimList;
    typedef typename std::vector< std::pair< Point, int > > PointDimList;


    struct PointDimCmp
    {
      bool operator() ( const std::pair< Point, int> & oLHS,
                        const std::pair< Point, int> & oRHS ) const
      {
        return oLHS.second > oRHS.second;
      }
    };

    struct NotQuadPoint
    {
      bool operator() ( const std::pair< Point, int> & oObj ) const
      {
        return oObj.second < 4;
      }
    };
    
  private:
    Tree oBndPointTree;
    WeightedPointList oNewBndPoints;
    
  public:
    
    BoundaryPointProcessor( PointIter pFirst, PointIter pEnd )
      : oBndPointTree( pFirst, pEnd )
    {  }

    BoundaryPointProcessor(){}

    //---------------------------------------------
    //  Initialize
    //---------------------------------------------
    void Initialize( WPointDimList & oWpList,
                     float fSmoothingLength,
                     float fPointSpacing,
                     float fMinEdgeLength )
    {
      std::vector< Point > oBndPointList;
      for( int i = 0; i < oWpList.size(); i ++ )
      {
        const WeightedPoint & p = oWpList[i].first;
        oBndPointList.push_back( Point( p.x, p.y, p.z ) );
      }
      oBndPointTree = Tree( oBndPointList.begin(),
                            oBndPointList.end() );
      for( int i = 0; i < oWpList.size(); i ++ )
        oNewBndPoints.push_back( oWpList[i].first );

      //      PointDimList oTestSmoothedPoints = SmoothPoints( oWpList, 30, 40 ;)
      PointDimList oTestSmoothedPoints = SmoothPoints( oWpList,
                                                       fSmoothingLength,
                                                       fPointSpacing );
      //      oNewBndPoints = CalculateWeight( oTestSmoothedPoints, 10 );
      oNewBndPoints = CalculateWeight( oTestSmoothedPoints, fMinEdgeLength );
    }
    
    //---------------------------------------------
    //  SmoothPoints
    //
    //  A dumb way to smooth points
    //
    //---------------------------------------------
    PointDimList SmoothPoints( const WPointDimList & oWpList,
                               float fSmoothingLength,
                               float fPointSpacing)  
    {
      PointDimList oNewBndPtrList;
      for( int i = 0; i < oWpList.size(); i ++ )
      {
        // only smooth non-quad points
        if( oWpList[i].second < 4 )
        {
          WeightedPoint p = oWpList[i].first;
          Point Query(p.x, p.y, p.z);
          Neighbor_search search( oBndPointTree, Query, 20 );
          int nNgbCounted = 1;
          WeightedPoint pNewPoint = oWpList[i].first;

          float TotalWeight = 1;
          for( Search_iterator it = search.begin();
               it != search.end(); ++ it )
          {
            if( it->second < fSmoothingLength )
            {
              Point pNew = it->first;
              float fDx2 = ( pNew.x() - p.x ) * ( pNew.x() - p.x );
              float fDy2 = ( pNew.y() - p.y ) * ( pNew.y() - p.y );
              float fDz2 = ( pNew.z() - p.z ) * ( pNew.z() - p.z );
              float r2 = ( fDx2 + fDy2 + fDz2 );

              if( r2 > 1e-4 )  // larger than some epsilon
              {
                float ri = 1./ std::sqrt( r2 );
                pNewPoint.x += ( ri * pNew.x() );
                pNewPoint.y += ( ri * pNew.y() );
                pNewPoint.z += ( ri * pNew.z() );
                nNgbCounted ++;
                TotalWeight += ri;
              }
              else
              {
                pNewPoint.x += ( pNew.x() );
                pNewPoint.y += ( pNew.y() );
                pNewPoint.z += ( pNew.z() );
                TotalWeight += 1;
              }
            }
          }
//           pNewPoint.x /= float( nNgbCounted );
//           pNewPoint.y /= float( nNgbCounted );
//           pNewPoint.z /= float( nNgbCounted );
          pNewPoint.x /= TotalWeight;
          pNewPoint.y /= TotalWeight;
          pNewPoint.z /= TotalWeight;
          
          Point ToInsert( pNewPoint.x, pNewPoint.y, pNewPoint.z );
          oNewBndPtrList.push_back( std::make_pair( ToInsert, oWpList[i].second ) );
        }
        else
        {
          Point p =  Point( oWpList[i].first.x, oWpList[i].first.y, oWpList[i].first.z );
          oNewBndPtrList.push_back( std::make_pair( p, oWpList[i].second ) );
        }
        // else do nothing
      }

      oNewBndPtrList.erase( CGAL::XDM_PointProcessing::
                            grid_simplify_point_set( oNewBndPtrList.begin(),
                                                     oNewBndPtrList.end(),
                                                     CGAL::First_of_pair_property_map< std::pair<Point, int> >(),
                                                     CGAL::Second_of_pair_property_map< std::pair<Point, int> >(),
                                                     fPointSpacing ),
                           oNewBndPtrList.end() );

      PointDimList (oNewBndPtrList).swap( oNewBndPtrList ); // swap trick
      
      return oNewBndPtrList;
    }
    
    //---------------------------------------------
    //  Calculate weight of each point
    //---------------------------------------------
    void CalculateWeight( WPointDimList & oWpList )
    {
      for( int i = 0; i < oWpList.size(); i ++ )
      {
        WeightedPoint p = oWpList[i].first;
        Point Query(p.x, p.y, p.z);
        Neighbor_search search( oBndPointTree, Query, 1 ); // find 1 neighbor
        oWpList[i].first.w = std::sqrt( CGAL::squared_distance( search.begin()->first,
                                                                Query ) )* float( 4 ) / float(3);
      }
    }


    //---------------------------------------------
    //  Calculate weight of each point
    //---------------------------------------------
    WeightedPointList CalculateWeight( PointDimList & oPointList,
                                       float fEdgeLength )
    {
      WeightedPointList oWpList;
      for( int i = 0; i < oPointList.size(); i ++ )
      {
        Point Query = oPointList[i].first;
        Neighbor_search search( oBndPointTree, Query, 1 ); // find 1 neighbor
        WeightedPoint p( oPointList[i].first.x(),
                         oPointList[i].first.y(),
                         oPointList[i].first.z(), 0 );

        if( oPointList[i].second >=4 )   // quad points get special treatment
        {
          p.w = std::min( double(5)/double(3) * fEdgeLength,
                          std::sqrt( CGAL::squared_distance( search.begin()->first,
                                                             Query ) ) );
        }
        else
        {
          p.w = std::min( double(4)/double(3) * fEdgeLength,
                          std::sqrt( CGAL::squared_distance( search.begin()->first,
                                                             Query ) ) ) * double( 3 ) / double( 4 );
        }
        oWpList.push_back( p );
      }
      return oWpList;
    }
    //---------------------------------------------
    //  GetWeightedPoints
    //---------------------------------------------
    WeightedPointList GetWeightedPoints() const
    {
      return oNewBndPoints;
    }

  };
  
}


#endif
