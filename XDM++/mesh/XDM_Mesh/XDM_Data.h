//--------------------------------------------
//
//  Test data structure for meshing.
//  THIS IS ONLY A TEST!  DO NOT USE THIS FOR PRODUCTION
//  CODE!!!!
//--------------------------------------------


#ifndef XDM_DATA_H_
#define XDM_DATA_H_

#include <boost/multi_array.hpp>
#include <CGAL/Bbox_3.h>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/Memory_sizer.h>
#include "BoundaryPointProcessor.h"
#include "3dMath.h"
namespace CGAL
{
  namespace XDM_test
  {
    
    class FIndicator : public std::unary_function<int, double>
    {
      const int label;
    public:
      FIndicator(int i) : label(i) {};
      double operator()(int x) const { return (x == label) ? 1. : 0.; }
    };
    
    static const int nVertices = 8;
    static const int nCubeDimension = 3;
    static const int VertexOrderGrayCode[ nVertices ][ nCubeDimension ] =
        {
          { 0, 0, 0 },  // 0
          { 1, 1, 1 },  // 1
          { 0, 0, 1 },  // 2
          { 1, 1, 0 },  // 3
          { 0, 1, 0 },  // 4
          { 1, 0, 1 },  // 5
          { 0, 1, 1 },  // 6
          { 1, 0, 0 }   // 7
        };

    template< class Kernel >
    class XDM_Data
    {
    public:
      typedef CGAL::Bbox_3 Bbox_3;
      
    public:
      
      Bbox_3 compute_bounding_box() const
      {
        
        float dI = fLengthI / float( 100 );
        float dJ = fLengthJ / float( 100 );
        float dK = fLengthK / float( 100 );
        
        return Bbox_3( -dI, -dJ, -dK,
                       fLengthI * nPixelI + dI,
                       fLengthJ * nPixelJ + dJ,
                       fLengthK * nPixelK + dK );
       
      }
      
      int trilinear_interpolation( float x, float y, float z ) const
      {
        return TrilinearLabelInterpolation( x,  y, z );
      }
    
      int   operator()( float x, float y, float z ) const
      {
        int nI = static_cast<int>( x / fLengthI );
        int nJ = static_cast<int>( y / fLengthJ );
        int nK = static_cast<int>( z / fLengthK );
        if( InBound( nI, nJ, nK ) )
          return oData[ nI ][ nJ ][ nK ];     
        else
          return 0;
      }

      void Set( int nI, int nJ, int nK, int nValue )
      {
        oData[ nI ][ nJ ][ nK ] = nValue;
      }
      void  Resize( int nI, int nJ, int nK )
      {
        nPixelI = nI;
        nPixelJ = nJ;
        nPixelK = nK;
        oData.resize( boost::extents[nI][nJ][nK] );
      }
      void  SetPixelSize( float fI, float fJ, float fK )
      {
        fLengthI = fI;
        fLengthJ = fJ;
        fLengthK = fK;
      }

      //-----------------------------------------
      //  CheckNgb
      //-----------------------------------------
      void CheckNgb( std::set<int> & oIDMap,
                     std::map< int, TestPoint > & oTestPointMap,
                     const int i, const int j, const int k,
                     const int dx, const int dy, const int dz ) const
      {
        int nCurIndex = oData[i + dx ][ j + dy ][k + dz ];
        
        oIDMap.insert( nCurIndex );
        if( oTestPointMap.find( nCurIndex ) == oTestPointMap.end() )
        {
          oTestPointMap.insert( std::make_pair( nCurIndex,
                                                TestPoint( i + dx,
                                                           j + dy,
                                                           k + dz,
                                                           1 ) ) );
        }
        else
        {

          oTestPointMap[ nCurIndex ].x += (i + dx);
          oTestPointMap[ nCurIndex ].y += (j + dy);
          oTestPointMap[ nCurIndex ].z += (k + dz);
          oTestPointMap[ nCurIndex ].w += 1;
        }
        
      }
      
      //-----------------------------------------
      //  GetBndPoints
      //-----------------------------------------
      std::vector< TestPoint > GetBndPoints() const
      {
        return oBndPoints_;        
      }

      std::vector< TestPoint > GetTopologyPoints() const
      {
        return oTopologyPoints_;        
      }

      //-----------------------------------------
      //  MakeTopologyPreservingPoints
      //-----------------------------------------
      void MakeTopologyPreservingPoints()
      {
        typedef GeneralLib::SVector3 SVector3;
        std::map<int, SVector3> GrainToCenterMap;
        std::map<int, int>      GrainToSizeMap;
        for( int i = 1; i < nPixelI - 1; i ++ )
        {
          for( int j = 1; j < nPixelJ -1 ; j ++ )
          {
            for( int k = 1; k < nPixelK -1 ; k ++ )
            {
              int nCurIndex = oData[i][j][k];
              if( GrainToCenterMap.find( nCurIndex ) == GrainToCenterMap.end() )
              {
                GrainToCenterMap[ nCurIndex ] = SVector3(i, j, k); 
                GrainToSizeMap  [ nCurIndex ] = 1;
              }
              else   // else do running average
              {
                int nGrainSize = GrainToSizeMap[ nCurIndex ];
                SVector3 v  = static_cast<Float>(nGrainSize) * GrainToCenterMap[ nCurIndex ];
                nGrainSize ++;
                GrainToCenterMap[ nCurIndex ] = ( v + SVector3(i, j, k) ) / static_cast<Float>( nGrainSize ); 
                GrainToSizeMap  [ nCurIndex ] = nGrainSize;                
              }
            }
          }
        }

        for( std::map<int, SVector3>::iterator pIter = GrainToCenterMap.begin();
             pIter != GrainToCenterMap.end(); ++ pIter )
        {
          SVector3 v = pIter->second;
          oTopologyPoints_.push_back( TestPoint( v.m_fX, v.m_fY, v.m_fZ, 0 ) );
        }
        
      }
      
      //-----------------------------------------
      //  MakeBndPoints
      //-----------------------------------------
      void MakeBndPoints( float fSmoothingLength,
                          float fPointSpacing,
                          float fMinEdgeLength  )
      {
        std::vector< std::pair< TestPoint, int > > oBndPointList;
        std::ofstream os("test.csv");
        std::cerr << nPixelI << " " << nPixelJ << " " << nPixelK << std::endl;
        std::cout << "Pixel Size: " <<  fLengthI << " " << fLengthJ << " " << fLengthK << std::endl;
        for( int i = 0; i < nPixelI - 1; i ++ )
        {
          for( int j = 0; j < nPixelJ -1 ; j ++ )
          {
            for( int k = 0; k < nPixelK -1 ; k ++ )
            {
              int nJunction = 0;
              std::set<int> oIDMap;
              std::map< int, TestPoint > oTestPointMap;

              for( int n = 0; n < nVertices; n ++ )
                CheckNgb( oIDMap, oTestPointMap, i, j, k,
                          VertexOrderGrayCode[n][0],
                          VertexOrderGrayCode[n][1],
                          VertexOrderGrayCode[n][2] );
              
              
              std::map<int, TestPoint>::iterator pCur;

              TestPoint oBndPoint(0, 0, 0, 0);
              for( pCur = oTestPointMap.begin();
                   pCur != oTestPointMap.end(); ++ pCur )
              {
                oBndPoint.x += pCur->second.x;
                oBndPoint.y += pCur->second.y;
                oBndPoint.z += pCur->second.z;
                oBndPoint.w += pCur->second.w;
              }
              
              oBndPoint.x /= oBndPoint.w;
              oBndPoint.y /= oBndPoint.w;
              oBndPoint.z /= oBndPoint.w;
              
              float fX = oBndPoint.x * fLengthI;
              float fY = oBndPoint.y * fLengthJ;
              float fZ = oBndPoint.z * fLengthK;
              
              float fScale = sqrt( ( fLengthI * fLengthI +  fLengthJ * fLengthJ + fLengthK * fLengthK ) / float(3) );
              if( oIDMap.size() == 3 )
              {
                TestPoint t( fX, fY, fZ, 0.3 * fScale );
                oBndPoints_.push_back( t );
                oBndPointList.push_back( std::make_pair( t, oIDMap.size()  ) );
              }
              else if( oIDMap.size() > 3 )
              {
                TestPoint t( fX, fY , fZ, 0.3 * fScale );
                oBndPoints_.push_back( t );
                oBndPointList.push_back( std::make_pair( t, oIDMap.size()  ) );
              }

            }
          }
        }

        XDM::BoundaryPointProcessor<Kernel> oBndPointProcess;
        
        oBndPointProcess.Initialize( oBndPointList, fSmoothingLength,
                                     fPointSpacing, fMinEdgeLength  );
        oBndPoints_ = oBndPointProcess.GetWeightedPoints();
        
        for( int i = 0; i < oBndPoints_.size(); i ++ )
        {
          os << oBndPoints_[i].x << " "
             << oBndPoints_[i].y << " "
             << oBndPoints_[i].z << " "
             << oBndPoints_[i].w << std::endl;
        }
        os.close();
        
      }
      
      void ReadBndPoints( const std::string & oFilename )
      {
        std::ifstream is;
        is.open( oFilename.c_str() );


        int count = 0;
        int tmp = 0;
        is >> std::skipws
             >> count >> std::skipws
             >> tmp >> std::skipws
             >> tmp >> std::skipws
             >> tmp;
          
        for( int i = 0; i < count; i ++ )
        {
          TestPoint t(0, 0, 0, 0);
          is >> std::skipws
             >> t.x >> std::skipws
             >> t.y >> std::skipws
             >> t.z >> std::skipws
             >> t.w;
          oBndPoints_.push_back( t );
          std::cout << t.x << " " << t.y << " " << t.z << " " << t.w << std::endl;
        }
        std::cout << oBndPoints_.size() << std::endl;
        is.close();
      }
      
    private:
      typedef boost::multi_array< int, 3 > GridType;
      GridType oData;

      float fLengthI;
      float fLengthJ;
      float fLengthK;
      int nPixelI;
      int nPixelJ;
      int nPixelK;
      std::vector< TestPoint > oBndPoints_;
      std::vector< TestPoint > oTopologyPoints_;
     
      
      bool InBound( int nI, int nJ, int nK) const
      {
        return nI >= 0 && nI < nPixelI
          && nJ >= 0 && nJ < nPixelJ
          && nK >= 0 && nK < nPixelK ;
      }

      
      //
      //  linear interpolate f(x) at x between end points x0 and x1
      //  Note that the interval is [x0, x1), x0 <= x < x1, 
      //  x0 < x1
      float Interpolate( float x0, float x1,
                         float f0, float f1, float x) const
      {
        float f = ( x - x0 ) * f1 + ( x1 - x ) * f0  ;
        return f;
      }

      //-------------------------------------------------------------------
      //  interpolate given a label, find the indicator function
      //  value at x, y, z
      //
      //-------------------------------------------------------------------
      float InterpolateLabelIndicator( float x, float y, float z, int nLabel ) const
      {
        double pValues[nVertices];
        float cX = x / fLengthI;
        float cY = y / fLengthJ;
        float cZ = z / fLengthK;
        
        int nX = static_cast<int>( cX );
        int nY = static_cast<int>( cY );
        int nZ = static_cast<int>( cZ );

        FIndicator FLabelIndicatorFn( nLabel );
        
        for( int i = 0; i < nVertices; i ++)
        {
          int nNgbX = VertexOrderGrayCode[i][0] + nX;
          int nNgbY = VertexOrderGrayCode[i][1] + nY;
          int nNgbZ = VertexOrderGrayCode[i][2] + nZ;

          if( InBound( nNgbX, nNgbY, nNgbZ ) )
            pValues[i] = FLabelIndicatorFn( oData[ nNgbX ][ nNgbY ][ nNgbZ ] );
     
        }

   
        // numbering comes from gray code
        //  (look at the definition of VertexOrderGrayCode)
        float fx00 = Interpolate( nX, nX + 1, pValues[0], pValues[7], cX );
        float fx01 = Interpolate( nX, nX + 1, pValues[2], pValues[5], cX );
        float fx10 = Interpolate( nX, nX + 1, pValues[4], pValues[3], cX );
        float fx11 = Interpolate( nX, nX + 1, pValues[6], pValues[1], cX );
        float fxy0 = Interpolate( nY, nY + 1, fx00, fx10, cY );
        float fxy1 = Interpolate( nY, nY + 1, fx01, fx11, cY ); 

        return Interpolate( nZ, nZ + 1, fxy0, fxy1, cZ ); 
      }

      //--------------------------------------------
      //  Trilinear interpolation of labels using
      //   an indicator function
      //--------------------------------------------
      int TrilinearLabelInterpolation( float x, float y, float z ) const
      {

        if( x < 0 || y < 0 || z < 0 )
          return 0;
        
        int nX = static_cast<int>( x / fLengthI );
        int nY = static_cast<int>( y / fLengthJ );
        int nZ = static_cast<int>( z / fLengthK );
        std::set<int> oLabels;
        for( int i = 0; i < nVertices ; i ++ )
        {
          int nNewX = VertexOrderGrayCode[i][0] + nX;
          int nNewY = VertexOrderGrayCode[i][1] + nY;
          int nNewZ = VertexOrderGrayCode[i][2] + nZ;

          if( InBound( nNewX, nNewY, nNewZ ) )
            oLabels.insert( oData[ nNewX ][ nNewY ][ nNewZ ] );
          else
            return 0;
        }
        if( oLabels.size() == 1)
          return *oLabels.begin();
        

        float fBestValue = 0;
        int nBestLabel   = 0;
        for( std::set<int>::const_iterator pIter = oLabels.begin();
             pIter != oLabels.end(); ++ pIter )
        {
          float fCurVal = InterpolateLabelIndicator( x, y, z, *pIter );
          if( fCurVal > fBestValue )
          {
            nBestLabel = *pIter;
            fBestValue = fCurVal;
          }
        }

        return nBestLabel;
      }
      
    };
    
  }
  
}







#endif
