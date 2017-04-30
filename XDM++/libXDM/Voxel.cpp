///////////////////////////////////////////////////////////////////
//  File:    Voxel.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Implementation of Voxel.h
//
//
/////////////////////////////////////////////////////////////////
#include "Voxel.h"


//-----------------------------------
//
// Useful Voxel Operations
//
//-----------------------------------
// namespace MeshGeometry
// {
//   //-----------------------------------------------------------------
//   // VoxelDifference  - returns the sum of the square sum of the difference between the distances of
//   //                    the three vertices
//   //---------------------------------------------------------------
//   Float VoxelDifference( const CXDMVoxel &v1,  const CXDMVoxel &v2 )
//   {
//     SVector3 tmp;
//     Float fDiff;
//     Float fMinDiff = MAX_FLOAT;

//     // vertices may be cyclic permutations of each other
//     for ( int i = 0; i < 3; i ++ )
//     {
//       fDiff = 0;
//       for (int j = 0; j < 3; j ++ )
//       {
//         tmp = v1.pVertex[j] - v2.pVertex[ (i + j)% 3];
//         fDiff += tmp.GetLength();
//       }
//       if ( fDiff < fMinDiff )
//         fMinDiff = fDiff;
//     }

//     return fMinDiff;
//   }
  
  
// //   //---------------------------------------------------------------
// //   // GetBoundingVoxel   - Returns a minimum voxel that bounds both of v1 and v2.
// //   //  
// //   // Note: the use of oBoundingVoxel is to prepare the compiler for return value optimization
// //   //---------------------------------------------------------------
// //   CXDMVoxel GetBoundingVoxel( const CXDMVoxel & v1, const CXDMVoxel & v2)
// //   {
    
// //     CXDMVoxel oBoundingVoxel;
    
// //     // if the two are overlapping, return one of them
// //     if ( v1.Contains( v2 ) && v2.Contains( v1 ) )
// //     {
// //       oBoundingVoxel = v1;
// //     }
// //     else if ( v1.Contains( v2 ) && ! v2.Contains( v1) ) // if v2 in v1
// //     {
// //       oBoundingVoxel = v1;
// //     }
// //     else if ( !v1.Contains( v2 ) &&  v2.Contains( v1) ) // if v1 in v2
// //     {
// //       oBoundingVoxel = v2;
// //     }
// //     else   // there is an intersection - this is tricky.
// //     {
// //       BBox2D oVoxelBBox( Point(0, 0), Point(0, 0));
      
// //       for ( int i = 0; i < 3; i ++ ) 
// //       {
// //         oVoxelBBox = Union( oVoxelBBox, Point( v1.pVertex[i].m_fX, v1.pVertex[i].m_fY ) ); 
// //         oVoxelBBox = Union( oVoxelBBox, Point( v2.pVertex[i].m_fX, v2.pVertex[i].m_fY ) ); 
// //       }
      
// //      // Find max sidelength
// //      Float fMaxX = oVoxelBBox.pMax.x - oVoxelBBox.pMin.x;
// //      Float fMaxY = oVoxelBBox.pMax.y - oVoxelBBox.pMin.y;
// //      Float fMaxEdge = max( fMaxX, fMaxY );
     
// //      // Expand the triangle about its center
// //      oBoundingVoxel = Larger( v1, v2 );             // find the larger voxel
// //      oBoundingVoxel.Resize( fMaxEdge );
     
     
// //     }
    
// //     return oBoundingVoxel;
// //   }
  
//   //---------------------------------------------------------------
//   //  Larger - Return the larger of the two voxels
//   //---------------------------------------------------------------
//   CXDMVoxel Larger( const CXDMVoxel &v1, const CXDMVoxel & v2)
//   {
//     CXDMVoxel oRet;
    
//     if( v1.fSideLength > v2.fSideLength )
//     {
//       oRet = v1;
//     }
//     else
//     {
//       oRet = v2;
//     }
    
//     return oRet;
//   }

//   //---------------------------------------------------------------
//   //  HasOverlappingEdge
//   //
//   //  Return true if there exist an edge that's overlapping
//   //
//   //---------------------------------------------------------------
//   Bool HasOverlappingEdge( const CXDMVoxel & v1, const CXDMVoxel & v2)
//   {
//     Float fMinSideLength = std::min( v1.fSideLength, v2.fSideLength );
//     Float fMaxSideLength = std::max( v1.fSideLength, v2.fSideLength );
//     Float fAngularTolerance = std::min( 0.025, Float( 0.25 ) * fMinSideLength / fMaxSideLength );
//     Float fPointTolerance = Float( 0.05 ) * fMinSideLength;
    
//     for( Int i = 0; i < 3; i ++ )
//     {
//       for ( Int j = 0; j < 3; j ++)
//       {
//         Bool bEdgeOverlap = MeshGeometry::EdgesOverlap( v1.pVertex[ i ], v1.pVertex[ (i + 1) % 3 ],
//                                                         v2.pVertex[ j ], v2.pVertex[ (j + 1) % 3 ],
//                                                         fAngularTolerance, fPointTolerance );
//         if ( bEdgeOverlap )
//           return true;
        
//       }
//     }
//     return false;
//   }

//   //---------------------------------------------------------------
//   //  GetOverlappingEdge
//   //
//   //  Return true if there exist an edge that's overlapping
//   //
//   //---------------------------------------------------------------
//   Bool GetOverlappingEdge( std::pair<SVector3, SVector3> &oResEdge,
//                            const CXDMVoxel &v1, const CXDMVoxel &v2 )
//   {
//     Float fMinSideLength = std::min( v1.fSideLength, v2.fSideLength );
//     Float fMaxSideLength = std::max( v1.fSideLength, v2.fSideLength );
//     Float fAngularTolerance = std::min( 0.025, Float( 0.25 ) * fMinSideLength / fMaxSideLength );
//     Float fPointTolerance = Float( 0.05 ) * fMinSideLength;
    
//     for( Int i = 0; i < 3; i ++ )
//     {
//       for ( Int j = 0; j < 3; j ++)
//       {
//         Bool bEdgeOverlap = MeshGeometry::EdgesOverlap( v1.pVertex[ i ], v1.pVertex[ (i + 1) % 3 ],
//                                                         v2.pVertex[ j ], v2.pVertex[ (j + 1) % 3 ],
//                                                         fAngularTolerance, fPointTolerance );
//         if ( bEdgeOverlap )
//         {
//           oResEdge.first  = v1.pVertex[ i ];
//           oResEdge.second = v1.pVertex[ (i + 1) % 3 ];
//           return true;
//         }
//       }
//     }
//     return false;
//   }

//   //-----------------------------------
//   //
//   //  GetConnectingEdge
//   //
//   //  Return a pair of vertices in SVector3
//   //-----------------------------------
//   void GetConnectingEdge( std::pair<SVector3, SVector3> &oResEdge,
//                           const CXDMVoxel &v1, const CXDMVoxel &v2 )
//   {
//     oResEdge.first  = v1.GetCenter();
//     oResEdge.second = v2.GetCenter();
//   }


  
// };
