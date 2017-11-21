//----------------------------------------------------
//  MeshAlignment.h
//
//  Author:   Frankie Li
//  Purpose:  Utilities for mesh alignment
//
//     --------------  This is buggy and mostly does not work --------------------
//
//
//----------------------------------------------------

#ifndef MESH_ALIGNMENT_
#define MESH_ALIGNMENT_

#include <set>
#include <map>
#include <CGAL/IO/File_medit.h>
#include <CGAL/random_selection.h>
#include "3dMath.h"
#include "Sampling.h"
#include "Quaternion.h"

namespace Pandora
{
  template< class C3T3 >
  class MeshAligner
  {

  public:
    typedef typename C3T3::Triangulation::Vertex_iterator Vertex_iterator;
    typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertex_iterator;
    typedef typename C3T3::Triangulation::Vertex          Vertex;
    typedef typename C3T3::Triangulation::Geom_traits::Bare_point Point3;
    typedef typename C3T3::Triangulation::Geom_traits::RT RT;
    typedef typename C3T3::Cell_handle                    Cell_handle;
    typedef GeneralLib::SMatrix3x3 SMatrix3x3;
    typedef GeneralLib::SVector3   SVector3;
    mutable GeneralLib::CUniformRandomReal RandGen;

    //-------------------------------------------------
    // SelectSampleSurfacePoint
    // Purpose:  Select points from the surface of the sample
    //           for alignment purpose.
    //-------------------------------------------------
    void SelectSampleSurfacePoint()
    {
      
    }

    //-------------------------------------------------
    //  DistanceCost
    //-------------------------------------------------
    template< class Iter >
    float DistanceCost( const C3T3 & c3t3_ref,
                        Iter pCur, Iter pEnd,
                        const SMatrix3x3 & R, const SVector3 & T ) 
    {
      float cost = 0;
      int nCounted = 0;
      int nStart = pEnd - pCur;
      for(; pCur != pEnd; ++ pCur )
      {
        Point3 p = pCur->point();
        TransformPoint( p, R, T );
        Cell_handle ch = c3t3_ref.triangulation().locate( p );
        if( ! c3t3_ref.triangulation().is_infinite( ch ) )
        {
          cost += ( ch->vertex(1)->point() - p ).squared_length();
          nCounted ++;
        }
      }
      nCounted = std::max( nCounted, 1 );
        
      return cost / float( nCounted );
    }
    
    //-------------------------------------------------
    //  PointAlignMesh
    //-------------------------------------------------
    void PointAlignMesh( const C3T3 & c3t3_ref,
                         const C3T3 & c3t3_,
                         SMatrix3x3 & Rotation,
                         SVector3   & Translation,
                         const SVector3   & TransStepSize,
                         int nTestPoints, int nSteps,
                         float SO3_Radius )
    {
      std::vector< Vertex > VertexList;  // Thank you C++ for making this necessary

      //---------------------------------
      //  Need to change to select boundar vertices only
      //---------------------------------
      for( Finite_vertex_iterator pIter = c3t3_.triangulation().finite_vertices_begin();
           pIter != c3t3_.triangulation().finite_vertices_end(); ++ pIter )
        VertexList.push_back( *pIter );
      
      std::random_shuffle( VertexList.begin(), VertexList.end() );
      assert( (VertexList.end() - VertexList.begin() ) > nTestPoints );
      using namespace GeneralLib;
      UniformGrid::CQuaternionGrid RandSO3Gen;
      std::vector<SQuaternion> RandomLocalPoints
        = RandSO3Gen.GetRandomLocalGrid( SO3_Radius, nSteps );

      float fBestCost =  DistanceCost( c3t3_ref, VertexList.begin(),
                                       VertexList.begin() + nTestPoints,
                                       Rotation, Translation );
      for( int i = 0; i < nSteps; i ++ )
      {
        float fDx = RandGen.GetRandomVariable( -TransStepSize.m_fX, TransStepSize.m_fX );
        float fDy = RandGen.GetRandomVariable( -TransStepSize.m_fY, TransStepSize.m_fY );
        float fDz = RandGen.GetRandomVariable( -TransStepSize.m_fZ, TransStepSize.m_fZ );
        SMatrix3x3 CurrentRot = RandomLocalPoints[i].GetRotationMatrix3x3() * Rotation ;
        SVector3 CurrentPos   = Translation + SVector3( fDx, fDy, fDz );
        float fCost = DistanceCost( c3t3_ref, VertexList.begin(),
                                    VertexList.begin() + nTestPoints,
                                    CurrentRot, CurrentPos ); 
        if( fCost < fBestCost )
        {
          std::cout << "Best Cost " << fBestCost <<" New Cost " << fCost << std::endl;
          fBestCost = fCost;
          Translation = CurrentPos;
          Rotation    = CurrentRot;

          std::cout << " Trans : " << Translation << std::endl;
          std::cout << " Rotation : " << std::endl << Rotation << std::endl;        
        }
      }


    }

  };  // end MeshAlignment
}

#endif
