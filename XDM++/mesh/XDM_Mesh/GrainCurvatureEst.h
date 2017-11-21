//----------------------------------------------------
//  GrainCurvatureEst.h
//
//  Author:   Frankie Li
//  Purpose:  Estimate Grain curvature given a Complex_3_in_triangulation_3
//
//----------------------------------------------------

#ifndef GRAIN_CURVE_H
#define GRAIN_CURVE_H

#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Lapack/Linear_algebra_lapack.h>

#include <CGAL/IO/File_medit.h>
#include "SimpleMeshVTKWriter.h"
#include <map>
#include <vector>
#include "MeshUtilities.h"

namespace Pandora
{
  namespace Details
  {
    struct SCurvatureData
    {
      SCurvatureData(){}
      SCurvatureData( float fK1_, float fK2_, float fCondition_ ):
        fK1( fK1_ ), fK2( fK2_ ), fConditionNum( fCondition_ ) {}
      
      float fK1;
      float fK2;
      float fConditionNum;
    };
  }

  template< class C3T3 >
  void PrintBndIDToGrainID( std::ostream & os,
                            const C3T3 & c3t3 );

  //----------------------------------------------------
  //  CurveEstimator
  //
  //----------------------------------------------------
  template< class C3T3 >
  class CurveEstimator
  {
  public:
    typedef typename C3T3::Subdomain_index Subdomain_index;
    typedef typename C3T3::Surface_index   Surface_index;
    typedef typename C3T3::Triangulation::Vertex_handle Vertex_handle;
    typedef typename C3T3::Facet           Facet;
    typedef typename C3T3::Cell_handle     Cell_handle;

    typedef typename C3T3::Facet_iterator  Facet_iterator;
    typedef typename C3T3::Cell_iterator   Cell_iterator;
    
    typedef std::multimap< Subdomain_index, Cell_handle > GrainTetMap;
    typedef std::multimap< Subdomain_index, Facet >       GrainFacetMap;

    // for fitting
    struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
    typedef typename CGAL::Monge_via_jet_fitting<K>   Monge_via_jet_fitting;
    
  private:
    
    typedef typename Details::SCurvatureData                   SCurvatureData;
    typedef typename std::map< Vertex_handle, SCurvatureData > VertexCurveMap;
    typedef typename std::pair< SCurvatureData, SCurvatureData> DoubleSideCurve;
    typedef typename std::map< Vertex_handle, DoubleSideCurve > VertexDoubleCurveMap;
    
    typedef typename std::map< Facet, SCurvatureData >         FacetCurveMap; 
    typedef typename std::set< Vertex_handle >                 Vertex_set;

    typedef typename std::map< Subdomain_index, FacetCurveMap > SubdomainIDToCurveMap;
    
    const C3T3            & c3t3_;

    GrainTetMap           GrainIDToTetMap;
    GrainFacetMap         GrainIDToFacetMap;
    SubdomainIDToCurveMap GrainIDToCurveMap;
    std::set< Subdomain_index > oDomainIDSet;
    
    //-----------------------------------
    //  Property Maps
    //-----------------------------------
    struct MeanCurvatureMap
    {
      float operator() ( const SCurvatureData & o ) const
      {  return ( o.fK1 + o.fK2 ) / 2.; }
    };
    
    struct AbsMeanCurvatureMap
    {
      float operator() ( const SCurvatureData & o ) const
      {  return ( fabs( o.fK1 ) + fabs( o.fK2 ) ) / 2.; }
    };

    struct ConditionNumberMap
    {
      float operator() ( const SCurvatureData & o ) const
      {  return o.fConditionNum; }
    };
    
    //-----------------------------------
    // BuildGrainTetMap
    //-----------------------------------
//     void BuildGrainTetMap( )
//     {
//       CGAL::Default_cell_index_pmap<C3T3> oCellToID( c3t3_ );
//       for( Cell_iterator pCur = c3t3_.cells_begin();
//            pCur != c3t3_.cells_end(); ++pCur )
//       {
//         Subdomain_index nID = c3t3_.subdomain_index( pCur );
//         Cell_handle ch = pCur;
//         GrainIDToTetMap.insert( std::make_pair( nID, ch ) );
//         oDomainIDSet.insert( nID );
//       }
//     }

    //-----------------------------------
    // BuildGrainFacetMap
    //-----------------------------------
//     void BuildGrainFacetMap( )
//     {
//       for( Facet_iterator pCur = c3t3_.facets_begin();
//            pCur != c3t3_.facets_end(); ++pCur )
//       {
//         Surface_index facet_index = c3t3_.surface_index( *pCur );
//         GrainIDToFacetMap.insert( std::make_pair( facet_index.first,  *pCur ) );
//         GrainIDToFacetMap.insert( std::make_pair( facet_index.second, *pCur ) );
//       }
//     }

    //-----------------------------------------------------------------------------
    //  WriteFacetsToPolydataVtk
    //-----------------------------------------------------------------------------
    template< class PropertyMap >
    void WriteFacetData( std::ostream & os,
                         const PropertyMap & PMap )
    {
      typedef typename C3T3::Triangulation Tr;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename C3T3::Facet_iterator Facet_iterator;

      typedef typename SubdomainIDToCurveMap::iterator GrainIter;
      typedef typename FacetCurveMap::iterator FacetCurveIter;

      for( GrainIter pGrain = GrainIDToCurveMap.begin();
           pGrain != GrainIDToCurveMap.end(); ++ pGrain )
      {
        for( FacetCurveIter pFacetCurve = pGrain->second.begin();
             pFacetCurve != pGrain->second.end(); ++ pFacetCurve )
        {
          Facet f = pFacetCurve->first;
          os << PMap( pFacetCurve->second ) << std::endl;
        }
      }
    }

    //-----------------------------------------------------------------------------
    //  WriteFacetsToPolydataVtk
    //-----------------------------------------------------------------------------
    template< class VertexToIDMap >
    void WriteFacetsToPolydataVtk( std::ostream & os,
                                   const VertexToIDMap & VertexIDPMap )
    {
      typedef typename C3T3::Triangulation Tr;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename C3T3::Facet_iterator Facet_iterator;

      typedef typename SubdomainIDToCurveMap::iterator GrainIter;
      typedef typename FacetCurveMap::iterator FacetCurveIter;

      for( GrainIter pGrain = GrainIDToCurveMap.begin();
           pGrain != GrainIDToCurveMap.end(); ++ pGrain )
      {
        for( FacetCurveIter pFacetCurve = pGrain->second.begin();
             pFacetCurve != pGrain->second.end(); ++ pFacetCurve )
        {
          Facet f = pFacetCurve->first;
          os << 3 << " ";
          for (int i=0; i<4; i++)
          {
            if (i != f.second)   // facet is represeted by a cell and an index opposite the facet
                                 // Therefore, f.second (opp index) must not be included
            {
              const Vertex_handle& vh = f.first->vertex(i);
              typename VertexToIDMap::const_iterator pFound = VertexIDPMap.find( vh );
              if( pFound != VertexIDPMap.end() )
                os << pFound->second << " ";
              else
                os << -1 << " ";  // error condition
            }
          }
          os << std::endl;
        }
      }
      
    }
    
    //-----------------------------------------------------------------------------
    //  EstimateGrainFacetCurvature
    //-----------------------------------------------------------------------------
    void EstimateGrainFacetCurvature( FacetCurveMap & FacetToCurveMap,
                                      Subdomain_index nID ) const;
    
    //-----------------------------------------------------------------------------
    //  CalculateGrainCurvature
    //-----------------------------------------------------------------------------
    void CalculateGrainCurvature   ( VertexCurveMap & VertexToCurveMap,
                                     Subdomain_index nID ) const;

    void GetNthNeighboringVertices( Vertex_set & oNgbSet,
                                    Vertex_handle oVertex,
                                    int nMinElements,
                                    int nMaxDepth ) const;

    void GetNSurfaceNeighbors( Vertex_set & oNgbSet,
                               Vertex_handle oVertex,
                               Subdomain_index nDomainID,
                               int nMinElements,
                               int nMaxDepth ) const;
    //-----------------
    //  Default ctor
    //-----------------
    CurveEstimator(); // not allowed

  public:

    CurveEstimator( const C3T3 & obj ) : c3t3_( obj ) {}
    
    void WriteGrainCurveVtkFile( const std::string & filename );
  };
  
}


//----------------------------------------------------
//
//  I M P L E M E N T A T I O N
//
//----------------------------------------------------
namespace Pandora
{
  template< class C3T3 >
  void PrintBndIDToGrainID( std::ostream & os,
                            const C3T3 & c3t3 )
  {
    typedef typename C3T3::Facet_iterator Facet_iterator;
    typedef typename C3T3::Surface_index  Surface_index;
    
    CGAL::Default_facet_index_pmap< C3T3 > oBndToIDMap( c3t3 );
    for( Facet_iterator pCur = c3t3.facets_begin();
         pCur != c3t3.facets_end(); ++pCur )
    {
      Surface_index facet_index = c3t3.surface_index( *pCur );
      os << oBndToIDMap.surface_index( facet_index ) << " "
         << facet_index.first << " "
         << facet_index.second << std::endl;
    }
  }

  
  //-----------------------------------------------------------------------------
  //  EstimateGrainFacetCurvature
  //-----------------------------------------------------------------------------
  template< class C3T3 >
  void CurveEstimator<C3T3>::EstimateGrainFacetCurvature( FacetCurveMap & FacetToCurveMap,
                                                          Subdomain_index nID ) const
  {
    typedef typename GrainFacetMap::const_iterator GrainFacetIter;
    VertexCurveMap oVertexCurvature;
    CalculateGrainCurvature( oVertexCurvature, nID );
    
    std::pair<GrainFacetIter, GrainFacetIter> pRange;
    pRange = GrainIDToFacetMap.equal_range( nID );
    
    for( GrainFacetIter pCur = pRange.first; pCur != pRange.second; ++ pCur )
    {
      SCurvatureData oFacetCurve(0, 0, 0);
      
      for( int i = 0; i < 4; i ++ )
      {
        if( i != pCur->second.second )
        {
          Cell_handle  c = pCur->second.first;
          Vertex_handle v = c->vertex( i );
          SCurvatureData & oVertexCurve = oVertexCurvature[v]; 
          oFacetCurve.fK1 += oVertexCurve.fK1;
          oFacetCurve.fK2 += oVertexCurve.fK2;
          oFacetCurve.fConditionNum += oVertexCurve.fConditionNum;
        }
      }
      oFacetCurve.fK1 /= float( 3 );
      oFacetCurve.fK2 /= float( 3 );
      oFacetCurve.fConditionNum /= float( 3 );
      FacetToCurveMap.insert( std::make_pair( pCur->second, oFacetCurve ) );
    }
  }
  
  //-----------------------------------------------------------------------------
  //  CalculateGrainCurvature
  //-----------------------------------------------------------------------------
  template< class C3T3>
  void CurveEstimator<C3T3>::CalculateGrainCurvature ( VertexCurveMap & VertexToCurveMap,
                                                       Subdomain_index nID ) const
  {
    typedef typename GrainTetMap::const_iterator GrainTetIter;
    std::pair< GrainTetIter, GrainTetIter > pRange;
    pRange = GrainIDToTetMap.equal_range( nID );
    
    Vertex_set UniqueVertices;
    for( GrainTetIter pCur = pRange.first; pCur != pRange.second; ++ pCur )
    {
      for( int i = 0; i <= 3; i ++ )
      {
        Cell_handle  CH = pCur->second; // implicit cast
        if( c3t3_.in_dimension( CH->vertex(i) ) <= 2 )  // if Facet
        {
          // find all neighbors on surface
          // calculate curvature
          UniqueVertices.insert( CH->vertex(i) );
        }
      }
    }
    VertexToCurveMap.clear(); // clearing existing solution
    typedef typename Vertex_set::iterator SetIter;
    typedef typename Monge_via_jet_fitting::Monge_form   Monge_form;
    typedef typename C3T3::Triangulation::Vertex::Point  Point;
    for( SetIter pCur = UniqueVertices.begin(); pCur != UniqueVertices.end();
         ++ pCur )
    {
      Vertex_set oNgbSet;
      GetNSurfaceNeighbors( oNgbSet, *pCur, nID, 30,  2 );
      
      std::vector<Point> oPointList;
      for( SetIter it = oNgbSet.begin(); it != oNgbSet.end(); ++ it )
      {
        oPointList.push_back( (*it)->point() );
      }
      
      SCurvatureData oCurve;
      if( oPointList.size() > 15 )
      {
        Monge_form FittedMorgeForm;
        Monge_via_jet_fitting MongeJetFit;
        FittedMorgeForm = MongeJetFit( oPointList.begin(), oPointList.end(), 2, 2 );
        

        oCurve.fConditionNum = MongeJetFit.condition_number();
        oCurve.fK1           = FittedMorgeForm.coefficients()[0];
        oCurve.fK2           = FittedMorgeForm.coefficients()[1];
      }
      VertexToCurveMap.insert( std::make_pair( *pCur, oCurve ) );
    }
  
  }
  
  //-----------------------------------------------------------------------------
  // GetNthNeighboringVertices
  //-----------------------------------------------------------------------------
  template< class C3T3 >
  void CurveEstimator<C3T3>::GetNthNeighboringVertices( Vertex_set & oNgbSet,
                                                        Vertex_handle oVertex,
                                                        int nMinElements,
                                                        int nMaxDepth ) const
  {
    if( nMinElements < oNgbSet.size() )
      return;
    
    if( nMaxDepth <= 0 )
      return;
    
    std::vector< Vertex_handle > oNearestNgb;
    c3t3_.triangulation().finite_adjacent_vertices( oVertex, std::back_inserter( oNearestNgb ) );
    typedef typename std::vector< Vertex_handle >::iterator VS_iter;
    for( VS_iter pCur = oNearestNgb.begin();
         pCur != oNearestNgb.end();
         ++ pCur )
      GetNthNeighboringVertices( oNgbSet, *pCur, nMinElements, nMaxDepth - 1 );
    oNgbSet.insert( oNearestNgb.begin(), oNearestNgb.end() );
  }
  
  //-----------------------------------------------------------------------------
  // GetNSurfaceNeighbors
  //-----------------------------------------------------------------------------
  template< class C3T3 >
  void CurveEstimator<C3T3>::GetNSurfaceNeighbors( Vertex_set & oNgbSet,
                                                   Vertex_handle oVertex,
                                                   Subdomain_index nDomainID,
                                                   int nMinElements,
                                                   int nMaxDepth ) const
  {
    if( nMaxDepth <= 0 )
      return;

    std::vector< Cell_handle > oIncidentCells;
    c3t3_.triangulation().incident_cells( oVertex, std::back_inserter( oIncidentCells ) );
    
    typedef typename std::vector< Cell_handle >::iterator CS_iter;
    Vertex_set oNearestNgb;
    for( CS_iter pCur = oIncidentCells.begin(); pCur != oIncidentCells.end(); ++ pCur )
    {
      if( c3t3_.subdomain_index( *pCur ) == nDomainID )
      {
        for( int i = 0; i <=3; i ++ )  // for each vertex
        {
          Cell_handle CH = *pCur;
          if( c3t3_.in_dimension( CH->vertex( i ) ) <= 2 ) // is facet
            oNearestNgb.insert( CH->vertex(i) );
        }
      }
    }

    oNgbSet.insert( oNearestNgb.begin(), oNearestNgb.end() );
  
    typedef typename Vertex_set::iterator VS_iter; 
    if( oNgbSet.size() < nMinElements )
      for( VS_iter pCur = oNearestNgb.begin(); pCur != oNearestNgb.end(); ++ pCur )
        GetNSurfaceNeighbors( oNgbSet, *pCur, nDomainID, nMinElements, nMaxDepth - 1 );
  }

  //-----------------------------------------------------------------------------
  // WriteGrainCurveVtkFile
  //-----------------------------------------------------------------------------
  template< class C3T3 >
  void CurveEstimator<C3T3>::WriteGrainCurveVtkFile( const std::string & filename )
  {
    Pandora::Details::BuildGrainTetMap  ( GrainIDToTetMap,   c3t3_ );
    Pandora::Details::BuildGrainFacetMap( GrainIDToFacetMap, c3t3_ );

    for( Cell_iterator pCur = c3t3_.cells_begin();
         pCur != c3t3_.cells_end(); ++pCur )
      oDomainIDSet.insert( c3t3_.subdomain_index( pCur ) );
    
    int Number_of_facets = 0;
    
    typedef typename std::set< Subdomain_index >::iterator Iter;
    for( Iter pCur = oDomainIDSet.begin(); pCur != oDomainIDSet.end();
         ++ pCur )
    {
      FacetCurveMap FacetToCurveMap;
      EstimateGrainFacetCurvature( FacetToCurveMap, *pCur );  
      GrainIDToCurveMap.insert( std::make_pair( *pCur, FacetToCurveMap ) );
      Number_of_facets += FacetToCurveMap.size();
    }

    std::ofstream os( filename.c_str() );
    //-----------------------------------------------------------
    //  This is really a special case...
    //
    // --- header information
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << "Polydata file from Pandora" << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET POLYDATA" << std::endl;
    // --- end header information
    
    typedef typename C3T3::Facet_iterator Facet_iterator;
    typedef typename C3T3::Facet Facet;
    typedef typename C3T3::Surface_index  Surface_index;
    typedef typename C3T3::Triangulation::Vertex_handle Vertex_handle;
    CGAL::Default_facet_index_pmap<C3T3>   FacetToIDMap( c3t3_  );
    
    typedef typename C3T3::Triangulation Tr;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename std::map<Vertex_handle, int> VertexPMap;
    VertexPMap VertexToIDMap;
    
    // section header information
    os << "POINTS " << c3t3_.triangulation().number_of_vertices() << " float" << std::endl;
    Details::WriteVerticesToPolydataVtk( os, c3t3_, VertexToIDMap );
    os << "POLYGONS " << Number_of_facets << " " << Number_of_facets * 4 << std::endl;
    WriteFacetsToPolydataVtk( os, VertexToIDMap );

    //  Grain IDs
    os << "CELL_DATA " << Number_of_facets << std::endl;
    os << "SCALARS GrainID int 1 " << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    typedef typename SubdomainIDToCurveMap::iterator GrainIter;
    typedef typename FacetCurveMap::iterator FacetCurveIter;
    for( GrainIter pGrain = GrainIDToCurveMap.begin();
         pGrain != GrainIDToCurveMap.end(); ++ pGrain )
      for( FacetCurveIter pFacetCurve = pGrain->second.begin();
           pFacetCurve != pGrain->second.end(); ++ pFacetCurve )
        os << pGrain->first << std::endl;
    
    os << "SCALARS Mean_Curvature float 1 " << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    WriteFacetData( os, MeanCurvatureMap() );

    
    os << "SCALARS Abs_Mean_Curvature float 1 " << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    WriteFacetData( os, AbsMeanCurvatureMap() );

    os << "SCALARS Condition_Number float 1 " << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    WriteFacetData( os, ConditionNumberMap() );
    
    os.close();
    
    //-----------------------------------------------------------
  }
}

#endif
