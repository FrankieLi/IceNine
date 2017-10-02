////////////////////////////////////////////////////////
//
//  File:  UltraLightWeightGraph.h
//  Purpose:  An extremely barebone graph designed
//            for minimal memory footprint.
//
////////////////////////////////////////////////////////

#ifndef _ULW_GRAPH_H_
#define _ULW_GRAPH_H_

#include "Types.h"
#include <map>
#include <vector>

namespace ULWGraph
{

  //-----------------------------------
  //  Ultra light weight adjacency
  //  graph.  This does not conform to Boost Graph Concept
  //  at all. (hence the ultra light weight)
  //
  // -- if performance permits, conform to BGL concept,
  //    then this can be used with BGL libraries.
  //
  //-----------------------------------
  template< typename VertexT, class LessCmp = less<VertexT> >
  class CAdjacencyGraph
  {
  public:
    typedef bool                                           DataT;    
    typedef const VertexT &                                CFVertexT;
    typedef typename std::map<VertexT, DataT, LessCmp>     NgbListT; 
    typedef typename std::pair<VertexT, DataT>             VertexDataPair;
    typedef typename std::map<VertexT, NgbListT, LessCmp>  AdjGraphT;
    typedef typename std::pair<VertexT, NgbListT>          VertexNgbListPair;
    typedef typename NgbListT::const_iterator              ConstVertexIterT;
    typedef typename NgbListT::iterator                    VertexIterT;
    typedef typename AdjGraphT::const_iterator             ConstNgbListIterT;
    typedef typename AdjGraphT::iterator                   NgbListIterT;
    
  private:
    AdjGraphT oAdjList;

  public:

    //-------------------------------------------
    //  AddVertex
    //   Add vertex to graph if it isn't already there.
    //-------------------------------------------
    void AddVertex   ( CFVertexT u )
    {
      NgbListIterT pIter = oAdjList.find( u );
      if( pIter == oAdjList.end() )
        oAdjList.insert( VertexNgbListPair( u, NgbListT( ) ) );
    }

    //-------------------------------------------
    //  DeleteVertex
    //    Delete vertex specified and all of its associated
    //    edges.  If the vertex is not found, nothing is done.
    //-------------------------------------------
    void DeleteVertex( CFVertexT u )
    {
      ConstVertexIterT pFirst, pEnd;
      boost::tie( pFirst, pEnd ) = GetNeighbors( u );

      for ( ; pFirst != pEnd; ++ pFirst )
        DeleteEdge( u, pFirst->first );

      oAdjList.erase( u );
    }

    //-------------------------------------------
    //  AddEdge
    //    Add edge to graph if it doesn't already exist.
    //    If either u or v are not already in the graph,
    //    they too will be added.
    //-------------------------------------------
    void AddEdge   ( CFVertexT u, CFVertexT v )
    {
      AddVertex( u );
      AddVertex( v );

      NgbListIterT     pIter = oAdjList.find( u );
      ConstVertexIterT pNgb = pIter->second.find( v );
      if( pNgb == pIter->second.end() )
        pIter->second.insert( VertexDataPair( v, true ) );

      pIter = oAdjList.find( v );
      pNgb  = pIter->second.find( u );
      if( pNgb == pIter->second.end() )
        pIter->second.insert( VertexDataPair( u, true ) );
    }

    //-------------------------------------------
    //  DeleteEdge
    //    Delete the edge specified by (u, v).  If
    //    either of the vertex is not found, nothing
    //    is done.
    //-------------------------------------------
    void DeleteEdge( CFVertexT u, CFVertexT v)
    {
      NgbListIterT pUList = oAdjList.find( u );
      NgbListIterT pVList = oAdjList.find( v );

      if( pUList == oAdjList.end() || pVList == oAdjList.end() )
        return;

      pUList->second.erase( v );
      pVList->second.erase( u );
    }

    //-------------------------------------------
    //  GetNeighbors
    //  --  Find all neighbors of u, where ngb is defined
    //      to be anyone with an edge with u.  If u
    //      does not exist, it will be added to the graph.
    //      This is done purely so that a pair of iterators
    //      maybe returned.  (Pretty crappy, I know, but it's
    //      much lighter weight than the alternative.)
    //-------------------------------------------
    std::pair< ConstVertexIterT, ConstVertexIterT >
    GetNeighbors( CFVertexT u )
    {
      AddVertex( u );
      NgbListIterT pIter = oAdjList.find( u );

      typedef std::pair< ConstVertexIterT, ConstVertexIterT > RetT; 
      return RetT( pIter->second.begin(), pIter->second.end() ); 
    }


    ConstNgbListIterT Begin() const { return oAdjList.begin(); }
    ConstNgbListIterT End() const { return oAdjList.end(); }

    
    //-------------------------------------------
    // DEBUG
    //-------------------------------------------
    void PrintGraph()
    {
      NgbListIterT pFirst = oAdjList.begin();
      std::cout << "Num Vertices " << oAdjList.size() << std::endl;
      for( ; pFirst != oAdjList.end(); ++ pFirst )
      {
        VertexIterT pNgbList = pFirst->second.begin();
        for( ; pNgbList != pFirst->second.end(); ++ pNgbList )
          std::cout << "( " << pFirst->first << ", " << pNgbList->first << ")" << std::endl;
      }
    }
    
  };
  
}

#endif
