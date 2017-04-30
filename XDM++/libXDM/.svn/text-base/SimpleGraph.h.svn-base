/////////////////////////////////////////////////////////////////
//
//  File:    SimpleGraph.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Simple wrapper around Boost Graph Library (BGL) that builds
//  a graph out of any iterator.  This is to reduce the learning curve
//  of XDM++ library.  The detail implementation of the graph as well
//  as the mapping of vertices to BGL's internal vertex representation
//  are completely hidden from the user.  Note that knowledge of implemented
//  algorithms are still very important, as they affects the runtime performance
//  significantly.
//
//  Requirements:  Boost Graph Library
//
////////////////////////////////////////////////////////////////
#ifndef _SIMPLE_GRAPH_H_
#define _SIMPLE_GRAPH_H_

#include "Debug.h"
#include <functional>
#include <algorithm>
#include <vector>
#include "Types.h"

#define BOOST_NO_HASH   // for boost graph
// Boost Graphics Library
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using std::vector;
using std::less;

//-------------------------------------------------------------------------------
//
//   Template class CSimpleUGraph -- Simple Undirected Graph
//
//   Purpose:  This class generates a graph
//             using functions from Boost Graph Library.  The lookup table
//             used here to map between vertex and vertex type is std::map
//             The graph used here is an adjacency list.  Data is held at each
//             vertex, and they are stored at vertex property of the vertex list.
//
//
//   Requirement:  1.  ObjT has a well defined operator< 
//                 2.  ObjT has a well defined copy constructor 
//                 
//
//
//
//-------------------------------------------------------------------------------
template< class ObjT, class Compare = less<ObjT> >
class CSimpleUGraph
{
  //-----------------------------
  // Type definitions only  for graph algorithms
  //------------------------------
public:
  typedef ObjT   ElementType;
  typedef ObjT & ObjT_Ref; 
  typedef typename boost::property < boost::vertex_name_t, ObjT > VertexPropertyMapT;
  typedef typename boost::adjacency_list < boost::setS, boost::vecS, boost::undirectedS, VertexPropertyMapT > Graph;
  typedef Graph &Graph_Ref;
  typedef typename boost::property_map < Graph, boost::vertex_name_t >::type VertexNameT;
  typedef typename boost::graph_traits<Graph>::edge_iterator EdgeIteratorT;
  typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIteratorT;
  typedef typename boost::graph_traits<Graph>::vertices_size_type VertexIndexT;
  
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef std::map< ObjT, Vertex, Compare > DataVertexMapT;

  // Iterators
  typedef typename boost::graph_traits< Graph >::adjacency_iterator AdjVertexIterT;
  
protected:

  //------------------------------
  //  Note about Maps:
  //
  //  Map gaurantees that the operator[]
  //  is implemented.  Therefore oVertexToDataMap[v]
  //  returns the ObjT object (data) that corrosponds to the
  //  vertex, v
  //------------------------------
  
  VertexNameT     oVertexToDataMap;
  Graph           oBoostGraph;      // BGL graph
  DataVertexMapT  oDataToVertexMap;

public:
  
  //------------------------------
  //  Default Constructor
  //------------------------------
  CSimpleUGraph();

  //------------------------------
  //  Copy Constructor
  //------------------------------
  CSimpleUGraph( const CSimpleUGraph & oRhs );
  
  //------------------------------
  //  Accessors
  //------------------------------

  //------------------------------
  //
  //  GetNeighbors
  //
  //  Given a vertex, return the neighboring vertices as a pair of iterators
  //  naming the first and the last element
  //------------------------------
  vector< ObjT > GetNeighbors( const ObjT & oVertex ) const;
  
  //------------------------------
  //
  //  AddVertex( ObjT oNewVertex )
  //
  //  Note:  oNewVertex is actually the data.  Vertices are labled by data in
  //         CSimpleUGraph 
  //------------------------------
  void AddVertex( const ObjT &oNewVertex );

  //------------------------------
  //
  //  AddEdge
  //
  //  oV1 and oV2 will be added if they do not exist.  An edge will also be added
  //  between the two vertices
  //
  //------------------------------
  void AddEdge( const ObjT &oV1, const ObjT & oV2 );

  //------------------------------
  //  RemoveVertex
  //------------------------------
  bool RemoveVertex( const ObjT & oV1 );
  
  //------------------------------
  //  RemoveEdge
  // 
  //  The edge (oV1, oV2) will be removed.  True is returned
  //  upon successful removal.  False otherwise.
  //------------------------------
  bool RemoveEdge( const ObjT & oV1, const ObjT & oV2 );
  
  //------------------------------
  //  Algorithms
  //------------------------------

  //------------------------------
  //
  //  Get connected components of the graph
  //
  //------------------------------
  void GetConnectedComponents( vector< vector< ObjT > > &vObjTComponents );

  //------------------------------
  //  Accessor
  //------------------------------

  //------------------------------
  //  GetVertices
  //  Return a list of of vertices in the graph.
  //  WARNING:  Clearly this is not an efficient implementation.
  // 
  //------------------------------
  vector< ObjT > GetVertices() const;
  
  //------------------------------
  //  GetBoostGraph
  //  Somethings are just not covered by SimpleGraph
  //------------------------------
  Graph_Ref GetBoostGraph() const
  {
    return oBoostGraph;
  }

  //------------------------------
  // operator=
  //------------------------------
  CSimpleUGraph& operator=( const CSimpleUGraph & oRHS );
  
  //------------------------------
  //  [DEBUG]
  //  PrintGraph
  //
  //------------------------------
  void PrintGraph();
  
  //------------------------------
  //  GetDataVertexMap
  //------------------------------
  DataVertexMapT GetDataVertexMap() const
  {
    return oDataToVertexMap;
  };
  
};


#include "SimpleGraph.tmpl.cpp"

#endif
