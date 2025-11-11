/////////////////////////////////////////////////////////////////
//
//  File:    SimpleGraph.tmpl.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Implementation of SimpleGraph.h
//
///////////////////////////////////////////////////////////////


//-------------------------------------------------------------------------------
//  Default Constructor
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
CSimpleUGraph< ObjT, Compare >::CSimpleUGraph()
{
  oVertexToDataMap = get( boost::vertex_name, oBoostGraph);      // get the iterator of the name list
}

//-------------------------------------------------------------------------------
//  Copy Constructor
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
CSimpleUGraph< ObjT, Compare >::CSimpleUGraph( const CSimpleUGraph & oRhs )
{
  oBoostGraph      = oRhs.oBoostGraph;
  oDataToVertexMap = oRhs.oDataToVertexMap;
  oVertexToDataMap = std::get( boost::vertex_name, oBoostGraph );
}

//-------------------------------------------------------------------------------
//  AddVertex
//
//  Note:  oNewVertex is actually the data.  Vertices are labled by data in
//         CSimpleUGraph 
// 
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
void CSimpleUGraph< ObjT, Compare >::AddVertex( const ObjT &oNewVertex )
{
  Vertex u;
  typename DataVertexMapT::iterator vMapPos;
  bool bInserted;

  if( oDataToVertexMap.find( oNewVertex ) != oDataToVertexMap.end() )
    return;
  
  // note that Vertex() is just a place holder
  std::tie( vMapPos, bInserted ) = oDataToVertexMap.insert( std::make_pair( oNewVertex, Vertex() ) );

  if ( bInserted )                                    // If data is not already added in the map
  {
    u = add_vertex( oBoostGraph );                    // add new vertex to graph
    oVertexToDataMap[ u ] = oNewVertex;               // name the vertex with the current data
    vMapPos->second = u;                              // assign the correct vertex to the map
  }
  // If data (NewVertex) already exist in the map, do nothing

}

//-------------------------------------------------------------------------------
//  AddEdge
//
//  oV1 and oV2 will be added if they do not exist.  An edge will also be added
//  between the two vertices
//
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
void CSimpleUGraph< ObjT, Compare >::AddEdge( const ObjT &oV1, const ObjT & oV2 )
{
  typename boost::graph_traits<Graph>::edge_descriptor e;
  bool bAdded;
  AddVertex( oV1 );   // should lazy these two - check to see if vertex is there first
  AddVertex( oV2 );
  Vertex u = oDataToVertexMap[ oV1 ];
  Vertex v = oDataToVertexMap[ oV2 ];
  std::tie( e, bAdded ) = add_edge( v, u, oBoostGraph );    // add edge
}

//-------------------------------------------------------------------------------
//  RemoveVertex
//
//  Edges coming in and out of vertex oV will be removed.  The vertex is then
//  removed from the internal representations as well.
//
//  Return true if the vertex to remove is found and removed.  False is returned
//  otherwise.
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
bool CSimpleUGraph< ObjT, Compare >::RemoveVertex( const ObjT & oV )
{
  typename DataVertexMapT::iterator pMapPos = oDataToVertexMap.find( oV );

  if( pMapPos == oDataToVertexMap.end() )
    return false;
  Vertex u = pMapPos->second;
  clear_vertex ( u, oBoostGraph );
  remove_vertex( u, oBoostGraph );
  return true;
}

//-------------------------------------------------------------------------------
//  RemoveEdge
//
//  The edge between oV1 and oV2 will be removed.  If the edge or any of its
//  vertices is not found, false will be returned.  True is returned upon
//  successful edge removal.
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
bool CSimpleUGraph< ObjT, Compare >::RemoveEdge( const ObjT & oV1, const ObjT & oV2 )
{
  typename DataVertexMapT::iterator pMapPos1 = oDataToVertexMap.find( oV1 );
  typename DataVertexMapT::iterator pMapPos2 = oDataToVertexMap.find( oV2 );

  if( pMapPos1 == oDataToVertexMap.end() || pMapPos2 == oDataToVertexMap.end() )
    return false;

  remove_edge( pMapPos1->second, pMapPos2->second, oBoostGraph );
}

//-------------------------------------------------------------------------------
//
//  GetNeighbors
//
//  Given a vertex, return the neighboring vertices as a pair of iterators
//  naming the first and the last element
//
//  This function is NOT efficient.  Need to change over to using iterators.
//
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
vector< ObjT > CSimpleUGraph<ObjT, Compare>::GetNeighbors( const ObjT & oVertex ) const
{
  vector<ObjT> vNeighbors(0);
  typename DataVertexMapT::const_iterator pVertexIt;
  pVertexIt = oDataToVertexMap.find( oVertex );
  if( pVertexIt == oDataToVertexMap.end() )
  {
    return vNeighbors;
  }

  AdjVertexIterT pNgb, pEnd;
  // vertices are unique - there can only be one element in the iterator
  std::tie( pNgb, pEnd ) = adjacent_vertices( pVertexIt->second, oBoostGraph );

 
  for(; pNgb != pEnd; pNgb ++ )
  {
    vNeighbors.push_back( oVertexToDataMap[ *pNgb ]  );
  }

  return vNeighbors;
}

//-------------------------------------------------------------------------------
//
//  GetConnectedComponents
//
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
void CSimpleUGraph<ObjT, Compare>::GetConnectedComponents( vector< vector< ObjT > > &vObjTComponents )
{
  // This is not pretty, but it works
  vector< VertexIndexT > vComponents( num_vertices( oBoostGraph ) );
  Int nNumComponents = connected_components( oBoostGraph, &vComponents[0] );

  vObjTComponents.clear();
  vObjTComponents.resize( nNumComponents );

  for ( Size_Type nVertexIndex = 0; nVertexIndex < vComponents.size(); nVertexIndex ++ )
  {
    VertexIndexT nComponentIndex;

    nComponentIndex = vComponents[ nVertexIndex ];
    vObjTComponents[ nComponentIndex ].push_back( oVertexToDataMap[ nVertexIndex ] ); 
  }
  
}

//-------------------------------------------------------------------------------
//  GetVertices
//  Return a list of of vertices in the graph.
//  WARNING:  Clearly this is not an efficient implementation.
// 
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
vector< ObjT > CSimpleUGraph< ObjT, Compare >::GetVertices( ) const
{
  VertexIteratorT pCur, pEnd;
  std::tie( pCur, pEnd )  = boost::vertices( oBoostGraph );
  vector< ObjT > vAllVertices(0);

  for( ; pCur != pEnd; ++ pCur )
  {
    ObjT oCurObj = oVertexToDataMap[ *pCur ];
    vAllVertices.push_back( oCurObj );
  }

  return vAllVertices;
}

//-------------------------------------------------------------------------------
// operator=
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
CSimpleUGraph< ObjT, Compare> & CSimpleUGraph< ObjT, Compare >::operator=( const CSimpleUGraph & oRhs )
{
  oBoostGraph      = oRhs.oBoostGraph;
  oDataToVertexMap = oRhs.oDataToVertexMap;
  oVertexToDataMap = std::get( boost::vertex_name, oBoostGraph );

  return *this;
}

//-------------------------------------------------------------------------------
//
//  PrintGraph  [ DEBUG ]
//
//-------------------------------------------------------------------------------
template< class ObjT, class Compare >
void CSimpleUGraph< ObjT, Compare >::PrintGraph()
{
  std::cout << "edges(g) = ";
  EdgeIteratorT  ei;
  EdgeIteratorT  ei_end;
  
  for ( std::tie(ei, ei_end) = edges( oBoostGraph ); ei != ei_end; ++ei)
  {
    std::cout << "(" << oVertexToDataMap[source( *ei, oBoostGraph )] 
              << "," << oVertexToDataMap[target( *ei, oBoostGraph )] << ") ";
  }
  std::cout << std::endl;
}
