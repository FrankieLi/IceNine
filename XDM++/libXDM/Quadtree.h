#ifndef _QUADTREE_H_
#define _QUADTREE_H_

#include <vector>
#include "Debug.h"
#include "BBox.h"
#include <iostream>

namespace RangeSearch
{

using std::vector;

template <class NodeData> 
class CQuadNode
{
public:

  CQuadNode *opChildren[4];
  vector<NodeData> vData;
  PBRMath::BBox2D Bound;  
  
  //----------------------
  //  Default constructor
  //----------------------
  CQuadNode():vData()
  {
    for( Int i = 0; i < 4; i++ )
      opChildren[i] = NULL;
  }
  
  //----------------------
  //  Copy Constructor
  //----------------------
  CQuadNode( const CQuadNode & oRHS ):vData( oRHS.vData ), Bound( oRHS.Bound )
  {
    for( Int i = 0; i < 4; i++ )
    {
      if( oRHS.opChildren[i] != NULL )
        opChildren[i] = new CQuadNode<NodeData>( *oRHS.opChildren[i] );
      else
        opChildren[i] = NULL;
    }
  }

  //----------------------
  //  operator=
  //----------------------
  CQuadNode & operator=( const CQuadNode & oRHS )
  {
    vData = oRHS.vData;
    Bound = oRHS.Bound;
    for( Int i = 0; i < 4; i++ )
    {
      if( oRHS.opChildren[i] != NULL )
        opChildren[i] = new CQuadNode<NodeData>( *oRHS.opChildren[i] );
      else
        opChildren[i] = NULL;
    }
    return *this;
  }
  
  //----------------------
  //  Destructor
  //----------------------
  ~CQuadNode()
  {
    for( Int i = 0; i < 4; i++ )
      delete opChildren[i];
  }
};


//----------------------------------------------------------------------
//
//  LookupProc is a functor that takes two parameters
//  and does the intended operations.  i.e., 
//  AddNode(P, node->data[i]) would add the data point
//  node->data[i] when ever it is called, thus associating
//  it to be a data point that's within the same space as
//  the leaf node that occupies point P
//----------------------------------------------------------------------
template < class NodeData > 
class CQuadtree
{

private:
  
  struct SAccumulator
  {
    vector< NodeData > oRes;
    void operator()( const PBRMath::BBox2D & oBox, const NodeData & oObj  ) { oRes.push_back( oObj ); }
  };

protected:
  PBRMath::BBox2D oTreeBound;
  Int nMaxDepth;
  CQuadNode< NodeData > oRoot;
  
  //--------------------------
  //  myAdd -- helper for adding element
  //--------------------------
  void myAdd( CQuadNode<NodeData> *node, const PBRMath::BBox2D &nodeBound,
              const NodeData &data, const PBRMath::BBox2D &dataBound, 
              Float diag2, Int depth = 0 );
  
  //--------------------------
  //  myLookup -- helper for point lookup
  //--------------------------
  template< typename LookupProc >
  void myLookup( const CQuadNode<NodeData> *node, const PBRMath::BBox2D &nodeBound,
                 const PBRMath::Point &P, LookupProc &process ) const;

  template< typename EqualCmp >
  void myDelete( CQuadNode<NodeData> *node, 
                 const PBRMath::BBox2D &nodeBound,
                 const PBRMath::Point &p,
                 const EqualCmp & FEqualCmp );
  
  //--------------------------
  //  myRangeSearch -- helper for range lookup
  //
  //  Requirements:  LookupProc, RetType
  //
  //    Functions defined:
  //     void LookupProc( PBRMath::BBox2D, NodeData )
  //
  //--------------------------
  template< typename LookupProc >
  void myRangeSearch( const CQuadNode<NodeData> *node, const PBRMath::BBox2D &nodeBound,
                      const PBRMath::BBox2D & oBoxRange, LookupProc &process ) const;

  //--------------------------
  //  myPrintTree
  //--------------------------
  void myPrintTree( const CQuadNode<NodeData> *node, Int currentDepth ) const;

  //--------------------------
  //
  //--------------------------
  
public:

  //------------------------------------------------------------------------
  //                    PUBLIC FUNCTIONS
  //------------------------------------------------------------------------
  
  //--------------------------
  //  CQuadtree
  //-------------------------
  CQuadtree( const PBRMath::BBox2D &b, Int md )
    : oTreeBound(b), nMaxDepth( md ), oRoot()
  {
    oRoot.Bound = b;  // should check bound
  }
  
  //--------------------------
  //  CQuadtree
  //-------------------------
  CQuadtree()
    :oTreeBound( PBRMath::Point( MAX_FLOAT, MAX_FLOAT ),
                 PBRMath::Point( MIN_FLOAT, MIN_FLOAT ) ),
     nMaxDepth( 0 ), oRoot()
  { }

  //--------------------------
  //  SetMaxDepth
  //--------------------------
  void SetMaxDepth( Int nNewMaxDepth )
  {
    nMaxDepth = nNewMaxDepth;
  };

  //--------------------------
  //  SetTreeBound
  //--------------------------
  void SetTreeBound( const PBRMath::BBox2D & oNewBound )
  {
    oTreeBound = oNewBound;
    oRoot.Bound = oNewBound;
  }

  //--------------------------
  //  Add -- insert element
  //--------------------------
  void Add( const NodeData &data, const PBRMath::BBox2D &dataBound )
  {
    myAdd( &oRoot, oTreeBound, data, dataBound, 
           DistanceSquared(dataBound.pMin, dataBound.pMax) ); 
  }

  //--------------------------
  //  Lookup - point location
  //--------------------------
  template< typename LookupProc >
  void Lookup( const PBRMath::Point &P, LookupProc &process ) const
  {
    myLookup( &oRoot, oTreeBound, P, process ); 
  }

  template< typename DeleteCmp >
  void Delete( const PBRMath::Point & p, DeleteCmp FToDelete )
  {
    if( oRoot.Bound.Inside( p ) )
      myDelete( &oRoot, oTreeBound, p, FToDelete ); 
  }
  
  //--------------------------
  //  RangeSearch - point location
  //--------------------------
  template< typename LookupProc >
  void RangeSearch( const PBRMath::BBox2D &oBoxRange, LookupProc &process ) const
  {
    myRangeSearch( &oRoot, oTreeBound, oBoxRange, process ); 
  }

  //--------------------------
  //  Non-unique points do occur
  //--------------------------
  vector<NodeData> GetAllData() const
  {
    SAccumulator oProcess;
    RangeSearch( oTreeBound, oProcess );
    return oProcess.oRes;
  }
  
  //--------------------------
  //  PrintTree
  //--------------------------
  void PrintTree() const
  {
    myPrintTree(&oRoot, 0);
  }
  
};


//////////////////////////////////////////////////////////////
//
//  Quad tree nodes look like this:
//
//
//    ___________
//    |    |    |
//    | 2  |  3 |
//    |____|____|
//    |    |    |
//    | 0  |  1 |             
//    |____|____|                       
//                    yx                             false    true
//    in binary:  0 - 00    least significant big is x_lower, x_upper         
//                1 - 01    most significant bit is  y_lower, y_upper
//                2 - 10
//                3 - 11
//
//    The above fact will be taken advantage of in the add function
// 
// Private Function:  myAdd
// 
// Action:  Add a new node to the tree
//
//
//////////////////////////////////////////////////////////////

template <class NodeData>
void	CQuadtree<NodeData>::myAdd(CQuadNode<NodeData> *node,
                                               const PBRMath::BBox2D &nodeBound,
                                               const NodeData &data,
                                               const PBRMath::BBox2D &dataBound, 
                                               Float itemDiag2, Int depth)
{
  //-------------------
  // leaf node condition - add to current node
  //
  //
  // Slight subtlety here:  One might worry that item's diagonal
  // will be larger than the current cell.  This will not be a problem
  // because the neighboring cells will also contain this item.
  //-------------------
  if(   ( depth == nMaxDepth )
     || DistanceSquared(nodeBound.pMin, nodeBound.pMax) < itemDiag2 )
  {
    node->vData.push_back(data);
    return;
  }

  // find the child to add this data to
	
  PBRMath::Point pMid = 0.5 * nodeBound.pMin  +  nodeBound.pMax * 0.5; 
  Bool overlap[4];
	
  // node 0 and 2 are the ones with the lower x values => x in [x_min, x_mid]
  // node 1 and 3 are the ones with the upper x values =? x in [x_mid, x_max]
  overlap[0] = overlap[2] = dataBound.pMin.x < pMid.x;
  overlap[1] = overlap[3] = dataBound.pMax.x >= pMid.x;
	
  // y lower
  overlap[0] &= (dataBound.pMin.y < pMid.y);  
  overlap[1] &= (dataBound.pMin.y < pMid.y);

  // y upper
  overlap[2] &= (dataBound.pMax.y >= pMid.y);
  overlap[3] &= (dataBound.pMax.y >= pMid.y);

  for( Int child = 0; child < 4; child ++ )
  {
    if( overlap[child] )
    {
      // create new nodes if it isn't there already
      if( node->opChildren[child] == NULL )
      {
        node->opChildren[child] = new CQuadNode<NodeData>;
      }
      
      // recursively search, create the bound of the children
      PBRMath::BBox2D childBound;
      
      // x direction  ture => upper, false => lower
      childBound.pMin.x = ( 1 & child) ? pMid.x : nodeBound.pMin.x;
      childBound.pMax.x = ( 1 & child) ? nodeBound.pMax.x : pMid.x;
      
      // y direction: true => lower, false => upper
      childBound.pMin.y = ( 2 & child) ? pMid.y : nodeBound.pMin.y;
      childBound.pMax.y = ( 2 & child) ? nodeBound.pMax.y : pMid.y;

      QUADTREE_DEBUG( node->opChildren[child]->Bound = childBound );
      
      myAdd( node->opChildren[child], childBound, data, dataBound, 
             itemDiag2, depth +1 );
    }
  }
}

//------------------------------------------------------------------------------
//
//  myLookup
//  Generic Point Lookup
//  
//------------------------------------------------------------------------------
template < class NodeData >
template < typename LookupProc >
void CQuadtree< NodeData >::myLookup( const CQuadNode<NodeData> *node, 
                                      const PBRMath::BBox2D &nodeBound,
                                      const PBRMath::Point &p, 
                                      LookupProc &process ) const
{
  for( Size_Type i = 0; i < node->vData.size(); ++i )
    process( p, node->vData[i] );   // not checking local bbounding box for performance reasons
	
  // Determine which quad child node _p_ is inside
  PBRMath::Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;

  Int child = (p.x > pMid.x ? 1 : 0) + (p.y > pMid.y ? 2 : 0);


  if (node->opChildren[child]) {
    // Compute _childBound_ for octree child _child_
    PBRMath::BBox2D childBound;

    // x direction  ture => upper, false => lower
    childBound.pMin.x = ( 1 & child) ? pMid.x : nodeBound.pMin.x;
    childBound.pMax.x = ( 1 & child) ? nodeBound.pMax.x : pMid.x;
		
    // y direction: true => lower, false => upper
    childBound.pMin.y = ( 2 & child) ? pMid.y : nodeBound.pMin.y;
    childBound.pMax.y = ( 2 & child) ? nodeBound.pMax.y : pMid.y;
		
    myLookup( node->opChildren[child], childBound, p, process);
  }
}


//------------------------------------------------------------------------------
//
//  myDelete  -- delete all object that intersects the point p
//  
//------------------------------------------------------------------------------
template < class NodeData >
template < class DeleteCmp >
void CQuadtree< NodeData >::myDelete( CQuadNode<NodeData> *node, 
                                      const PBRMath::BBox2D &nodeBound,
                                      const PBRMath::Point &p,
                                      const DeleteCmp & FToDelete )
{
  typename vector<NodeData>::iterator pToRemove;
  //pToRemove
  std::remove_if( node->vData.begin(), node->vData.end(), FToDelete );
  node->vData.erase( pToRemove, node->vData.end() );
  
  // Determine which quad child node _p_ is inside
  PBRMath::Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;

  Int child = (p.x > pMid.x ? 1 : 0) + (p.y > pMid.y ? 2 : 0);
  
  if (node->opChildren[child]) {
    // Compute _childBound_ for octree child _child_
    PBRMath::BBox2D childBound;

    // x direction  ture => upper, false => lower
    childBound.pMin.x = ( 1 & child) ? pMid.x : nodeBound.pMin.x;
    childBound.pMax.x = ( 1 & child) ? nodeBound.pMax.x : pMid.x;
		
    // y direction: true => lower, false => upper
    childBound.pMin.y = ( 2 & child) ? pMid.y : nodeBound.pMin.y;
    childBound.pMax.y = ( 2 & child) ? nodeBound.pMax.y : pMid.y;
		
    myDelete( node->opChildren[child], childBound, p, FToDelete );
  }
}


//------------------------------------------------------------------------------
//  myRangeSearch
//  Generic Range Lookup
//------------------------------------------------------------------------------
template <class NodeData>
template< typename LookupProc >
void CQuadtree<NodeData>::myRangeSearch( const CQuadNode<NodeData> *node, 
                                         const PBRMath::BBox2D &nodeBound,
                                         const PBRMath::BBox2D &oSearchRange, 
                                         LookupProc &process ) const
{
  for( Size_Type i = 0; i < node->vData.size(); ++i )
    process( oSearchRange, node->vData[i] );

  Bool overlap[4];

  PBRMath::Point pMid = 0.5 * nodeBound.pMin  +  nodeBound.pMax * 0.5;
  // node 0 and 2 are the ones with the lower x values => x in [x_min, x_mid]
  // node 1 and 3 are the ones with the upper x values =? x in [x_mid, x_max]
  overlap[0] = overlap[2] = oSearchRange.pMin.x < pMid.x;
  overlap[1] = overlap[3] = oSearchRange.pMax.x >= pMid.x;
	
  // y lower
  overlap[0] &= ( oSearchRange.pMin.y < pMid.y );  
  overlap[1] &= ( oSearchRange.pMin.y < pMid.y );

  // y upper
  overlap[2] &= ( oSearchRange.pMax.y >= pMid.y );
  overlap[3] &= ( oSearchRange.pMax.y >= pMid.y );
  
  for( Int nChild = 0; nChild < 4; ++ nChild ) 
  {
    if ( overlap[ nChild ] && node->opChildren[ nChild ] )
    {
      PBRMath::BBox2D childBound;
      
      // x direction  ture => upper, false => lower
      childBound.pMin.x = ( 1 & nChild ) ? pMid.x : nodeBound.pMin.x;
      childBound.pMax.x = ( 1 & nChild ) ? nodeBound.pMax.x : pMid.x;
      
      // y direction: true => lower, false => upper
      childBound.pMin.y = ( 2 & nChild ) ? pMid.y : nodeBound.pMin.y;
      childBound.pMax.y = ( 2 & nChild ) ? nodeBound.pMax.y : pMid.y;
      
      myRangeSearch( node->opChildren[ nChild ], childBound, oSearchRange, process );
    }
  }
}
//------------------------------------------------------------------------------
//
//  myPrintTree
//  Generic Lookup
//  
//------------------------------------------------------------------------------
template <class NodeData>
void CQuadtree< NodeData >::myPrintTree( const CQuadNode<NodeData> *node,
                                         Int currentDepth ) const
{
  // print myself
	
  if(nMaxDepth == currentDepth){
    return;
  }

  if( node->vData.size() > 0 )
  {
    std::cout << currentDepth << ", "  
              << node->Bound.pMin.x << ", " << node->Bound.pMin.y << ", " 
              << node->Bound.pMax.x << ", " << node->Bound.pMax.y << ", "
              << node->vData[0]->oPos << std::endl;
  }
  
  // print children
  for( Int i = 0; i < 4; i ++ ){
    if(node->opChildren[i] != NULL)
      myPrintTree(node->opChildren[i], currentDepth +1);
  }
}



};//  end of RangeSearch namespace

#endif
