//-----------------------------------------------------
//
//  Filename:  RangeSearch.h
//
//  Author:    Frankie Li
//  e-mail:    sfli@cmu.edu
//
//  Implementation of range search
//
//-----------------------------------------------------

#ifndef RANGE_SEARCH_H_
#define RANGE_SEARCH_H_
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <map>
#include <vector>
#include <boost/function.hpp>

#include "Types.h"
#include "Debug.h"


namespace RangeSearch
{

  using boost::tuple;
  using boost::shared_ptr;
  using std::map;
  using std::vector;
  using namespace boost::lambda;
  using boost::get;
  
  template< typename CoordinateT, typename DataT >
  struct SYNode;    // Forward declaration

  template<typename CoordinateT, typename DataT>
  struct SXNode
  {
    SXNode(): pLeft(), pRight() {};
    CoordinateT XValue;
    vector< SYNode<CoordinateT, DataT> > vYCoordList;
    shared_ptr< SXNode > pLeft;
    shared_ptr< SXNode > pRight;
    Bool bIsLeaf;
  };

  template< typename CoordinateT, typename DataT >
  struct SYNode
  {
    typedef typename vector<SYNode>::iterator NextPtrT;
    typedef typename vector<SYNode>::difference_type PtrOffsetT;
    static const PtrOffsetT NIL = -1;    

    SYNode( const CoordinateT & _x,
            const CoordinateT & _y,
            const DataT & _d ): XValue( _x ),
                                YValue( _y ),
                                oDataPoint( _d ),
                                nLeftNextOffset ( NIL ),
                                nRightNextOffset( NIL )  { };
    CoordinateT XValue;
    CoordinateT YValue;
    DataT oDataPoint;
    PtrOffsetT nLeftNextOffset;
    PtrOffsetT nRightNextOffset;
  };

  //-----------------------------------------------------
  //
  //  CLayeredRangeTree2D
  //
  //  Purpose:  A simple 2D range tree for orthogonal range
  //            search.  Range is defined over CoordinateT.
  //            (i.e., Float, Int )  Note that this is a
  //            static data structure, implemented with
  //            fractional cascading.
  //-----------------------------------------------------
  template< class CoordinateT, class DataT >
  class CLayeredRangeTree2D
  {
    
  public:

    //----------------------------------------------
    //  Typedefs only  ( please don't shorten AssociatedNode to AssNode ) 
    //----------------------------------------------
    typedef CoordinateT FieldT;

    typedef SXNode<FieldT, DataT> PrimaryNodeT;
    typedef SYNode<FieldT, DataT> AssociatedNodeT;
    typedef typename boost::tuple< FieldT, FieldT, DataT > InputT; 
    typedef typename boost::shared_ptr< SXNode< FieldT, DataT > > TreeNodePtrT;
    typedef typename std::vector< InputT > PrimaryCoordMapT;
    typedef typename PrimaryCoordMapT::iterator PrimaryCoordMapIterT;
    typedef typename AssociatedNodeT::NextPtrT AssociatedNodePtrT;
    typedef typename AssociatedNodeT::PtrOffsetT AssociatedNodeOffsetT;
    
    //----------------------------------------------
    //  First coordinate comparison
    //----------------------------------------------
    template < int N >
    struct FTupleLessCmp
    {
      bool operator()( const InputT  &t1,  const InputT  &t2 ) const
      {
        return ( t1.get<N>() < t2.get<N>() );
      };
    };
    
  private:

    //----------------------------------------------
    //  ConnectAssociatedVector
    //  Purpose:  Connect the pointers for each of the nodes in the
    //            associated vectors.  This is the construction
    //            of the fractional cascading structure.
    //
    //  Parameters:  pFirst -  first element of the new associated vector
    //                         of the child node
    //               pEnd   -  past last element
    //  In another words, the child associated vector has elements,
    //  [ pFirst, pEnd ).
    //
    //  Note that pLeftFirst == pLeftEnd will result in pLeftEnd being
    //  assigned to each of pLeftNext of pParentNode's associated vector.
    //  (i.e., it satisfies the condition of nodes to be ignored.)
    //----------------------------------------------
    template< typename SelectFn >
    void ConnectAssociatedVector( const TreeNodePtrT & pChildNode,
                                  const TreeNodePtrT & pParentNode,
                                  SelectFn MemberSelector );
    
    //----------------------------------------------
    //  FillAssociatedVectors
    //
    //  Purpose:  Fillts the tree node's associated vector 
    //
    //----------------------------------------------
    void FillAssociatedVector( const PrimaryCoordMapIterT & pFirst,
                               const PrimaryCoordMapIterT & pEnd,
                               const TreeNodePtrT & pParentNode );

    //----------------------------------------------
    //  FindSplitNode
    //----------------------------------------------
    TreeNodePtrT FindSplitNode( const TreeNodePtrT pTree,
                                FieldT MinX, FieldT MaxX ) const;

    //----------------------------------------------
    //  SplitPrimary
    //----------------------------------------------
    PrimaryCoordMapIterT SplitPrimary( const PrimaryCoordMapIterT & pFirst,
                                       const PrimaryCoordMapIterT & pEnd );
    
    //----------------------------------------------
    //  IsLeaf
    //----------------------------------------------
    Bool IsLeaf( const TreeNodePtrT &pTree ) const
    {
      return ( pTree->pLeft.get() == NULL && pTree->pRight.get() == NULL );
    };

    //----------------------------------------------
    //  ReportQuery
    //----------------------------------------------
    void ReportQuery( vector<DataT> & oResult,
                      const AssociatedNodePtrT &pFirst,
                      const AssociatedNodePtrT &pEnd,
                      FieldT MaxY ) const;
    
    void ReportLeaf( vector<DataT> & vResult,
                     const AssociatedNodePtrT &pFirst,
                     const AssociatedNodePtrT &pEnd,
                     FieldT MinX, FieldT MaxX,
                     FieldT MaxY ) const;
      
    //----------------------------------------------
    //  TreeWalkLeft
    //----------------------------------------------
    void TreeWalkLeft( vector<DataT> & vResult,
                       const TreeNodePtrT & pSplitNode,
                       AssociatedNodePtrT pMin,
                       FieldT MinX, FieldT MaxX,
                       FieldT MinY, FieldT MaxY ) const;

    //----------------------------------------------
    //  TreeWalkRight
    //----------------------------------------------
    void TreeWalkRight( vector<DataT> & vResult,
                        const TreeNodePtrT & pSplitNode,
                        AssociatedNodePtrT pMin,
                        FieldT MinX, FieldT MaxX,
                        FieldT MinY, FieldT MaxY ) const;

    //----------------------------------------------
    //  MoveAssociatedPtr
    //----------------------------------------------
    std::pair< bool, AssociatedNodePtrT >
    MoveAssociatedPtr( AssociatedNodePtrT    pNode,
                       AssociatedNodeOffsetT nOffset) const
    {
      std::pair<bool, AssociatedNodePtrT> oRes;
      if( nOffset == AssociatedNodeT::NIL )
      {
        oRes.first  = false;
        oRes.second = pNode;
      }
      else
      {
        oRes.first  = true;
        oRes.second = pNode + nOffset;
      }
      return oRes;
    }
    
    //----------------------------------------------
    //  BuildTree
    //
    //  Purpose:  Construct the binary tree and the associated structure
    //            of fractional cascading.
    //
    //  Precondition:  pFirst is an iterator pointing to a list of
    //                 elements sorted by the primary coordinate.
    //                 The list has a total of pEnd - pFirst points.
    //                 (i.e., the list goes from [pFirst, pEnd ) ) 
    //
    //----------------------------------------------
    TreeNodePtrT BuildTree( const PrimaryCoordMapIterT & pFirst,
                            const PrimaryCoordMapIterT & pEnd,
                            const PrimaryCoordMapIterT & pPointSetFirst,
                            const PrimaryCoordMapIterT & pPointSetEnd,
                            const TreeNodePtrT & pParentNode );
    
    //----------------------------------------------
    //  Data Members
    //----------------------------------------------
    TreeNodePtrT pRoot;
    PrimaryCoordMapT oInputMap;
    
  public:

    //----------------------------------------------
    //  Default constructor
    //----------------------------------------------
    CLayeredRangeTree2D(): pRoot(){};
    
    //----------------------------------------------
    //  InsertElement
    //----------------------------------------------
    void InsertElement( const FieldT & oFirstCoord,
                        const FieldT & oSecondCoord,
                        const DataT  & oDataPoint );
    
    //----------------------------------------------
    //  BuildTree
    //----------------------------------------------
    void BuildTree( )
    {
      sort( oInputMap.begin(), oInputMap.end(), FTupleLessCmp<0>() );
      pRoot = BuildTree( oInputMap.begin(), oInputMap.end(),
                         oInputMap.begin(), oInputMap.end(),pRoot );
    };
    
    //----------------------------------------------
    //  Query
    //----------------------------------------------
    template< class ProcessingFn >
    Bool Query( ProcessingFn & Fn ) const;    

    //----------------------------------------------
    //  Query
    //----------------------------------------------
    vector<DataT> Query( FieldT MinX, FieldT MaxX,
                         FieldT MinY, FieldT MaxY ) const;
    
    //----------------------------------------------
    //  TestTree
    //----------------------------------------------
    void TestPrimaryTree( const TreeNodePtrT & pTree );
    void TestPrimaryTree( )
    {
      TestPrimaryTree( pRoot );
    }
  };

};

#include "RangeSearch.tmpl.cpp"

#endif
