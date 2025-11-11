//-----------------------------------------------------
//
//  Filename:  RangeSearch.tmpl.cpp
//
//  Author:    Frankie Li
//  e-mail:    sfli@cmu.edu
//
//  Implementation of range search
//
//-----------------------------------------------------


namespace RangeSearch
{
  using std::tuples::make_tuple;
  using std::shared_ptr;
  using std::multimap;
  using std::map;
  using std::make_pair;
  using std::pair;
  using std::vector;
  using std::get;
  
  //----------------------------------------------
  //  InsertElement
  //----------------------------------------------
  template< class FieldT, class DataT >
  void
  CLayeredRangeTree2D<FieldT, DataT>::InsertElement( const FieldT & oFirstCoord,
                                                     const FieldT & oSecondCoord,
                                                     const DataT  & oDataPoint )
  {
    oInputMap.push_back( InputT( oFirstCoord, oSecondCoord, oDataPoint ) );
  }

  //----------------------------------------------
  //  FindSplitNode
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  typename CLayeredRangeTree2D<CoordinateT, DataT>::TreeNodePtrT
  CLayeredRangeTree2D<CoordinateT, DataT>::FindSplitNode( const TreeNodePtrT pTree,
                                                          FieldT MinX, FieldT MaxX ) const
  {
    
    DEBUG_ASSERT( pTree.get() != NULL, "FindSplitNode: Broken precondition!\n" );
    if( IsLeaf( pTree )  )     // reached leaf node
    {
      return pTree;
    }
    else if( MinX > pTree->XValue ) 
    {
      return FindSplitNode( pTree->pRight, MinX, MaxX );
    }
    else if( MaxX <= pTree->XValue )
    {
      return FindSplitNode( pTree->pLeft, MinX, MaxX );
    }
    else                                   //  ( MinX, MaxX ]
    {
      return pTree;
    }
    
  }

  //----------------------------------------------
  //  ReportQuery
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  void
  CLayeredRangeTree2D<CoordinateT, DataT>::ReportQuery( vector<DataT> & vResult,
                                                        const AssociatedNodePtrT &pFirst,
                                                        const AssociatedNodePtrT &pEnd,
                                                        FieldT MaxY ) const
  {
    AssociatedNodePtrT pCur = pFirst;
    for(; pCur != pEnd && pCur->YValue <= MaxY; ++ pCur )
      vResult.push_back( pCur->oDataPoint );
  }

  //----------------------------------------------
  //  ReportLeaf
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  void
  CLayeredRangeTree2D<CoordinateT, DataT>::ReportLeaf( vector<DataT> & vResult,
                                                       const AssociatedNodePtrT &pFirst,
                                                       const AssociatedNodePtrT &pEnd,
                                                       FieldT MinX, FieldT MaxX,
                                                       FieldT MaxY ) const
  {
    AssociatedNodePtrT pCur = pFirst;
    for(; pCur != pEnd && pCur->YValue <= MaxY; ++ pCur )
    {
      if( pCur->XValue >= MinX && pCur->XValue <= MaxX )
        vResult.push_back( pCur->oDataPoint );
    }
  }
 
  
  //----------------------------------------------
  //  ConnectAssociatedVector
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  template< typename SelectFn > 
  void CLayeredRangeTree2D<CoordinateT, DataT>::
  ConnectAssociatedVector( const TreeNodePtrT & pChildNode,
                           const TreeNodePtrT & pParentNode,
                           SelectFn MemberSelector  ) 
  {
    AssociatedNodePtrT pFirst;
    AssociatedNodePtrT pEnd;
    
    if( pChildNode.get() != NULL )
    {
      pFirst = pChildNode->vYCoordList.begin();
      pEnd   = pChildNode->vYCoordList.end();
     
      for( AssociatedNodePtrT pCur = pParentNode->vYCoordList.begin();
           pCur != pParentNode->vYCoordList.end(); ++ pCur )
      {
        AssociatedNodePtrT pLowerBound =
          find_if( pFirst, pEnd, 
                   bind( &AssociatedNodeT::YValue, _1 ) >= pCur->YValue );

        if( pLowerBound == pEnd )
        {
          MemberSelector( *pCur ) = AssociatedNodeT::NIL;
        }
        else
        {
          MemberSelector( *pCur ) = pLowerBound - pFirst;
        }
      }
    }

  }
  
  //----------------------------------------------
  //  FillAssociatedVector
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  void CLayeredRangeTree2D<CoordinateT, DataT>::
  FillAssociatedVector( const PrimaryCoordMapIterT & pFirst,
                        const PrimaryCoordMapIterT & pEnd,
                        const TreeNodePtrT & pParentNode )
  {
    PrimaryCoordMapIterT pCur = pFirst;
    for( ; pCur != pEnd; ++ pCur )
    {
      AssociatedNodeT oNewYNode( pCur->get<0>(), pCur->get<1>(), pCur->get<2>() );
      pParentNode->vYCoordList.push_back( oNewYNode );
    }

    //  sort by y coordinate
    sort( pParentNode->vYCoordList.begin(),
          pParentNode->vYCoordList.end(),
          bind( &AssociatedNodeT::YValue, _1) < bind( &AssociatedNodeT::YValue, _2) );
  }

  //----------------------------------------------
  //  SplitPrimary
  //  Purpose: split set by median so that
  //
  //  Start with the set [ pFirst, pEnd )
  //
  //  [ pFirst, pMid) contains all points <= median
  //  [ pMid, pEnd )  contains all points > median
  //
  //  Note that this is a split along the primary coordinate
  //
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  typename CLayeredRangeTree2D<CoordinateT, DataT>::PrimaryCoordMapIterT
  CLayeredRangeTree2D<CoordinateT, DataT>::SplitPrimary( const PrimaryCoordMapIterT & pFirst,
                                                         const PrimaryCoordMapIterT & pEnd )
  {
    typedef typename PrimaryCoordMapIterT::difference_type IterDistanceT;
    IterDistanceT nElements  = pEnd - pFirst;
    
    if( nElements == 1 )
      return pEnd;

    if( nElements == 2 )
      return ( pFirst + 1 );

    PrimaryCoordMapIterT pMidpoint = pFirst + (nElements ) / 2;
    pMidpoint = upper_bound( pFirst, pEnd, pMidpoint->get<0>(),
                             FTupleLessCmp<0>() );
    return pMidpoint;
  }

  
  //----------------------------------------------
  //  BuildTree
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  typename CLayeredRangeTree2D<CoordinateT, DataT>::TreeNodePtrT
  CLayeredRangeTree2D<CoordinateT, DataT>::BuildTree( const PrimaryCoordMapIterT & pFirst,
                                                      const PrimaryCoordMapIterT & pEnd,
                                                      const PrimaryCoordMapIterT & pPointSetFirst,
                                                      const PrimaryCoordMapIterT & pPointSetEnd,
                                                      const TreeNodePtrT & pParentNode )
  {
    
    typedef typename PrimaryCoordMapIterT::difference_type IterDistanceT;
    IterDistanceT nElements  = pEnd - pFirst;
    TreeNodePtrT pNewNode( new PrimaryNodeT() );
    
    FillAssociatedVector( pPointSetFirst, pPointSetEnd, pNewNode );

    if( pFirst == pEnd )     // empty tree  ( note that X value is not defined )
    {
      pNewNode->bIsLeaf = true;
      TreeNodePtrT pEmpty;
      pNewNode->pLeft  = pEmpty;
      pNewNode->pRight = pEmpty;
      return pNewNode;
    }

    if( nElements == 1 )  // one element
    {
      pNewNode->XValue = pFirst->get<0>();
      pNewNode->bIsLeaf = false;
      PrimaryCoordMapIterT pSetMidpoint = upper_bound( pPointSetFirst, pPointSetEnd,
                                                       pFirst->get<0>(),
                                                       FTupleLessCmp<0>() );
      
      pNewNode->pLeft  = BuildTree( pFirst, pFirst,
                                    pPointSetFirst, pSetMidpoint, pNewNode );  // less than or equal to
      pNewNode->pRight = BuildTree( pFirst, pFirst,
                                    pSetMidpoint, pPointSetEnd, pNewNode );    // strictly greater than
      
    }
    else
    {
      pNewNode->bIsLeaf = false;
      PrimaryCoordMapIterT pMidpoint = SplitPrimary( pFirst, pEnd );
      PrimaryCoordMapIterT pSetMidpoint = upper_bound( pPointSetFirst, pPointSetEnd, (pMidpoint -1)->get<0>(),
                                                       FTupleLessCmp<0>() );
      
      pNewNode->XValue = ( pMidpoint - 1 )->get<0>();
      pNewNode->pLeft  = BuildTree( pFirst, pMidpoint - 1 ,
                                    pPointSetFirst, pSetMidpoint, pNewNode );  // less than or equal to
      pNewNode->pRight = BuildTree( pMidpoint, pEnd,
                                    pSetMidpoint, pPointSetEnd, pNewNode ); // strictly greater than
    }
          
    ConnectAssociatedVector( pNewNode->pLeft,
                             pNewNode,
                             bind( &AssociatedNodeT::nLeftNextOffset, _1 ) );
    ConnectAssociatedVector( pNewNode->pRight,
                             pNewNode,
                             bind( &AssociatedNodeT::nRightNextOffset, _1 ) );
        
    return pNewNode;
  }
    
  //----------------------------------------------
  //  TestTree
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  void CLayeredRangeTree2D<CoordinateT, DataT>::TestPrimaryTree( const TreeNodePtrT & pTree )
  {
    
    if( pTree.get() == NULL )
    {
      return;
    }
    else
    {
      TestPrimaryTree( pTree->pLeft );
      if( pTree->XValue == 140 ) //   if( pTree->XValue >= 139 && pTree->XValue <= 141 )
      {
        std::cout << "[" << pTree->XValue << "] " << std::endl;

        if( IsLeaf( pTree ) )
          std::cout << "140 is leaf " << std::endl;
      for( AssociatedNodePtrT pCur =  pTree->vYCoordList.begin();
           pCur != pTree->vYCoordList.end(); ++ pCur )
      {
  
        DEBUG_ASSERT( pCur->nLeftNextOffset >= -1, "undefined next \n");
        DEBUG_ASSERT( pCur->nRightNextOffset >= -1, "undefined next \n");
        std::cout << pCur->YValue << " ";

      }
      std::cout << std::endl;
      }
      TestPrimaryTree( pTree->pRight );
    }
  }

  //----------------------------------------------
  //  TreeWalkLeft
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  void
  CLayeredRangeTree2D<CoordinateT, DataT>::TreeWalkLeft( vector<DataT> & vResult,
                                                         const TreeNodePtrT & pLeft,
                                                         AssociatedNodePtrT pMin,
                                                         FieldT MinX, FieldT MaxX,
                                                         FieldT MinY, FieldT MaxY ) const
  {
    TreeNodePtrT pCur = pLeft;
    Bool bInYRange = true;
    AssociatedNodeOffsetT pMinShift = pMin->nLeftNextOffset ;
    std::tie( bInYRange, pMin ) = MoveAssociatedPtr( pCur->vYCoordList.begin(),
                                                       pMinShift );
    while ( pCur.get() != NULL && !IsLeaf( pCur ) && bInYRange )
    {
      if( MinX <= pCur->XValue )
      {
        if ( pMin->nRightNextOffset != AssociatedNodeT::NIL && pCur->pRight.get() != NULL )
        {
          AssociatedNodePtrT pReportFirst = pMin->nRightNextOffset + pCur->pRight->vYCoordList.begin();
          ReportQuery( vResult, pReportFirst, pCur->pRight->vYCoordList.end(), MaxY );
        }
        pMinShift = pMin->nLeftNextOffset;
        pCur = pCur->pLeft;
      }
      else
      {
        pMinShift = pMin->nRightNextOffset;
        pCur = pCur->pRight; 
      }
      std::tie( bInYRange, pMin ) = MoveAssociatedPtr( pCur->vYCoordList.begin(),
                                                         pMinShift );
    }


    if( pCur.get() != NULL && bInYRange )
      ReportLeaf( vResult, pMin, pCur->vYCoordList.end(), MinX, MaxX, MaxY );
  }

  //----------------------------------------------
  //  TreeWalkRight
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  void
  CLayeredRangeTree2D<CoordinateT, DataT>::TreeWalkRight( vector<DataT> & vResult,
                                                          const TreeNodePtrT & pRight,
                                                          AssociatedNodePtrT pMin,
                                                          FieldT MinX, FieldT MaxX,
                                                          FieldT MinY, FieldT MaxY ) const
  {

    TreeNodePtrT pCur = pRight;
    DEBUG_ASSERT( !IsLeaf( pCur ), "TreeWalk unexpectedly applied to leaf\n");
    Bool bInYRange = true;
    AssociatedNodeOffsetT pMinShift = pMin->nRightNextOffset;
    std::tie( bInYRange, pMin ) = MoveAssociatedPtr( pCur->vYCoordList.begin(),
                                                       pMinShift );
    Bool bMovedLeft;
    while (  pCur.get() != NULL && !IsLeaf( pCur ) && bInYRange )
    {
      if( MaxX > pCur->XValue )
      {
        if ( pMin->nLeftNextOffset != AssociatedNodeT::NIL  && pCur->pLeft.get() != NULL )
        {
          AssociatedNodePtrT pReportFirst = pMin->nLeftNextOffset + pCur->pLeft->vYCoordList.begin();
          ReportQuery( vResult, pReportFirst, pCur->pLeft->vYCoordList.end(), MaxY );
        }
        
        pMinShift = pMin->nRightNextOffset;
        pCur = pCur->pRight;
        bMovedLeft = false;
      }
      else
      {
        pMinShift = pMin->nLeftNextOffset;
        pCur = pCur->pLeft;
        bMovedLeft = true;
      }
      std::tie( bInYRange, pMin ) = MoveAssociatedPtr( pCur->vYCoordList.begin(),
                                                         pMinShift );
    }

    if( pCur.get() != NULL && bInYRange )
      ReportLeaf( vResult, pMin, pCur->vYCoordList.end(), MinX, MaxX, MaxY );
    
  }
  
  //----------------------------------------------
  //  Query
  //----------------------------------------------
  template< class CoordinateT, class DataT >
  vector<DataT>
  CLayeredRangeTree2D<CoordinateT, DataT>::Query( FieldT MinX, FieldT MaxX,
                                                  FieldT MinY, FieldT MaxY ) const
  {
    vector<DataT> vData;
    if( pRoot.get() == NULL )
      return vData;
    
    TreeNodePtrT pSplitNode = FindSplitNode( pRoot, MinX, MaxX );

    if( pSplitNode.get() == NULL )
      return vData;
    
    AssociatedNodePtrT pMin   = find_if( pSplitNode->vYCoordList.begin(),
                                         pSplitNode->vYCoordList.end(),
                                         bind( &AssociatedNodeT::YValue, _1 ) >= MinY );
    if( IsLeaf( pSplitNode ) )
    {
      ReportLeaf( vData, pMin, pSplitNode->vYCoordList.end(), MinX, MaxX, MaxY );
      return vData;
    }
    else  // tree walk
    {
      if( pMin == pSplitNode->vYCoordList.end() )
        return vData;
      
      if( pMin->nLeftNextOffset != AssociatedNodeT::NIL )
      {
        TreeWalkLeft ( vData, pSplitNode->pLeft, pMin, MinX, MaxX, MinY, MaxY );
      }

      if( pMin->nRightNextOffset != AssociatedNodeT::NIL )
      {
        TreeWalkRight ( vData, pSplitNode->pRight, pMin, MinX, MaxX, MinY, MaxY );
      }
      
      
    }
    return vData;
  }
};
