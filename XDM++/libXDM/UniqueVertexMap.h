/////////////////////////////////////////////////////////////////
//
//  File:    UniqueVertexMap
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Purpose: Specialized class to identify unique vertices
//
/////////////////////////////////////////////////////////////////

#ifndef _UNIQUE_VERTEX_MAP
#define _UNIQUE_VERTEX_MAP

#include "Types.h"
#include "3dMath.h"
#include "Quadtree.h"
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace XDMUtility
{
  //---------------------------------------------------
  //  The struct of simple vertex
  //---------------------------------------------------
  class SVertex
  {
  public:

    SVector3 oPos;
    Int nID;

    SVertex():oPos(), nID(){}
    SVertex( const SVector3 & v, Int n ): oPos( v ), nID( n ){}

    // Less than operator
    Bool operator<( const SVertex & oRHS ) const
    {
      return nID < oRHS.nID;
    }
  
  };


  
  typedef boost::shared_ptr< SVertex > CSVertexPtr;
  
  
  //---------------------------------------------------
  //  A simple unique vertex locator 
  //---------------------------------------------------
  template< class CVertexPtr = CSVertexPtr >
  class CUniqueVertexMap: public RangeSearch::CQuadtree< CVertexPtr > 
  {
    
  private:

    //---------------------------------------------------
    //   Default function for quadtree operations
    //
    //---------------------------------------------------
    struct SVertexOverlapLookup
    {
      Int nRepeat;    // multiple insertion detection.  i.e., detect when fTolerance is incorrect
      Bool bFound;
      Float fTolerance;
      CVertexPtr pFoundVertex;
    
      SVertexOverlapLookup( Float _fTol ): nRepeat(0), bFound( false ),
                                           fTolerance( _fTol ),
                                           pFoundVertex() {};
      
      void operator()( const PBRMath::Point & p, CVertexPtr pVertex )
      {
        if( ! bFound )  // if item already found from previous lookup
        {
          Float fDx = p.x - pVertex->oPos.m_fX;
          Float fDy = p.y - pVertex->oPos.m_fY;
          if( fabs( fDx ) < fTolerance && fabs( fDy ) < fTolerance )
          {
            DEBUG_ASSERT( nRepeat == 0, "ERROR!  Tolerance set too high for UniqueVertexMap.\n" );
            pFoundVertex = pVertex;
            bFound = true;
            nRepeat = 1;
          }
          else
          {
            bFound = false;
          }
        }
      }
    };
    
    //---------------------------------------------------
    //  Multi-vertex lookup
    //---------------------------------------------------
    struct SMultiVertexLookup
    {
      vector<CVertexPtr> oPointList;
      void operator()( const PBRMath::BBox2D & oSearchBox, CVertexPtr pVertex )
      {
        Point p( pVertex->oPos.m_fX, pVertex->oPos.m_fY );
        if( oSearchBox.Inside( p ) )
          oPointList.push_back( pVertex );          
      }
    };

    
    typedef typename RangeSearch::CQuadNode< CVertexPtr > VertexTreeNodeT;
    //  CUniqueVertexMap( const CUniqueVertexMap & );
    Float fTolerance;
    Float f2Tolerance;  // 2 times tolerance -- used for error checking
    
    //---------------------------------------------------
    //  UniqueAdd
    //---------------------------------------------------
    void UniqueAdd( VertexTreeNodeT *node,
                    const PBRMath::BBox2D &nodeBound,
                    const CVertexPtr &data,
                    const PBRMath::BBox2D &dataBound, 
                    Float itemDiag2, Int depth = 0 )
    {
      //-------------------
      // leaf node condition - add to current node
      //
      //
      // Slight subtlety here:  One might worry that item's diagonal
      // will be larger than the current cell.  This will not be a problem
      // because the neighboring cells will also contain this item.
      //-------------------
      if(   ( depth == this->nMaxDepth )
            || DistanceSquared(nodeBound.pMin, nodeBound.pMax) < itemDiag2 )
      {
        
        if ( node->vData.size() == 0 )
        {
          node->vData.push_back(data);
        }
        else
        {
          // using the same distance scheme as the rest of the unique vertex map.
          bool bUnique = true;
          for( Size_Type i = 0; i < node->vData.size(); i ++ )
          {
            Float fXDiff = node->vData[i]->oPos.m_fX - data->oPos.m_fX;
            Float fYDiff = node->vData[i]->oPos.m_fY - data->oPos.m_fY;
            
            if( fabs( fXDiff ) < f2Tolerance && fabs( fYDiff ) < f2Tolerance )
              bUnique = false;
          }
          if( bUnique )
            node->vData.push_back(data);
        }
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
            node->opChildren[child] = new RangeSearch::CQuadNode<CVertexPtr>;
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
      
          UniqueAdd( node->opChildren[child], childBound, data, dataBound, 
                     itemDiag2, depth +1 );
        }
      }
    }   
    
  public:

    //---------------------------------------------------
    // Default Ctor
    //---------------------------------------------------
    CUniqueVertexMap( ): RangeSearch::CQuadtree<CVertexPtr>(),
                         fTolerance( std::numeric_limits<Float>::epsilon() ),
                         f2Tolerance( fTolerance * Float(2) ) {}
    
    CUniqueVertexMap( const PBRMath::BBox2D & oDataBound, Int nMaxDepth ):
      RangeSearch::CQuadtree<CVertexPtr>( oDataBound, nMaxDepth ),
      fTolerance( std::numeric_limits<Float>::epsilon() ),
      f2Tolerance( fTolerance * Float(2) ){}

    CUniqueVertexMap( const PBRMath::BBox2D & oDataBound, Int nMaxDepth, Float fError ):
      RangeSearch::CQuadtree<CVertexPtr>( oDataBound, nMaxDepth ),
      fTolerance( fError ), f2Tolerance( fTolerance * Float(2) ){}
        
    //---------------------------------------------------
    //  returns true if the pVertex exists
    //          false otherwise
    //---------------------------------------------------
    std::pair< bool, CVertexPtr > find( const CVertexPtr &pVertex ) 
    {
      PBRMath::Point oTestPoint( pVertex->oPos.m_fX, pVertex->oPos.m_fY );
      SVertexOverlapLookup oProcess( fTolerance );
      Lookup( oTestPoint, oProcess );
      std::pair< bool, CVertexPtr > oRes( oProcess.bFound, oProcess.pFoundVertex );
      return oRes;
    }

    //---------------------------------------------------
    // Find all vertices that are within the bounding box specified
    //---------------------------------------------------
    vector< CVertexPtr > find( const PBRMath::BBox2D & oDataBound )
    {
      SMultiVertexLookup oLookup;
      this->RangeSearch( oDataBound, oLookup );
      return oLookup.oPointList;
    }
    
    //---------------------------------------------------
    // insert, return false if the value already exists
    //---------------------------------------------------
    std::pair< bool, CVertexPtr > insert( const CVertexPtr & pVertex )
    {
      std::pair< bool, CVertexPtr > oTest = find( pVertex );      
        
      if( oTest.first )
      {
        oTest.first = false;   // oTest.second will contain the existing vertex ptr.
      }
      else
      {
        PBRMath::BBox2D oDataBound( Point( pVertex->oPos.m_fX - fTolerance, pVertex->oPos.m_fY - fTolerance ),
                                    Point( pVertex->oPos.m_fX + fTolerance, pVertex->oPos.m_fY + fTolerance ) );
        
        UniqueAdd( &this->oRoot, this->oTreeBound, pVertex, oDataBound, 
                   DistanceSquared( oDataBound.pMin, oDataBound.pMax) ); 
        oTest.first  = true;
        oTest.second = pVertex;
      }
      return oTest;
    }
    
    //---------------------------------------------------
    //  SetTolerance
    //---------------------------------------------------
    void SetTolerance( Float fTol )
    {
      fTolerance = fTol;
      f2Tolerance = Float( 2 ) * fTolerance;
    }
    
    //---------------------------------------------------
    //  PrintUniqueVertex
    //---------------------------------------------------
    void PrintUniqueVertex() const
    {
      this->PrintTree( );
    }

  };
  
  
};  // end namespace


#endif
