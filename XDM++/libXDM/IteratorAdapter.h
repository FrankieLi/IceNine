/////////////////////////////////////////////////////////////////
//
//  File:    PointIteratorAdapter
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Purpose: An adapter of iterators of pointer type
//
//
//
//
/////////////////////////////////////////////////////////////////

#ifndef _INTERATOR_ADAPTER_H_
#define _INTERATOR_ADAPTER_H_
#include <boost/any.hpp>

namespace XDMUtility
{
  template< typename DataT >
  struct PointerToDataMap
  {
    template< typename DataPtr >
    DataT operator()( const DataPtr & p ) const
    {
      return *p;
    }

    template<  typename DataPtr >
    DataT & operator()( const DataPtr p )
    {
      return *p;
    }
  };
  
  template< typename DataT >
  struct SecondOfPtrMap
  {
    template< typename DataPtr >
    DataT operator()( const DataPtr & p ) const
    {
      return p->Second;
    }

    template<  typename DataPtr >
    DataT & operator()( const DataPtr p )
    {
      return p->second;
    }
  };
  
  //-----------------------------------
  //
  //  IteratorAdapter
  //
  //  Purpose:  A simple adaptor that converts iterators supplied by
  //            the data structures of MicMesh into an iterator across
  //            DataT.  In another words, this acts as a pointer to
  //            the linked-list of the data, even though the data is only
  //            stored once by the internal structure.  This is taylored
  //            to the MicMesh class and SDataPoint.
  //
  //  Status:   Not tested
  //
  //-----------------------------------
  template< class DataT, class IteratorT, class PropertyMap = PointerToDataMap<DataT> >
  class IteratorAdapter
  {
  private:
    IteratorT pCurrent;
    PropertyMap GetData;
  public:
    
    //-----------------------------------
    //  Default ctor
    //-----------------------------------
    IteratorAdapter( ) {}
    
    //-----------------------------------
    //  Constructor
    //-----------------------------------
    IteratorAdapter( IteratorT pIter ): pCurrent( pIter ) {};

    //-----------------------------------
    // operator==
    //-----------------------------------
    bool operator==( const IteratorAdapter & oRHS ) const
    {
      return ( pCurrent == oRHS.pCurrent );
    }

    //-----------------------------------
    //  operator!=
    //-----------------------------------
    bool operator!= ( const IteratorAdapter & oRHS ) const
    {
      return ( pCurrent != oRHS.pCurrent );
    }

    //-----------------------------------
    //  operator->
    //-----------------------------------
    const DataT * operator->() const
    {
      return & GetData( *pCurrent );
    }
    //-----------------------------------
    //  operator->
    //-----------------------------------
    DataT * operator->()
    {
      return & GetData( *pCurrent );
    }

    //-----------------------------------
    //  operator*
    //-----------------------------------
    const DataT & operator*() const
    {
      return GetData( *pCurrent );
    }

    //-----------------------------------
    //  operator*
    //-----------------------------------
    DataT & operator*()
    {
      return GetData( *pCurrent );
    }
    
    //-----------------------------------
    // operator ++ (prefix)
    //-----------------------------------
    IteratorAdapter<DataT, IteratorT> & operator++()
    {
      pCurrent ++;
      return *this;
    }
    
    //-----------------------------------
    // operator ++ (postfix)
    //-----------------------------------
    IteratorAdapter<DataT, IteratorT> operator++( int )
    {
      pCurrent ++;
      return *this;
    }
  };




  //---------------------------------------------------
  //
  //  SerializingIterator  - convert 2D matrix to an iterator
  //
  //---------------------------------------------------
  template< class MatrixDataT >
  class SerializingIterator
  {
  protected:

    typedef SerializingIterator<MatrixDataT>    Self;
    typedef typename MatrixDataT::element       element;
    typedef typename MatrixDataT::element &     reference;
    typedef const typename MatrixDataT::element & const_reference;
    typedef typename MatrixDataT::index         index;

    MatrixDataT & Data;
    index Idx1;
    index Idx2;
    index Size1;
    index Size2;

  public:
    
    SerializingIterator( MatrixDataT & D_,
                         index Size1_, index Size2_,
                         index idx1_, index idx2_ )
      : Data( D_ ),
        Idx1( idx1_ ), Idx2( idx2_ ),
        Size1( Size1_ ), Size2( Size2_ ) {}
    
    bool operator==( const Self & RHS ) const
    {
      return ( Idx1 == RHS.Idx1 ) && ( Idx2 == RHS.Idx2 );
    }
    
    bool operator!=( const Self & RHS ) const
    {
      return !(*this == RHS);
    }
    
    SerializingIterator & operator++( )
    {
      if( Idx1 < Size1 )
      {
        if( Idx2 < (Size2 - 1) )
        {
          Idx2 ++;
        }
        else
        {
          Idx1 ++;
          Idx2 = 0;
        }
      }
      if( Idx1 >= Size1 )
        Idx2 = Size2;
      
      return *this;
    }
    
    SerializingIterator & operator++( int )
    {
      if( Idx1 < Size1 )
      {
        if( Idx2 < (Size2 - 1) )
        {
          Idx2 ++;
        }
        else
        {
          Idx1 ++;
          Idx2 = 0;
        }
      }
      
      if( Idx1 >= Size1 )
        Idx2 = Size2;
      return *this;
    }
    
    reference       operator*()       { return Data[Idx1][Idx2]; } 
    const_reference operator*() const { return Data[Idx1][Idx2]; }

    index Index1() const { return Idx1;      }
    index Index2() const { return Idx2; }


    element * operator->()             { return &( Data[Idx1][Idx2] ); }
    const element * operator->() const { return &( Data[Idx1][Idx2] ); }
    
  };
  
  template< class MatrixDataT >
  class Const_SerializingIterator: public SerializingIterator< MatrixDataT >
  {
  public:
    typedef SerializingIterator<MatrixDataT>    Base;
    typedef typename MatrixDataT::element &     reference;
    typedef const typename MatrixDataT::element & const_reference;
    typedef typename MatrixDataT::element       element;
    typedef typename MatrixDataT::index         index;

  public:

    Const_SerializingIterator( MatrixDataT & D_,
                               index Size1_, index Size2_,
                               index idx1_, index idx2_ )
      : Base( D_, Size1_, Size2_, idx1_, idx2_ ) {}
    
    const_reference operator*()        { return   Base::Data[Base::Idx1][Base::Idx2]; } 
    const_reference operator*() const  { return   Base::Data[Base::Idx1][Base::Idx2]; }
    const element * operator->()       { return &(Base::Data[Base::Idx1][Base::Idx2]); }
    const element * operator->() const { return &(Base::Data[Base::Idx1][Base::Idx2]); }    
  };
  
}

#endif
