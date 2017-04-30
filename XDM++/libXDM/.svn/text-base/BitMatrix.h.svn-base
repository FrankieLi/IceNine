////////////////////////////////////////////////////////////
//
//  Filename:  Raster.h
//  Author:    Frankie Li
//
//
//  Purpose:   Detail implementation of rastering for convex
//             polygons and possibly generalized nonconvex ones.
//
//  NOTE:      This file was originally in the IceNine source tree.
//             It has been copied over to the XDM++ tree.  We'll have to
//             keep track of these two verions manually for now.
//
////////////////////////////////////////////////////////////
#ifndef BIT_MATRIX_H_
#define BIT_MATRIX_H_


#include <boost/dynamic_bitset.hpp>
#include <string>
#include <iostream>

namespace GeneralLib
{
  //-------------------------------------
  //  Forward declarations
  //-------------------------------------
  class BitMatrix;
  
  template< class MatrixType >
  class BitMatrixIterator1;

  template< class MatrixType >
  class BitMatrixIterator2;

  
  //-------------------------------------------------------
  //  BitMatrixIterator2
  //-------------------------------------------------------
  template< class MatrixType >
  class BitMatrixIterator2
  {
  public:
    typedef std::string::size_type      size_type;
    typedef typename MatrixType::reference  reference;

  protected:
    size_type CurrentIdx1;
    size_type CurrentIdx2;
    size_type Size2;

    MatrixType & Data;
  public:

    BitMatrixIterator2( size_type Idx1, size_type Idx2, size_type Size2_,
                        MatrixType & Data_)
      : CurrentIdx1( Idx1 ), CurrentIdx2( Idx2 ), Size2( Size2_ ), Data( Data_ ) {}

    bool operator==(const BitMatrixIterator2 & RHS ) const
    {
      return ( CurrentIdx2 == RHS.CurrentIdx2 ) && ( CurrentIdx1 == RHS.CurrentIdx1 );
    }

    bool operator!=(const BitMatrixIterator2 & RHS ) const
    {
      return !(*this == RHS);
    }
    
    BitMatrixIterator2 & operator++( )
    {
      if( CurrentIdx2 < Size2 )
        CurrentIdx2 ++;
      while( CurrentIdx2 < Size2
              && ! Data( CurrentIdx1, CurrentIdx2 ) )
        CurrentIdx2 ++;
      

      return *this;
    }

    BitMatrixIterator2 & operator++( int )
    {
      if( CurrentIdx2 < Size2 )
        CurrentIdx2 ++;
      
      while( CurrentIdx2 < Size2
             && ! Data( CurrentIdx1, CurrentIdx2 ) )
        CurrentIdx2 ++;
       
      return *this;
    }
    
    reference operator*()       { return Data( CurrentIdx1, CurrentIdx2 ); } 
    bool      operator*() const { return Data( CurrentIdx1, CurrentIdx2 ); }

    size_type index1() const { return CurrentIdx1; }
    size_type index2() const { return CurrentIdx2; }

  };

 

  template< class MatrixType >
  class Const_BitMatrixIterator2: public BitMatrixIterator2< MatrixType >
  {

  public:
    typedef std::string::size_type      size_type;
    typedef BitMatrixIterator2<MatrixType> Base;
    Const_BitMatrixIterator2( size_type Idx1, size_type Idx2, size_type Size2_,
                              MatrixType & Data_)
      :Base( Idx1, Idx2, Size2_, Data_ ) {}
    bool operator*()       { return Data( Base::CurrentIdx1, Base::CurrentIdx2 ); } 
    bool operator*() const { return Data( Base::CurrentIdx1, Base::CurrentIdx2 ); }
  };

  
  //-------------------------------------------------------
  //  BitMatrixIterator1
  //-------------------------------------------------------
  template< class MatrixType >
  class BitMatrixIterator1
  {
  public:
    typedef BitMatrixIterator2<MatrixType>  BitMatrixIter2;
    typedef BitMatrixIterator1              BitMatrixIter1; 
    
    typedef std::string::size_type size_type;

  protected:
    size_type CurrentIdx1;
    size_type Size1;
    size_type Size2;
    
    MatrixType & Data;
  public:

    BitMatrixIterator1( size_type Idx1,
                        size_type Size1_, size_type Size2_,
                        MatrixType & Data_)
      : CurrentIdx1( Idx1 ), Size1( Size1_ ), Size2( Size2_ ), Data( Data_ ) {}



    BitMatrixIter1 & operator++( )
    {
      if( CurrentIdx1 < Size1 )
        CurrentIdx1 ++;
      return *this;
    }

    BitMatrixIter1 & operator++( int )
    {
      if( CurrentIdx1 < Size1 )
        CurrentIdx1 ++;
      return *this;
    }
    
    BitMatrixIter2 begin() const
    {
      return BitMatrixIter2( CurrentIdx1, 0, Size2, Data );
    }

    BitMatrixIter2 end() const
    {
      return BitMatrixIter2( CurrentIdx1, Size2, Size2, Data );
    }

    bool operator==(const BitMatrixIter1 & RHS ) const
    {
      return ( CurrentIdx1 == RHS.CurrentIdx1 );
    }
    
    bool operator!=(const BitMatrixIter1 & RHS ) const
    {
      return !(*this == RHS);
    }
    
  };

  template< class MatrixType >
  class Const_BitMatrixIterator1: public BitMatrixIterator1< MatrixType >
  {
  public:
    typedef BitMatrixIterator1<MatrixType> Base;
    typedef std::string::size_type      size_type;
    
    typedef Const_BitMatrixIterator2<MatrixType>  Const_BitMatrixIter2;
    Const_BitMatrixIterator1( size_type Idx1,
                              size_type Size1_, size_type Size2_,
                              MatrixType & Data_)
      :Base(Idx1, Size1_, Size2_, Data_ ) {}

    Const_BitMatrixIter2 begin() const
    {
      return Const_BitMatrixIter2( Base::CurrentIdx1, 0, Base::Size2, Base::Data );
    }

    Const_BitMatrixIter2 end() const
    {
      return Const_BitMatrixIter2( Base::CurrentIdx1, Base::Size2, Base::Size2, Base::Data );
    }

  };
  
  //-------------------------------------------------------
  //  BitMatrix Class
  //-------------------------------------------------------
  class BitMatrix
  {
  public:
    typedef boost::dynamic_bitset<> BitArrayT;
    typedef BitArrayT::size_type    size_type;

    typedef BitMatrixIterator1<BitMatrix>  iterator1;
    typedef BitMatrixIterator2<BitMatrix>  iterator2;
    
    typedef Const_BitMatrixIterator1<const BitMatrix>  const_iterator1;
    typedef Const_BitMatrixIterator2<const BitMatrix>  const_iterator2;

    typedef bool   value_type;
    typedef bool   const_reference;
    typedef BitArrayT::reference  reference;    // forcing assignment to not work in a normal way
    
    
  private:

    BitArrayT        Data;
    size_type        Size1;    // first index (slower moving one)
    size_type        Size2;
    
  public:

    //----------------------------------
    //  default ctor
    //----------------------------------
    BitMatrix(): Data(), Size1( 0 ), Size2( 0 )  {}

    //----------------------------------
    //  ctor  
    //----------------------------------
    BitMatrix( size_type Size1_, size_type Size2_):
      Size1( Size1_ ), Size2( Size2_)
    {
      Data.resize( Size1 * Size2, false );
    }

    //----------------------------------
    //  Copy Ctor and operator == are defaults
    //----------------------------------

    //----------------------------------
    //  Resize  (preserve is there as a compatiblity for space matrix)
    //----------------------------------
    void resize( size_type Size1_, size_type Size2_, bool preserve = false )
    {
      Size1 = Size1_;
      Size2 = Size2_;
      Data.resize( Size1 * Size2_, false);
    }

    void clear()   {  Data.clear(); }

    bool  operator()( size_type i, size_type j ) const  {  return Data[i * Size1 + j]; }
    //----------------------------------
    //  Accessors
    //----------------------------------
    reference operator()( size_type i, size_type j ) {  return  Data[i * Size1 + j]; }

    size_type size1() const { return Size1; }
    size_type size2() const { return Size2; }

        

    iterator1 begin1()  { return iterator1( 0,     Size1, Size2, *this ); }
    iterator1 end1()    { return iterator1( Size1, Size1, Size2, *this ); }
    
    // ----- const_iterator is actually not implemented correctly...
    const_iterator1 begin1() const { return const_iterator1( 0,     Size1, Size2, *this ); }
    const_iterator1 end1()   const { return const_iterator1( Size1, Size1, Size2, *this ); }
  };


  //---------------------------------------
  //  Define += operator for reference
  //---------------------------------------
  template< class T >
  BitMatrix::reference operator+=( BitMatrix::reference LHS,  const T & RHS  )
  {
    LHS = bool( LHS ) || (RHS > 0);
    return LHS;
  }

  template<  >
  BitMatrix::reference operator+=( BitMatrix::reference LHS,  const bool & RHS  );


  
  
  
}// end namespace GeneralLib


#endif
