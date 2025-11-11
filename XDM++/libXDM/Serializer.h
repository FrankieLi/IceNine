////////////////////////////////////////////////////////////
//
//  Serializer.h
//   Author:   Frankie Li ( sfli@cmu.edu )
//
//
//   Purpose:  A simple Serializer and Deserializer that
//             packs and unpacks arbitary data into
//             a byte array.
//
//
//   WARNING:  Current implementation of the Serializer is
//             a VERY THIN abstraction over type casting.
//             All functions with "Compact" are intended for
//             use with what C++ defines as P.O.D, or Plain Old
//             Data.  This is a very well defined concept under
//             C++ standard.  (Check google for detail.)
//
//             When in doubt, use the "Complex" version of the
//             serialization and deserialization functions.
//             This would require an implementation of Save and
//             Restore for the serializable class object.
//
//
////////////////////////////////////////////////////////////


#ifndef _SERIALIZER_H_
#define _SERIALIZER_H_

#include <string>
#include <vector>
#include "Debug.h"
#include "Types.h"
#include <tuple>

using std::vector;
using std::string;

namespace GeneralLib
{

  //----------------------------------
  //  DefaultCompactT
  //
  //  Default compact type.  The purpose of
  //  this is to make the insertion code uniform.
  //
  //----------------------------------
  template< typename ElementT >
  class DefaultCompactT
  {
  public:
    ElementT e;
    
    DefaultCompactT(){};
    DefaultCompactT( const ElementT & _e ): e(_e){};
  };

  
//----------------------------------
//  Class CSerializer
//----------------------------------
class CSerializer
{
private:
  string oBuf;  // using a string as storage
  Size_Type nLastInsertedPos;

public:


  
  CSerializer(): oBuf(), nLastInsertedPos( 0 )
  {};


  //----------------------------------
  //   InsertHeader
  //----------------------------------
  Size_Type InsertHeader( const Size_Type nElements,  const Size_Type nElementSize )
  {
    Size_Type nSize = sizeof( Size_Type ) / sizeof( char );

    oBuf.insert( nLastInsertedPos, reinterpret_cast< const char * >( &nElements ), nSize );
    nLastInsertedPos += nSize;
    oBuf.insert( nLastInsertedPos, reinterpret_cast< const char * >( &nElementSize ), nSize );
    nLastInsertedPos += nSize;

    return nLastInsertedPos;
  }

  //----------------------------------
  //  Generalized InsertCompactObj
  //
  //  Requirement:   oElement must be a continuous block, i.e., POD (Plain Old Data)
  //                 as defined by C++ standard. (see below)
  //
  //                 C++ compilers gaurantees that there will be no memory rearragement,
  //                 hidden virtual pointers, or other magical treatment in optimization
  //                 and processing of Plain Old Data (POD).  (Note that padding may still
  //                 happen.)  POD is defined as follows:
  //
  //                 3.9(10): "Arithmetic types (3.9.1), enumeration types,
  //                 pointer types, and pointer to member types (3.9.2) and
  //                 cv-qualified versions of these types (3.9.3) are collectively
  //                 caller scalar types. Scalar types, POD-struct types, POD-union
  //                 types (clause 9), arrays of such types and cv-qualified versions
  //                 of these types (3.9.3) are collectively called POD types"
  //
  //                 9(4): "A POD-struct is an aggregate class that has no non-static
  //                 data members of type non-POD-struct, non-POD-union (or array of
  //                 such types) or reference, and has no user-define copy operator and
  //                 no user-defined destructor. Similarly a POD-union is an aggregate
  //                 union that has no non-static data members of type non-POD-struct,
  //                 non-POD-union (or array of such types) or reference, and has no
  //                 user-define copy operator and no user-defined destructor.
  //
  //                 8.5.1(1): "An aggregate is an array or class (clause 9) with no
  //                 user-declared constructors (12.1), no private or protected
  //                 non-static data members (clause 11), no base classes (clause 10)
  //                 and no virtual functions (10.3)."
  //
  //
  //                 i.e., this is a shallow copy.  If oElement
  //                 contains any pointers, vectors, or complex
  //                 data structure that is not compact, InsertComplexObj
  //                 must be called.
  //
  //  Action:  Reinterpret oElement as a collection of characters.
  //           Note that this is possible because character is supposed
  //           to be a 1 byte, by C/C++ definition.
  //
  //  Rationale:  The reason to keep "char" instead of typedef-ing it to
  //              byte is to make sure that in case future changes of C/C++
  //              definition occurs, this wouldn't be become a hidden bug
  //              somehow.
  //
  //  Return:    Return true if operation is succesful.  False is returned
  //             otherwise. 
  //----------------------------------
  template< typename ElementT >
  Bool InsertCompactObj( const ElementT & oElement )
  {
    const char * p = reinterpret_cast< const char * > ( & oElement );
    Size_Type nCharSize = sizeof( oElement ) / sizeof( char );

    InsertHeader( 1, nCharSize );
    oBuf.insert( nLastInsertedPos, p, nCharSize );
    nLastInsertedPos += nCharSize;
    
    return ( oBuf.max_size() >= nLastInsertedPos );
  }
  
  //----------------------------------
  //  InsertComplexObj
  //
  //  Requirement:  oElement must implement the function:
  //                bool Save( CSerializer & oSerialBuf )
  //                bool Restore( CDeserialier & oSerialBuf )
  //              
  //
  //----------------------------------
  template< typename ElementT >
  Bool InsertComplexObj( const ElementT & oElement )
  {
    return oElement.Save( *this );
  }
  
  //----------------------------------
  //  InsertCompactVector
  //  
  //  Requirement:  ElementT must be a compact type.  i.e., it must
  //                be all located in one contiuous block of memory.
  //
  //
  //  Rationale:    Since vector insertion is used so often, this is
  //                a specialized version.  (reduced memory requirement
  //                easier to call)
  //
  //  Return:    Return true if operation is succesful.  False is returned
  //             otherwise. 
  //----------------------------------
  template< typename ElementT >
  Bool InsertCompactVector( const vector< ElementT > & vElements )
  {
    Size_Type nCharSize = sizeof( ElementT ) / sizeof( char );
    Size_Type nElements = vElements.size();
    InsertHeader( nElements, nCharSize );
    
    for( Size_Type i = 0; i < vElements.size(); i ++ )
    {
      const char * p = reinterpret_cast< const char * > ( & vElements[i] );
      oBuf.insert( nLastInsertedPos, p, nCharSize );
      nLastInsertedPos += nCharSize;
    }

    return ( oBuf.max_size() >= nLastInsertedPos );
  }


  //----------------------------------
  //
  //  InsertComplexVector
  //
  //  Requirement:  oElement must implement the function:
  //                bool Save( CSerializer & oSerialBuf )
  //                bool Restore( CDeserialier & oSerialBuf )
  //                static Size_Type GetCompactSize();
  //----------------------------------
  template< typename ElementT >
  Bool InsertComplexVector( const vector< ElementT > & vElements )
  {
    Size_Type nCharSize = 0; // this is not well defined
    Size_Type nElements = vElements.size();
    InsertHeader( nElements, nCharSize );
    
    for( Size_Type i = 0; i < vElements.size(); i ++ )
    {
      vElements[i].Save( *this );
    }
    return ( oBuf.max_size() >= nLastInsertedPos );
  }

                          
  //----------------------------------
  //  GetBuffer
  //
  //  Return a pointer to the buffer
  //
  //  Warning:  The pointer returned is only garanteed to
  //            remain unchanged until the next non-const function
  //            is called!!!!
  //
  //----------------------------------
  const char * GetBuffer() const
  {
    return oBuf.data();
  }

  //----------------------------------
  //  GetSize
  //----------------------------------
  Size_Type GetSize() const
  {
    return oBuf.size();
  }

};

//----------------------------------
//  Class CDeserializer
//----------------------------------
class CDeserializer
{
  
private:
  
  char *pBuf;
  Size_Type nBufSize;
  Size_Type nDeserializerPos;

  CDeserializer( const CDeserializer & oRHS ); 
  
public:

  //----------------------------------
  //  Default constructor
  //----------------------------------
  CDeserializer(): pBuf( NULL ), nBufSize(0), nDeserializerPos( 0 ){};

  //----------------------------------
  //  Constructor
  //----------------------------------
  CDeserializer( Size_Type nLength ): nBufSize( nLength ), nDeserializerPos( 0 )
  {
    pBuf = new char[ nLength ];
  };
  
  //----------------------------------
  //  Destructor
  //----------------------------------
  ~CDeserializer()
  {
    if( pBuf )
      delete [] pBuf;
    nBufSize = 0;
  }

  //----------------------------------
  //  Set
  //----------------------------------
  void Set( Size_Type nLength )
  {
    if( pBuf )
      delete [] pBuf;

    pBuf = new char[ nLength ];
    nBufSize = nLength;
  }

  //----------------------------------
  //   GetHeader
  //----------------------------------
  Size_Type GetHeader( Size_Type *pnElements,  Size_Type * pnElementSize )
  {
    Size_Type nSize = sizeof( Size_Type ) / sizeof( char );

    const char * p1 = pBuf + nDeserializerPos; 
    memcpy( reinterpret_cast< char* >( pnElements ), p1, nSize );
    nDeserializerPos += nSize;

    const char * p2 = pBuf + nDeserializerPos; 
    memcpy( reinterpret_cast< char* >( pnElementSize ), p2, nSize );
    nDeserializerPos += nSize;

    if ( nDeserializerPos >= nBufSize )
      return string::npos;
    else
      return nDeserializerPos;
  }
  
  
  //----------------------------------
  //
  //   GetCompactObj
  //
  //   Action:    Given a variable, unpack from the
  //              buffer to variable.  Note that ElementT
  //              must be a compact object.  In another words,
  //              it must reside in one continuous memory block.
  //
  //   Return:    Return true if object extraction is successful.
  //              False otherwise.
  //
  //   Rationale:  pRes is a pointer instead of a reference to emphasize that
  //               it will be changed.
  //
  //----------------------------------
  template< typename ElementT >
  Bool GetCompactObj( ElementT * pRes )
  {
    Size_Type nElements;
    Size_Type nElementSize;
    GetHeader( &nElements, &nElementSize );
    Size_Type nCharSize = sizeof( (*pRes) ) / sizeof( char );
    
    DEBUG_ASSERT( nElements == 1, "[Deserializer::Get] -- More than 1 element!  Use the iterator get!\n");
    DEBUG_ASSERT( nElementSize == nCharSize, "[Deserializer::Get] -- Type size mismatch!" );

    const char * p = pBuf + nDeserializerPos; 
    
    if( nCharSize + nDeserializerPos > nBufSize )
      return false;
    
    memcpy( reinterpret_cast< char * >(  pRes ), p, nCharSize );
    nDeserializerPos += nCharSize;

    return ( nBufSize >= nDeserializerPos );
  }

  //----------------------------------
  //   GetComplexObj
  //  Requirement:  oElement must implement the function:
  //                bool Save( CSerializer & oSerialBuf )
  //                bool Restore( CDeserialier & oSerialBuf )
  //                
  //----------------------------------
  template< typename ElementT >
  Bool GetComplexObj( ElementT * pRes )
  {
    return pRes->Restore( *this );
  }

  //----------------------------------
  //   GetCompactVector
  //
  //   Purpose:   A further specialized version
  //              for cases where ElementT is already compact.
  //              i.e., there shallow copy for elements in vList.
  //
  //   Note:      There is a lot of overlapping code between
  //              this and the generalized GetVector.
  //              This is so that we can remove the requirement
  //              of a default constructor for this function.
  //----------------------------------
  template< typename ElementT >
  Bool GetCompactVector( vector<ElementT> & vList )
  {

    Size_Type nElements;
    Size_Type nElementSize;
    GetHeader( &nElements, &nElementSize );
    Size_Type nCharSize = sizeof( ElementT ) / sizeof( char );
    
    Size_Type nArraySize = nElements * nElementSize;
    DEBUG_ASSERT( nElementSize == nCharSize, "[Deserializer::GetVector] -- Type size mismatch!" );
    DEBUG_ASSERT( nDeserializerPos + nArraySize <= nBufSize,
                  "[Deserializer::GetVector] -- Exceeding Buffer Length!\n" );

    if( nDeserializerPos + nArraySize > nBufSize )
    {
      return false;
    }
    
    char * p = pBuf + nDeserializerPos;
    ElementT *pCur   = reinterpret_cast< ElementT * >( p );
    ElementT *pLast  = reinterpret_cast< ElementT * >( p + nArraySize );
    vList.clear();
    for( ; pCur != pLast; pCur ++ )
    {
      vList.push_back( *pCur );
    }
    
    nDeserializerPos += nArraySize;
    
    return ( nBufSize >= nDeserializerPos );
  }

  //----------------------------------
  //
  //  GetComplexVector
  //  Requirement:  oElement must implement the function:
  //                bool Save( CSerializer & oSerialBuf )
  //                bool Restore( CDeserialier & oSerialBuf )
  //                
  //----------------------------------
  template< typename ElementT >
  Bool GetComplexVector( vector<ElementT> & vList )
  {
    Size_Type nElements;
    Size_Type nElementSize;   // nElementSize is a placeholder
    
    GetHeader( &nElements, &nElementSize );
            
    Bool bSuccess = true;
    vList.clear();
    vList.resize( nElements );
    for( Size_Type i = 0; i < nElements; i ++  )
    {
      bSuccess = bSuccess && vList[i].Restore( *this );
    }

    return bSuccess;
  }
  
  //----------------------------------
  //  ResetDeserializer
  //----------------------------------
  void Reset()
  {
    nDeserializerPos = 0; 
  }

  //----------------------------------
  //  GetBuffer
  //
  //  Return a pointer to the buffer
  //
  //----------------------------------
  char * GetBuffer()
  {
    return pBuf;
  }

  //----------------------------------
  //  GetSize
  //----------------------------------
  Size_Type GetSize() const
  {
    return nBufSize;
  }
  
  void Resize( Size_Type n )
  {
    if( nBufSize > 0 )
    {
      char * pNewBuf = new char[ n ];
      memcpy( pNewBuf, pBuf, nBufSize );
      delete pBuf;
      pBuf = pNewBuf;
      nBufSize = n;
    }
    else
    {
      pBuf = new char [ n ];
      nBufSize = n;
    }
  }
};
  
}  // namespace GeneralLib

#endif
