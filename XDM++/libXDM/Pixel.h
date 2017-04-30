#ifndef PIXEL_H_
#define PIXEL_H_

#include "Types.h"
#include "BBox.h"

using namespace PBRMath;

class Pixel
{
 public:

  Pixel() {};
  Pixel( Int _x, Int _y ): x(_x), y(_y) {};
  Pixel(const Point & p):  x(p.x), y(p.y) { };

  Int x;
  Int y;
  IntensityT fIntensity;
  
  bool operator<( const Pixel &rhs ) const
  {
    return ( fIntensity < rhs.fIntensity );
  } 
  
  bool operator>( const Pixel &rhs ) const
  {
    return ( fIntensity > rhs.fIntensity );
  }
  
};




//-------------------------------------------------------------------------------------
//
//  CComparablePixel:Pixel  - compares Pixel types
//
//-------------------------------------------------------------------------------------
class CPixelComparable //: public Pixel
{
private:
  
  Int nWidth;
  Int nHeight;

  
  // greater than is not defined
  bool operator > ( const CPixelComparable& oRHS ) const;
  bool operator <= ( const CPixelComparable& oRHS ) const;
  bool operator == ( const CPixelComparable& oRHS ) const;
  bool operator >= ( const CPixelComparable& oRHS ) const;
public:

    
  Int x;
  Int y;
  Float fIntensity;
  
  CPixelComparable(){};
  CPixelComparable( Int nW, Int nH):nWidth( nW ), nHeight( nH ){};
  CPixelComparable(const Point & p):  x(p.x), y(p.y) { };

  
  //-----------------------------------------
  // Strictly Weak Ordering.
  //
  // return true if oPixel1 comes before oPixel2 in the ordering specified
  //-----------------------------------------
  inline bool operator < (const CPixelComparable& oRHS ) const
  {
    return ( ( x * nWidth + y)  < ( oRHS.x * oRHS.nWidth + oRHS.y )  ); 
  }

//   inline bool operator () (const CPixelComparable& oLHS, const CPixelComparable& oRHS ) const
//   {
//     return ( ( oLHS.x * oLHS.nWidth + oLHS.y)  < ( oRHS.x * oRHS.nWidth + oRHS.y )  ); 
//   }

  inline CPixelComparable  operator=( const CPixelComparable & oLHS )
  {

    nWidth = oLHS.nWidth;
    nHeight = oLHS.nHeight;
    x = oLHS.x;
    y = oLHS.y;
    fIntensity = oLHS.fIntensity;
    
    return *this;
  }

  inline Int MapToInt() const
  {
    return x * nWidth + y;
  }
};





#endif
