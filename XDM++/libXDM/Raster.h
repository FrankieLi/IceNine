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
//  - reverted
////////////////////////////////////////////////////////////


#ifndef _RASTER_H_
#define _RASTER_H_


#define BOOST_NO_HASH   // for boost graph

#include "Debug.h"
#include <vector>
#include <string>
#include <fstream>
#include "Types.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Pixel.h"
#include "SutherlandHodgman.h"
#include "Error.h"


#include <boost/function.hpp>

#define BOOST_UBLAS_STRICT_MATRIX_SPARSE
using std::ofstream;
using std::vector;
using std::cout;
using std::endl;
using std::string;

#include <map>
using std::map;

using namespace PBRMath;

namespace GeneralLib
{

  
  //------------------------------------------------------------------------------------
  //
  //  Trait class used for Raster
  //  This gives us the trains of a 2D matrix concept (?)
  //
  //  TODO:  Complete this trait class to generalize Raster?  Or just switch out
  //         the data type and use boost?  (perhaps having an abstract in between would
  //         make life easier.)
  //
  //         Note that we want to do template specialization to different traits of different 
  //         kinds of CMatrix2D
  //------------------------------------------------------------------------------------
  template < class CMatrix2D >
  struct SMatrix2DTrait
  {
    typedef typename CMatrix2D::MatrixStorageT MatrixT;
    typedef typename CMatrix2D::ElementConstRefT MatrixElementConstRefT;
    typedef typename CMatrix2D::ElementRefT MatrixElementRefT;
      
    // NOTE:  Visual Studio may not work with this.  We're using
    //        nested iterators here.  Note also that the number '1'
    //        and '2' are here to name the dimensions.
    //
    typedef typename CMatrix2D::const_iterator1 Iterator1;
    typedef typename CMatrix2D::const_iterator2 Iterator2;
  };

//------------------------------------------------------------------------------------
//
//  class  CRaster
//
//  An abstraction of a raster to hide the detailed implementation, i.e., how pixels
//  are stored and shapes rasterized.
//
//------------------------------------------------------------------------------------
template< typename DataT, class MatrixT = boost::numeric::ublas::compressed_matrix<DataT> >
class CRaster
{
  
public:  // type declerations
  
  typedef MatrixT  SparseMatrixT; 
  typedef typename SparseMatrixT::const_reference MatrixElementT_ConstRef;
  typedef typename SparseMatrixT::value_type MatrixElementT;
  typedef typename SparseMatrixT::reference MatrixElementT_Reference;
  
  // NOTE:  Visual Studio may not work with this.  We're using
  //        nested iterators here.  Note also that the number '1'
  //        and '2' are here to name the dimensions.
  //
  typedef typename SparseMatrixT::const_iterator1 SparseIterator1;
  typedef typename SparseMatrixT::const_iterator2 SparseIterator2;


private:
  //------------------------------------------
  //  Scan line related data structure
  //
  //  A buffer used for rasterising of convex polygon
  //  It is a nDetectorHeight x 2 array of int, used to
  //  designate ends of the scanlines.
  //
  //  It will be initialized once at the construction of this object
  //
  //  The reason that this is not a local variable is because we do not
  //  know nDetectorHeight at compile time.  (That's something from the
  //  configuation file.)
  //
  //------------------------------------------
  mutable vector< vector<Int> > oEdgeTable;

  //------------------------------------------
  //   Allocate the memory required for pRasterEdgeBuffer based on the hight
  //   of the detector.  Return false if the buffer is already allocated
  //   or if there is not enough memory.
  //------------------------------------------
  void InitializeEdgeBuffer( Int nRasterHight );

  //------------------------------------------
  //  Safely remove pRasterEdgeBuffer and delete its memory
  //------------------------------------------
  bool DeleteEdgeBuffer( Int nRasterHeight );
  
  //------------------------------------------
  //  Reset the edge buffer between nYMin, and nYMax, inclusively 
  //------------------------------------------
  void ResetEdgeBuffer( Int nYMin, Int nYMax );

  //------------------------------------------
  //  Copy, make *this the same as oRhs
  //------------------------------------------
  void Copy( const CRaster< DataT, MatrixT > & oRHS );  
  
  //------------------------------------------
  //  Save endpoints of the scanline to pEdgeTable
  //  
  //  Precondition:  pEdgeTable initialized to the size of the screen.  oVertex0, oVertex1
  //  will not have pixels outside of pEdgeTable
  //------------------------------------------
  void CalculateScanline( const Pixel & oVectex0, const Pixel & oVertex1 );

  //------------------------------------------
  //  Delete current image buffer
  //------------------------------------------
  void DeleteImageBuffer();

  //------------------------------------------
  //  InitializeImageBuffer
  //------------------------------------------
  void InitializeImageBuffer( Int nWidth, Int nHeight );

 
  //-------------------------
  //  Calculation Helpers
  //-------------------------

  
  //------------------------------------------
  //  ClipPolygon
  //------------------------------------------
  inline vector<Pixel> ClipPolygon( const vector<Point> &oPolygon ) const
  {
    vector<Pixel> oClippedPolygon;
    SutherlandHodgman oClipper( RectF( 0, 0, nRasterWidth - 1 , nRasterHeight - 1) );
    oClipper.Clip( oClippedPolygon, oPolygon);
    return oClippedPolygon;
  }


protected:
  
  Int nRasterHeight;
  Int nRasterWidth;

 
  SparseMatrixT pImage;

  //------------------------------------------
  //  Copy constructor is forbidden
  //------------------------------------------
  // CRaster( const CRaster & oRaster); 
public:
  
  CRaster( Int nRasterWidth, Int nRasterHeight );

  //------------------------------------------
  //  Default constructor sets size to (0, 0)  A
  //  resize is necessary.  Use with caution
  //------------------------------------------
  CRaster() :  oEdgeTable(), nRasterHeight( 0 ), nRasterWidth( 0 ){}
  
  //------------------------------------------
  //  Copy Constructor
  //------------------------------------------
  CRaster( const CRaster< DataT, MatrixT > &  oRHS );
  
  //------------------------------------------
  //  operator=
  //------------------------------------------
  CRaster<DataT, MatrixT> & operator= ( const CRaster<DataT, MatrixT>  &oRHS );

  //------------------------------------------
  // Destructor
  //------------------------------------------
  ~CRaster();


  //------------------------------------------
  //  GeneralRasterizePolygon
  //
  //  TODO:  Figure out a way to remove the last parameter from
  //  the list.  I.e., somehow make it so that passing the type is
  //  all that's necessary.
  //------------------------------------------
  template< typename HelperFnT >
  void GeneralRasterizePolygon( const vector<Point> & oPolygon, DataT fFillValue,
                                HelperFnT  &HelperFn );
  
  //------------------------------------------
  // Rasterize, or fill the polygon specified by vertices.
  // Vertices must be in winding order. (i.e, right hand rule).
  //------------------------------------------
  void RasterizePolygon( const vector<Point> & oVertices, DataT fFillValue );

  //------------------------------------------
  // Rasterize, or fill the triangle specified by vertices.
  // Vertices must be in winding order. (i.e, right hand rule).
  //------------------------------------------
  void RasterizeTriangle( const Point & v0,
                          const Point & v1, const Point &v2,
                          DataT fFillValue );

  //------------------------------------------
  //  Modifier functions 
  //------------------------------------------
  void Resize( Int nRasterWidth, Int nRasterHeight );
  void Fill( DataT value );
  void SetPixel( Int nI, Int nJ, DataT fValue );
  void AddToPixel( Int nI, Int nJ, DataT fValue );
  
  //------------------------------------------
  //   Clears the image
  //------------------------------------------
  void ClearImage();

  //------------------------------------------
  //   Print raster to the file specified by the filename
  //------------------------------------------
  void PrintRaster( string sFilename ) const;

  //------------------------------------------
  // Accessor:
  //------------------------------------------
  
  Int GetHeight() const;
  Int GetWidth() const;
  SparseIterator1 GetBegin1( ) const;
  SparseIterator1 GetEnd1( ) const;
  typename SparseMatrixT::const_reference operator() ( Int nI, Int nJ ) const;
  typename SparseMatrixT::reference       operator() ( Int nI, Int nJ );

  //------------------------------------------
  //  Check to see if nI, nJ is in the bound of the image
  //------------------------------------------
  Bool IsInBound( Int nI, Int nJ ) const;

  //------------------------------------------------------------------------------------
  //
  //  struct SPixelProcessor
  //
  //
  //  Templated objects with state memory (structs) for generalized usage of the
  //  rasterize function ( GenRasterizePolygon ).  Note that an empty templated class
  //  is presented for the general case.  Specific implementation (for example, pixel by
  //  pixel Lorentz factor) may be implemented using template specialization, as exemplified
  //  below.  
  //
  //  Note:  The use of template is to prevent dynamic polymorphism.  This is essentially
  //         a selector, so there is little reason that we need anything more than static
  //         polymorphic types.  Since rastering is something that's computationally intensitve,
  //         and the pixel processing function is called so many times, the idea of using static
  //         polymorphic type seems more appealing.
  //
  //  TODO:  Make generalized rastering functions availble to the public.
  //
  //  TODO:  Make this something that the Trait Class specifies.  This way no instantiation
  //         of function object is required by the user.
  //
  //
  //------------------------------------------------------------------------------------

  
  //------------------------------------------
  //  Specialized Template for accumulate the
  //  intensity of the projected peaks.
  //------------------------------------------
  struct AccumulatorT 
  {
    inline void operator() ( MatrixElementT_Reference p, DataT f)
    {
      p += f;
    }
  };

  //------------------------------------------
  //  Specialized for assignment of value
  //------------------------------------------
  struct AssignerT
  {
    inline void operator() ( MatrixElementT_Reference p, DataT f)
    {
      p = f;
    }
  };
    
  //------------------------------------------
  //  Specialized Template to count the number
  //  of pixel overlap.  
  // 
  //  Note that nPixelOverlap must be initialized
  //  before use.
  // 
  //
  //------------------------------------------
  struct OverlapCounterT
  {
    Int nPixelOverlap;
    inline void operator() ( MatrixElementT_ConstRef p, DataT f)
    {
      if ( p > 0 )
        nPixelOverlap ++;
    }
  };

  //------------------------------------------
  //  Specialized Template to count the number
  //  of pixels to be rastered on the detector
  //------------------------------------------
  struct HitCounterT
  {
    Int nPixelOnDetector;
    inline void operator() ( MatrixElementT_ConstRef p, DataT f)
    {
      nPixelOnDetector ++;
    }
  };
  

  
};  
  

} // end of GeneralLib namespace

#include "Raster.tmpl.cpp"


#endif
