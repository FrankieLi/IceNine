/////////////////////////////////////////////////////////////////
//
//  File:    PeakFile.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Purpose: Binary peak file format.
//
//
//
//
/////////////////////////////////////////////////////////////////

#include "Types.h"
#include "UFF.h"
#include "Pixel.h"
#include <set>

namespace GeneralLib
{

  //----------------------------------------------------------------
  // class CPeakFileInput  (for input only)
  //
  // Details:  Peak File format in binary for n pixels: 
  //           [float(32) version header]
  //           [SBlockHeader]
  //           [SubHeader for first coordinate]
  //           [ n uint(16) pixel first coordinate ]
  //           [SubHeader for second coordinate]
  //           [ n uint(16) pixel second coordinate ]
  //           [SubHeader for intensity ]
  //           [ n float(32) pixel intensity ]
  //           [SubHeader for peakID ]
  //           [ n uint(32) for peakID ]
  //
  //
  //----------------------------------------------------------------
  class CInputPeakFile
  {

  private:
    SBlock oPeakDataBlock;
    U32 nCurPixel;
    U32 nNumPixels;
    U32 nNumPeaks;
    CUFF oDataFile;
    
  public:

    //----------------
    // ctor
    //----------------
    CInputPeakFile(): nCurPixel( 0 ), nNumPixels( 0 )
    {}

    ~CInputPeakFile() { oDataFile.DeleteData(); }
    //----------------
    // Read
    //----------------
    bool  Read( const string & sFilename );
    
    //----------------
    //  ACCESSORS
    //----------------
    Int  NumPixels      ()  const { return static_cast<Int>( nNumPixels ); }
    Int  NumPeaks       ()  const { return static_cast<Int>( nNumPeaks  ); }
    bool GetNextPixel   ( Pixel & oRes, Int & nPeakID );
    void ResetPixelPos  ()        { nCurPixel = 0;  } 
    bool eof            ()        { return ( nCurPixel >= nNumPixels ); }
  };


  //----------------------------------------------------------------
  //  COutputPeakFile
  //----------------------------------------------------------------
  class COutputPeakFile
  {

  private:
    SBlock oPeakDataBlock;
    U32 nNumPixels;
    U32 nNumPeaks;
    U16* pPixelX;
    U16* pPixelY;
    F32* pIntensity;
    U16* pPixelID;
    CUFF oDataFile;
        
    string pBlock0Name;
    string pBlock1Name;
    string pBlock2Name;
    string pBlock3Name;
    string pBlock4Name;
  public:

    //----------------
    // constructor
    //----------------
    COutputPeakFile(): nNumPixels( 0 ), nNumPeaks(0),
                       pBlock0Name( "PeakFile\0" ),
                       pBlock1Name( "PixelCoord0\0" ),
                       pBlock2Name( "PixelCoord1\0" ),
                       pBlock3Name( "Intensity\0" ),
                       pBlock4Name( "PeakID\0" ) {}

    //----------------
    // destructor
    //----------------
    ~COutputPeakFile()
    {
      oDataFile.DeleteData();
      delete [] pPixelX;
      delete [] pPixelY;
      delete [] pIntensity;
      delete [] pPixelID;
      
    }
      
    //----------------
    // Write
    //----------------
    bool Write( const string & sFilename );
    
    //----------------
    //  InsertPixel
    //----------------
    bool InsertPixels( const vector< vector<Pixel> > & oPeakList,
                       const vector<Int> & oIDMap );

    //----------------
    //  InsertPixel
    //----------------
    bool InsertPixels( const vector< vector<Pixel> > & oPeakList );

    //----------------
    //  Mutators
    //----------------    
    //----------------
    //  SetNumPeaks
    //----------------
    void SetNumPeaks( Int n ) { nNumPeaks = static_cast<U32>( n ); }
    
  };
  
}
