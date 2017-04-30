/////////////////////////////
//  
//  File:    TiffReader.h
//  Author:  Frankie Li
//  Contact: sfli@cmu.edu
//
//  Purpose: A simple wrapper to facilitate typical TIFF input/output needs
//
/////////////////////////////
#ifndef TIFF_READER_H_
#define TIFF_READER_H_
#include "Types.h"
#include <string>
#include <iostream>
#include <tiffio.h>
#include <stdlib.h>
#include "Error.h"
#include <algorithm>
#include <fstream>

using std::ofstream;
using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::max;

//
//  CTiffTags
//  Interfacing with tiffio.h
//
//  These are minimal tags required to read and write a Tiff file
//
class CTiffTags{

public:
  uint32 nImageWidth;
  uint32 nImageLength;
  uint16 nBitsPerSample;
  uint16 nSamplesPerPixel;
  uint32 nRowsPerStrip;
  uint16 nCompression;
  uint16 nPhotometric;
  uint16 nFillOrder;
  uint16 nPlanarConfig;
  float  fXResolution;
  float  fYResolution;
  uint16 nResolutionUnit;

  int32 nStripSize;
  int32 nNumStrips;
};

class CImage{

private:
  bool IsBigEndian(const string & filename);

  
  // Tiff information
  CTiffTags oTags;

  //
  //  Copy values of images
  //
  void DeepCopy( const CImage &oImage );
  
  void Allocate(UInt** & pImage, UInt nRow, UInt nCol);

  void Delete(UInt** & pImage, UInt nRow);
  
public:
  UInt **pPixelValue;
  UInt nNumRow;
  UInt nNumCol;
  

  CImage();
  
  ~CImage();

  /////////////////////////
  //
  // Copy Constructor
  //
  /////////////////////////
  CImage( const CImage & oImage );

  //
  //  Assignment
  //
  CImage &operator=( const CImage  &oRHS);

  UInt &       operator()( UInt nRow, UInt nCol );
  const UInt operator()( UInt nRow, UInt nCol ) const;
  
  /////////////////////////
  //
  //  SetSize
  //   
  /////////////////////////
  void SetSize(UInt row, UInt col);

  /////////////////////////
  //
  //  Resize 
  //   
  /////////////////////////
  void Resize(UInt row, UInt col);
  
  /////////////////////////
  //
  //   DumpToFile
  //
  /////////////////////////
  void DumpToFile(const string & filename) const;

  /////////////////////////
  //
  //  ReadTiff
  //
  /////////////////////////
  bool ReadTiff(const string &filename);

  /////////////////////////
  //
  //  WriteTiff
  //
  /////////////////////////
  bool WriteTiff(const string &filename) const;
  
};

 
#endif
