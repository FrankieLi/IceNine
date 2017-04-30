/////////////////////
//
//   TiffReader.cpp
//
//   Modified from Micheal Still (mikal@stillh.com) 
//   http://www.ibm.com/developerworks/linux/library/l-libtiff
//
//
///////////////////
#include "TiffReader.h"



///////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
//  
//  Private:  CImage::DeepCopy
//
/////////////////////////////////////////////////////////////
void CImage::DeepCopy( const CImage & oImage)
{
  Resize( oImage.nNumRow, oImage.nNumCol );
  
  for ( UInt i = 0; i < nNumRow; i ++ )
    for( UInt j = 0; j < nNumCol; j ++ )
      pPixelValue[i][j] = oImage.pPixelValue[i][j];
  
  
}

/////////////////////////////////////////////////////////////
//  
//  Private:  CImage::Allocate
//
/////////////////////////////////////////////////////////////
void CImage::Allocate(UInt** &pImage, UInt nRow, UInt nCol)
{

  pImage = new UInt*[nRow];
  for(UInt i = 0; i < nRow; i ++)
  {
    pImage[i] = new UInt[nCol];
    for(UInt j = 0; j < nCol; j ++)
      pImage[i][j] = 0;    
  }
  
}
/////////////////////////////////////////////////////////////
//  
//  Private:  CImage::Allocate
//
/////////////////////////////////////////////////////////////
void CImage::Delete(UInt** &pImage, UInt nRow)
{
  if(pImage != NULL){
 
    for(UInt i = 0; i < nRow; i++)
      delete [] pImage[i];
    
    delete[] pImage;
    pImage = NULL;

  }

}

/////////////////////////////////////////////////////////////
//  
//  Public:  CImage::CImage - copy constructor
//
/////////////////////////////////////////////////////////////
CImage::CImage( const CImage & oImage ): pPixelValue(NULL), nNumRow( 0 ) , nNumCol( 0 )

{
  DeepCopy( oImage );
}

/////////////////////////////////////////////////////////////
//  
//  Public:  CImage::operator=
//
/////////////////////////////////////////////////////////////
CImage & CImage::operator=( const CImage  &oRHS)
{
  RUNTIME_WARNING( &oRHS != this, "ERROR!  Self assignment of CImage\n" );
  DeepCopy( oRHS );
  return *this;
}

/////////////////////////////////////////////////////////////
//  
//  Public:  CImage::CImage
//
/////////////////////////////////////////////////////////////
CImage::CImage():pPixelValue(NULL), nNumRow(0), nNumCol(0)
{
	
}


/////////////////////////////////////////////////////////////
//
//  Public:  CImage::~CImage
//
/////////////////////////////////////////////////////////////
CImage::~CImage()
{
  Delete( pPixelValue, nNumRow );
}



/////////////////////////////////////////////////////////////
//
//  Public:  CImage::SetSize
//
/////////////////////////////////////////////////////////////
void CImage::SetSize(UInt row, UInt col)
{
  
  nNumRow = row;
  nNumCol = col;

  pPixelValue = new UInt*[nNumRow];
  for(UInt i = 0; i < nNumRow; i ++){
    pPixelValue[i] = new UInt[nNumCol];
    for(UInt j = 0; j < nNumCol; j ++){
      pPixelValue[i][j] = 0;
    }
  }
	
	
}



/////////////////////////////////////////////////////////////
//
//  Public:  CImage::SetSize
//
/////////////////////////////////////////////////////////////
void CImage::Resize(UInt row, UInt col)
{
  if( nNumRow == row && nNumCol == col )
    return;

  
  UInt **pNewPixelValue;
  Allocate( pNewPixelValue, row, col );

  UInt nMaxRow = max( row, nNumRow );
  UInt nMaxCol = max( col, nNumCol );

  if( nNumRow > 0 && nNumCol > 0)  // no copy if originally there was nothing
    for ( UInt i = 0; i < nMaxRow; i ++ )
      for( UInt j = 0; j < nMaxCol; j ++ )
        pNewPixelValue[i][j] = pPixelValue[i][j];
  
  Delete( pPixelValue, nNumRow );
  pPixelValue = pNewPixelValue;

  nNumRow = row;
  nNumCol = col;
  
}


/////////////////////////////////////////////////////////////
//
//  Public:  CImage::DumpToFile
//
/////////////////////////////////////////////////////////////
void CImage::DumpToFile(const string & filename) const
{
  ofstream os;
  os.open(filename.c_str());

  for(UInt i = 0; i < nNumRow; i ++){
    for(UInt j = 0; j < nNumCol; j ++){
      os << pPixelValue[i][j] << " ";
    }
    os << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////////////
//
//  Public:  CImage::ReadTiff
//
//  Adopted from read.c from http://www-128.ibm.com/developerworks/linux/library/l-libtiff/
//
//  Note that this version of ReadTiff is not all that robust.  Can use some improvement
//
////////////////////////////////////////////////////////////////////////////////////////
bool CImage::ReadTiff(const string & filename)
{
  TIFF *tif;

  
  unsigned long int bufferSize;
  U8 *buffer;
	
	
  if( (tif = TIFFOpen(filename.c_str(), "r")) == NULL){
    cerr << "Cannot open " << filename << endl;
    return false;
  }

  oTags.nNumStrips = TIFFNumberOfStrips(tif);
  oTags.nStripSize = TIFFStripSize(tif);
  bufferSize = oTags.nNumStrips * oTags.nStripSize;
	

  // allocate memory for read (using malloc)
  buffer = new U8[ oTags.nNumStrips * oTags.nStripSize ];

  /*
  TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &(oTags.nImageWidth) ) ;
  TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &(oTags.nImageLength) ) ;
  TIFFGetField( tif, TIFFTAG_BITSPERSAMPLE, &(oTags.nBitsPerSample) ) ;
  TIFFGetField( tif, TIFFTAG_SAMPLESPERPIXEL, &(oTags.nSamplesPerPixel) );
  TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &(oTags.nRowsPerStrip) );
  
  TIFFGetField( tif, TIFFTAG_COMPRESSION, &(oTags.nCompression) );
  TIFFGetField( tif, TIFFTAG_PHOTOMETRIC, &(oTags.nPhotometric) );
  TIFFGetField( tif, TIFFTAG_FILLORDER, &(oTags.nFillOrder) );
  TIFFGetField( tif, TIFFTAG_PLANARCONFIG, &(oTags.nPlanarConfig) );
  
  //  TIFFGetField( tif, TIFFTAG_XRESOLUTION, &(oTags.fXResolution) ) ;
   TIFFGetField( tif, TIFFTAG_YRESOLUTION, &(oTags.fYResolution) );
   //TIFFGetField( tif, TIFFTAG_RESOLUTIONUNIT, &(oTags.nResolutionUnit) );

  */
  // some of the lines are commented out because error is returned for no good reason...
   
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &(oTags.nImageWidth) ) != 0, "ImageWidth Error\n" ) ;
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &(oTags.nImageLength) ) != 0, "ImageLength Error\n" );
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_BITSPERSAMPLE, &(oTags.nBitsPerSample) ) != 0, "Bits Per Sample Error\n" );
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_SAMPLESPERPIXEL, &(oTags.nSamplesPerPixel) ) != 0, "Samples Per Pixel Error\n");
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &(oTags.nRowsPerStrip) ) != 0, "Rows Per Pixel Error\n" );
  
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_COMPRESSION, &(oTags.nCompression) ) != 0, "Compression Error \n" );
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_PHOTOMETRIC, &(oTags.nPhotometric) ) != 0, "Photmetric Error\n" );
  // RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_FILLORDER, &(oTags.nFillOrder) ) != 0, "Fill Order Error \n" );
  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_PLANARCONFIG, &(oTags.nPlanarConfig) ) != 0, "Planar Configuration Error \n" );

  //  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_XRESOLUTION, &(oTags.fXResolution) ) != 0, "X Resolution Error \n" );
  //  RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_YRESOLUTION, &(oTags.fYResolution) ) != 0, "Y Resolution Error \n " );
  // RUNTIME_ASSERT( TIFFGetField( tif, TIFFTAG_RESOLUTIONUNIT, &(oTags.nResolutionUnit) ) != 0, "Resolution Unit Error\n" );
  
  nNumCol = oTags.nImageWidth;
  nNumRow = oTags.nImageLength;

  
  tsize_t imageOffset = 0;
  tsize_t byteRead;
  for (Int i = 0; i < oTags.nNumStrips; i++){
    if( (byteRead = TIFFReadEncodedStrip(tif, i, buffer + imageOffset,
                                         oTags.nStripSize)) == -1){
      cerr << "Read error on input strip number " << i << endl;
      return false;
    }
		
    imageOffset += byteRead;
  }
	
  
	
  // expect 2 byte (16 bit images), need to change this in the future
  if( ((oTags.nNumStrips * oTags.nStripSize) / (nNumRow * nNumCol)) != 2){
    cerr << "Error:  don't know how many bits are are in a pixel." << endl;
  }
	

  
  UInt bufPos = 0;
  U16 *uBuf = (U16*)buffer;
  SetSize(nNumRow, nNumCol);
  for(UInt i = 0; i < nNumRow; i ++){
    for(UInt j = 0; j < nNumCol; j ++){
      pPixelValue[i][j] = uBuf[bufPos];
      bufPos++;
    }
  }

  

  TIFFClose(tif);
  delete[] buffer;
	
  return true;
}

////////////////////////////////////////////////////////////
//
//  Public:  CIMage:: WriteTiff
//
//  Adopted from write.c of http://www-128.ibm.com/developerworks/linux/library/l-libtiff/
//
////////////////////////////////////////////////////////////
bool CImage::WriteTiff(const string &filename) const
{

  TIFF *tif;
	
  if( (tif = TIFFOpen(filename.c_str(), "w")) == NULL ){
    cerr << "Cannot open " << filename << endl;
    return false;
  }


  std::cout << "TIFF INFO: " << std::endl
            << "NumStrips " << oTags.nNumStrips << endl
            << "ImageWidth " << oTags.nImageWidth << endl 
            << "ImageLength " << oTags.nImageLength << endl
            << "SamplePerPixel " << oTags.nBitsPerSample << endl
            << "RowsPerStrip " << oTags.nRowsPerStrip << endl
            << "Compression " << oTags.nCompression << endl
            << "Photometic " << oTags.nPhotometric << endl
            << "FillOrder " << oTags.nFillOrder << endl
            << "PlanarConfig " << oTags.nPlanarConfig << endl
            << "XResolution " << oTags.fXResolution << endl
            << "YResolution " << oTags.fYResolution << endl
            << "ResolutionUnit " << oTags.nResolutionUnit << endl;
	
  // Basic tags - specific to our uses only
  TIFFSetField( tif, TIFFTAG_IMAGEWIDTH, oTags.nImageWidth );
  TIFFSetField( tif, TIFFTAG_IMAGELENGTH, oTags.nImageLength );
  TIFFSetField( tif, TIFFTAG_BITSPERSAMPLE, oTags.nBitsPerSample );
  TIFFSetField( tif, TIFFTAG_SAMPLESPERPIXEL, oTags.nSamplesPerPixel );
  TIFFSetField( tif, TIFFTAG_ROWSPERSTRIP, oTags.nRowsPerStrip );

  TIFFSetField( tif, TIFFTAG_COMPRESSION, oTags.nCompression );
  TIFFSetField( tif, TIFFTAG_PHOTOMETRIC, oTags.nPhotometric );
  TIFFSetField( tif, TIFFTAG_FILLORDER,  oTags.nFillOrder );
  TIFFSetField( tif, TIFFTAG_FILLORDER,  FILLORDER_LSB2MSB );
  TIFFSetField( tif, TIFFTAG_PLANARCONFIG, oTags.nPlanarConfig );

  TIFFSetField( tif, TIFFTAG_XRESOLUTION, oTags.fXResolution );
  TIFFSetField( tif, TIFFTAG_YRESOLUTION, oTags.fYResolution );
  TIFFSetField( tif, TIFFTAG_RESOLUTIONUNIT, oTags.nResolutionUnit );

  // Write the information to the file
  U16 *buffer;
  UInt bufPos = 0;

  
  // allocate memory for read (using malloc)
  buffer = new U16[ nNumRow * nNumCol ];

  
  // could really make this faster based on the way that pPixelValue is allocated
  for(UInt i = 0; i < nNumRow; i ++){
    for(UInt j = 0; j < nNumCol; j ++){
      buffer[bufPos] =  pPixelValue[i][j];
      bufPos++;
    }
  }



    
  tsize_t imageOffset = 0;
  tsize_t byteRead;
  U8* outBuf = (U8*) buffer;
  for ( Int i = 0; i < oTags.nNumStrips; i ++ )
  {    
    if ( (byteRead = TIFFWriteEncodedStrip( tif, i, outBuf + imageOffset , oTags.nStripSize)) == -1 )
    {
      cerr << "Write error on input strip number " << i << endl;
      return false;
    }
    
    imageOffset += byteRead;
  }
  
  
  TIFFClose(tif);
  delete[] buffer;

  
  return true;
}
  


/////////////////////////////////////////////////////////////
//
//  Private: CImage::IsBigEndian
//
/////////////////////////////////////////////////////////////
bool CImage::IsBigEndian(const string & filename)
{

  U8 buf[4];
  FILE* fd;

  if( (fd = fopen(filename.c_str(), "rb")) == NULL){
    cerr << "cannot open file " << filename.c_str() << endl;
    return false;
  }
  fread(buf, sizeof(U8), 4, fd);
  fclose(fd);	

  if(buf[3] == 42){
    return true;
  }else{
    return false;
  }
}
