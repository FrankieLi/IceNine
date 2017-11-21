//----------------------------------------------------
//  XDM_Mesh_IO.h
//  Purpose:  Read in "dx" type file, or an image file
//            to the internal data format.
//
//----------------------------------------------------


#ifndef XDM_Mesh_IO_H
#define XDM_Mesh_IO_H
#include <string>
#include <fstream>
namespace XDMCGal
{
  //-----------------------------------------------------------------------------
  //
  //  ReadDxFile  
  //
  //-----------------------------------------------------------------------------
  template< class ImageData >
  void ReadDxFile( ImageData & oData, std::string oFilename,
                   int nX, int nY, int nZ,
                   float fX, float fY, float fZ )
  {
    using std::ifstream;
    std::ifstream is;
    
    is.open( oFilename.c_str() );
    if( !is.good() )
    {
      std::cerr << "File not opened " << std::endl;
      std::cerr << oFilename << std::endl;
      exit(0);
    }

    oData.SetPixelSize( fX, fY, fZ );
    oData.Resize( nX, nY, nZ );
    for( int i = 0; i < nX; i ++ )
      for( int j = 0; j < nY; j ++ )
        for( int k = 0; k < nZ; k ++ )
          oData.Set( i, j, k, -1 );
  
    for( int k = 1; k < nZ-1; k ++ )
    {
      for( int j = 1; j < nY-1; j ++ )
      {
      
        for( int i = 1; i < nX-1; i ++ )     
        {
          int nTmp;
          is >> std::skipws >> nTmp;
          oData.Set( i, j , k, nTmp );
        }
      }
    }
    is.close();
  }
}

#endif
