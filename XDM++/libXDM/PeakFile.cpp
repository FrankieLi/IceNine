/////////////////////////////////////////////////////////////////
//
//  File:    PeakFile.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Purpose: Implementation of binary file format for peak files
//
//
//
//
/////////////////////////////////////////////////////////////////
#include "PeakFile.h"



namespace GeneralLib
{

  //----------------------------------------------------------------
  //
  //  Read
  //
  //
  //----------------------------------------------------------------
  bool CInputPeakFile::Read( const string & sFilename )
  {
    bool bSuccess = oDataFile.OpenForRead( sFilename );

    if( !bSuccess )
      return false;
    
    bSuccess = oDataFile.ReadBlock( oPeakDataBlock );
   
    if( !bSuccess )              // error checking for block data
      return false;
    
    bSuccess = oDataFile.CloseFile();
    nNumPeaks  = ( (U32*)oPeakDataBlock.m_pData )[0];
    nNumPixels = oPeakDataBlock.m_pChildren[0].m_nDataSize / sizeof( U16 );
    return bSuccess;
  }

  //----------------------------------------------------------------
  //  Return true if getting next pixel is possible
  //  Return false if this is past last pixel already (no more pixel)
  //----------------------------------------------------------------
  bool CInputPeakFile::GetNextPixel( Pixel & oRes, Int & nPeakID )
  {
    if( nCurPixel >= nNumPixels )
      return false;
    
    const vector<SBlock> &pDataBlocks = oPeakDataBlock.m_pChildren;
    oRes.x          = static_cast<Int>  ( ( (U16*)pDataBlocks[0].m_pData )[ nCurPixel ] );
    oRes.y          = static_cast<Int>  ( ( (U16*)pDataBlocks[1].m_pData )[ nCurPixel ] );
    oRes.fIntensity = static_cast<Float>( ( (F32*)pDataBlocks[2].m_pData )[ nCurPixel ] );
    nPeakID         = static_cast<Int>  ( ( (U16*)pDataBlocks[3].m_pData )[ nCurPixel ] );
    nCurPixel ++;
    return true;
  }

  //----------------
  //  InsertPixel
  //----------------
  bool COutputPeakFile::InsertPixels( const vector< vector<Pixel> > & oPeakList )
  {
    vector<int> IDMap;
    IDMap.resize( oPeakList.size() );
    for( int i = 0; i < oPeakList.size(); i ++ )
      IDMap[i] = i;
    return InsertPixels( oPeakList, IDMap );
  }
  
  //----------------------------------------------------------------
  //  InsertPixel
  //
  //----------------------------------------------------------------
  bool COutputPeakFile::InsertPixels( const vector< vector<Pixel> > & oPeakList,
                                      const vector<Int> & oIDMap )
  {
    DEBUG_ASSERT( nNumPixels == 0, "InsertPixels: Multiple insert not yet supported!\n" );

    if( nNumPixels != 0 )
      return false;
    
    for( Size_Type i = 0; i < oPeakList.size(); i ++ )
      nNumPixels += oPeakList[i].size();
    
    pPixelX    = new U16[ nNumPixels ];
    pPixelY    = new U16[ nNumPixels ];
    pIntensity = new F32[ nNumPixels ];
    pPixelID   = new U16[ nNumPixels ];    // this is really the ID of the pixel, which is peak ID

    std::set<Int> UniqueIDSet;
    Int n = 0;
    for( Size_Type i = 0; i < oPeakList.size(); i ++ )
    {
      for( Size_Type j = 0; j < oPeakList[i].size(); j ++ )
      {
        pPixelX[ n ]    = oPeakList[i][j].x;
        pPixelY[ n ]    = oPeakList[i][j].y;
        pIntensity[ n ] = oPeakList[i][j].fIntensity;
        pPixelID[ n ]   = oIDMap[ i ] ;
        UniqueIDSet.insert( pPixelID[n] );
        n ++;
      }
    }
    nNumPeaks = UniqueIDSet.size();
    return true;
  }

  //----------------------------------------------------------------
  // Write
  //----------------------------------------------------------------
  bool COutputPeakFile::Write( const string & sFilename )
  {
    
    static const char *pBlockNames[5] = {
      "PixelCoord0\0",
      "PixelCoord1\0",
      "Intensity\0",
      "PeakID\0",
      "PeakFile\0"
    };
        
    bool bSuccess = oDataFile.OpenForWrite( sFilename, true );
    
    if( !bSuccess )
      return false;
    
    oPeakDataBlock.m_nBlockType  = 1;   // TODO -- use enums to name these guys, specify endian-ness
    oPeakDataBlock.m_nDataFormat = 1;
    oPeakDataBlock.m_sName = pBlockNames[0];
    oPeakDataBlock.m_pData = NULL;

    vector<SBlock> & vChildrenBlocks = oPeakDataBlock.m_pChildren;
    vChildrenBlocks.resize(4);
    for( Size_Type i = 0; i < 4; i ++ )
    {
      vChildrenBlocks[i].m_nBlockType  = oPeakDataBlock.m_nBlockType;
      vChildrenBlocks[i].m_nDataFormat = oPeakDataBlock.m_nDataFormat;
      vChildrenBlocks[i].m_sName = pBlockNames[i + 1];
    }
    
    vChildrenBlocks[0].m_nDataSize = nNumPixels * sizeof( U16 );
    vChildrenBlocks[1].m_nDataSize = nNumPixels * sizeof( U16 );
    vChildrenBlocks[2].m_nDataSize = nNumPixels * sizeof( F32 );
    vChildrenBlocks[3].m_nDataSize = nNumPixels * sizeof( U16 );

    vChildrenBlocks[0].m_pData     = reinterpret_cast< U8* >( pPixelX );
    vChildrenBlocks[1].m_pData     = reinterpret_cast< U8* >( pPixelY );
    vChildrenBlocks[2].m_pData     = reinterpret_cast< U8* >( pIntensity );
    vChildrenBlocks[3].m_pData     = reinterpret_cast< U8* >( pPixelID );
    
    Int nChildrenSize = 0;
    for( Size_Type i = 0; i < 4; i ++ )
      nChildrenSize += oDataFile.GetBlockSize( vChildrenBlocks[i] );

    oPeakDataBlock.m_nDataSize = sizeof( U32 );
    oPeakDataBlock.m_pData = (U8*)(&nNumPeaks);        // saving number of peaks
    bSuccess = oDataFile.WriteBlock( oPeakDataBlock );
   
    if( !bSuccess )              // error checking for block data
      return false;
    
    bSuccess = oDataFile.CloseFile();
    return bSuccess;
  }
  
}
