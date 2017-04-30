//#include "GeneralPCH.h"
#include "UFF.h"

namespace GeneralLib
{
  //----------------------------------------------------------------------------------------------
  // Public : OpenForWrite
  //----------------------------------------------------------------------------------------------
  bool CUFF::OpenForWrite( const string &sFileName, bool bBinary )
  {		
    m_pFile = fopen( sFileName.c_str(), "wb" );
    m_bBinary = bBinary;

    // Write file header
    F32 fVersion = 1;
    if ( m_bBinary )
    {
      if ( fwrite( (void *)&fVersion, sizeof(F32), 1, m_pFile ) <= 0 )
        return false;
    }
    else
    {
      char sBlockString[32];
      sprintf( sBlockString, "Version: %.2f\n\n", fVersion );

      if ( fwrite( (void *)sBlockString, strlen( sBlockString ), 1, m_pFile ) <= 0 )
        return false;
    }

    DEBUG_ASSERT( m_pFile != NULL, "Unable to open file!" );
    return (m_pFile != NULL);
  }

  //----------------------------------------------------------------------------------------------
  // Public : OpenForRead
  //----------------------------------------------------------------------------------------------
  bool CUFF::OpenForRead( const string &sFileName )
  {
    // Reading a binary file
    if ( ( m_pFile = fopen( sFileName.c_str(), "rb" )) == NULL )
      return false;

    // Get the size of the file
    fseek( m_pFile, 0L, SEEK_END ); // Position to end of file
    m_nBufferSize = ftell( m_pFile );		// Get file length 
    rewind( m_pFile );				// Back to start of file

    // Read in the entire file and close the file handle
    m_pBuffer = new U8[m_nBufferSize + 1];
    m_nBufferOffset = 0;
    if ( fread( m_pBuffer, m_nBufferSize, 1, m_pFile ) <= 0 )
    {
      fclose( m_pFile );
      return false;
    }

    fclose( m_pFile );
    m_pFile = NULL;

    // Read version #
    //    F32 fVersion = *((F32 *)&m_pBuffer[m_nBufferOffset]);
    m_nBufferOffset += sizeof(F32);

    return true;
  }

  //----------------------------------------------------------------------------------------------
  // Public : CloseFile
  //----------------------------------------------------------------------------------------------
  bool CUFF::CloseFile()
  {
    if ( m_pFile )
    {
      int nReturnCode = fclose( m_pFile );
      DEBUG_ASSERT( nReturnCode == 0 );
      return ( nReturnCode == 0 );
    }

    return true;
  }

  //----------------------------------------------------------------------------------------------
  // Public : DeleteData
  //----------------------------------------------------------------------------------------------
  void CUFF::DeleteData()
  {
    if ( m_pBuffer )
    {
      delete [] m_pBuffer;
      m_pBuffer = NULL;
    }
  }

  //----------------------------------------------------------------------------------------------
  // Public : ReadBlock
  //----------------------------------------------------------------------------------------------
  bool CUFF::ReadBlock( SBlock &oBlock )
  {
    // Read the dummy data at the beginning of every block
    U32 nDummy = *((U32 *)&m_pBuffer[m_nBufferOffset]);
    m_nBufferOffset += sizeof(U32);

    if ( nDummy != 0xFEEDBEEF )
    {
      DEBUG_ASSERT( false, "Read error" );
      return false;
    }

    // Typecast the header
    SBlockHeader *pBlockHeader = (SBlockHeader *)&m_pBuffer[m_nBufferOffset];
    m_nBufferOffset += sizeof(SBlockHeader);

    // Fill in some block data from the header
    oBlock.m_nBlockType = pBlockHeader->m_nBlockType;
    oBlock.m_nDataFormat = pBlockHeader->m_nDataFormat;

    // pBlockHeader->m_nDataSize == (name size + child pointers size + block data size + child blocks sizes)
    // Thus, to get block data size, do some subtraction
    oBlock.m_nDataSize = pBlockHeader->m_nDataSize;
    oBlock.m_nDataSize -= pBlockHeader->m_nNameSize;
    oBlock.m_nDataSize -= pBlockHeader->m_nNumChildren * sizeof(U32);
	
    // Read the name
    if ( pBlockHeader->m_nNameSize > 0 )
    {
      oBlock.m_sName = string( (char *)&m_pBuffer[m_nBufferOffset], pBlockHeader->m_nNameSize - 1 );
      m_nBufferOffset += pBlockHeader->m_nNameSize;
    }

		
    // Determine where the data block starts (after the children pointers)
    U32 nDataRegion = m_nBufferOffset + pBlockHeader->m_nNumChildren * sizeof(U32);

    // Read the child pointers & and child data
    oBlock.m_pChildren.resize( pBlockHeader->m_nNumChildren );
    for ( U32 nChild = 0; nChild < pBlockHeader->m_nNumChildren; nChild++ )
    {
      U32 nChildOffset = *((U32 *)&m_pBuffer[m_nBufferOffset]);
		
      // Temporarily set the buffer pointer to where the child is
      U32 nBackup = m_nBufferOffset + sizeof(U32);
      m_nBufferOffset = nDataRegion + nChildOffset;

      if ( !ReadBlock( oBlock.m_pChildren[nChild] ) )
        return false;

      // Restore buffer pointer
      m_nBufferOffset = nBackup;
    }

    // By this point, the buffer pointer should point to beginning of data block
    DEBUG_ASSERT( nDataRegion == m_nBufferOffset );

    // Set the data pointer
    oBlock.m_pData = (U8 *)&m_pBuffer[nDataRegion];

    return true;
  }

  //----------------------------------------------------------------------------------------------
  // Public : WriteBlock
  //----------------------------------------------------------------------------------------------
  bool CUFF::WriteBlock( const SBlock &oBlock )
  {
    if ( m_bBinary )
      WriteBinaryBlock( oBlock );
    else
      WriteAsciiBlock( oBlock );

    return true; 
  }

  //----------------------------------------------------------------------------------------------
  // Public : GetBlockSize
  //----------------------------------------------------------------------------------------------
  U32 CUFF::GetBlockSize( const SBlock &oBlock )
  {
    U32 nSize = 0;

    // Dummy 0xFEEDBEEF data at the front
    nSize += sizeof(U32);

    // Header
    nSize += sizeof( SBlockHeader );

    nSize += GetDataSize( oBlock );

    return nSize;
  }

  //----------------------------------------------------------------------------------------------
  // Private : GetDataSize
  //----------------------------------------------------------------------------------------------
  U32 CUFF::GetDataSize( const SBlock &oBlock )
  {
    U32 nSize = 0;

    // Name string
    nSize += (U32)GetStringLength(oBlock.m_sName);

    // Child pointers
    nSize += (U32)oBlock.m_pChildren.size() * sizeof(U32);

    // Data for block
    nSize += oBlock.m_nDataSize;

    // Size of each child block
    for ( U32 nChild = 0; nChild < (U32)oBlock.m_pChildren.size(); nChild++ )
      nSize += GetBlockSize( oBlock.m_pChildren[nChild] );

    return nSize;
  }

  //----------------------------------------------------------------------------------------------
  // Private : WriteBinaryBlock
  //----------------------------------------------------------------------------------------------
  bool CUFF::WriteBinaryBlock( const SBlock &oBlock )
  {
    // Write dummy data at beginning of block for error checking purposes
    U32 nDummy = 0xFEEDBEEF;
    if ( fwrite( (void *)&nDummy, sizeof(U32), 1, m_pFile ) <= 0 )
      return false;

    SBlockHeader oBlockHeader;
    oBlockHeader.m_nBlockType = oBlock.m_nBlockType;
    oBlockHeader.m_nDataFormat = oBlock.m_nDataFormat;
    oBlockHeader.m_nNumChildren = (U32)oBlock.m_pChildren.size();
    DEBUG_ASSERT( (size_t)oBlockHeader.m_nNumChildren == oBlock.m_pChildren.size() );

    if ( oBlock.m_sName.size() > 0 )
    {
      // Include the terminating \0 character 
      oBlockHeader.m_nNameSize = (U32)GetStringLength(oBlock.m_sName);		
      DEBUG_ASSERT( (size_t)oBlockHeader.m_nNameSize == GetStringLength(oBlock.m_sName) );
    }
    else
      oBlockHeader.m_nNameSize = 0;

    // Size of data = name size + child pointers size + block data size + child blocks sizes
    oBlockHeader.m_nDataSize = GetDataSize( oBlock );

    // Write the header
    if ( fwrite( (void *)&oBlockHeader, sizeof(SBlockHeader), 1, m_pFile ) <= 0 )
      return false;

    // Write the name
    if ( oBlock.m_sName.size() > 0 )
    {
      if ( fwrite( (void *)oBlock.m_sName.c_str(), oBlockHeader.m_nNameSize, 1, m_pFile ) <= 0 )
        return false;
    }

    // Write the child pointers
    // Child pointers point relative to the beginning of the data block
    U32 nChildBlockStart = oBlock.m_nDataSize;		// Child blocks start after the data chunk
    for ( U32 nChild = 0; nChild < (U32)oBlock.m_pChildren.size(); nChild++ )
    {
      //nChildBlockStart++;
      if ( fwrite( (void *)&nChildBlockStart, sizeof(U32), 1, m_pFile ) <= 0 )
        return false;

      nChildBlockStart += GetBlockSize( oBlock.m_pChildren[nChild] );
    }

    // Write the data for this block
    if ( oBlock.m_nDataSize > 0 )
    {
      if ( fwrite( (void *)oBlock.m_pData, oBlock.m_nDataSize, 1, m_pFile ) <= 0 )
        return false;
    }

    // Write out each child block
    for ( U32 nChild = 0; nChild < (U32)oBlock.m_pChildren.size(); nChild++ )
    {
      if ( !WriteBinaryBlock( oBlock.m_pChildren[nChild] ) )
        return false;
    }

    return true;
  }

  //----------------------------------------------------------------------------------------------
  // Private : WriteAsciiBlock
  //----------------------------------------------------------------------------------------------
  bool CUFF::WriteAsciiBlock( const SBlock &oBlock, const string &sPerLineHeader )
  {
    // Write out header
    string sPrePend = sPerLineHeader;

    char sBlockString[256];

    if ( oBlock.m_sName.size() > 0 )
      sprintf( sBlockString, "\n%s%s (Type: %d) {\n", sPrePend.c_str(), oBlock.m_sName.c_str(), oBlock.m_nBlockType );
    else
      sprintf( sBlockString, "\n%s(Type: %d) {\n", sPrePend.c_str(), oBlock.m_nBlockType );

    if ( fwrite( (void *)sBlockString, GetStringLength( sBlockString ), 1, m_pFile ) <= 0 )
      return false;

    // Indent
    sPrePend += '\t';

    sprintf( sBlockString, "%sData Format = %d\n", sPrePend.c_str(), oBlock.m_nDataFormat );
    if ( fwrite( (void *)sBlockString, strlen( sBlockString ), 1, m_pFile ) <= 0 )
      return false;

    sprintf( sBlockString, "%sNumber of Children = %d\n", sPrePend.c_str(), static_cast<int>( oBlock.m_pChildren.size() ) );
    if ( fwrite( (void *)sBlockString, strlen( sBlockString ), 1, m_pFile ) <= 0 )
      return false;

    sprintf( sBlockString, "%sName Size = %d\n", sPrePend.c_str(), GetStringLength(oBlock.m_sName) );
    if ( fwrite( (void *)sBlockString, strlen( sBlockString ), 1, m_pFile ) <= 0 )
      return false;

    sprintf( sBlockString, "%sData Size = %d\n", sPrePend.c_str(), GetDataSize( oBlock ) );
    if ( fwrite( (void *)sBlockString, strlen( sBlockString ), 1, m_pFile ) <= 0 )
      return false;

    for ( U32 nChild = 0; nChild < (U32)oBlock.m_pChildren.size(); nChild++ )
    {
      if ( !WriteAsciiBlock( oBlock.m_pChildren[nChild], sPrePend ) )
        return false;
    }

    sprintf( sBlockString, "%s}\n", sPerLineHeader.c_str() );
    if ( fwrite( (void *)sBlockString, strlen( sBlockString ), 1, m_pFile ) <= 0 )
      return false;

    return true;
  }

  //----------------------------------------------------------------------------------------------
  // Private : GetStringLength
  //----------------------------------------------------------------------------------------------
  U32 CUFF::GetStringLength( const string & sChar )
  {
    if ( sChar.size() == 0 )
      return 0;
    
    return (U32)( sChar.size() + 1);
  }

}
