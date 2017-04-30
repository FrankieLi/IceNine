/********************************************************************
	created:	2007/03/08
	author:		Shan-Min Chao
	
	purpose:	Code that provides a binary format that can be used to
				store mesh data, texture data, configuration data... 
				any sort of data.
				The file format stores data in "blocks", where each
				block can contain children blocks.  Each block contains
				a header followed by binary data.
*********************************************************************/

#ifndef _UFF_H_
#define _UFF_H_
#include "Types.h"
#include "Debug.h"
#include <vector>
#include <string>
#include <cstdio>

using std::vector;
using std::string;

namespace GeneralLib
{
  ////////////////////////////////////////////////////////////////////////////////////////////////

  // The application reads/writes from/to this block structure, but this is not the exact format
  // of the binary data
  struct SBlock
  {
    U32 m_nBlockType;
    U32 m_nDataFormat;
    string m_sName;
    U8 *m_pData;
    U32 m_nDataSize;

    vector<SBlock> m_pChildren;

    // Constructor
    SBlock() : m_sName(), m_pData( NULL ) 
    {
      m_nBlockType = m_nDataFormat = m_nDataSize = (U32)-1;
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////

  class CUFF
  {
  public:

    CUFF(): m_pBuffer( NULL ), m_nBufferOffset(0), m_nBufferSize(0), m_pFile( NULL ){}
    ~CUFF() { DeleteData(); }
    // Open a file for writing.  If bBinary is set to false, then writes out a ascii summary
    // of the blocks for debugging purposes
    bool OpenForWrite( const string &sFileName, bool bBinary );

    // Only binary reading supported
    bool OpenForRead( const string &sFileName );

    // Finish reading/writing
    bool CloseFile();

    // Delete read-in data
    void DeleteData();

    // Read a block (which includes all its data and its children's data)
    bool ReadBlock( SBlock &oBlock );
    bool WriteBlock( const SBlock &oBlock );

    // Get the total block size as it will appear on disk (includes headers and all children)
    U32 GetBlockSize( const SBlock &oBlock );

  private:

    // Each binary chunk is preceded by this block header
    struct SBlockHeader
    {
      U16 m_nBlockType;
      U16 m_nDataFormat;

      U16 m_nNumChildren;
      U16 m_nNameSize;

      U32 m_nDataSize;

      // Data can be composed of multiple blocks if it is too big.
      // Each chunk can contain only a set amount of bytes.  This is so that loading can be done in chunks.
      U16 m_nChunkNum;
      U16 m_nTotalChunks;
    };

    // Get the data block size used for m_nDataSize
    U32 GetDataSize( const SBlock &oBlock );

    // Write out the block in binary (includes header, binary data, and all children)
    bool WriteBinaryBlock( const SBlock &oBlock );

    // Writes out the block header for debugging purposes.  sPerLineHeader is preprended to each line and
    // mainly used to prepend \t characters for formatting purposes
    bool WriteAsciiBlock( const SBlock &oBlock, const string &sPerLineHeader = "" );

    // Same as "strlen() + 1" (to include the terminating string character) except checks for NULL pointer,
    // in which case it returns 0
    U32 GetStringLength( const string & sChar );

    // Whether binary mode read is on
    bool m_bBinary;

    // For read operations, read the entire file at once for speed purposes
    U8 *m_pBuffer;
    U32 m_nBufferOffset;	// Current read pointer into m_pBuffer
    U32 m_nBufferSize;		// Size of buffer read in

    FILE *m_pFile;
  };
}

#endif
