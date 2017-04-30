/////////////////////////////////////////////////////////////////
//
//  File:    XDMFileTypes.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Purpose: Definition of tags and names for all XDM file types
//
/////////////////////////////////////////////////////////////////


namespace GeneralLib
{

  // Binary file types defined used
  // in XDM++
  //
  enum EBinaryFileTypes{

    eBinaryDetectorPeaks,     // peaks and pixel information
    eBinaryDetectorPixels,    // purely pixels, no peak information

    eNumBlockTypes
  };

  enum EDataFormat
	{
      eU16_1,		// 16-bit unsigned integer
      eU16_2,		// 2 16-bit unsigned integers
      eU32_1,		// 32-bit unsigned integer
      eF32_1,		// 32-bit float
      eF32_2,		// 2 32-bit floats
      eF32_3,		// 3 32-bit floats
      eF32_4,		// 4 32-bit floats
      eColorARGB,	// 4D packed unsigned bytes mapped to 0. to 1. range
	};
};
