#ifndef VTK_FILE_IO_H
#define VTK_FILE_IO_H

#include <iostream>
#include <fstream>
#include <string>
namespace VtkFilesIO
{

  template< typename Iterator >
  void WriteCellPoints( Iterator pIter );
}

#endif
