// fine1: C++'s own split of a grid, for the grid validity check (no reconstruction is run).
//
//   c++ -O1 -std=gnu++14 -DBOOST_ALLOW_DEPRECATED_HEADERS -Wno-deprecated-declarations \
//       -I<repo>/XDM++/libXDM -I<boost include> -isystem <eigen3 include> \
//       cpp_split_check.cpp <repo>/CMakeFiles/IceNine.dir/Src/ReconstructionStrategies.cpp.o \
//       <repo>/XDM++/libXDM/libXDM++.a -o cpp_split_check
//   cpp_split_check <grid.mic> <MinSideLength> <out.mic>
//
// Reads the grid with CMic::Read, calls ReconstructionStrategies::Utilities::SetMinimumResolution
// (the very call the reconstruction makes, ReconstructionStrategies.tmpl.cpp:52) and writes the
// regridded voxel list with CMic::Write. fine1_cpp_split_check.py compares that file with the
// grid it made, within 1e-6 mm.
#include <iostream>
#include <string>
#include <cstdlib>
#include "MicIO.h"
#include "ReconstructionStrategies.h"

int main( int argc, char ** argv )
{
  if( argc != 4 )
  {
    std::cerr << "usage: cpp_split_check grid.mic MinSideLength out.mic" << std::endl;
    return 2;
  }
  CMic oMic;
  if( !oMic.Read( argv[1] ) )
  {
    std::cerr << "cannot read " << argv[1] << std::endl;
    return 3;
  }
  std::cout << "read " << oMic.GetNumVoxels() << " voxels, initial side "
            << oMic.GetInitialSideLength() << std::endl;
  ReconstructionStrategies::Utilities::SetMinimumResolution( oMic, static_cast<Float>( atof( argv[2] ) ) );
  std::cout << "after " << oMic.GetNumVoxels() << " voxels" << std::endl;
  oMic.Write( argv[3] );
  return 0;
}
