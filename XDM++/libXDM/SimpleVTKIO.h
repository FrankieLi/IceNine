//------------------------------------------------------------------------
//
//  File:     SimpleVTKIO.h
//  Purpose:  VERY simple VTK input/output operations
//  Author:   Frankie Li
//  email:    sfli@cmu.edu
//
//
//------------------------------------------------------------------------


#ifndef _SIMPLE_VTK_IO_
#define _SIMPLE_VTK_IO_

#include "Voxel.h"
#include "3dMath.h"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>



#include <tuple>

namespace Visualization
{
  using std::vector;
  using std::string;
  using namespace GeneralLib;
  //-----------------------------------------------------------------------------------
  //  WriteVTKFile
  //
  //-----------------------------------------------------------------------------------
  void WriteVTKMicFile( const vector<SVoxel> & vAllVoxels,
                        const string & oFilename,
                        const string & sScalarFieldName );
  
  //-----------------------------------------------------------------------------------
  //  Vtk Point File
  //-----------------------------------------------------------------------------------
  void WriteVTKPointFile( const vector< SVector3 > & vAllPoints,
                          const string & oFilename,
                          const string & sScalarFieldName,
                          const vector< Float > & oScalarField );
}



#endif
