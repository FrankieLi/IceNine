//------------------------------------------------------------------------
//
//  File:     SimpleVTKIO.cpp
//  Purpose:  VERY simple VTK input/output operations
//  Author:   Frankie Li
//  email:    sfli@cmu.edu
//
//
//------------------------------------------------------------------------



#include "SimpleVTKIO.h"


namespace Visualization
{
  //-----------------------------------------------------------------------------------
  //  Vtk Point File
  //-----------------------------------------------------------------------------------
  void WriteVTKPointFile( const vector< SVector3 > & vAllPoints,
                          const string & oFilename,
                          const string & sScalarFieldName,
                          const vector< Float > & oScalarField )
  {
    using std::ofstream;
    
    ofstream oOutFile;
    oOutFile.open( oFilename.c_str() );
    oOutFile << "# vtk DataFile Version 2.0" << std::endl;
    oOutFile << oFilename.c_str() << std::endl;
    oOutFile << "ASCII" << std::endl;
    oOutFile << "DATASET POLYDATA" << std::endl;
    oOutFile << "POINTS " << vAllPoints.size()  << " float" << std::endl; 
  
    for( Size_Type i = 0; i < vAllPoints.size(); i ++ )
      oOutFile << vAllPoints[i] << std::endl;
  
    oOutFile << "CELL_DATA " << vAllPoints.size() << std::endl;
    oOutFile << "SCALARS GrainID int " << 1 << std::endl;
    oOutFile << "LOOKUP_TABLE default" << std::endl;
  
    for( Size_Type i = 0; i < oScalarField.size(); i ++ )
      oOutFile << oScalarField[i] << std::endl;
    oOutFile.close();
  }

  //-----------------------------------------------------------------------------------
  //  WriteVTKFile
  //
  //-----------------------------------------------------------------------------------
  void WriteVTKMicFile( const vector<SVoxel> & vAllVoxels,
                        const string & oFilename,
                        const string & sScalarFieldName )
  {
    using std::ofstream;

    ofstream oOutFile;
    oOutFile.open( oFilename.c_str() );
    oOutFile << "# vtk DataFile Version 2.0" << std::endl;
    oOutFile << oFilename.c_str() << std::endl;
    oOutFile << "ASCII" << std::endl;
    oOutFile << "DATASET POLYDATA" << std::endl;
    oOutFile << "POINTS " << vAllVoxels.size() * 3 << " float" << std::endl; 

    std::map<Int, Int> oUniqueIDMap;
    Int nCurID = 0;
    for( Size_Type i = 0; i < vAllVoxels.size(); i ++ )
    {
      for( Int n = 0; n < 3; n ++ )
        oOutFile << vAllVoxels[i].pVertex[n] * Float( 1000 )  + SVector3( 500, 500, 0 )<< std::endl;  // DEBUG added * 1000 
      Bool bSuccess;
      std::map<Int, Int>::iterator pCur;
      boost::tie( pCur, bSuccess )
        = oUniqueIDMap.insert( std::pair<Int, Int>( vAllVoxels[i].nGrainID,  nCurID ) );
      if( bSuccess )
        nCurID ++;
    }

    vector<Int> vNewIDs( oUniqueIDMap.size() );
    std::map<Int, Int>::iterator pCur = oUniqueIDMap.begin();
  
    for( Size_Type i = 0; i < vNewIDs.size(); i ++ )
      vNewIDs[i] = i;
  
    std::random_shuffle( vNewIDs.begin(), vNewIDs.end() );
    std::cout << "CurrentID " << nCurID << std::endl;
    for( ; pCur != oUniqueIDMap.end(); ++ pCur )
      pCur->second = vNewIDs[ pCur->second ];
    
    std::cout << vAllVoxels.size() << std::endl;
  
    oOutFile << "POLYGONS " << vAllVoxels.size() << " " << 4 * vAllVoxels.size() << std::endl; 
    for( Size_Type i = 0; i < vAllVoxels.size(); i ++ )
    {
      oOutFile << 3 << " "
               << 3 * i     << " "
               << 3 * i + 1 << " "
               << 3 * i + 2 << " "
               << std::endl;
    }
    oOutFile << "CELL_DATA " << vAllVoxels.size() << std::endl;
    oOutFile << "SCALARS GrainID int " << 1 << std::endl;
    oOutFile << "LOOKUP_TABLE default" << std::endl;

    std::cout << vAllVoxels.size() << std::endl;
    for( Size_Type i = 0; i < vAllVoxels.size(); i ++ )
    {
      oOutFile << oUniqueIDMap[vAllVoxels[i].nGrainID] << std::endl;
    }
    std::cout << vAllVoxels.size() << std::endl;
    oOutFile << "SCALARS " << sScalarFieldName << " float " << 1 << std::endl;
    oOutFile << "LOOKUP_TABLE default" << std::endl;

    for( Size_Type i = 0; i < vAllVoxels.size(); i ++ )
      oOutFile << vAllVoxels[i].fConfidence << std::endl;
  
  }


  
}
