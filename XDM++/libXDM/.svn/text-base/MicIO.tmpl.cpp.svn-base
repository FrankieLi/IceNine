////////////////////////////////////////////////////////////////////////
//
//
//   File:  MicIO.cpp
//   E-mail:  sfli@cmu.edu
//   Author:  Frankie Li
//   Input output operations for mic files.
// 
////////////////////////////////////////////////////////////////////////



//---------------------------------------------------------
//
//  Private:
//
//
//
//---------------------------------------------------------
template< class VoxelType >
void MicFileBase<VoxelType>::ClearBuffer()
{
  if(pBuffer)
    delete[] pBuffer;
	
  pBuffer = NULL;
  nBufferSize = 0;
}

//---------------------------------------------------------
//
//  Private:
//  LeftMostVoxel
//
//
//---------------------------------------------------------
template< class VoxelType >
Int MicFileBase<VoxelType>::GetLeftMostVertex(const VoxelType & v ) const
{
  Float fMinX = MAX_FLOAT;
  Int nLeftMost = 0;
  for ( int i = 0; i < 3; i ++ )
  {
    if( v.pVertex[i].m_fX < fMinX)
    {
      fMinX = v.pVertex[i].m_fX;
      nLeftMost = i;
    }
  }
  return nLeftMost;
}

//---------------------------------------------------------
//
//  Private:
//  ReadFileToBuf
//
//
//---------------------------------------------------------
template< class VoxelType >
bool MicFileBase<VoxelType>::ReadFileToBuf(const string &filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  FILE *pFile;
	
	
  if ( (pFile = fopen(filename.c_str(), "r" )) == NULL ){
    cerr << "[CMic::ReadFileToBuf]Cannot Open File: " 
         << filename.c_str()  << endl;
    return false;
  }
	
  // Get the size of the file
  fseek( pFile, 0L, SEEK_END ); // Position to end of file
  nBufferSize = ftell( pFile );		// Get file length 
  rewind( pFile );				// Back to start of file
	
  // Read in the entire file and close the file handle
  pBuffer = new char[nBufferSize + 1];
	
  if ( fread( pBuffer, nBufferSize, 1, pFile ) <= 0 ){
    fclose( pFile );
    cerr << "[MicFileBase::ReadFileToBuf]: File read error " << filename << endl;
    return false;
  }
	
  fclose( pFile );	
		
  pBuffer[ nBufferSize ] = '\0'; // NULL termination of string
	
  return true;
}

//---------------------------------------------------------
//
//
//  Default Destructor
//
//
//---------------------------------------------------------
template< class VoxelType >
MicFileBase<VoxelType>::~MicFileBase()
{
  if(pBuffer)
    delete[] pBuffer;
}



//---------------------------------------------------------
//
//  Public:
//  SetNumCol
//
//
//---------------------------------------------------------
template< class VoxelType >
void MicFileBase<VoxelType>::SetNumCol( Int newNumCols )
{
  nNumCols = newNumCols;
}

//---------------------------------------------------------
//  SetVoxelList
//---------------------------------------------------------
template< class VoxelType>
void MicFileBase<VoxelType>::SetVoxelList(const vector<VoxelType> & oVoxels)
{
  oVoxelList = oVoxels;
}

//-----------------------------------
//  SetID
//  Set ID of voxel, starting from nStartID,
//  and return the next ID to be used.
//-----------------------------------
template< class VoxelType >
Int MicFileBase<VoxelType>::SetID( Int nStartID )
{
  Int nCurID = nStartID;
  for( Size_Type i = 0; i < oVoxelList.size(); i ++ )
  {
    oVoxelList[i].nID = nCurID;
    nCurID ++;
  }
  return nCurID;
}
//---------------------------------------------------------
//  AddVoxel
//---------------------------------------------------------
template< class VoxelType >
void MicFileBase<VoxelType>::AddVoxel(const VoxelType & oVoxel)
{
  oVoxelList.push_back( oVoxel );
}

//---------------------------------------------------------
//  SetInitialSideLength
//---------------------------------------------------------
template< class VoxelType >
void MicFileBase<VoxelType>::SetInitialSideLength(Float fSideLength)
{
  fInitialSideLength = fSideLength;
}

//---------------------------------------------------------
//  GetMinSideLength
//---------------------------------------------------------
template< class VoxelType >
Float MicFileBase<VoxelType>::GetMinSideLength() const
{
  Float fMinSideLength = std::numeric_limits<Float>::max();

  for( Size_Type i = 0; i < oVoxelList.size(); i ++ )
    fMinSideLength = std::min( fMinSideLength,  oVoxelList[i].fSideLength ) ;

  return fMinSideLength;
}


//---------------------------------------------------------
//
//  Accessors
//
//---------------------------------------------------------

//---------------------------------------------------------
//
//
//  GetInitialSideLength
//
//
//---------------------------------------------------------
template< class VoxelType >
Float MicFileBase<VoxelType>::GetInitialSideLength() const
{
  return fInitialSideLength;
}

//---------------------------------------------------------
//
//
// GetNuMVoxels
//
//
//---------------------------------------------------------
template< class VoxelType >
Float MicFileBase<VoxelType>::GetNumVoxels() const
{
  return oVoxelList.size();
}

//---------------------------------------------------------
// Clear
//---------------------------------------------------------
template< class VoxelType >
void MicFileBase<VoxelType>::Clear()
{
  oVoxelList.clear();
}


//---------------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------------



//---------------------------------------------------------
//
//  Public:
//  Read
//
//
//---------------------------------------------------------
//      MIC file format
//
// line 1: a0, the fundamental lattice constant for the triangular mesh
//
// others:
//
// columns 1 -3: xyz location of left vertex of a triangle
// column 4: 1/2 for triangle type -- 1 for upward pointing or 2 for downward
// pointing
// column 5: generation number -- triangle side length = a0/(2^generation),
// generation = 0,1,2,...
// column 6: 0/1 for phase -- 0 mean it wasn't (or hasn't yet) fitted to an
// orientation, 1 means it was
// columns 7 -9: Eulers in degrees
// column 10: confidence parameter: fraction of simulated Bragg peaks that hit
// experimental peaks
//---------------------------------------------------------
// template<>
// bool MicFile<SVoxel>::Read( const string &filename )
// {
//   vector< vector<string> > vsTokens;

//   if( !ReadFileToBuf(filename))
//   {
//     return false;
//   }

//   GeneralLib::Tokenize( vsTokens, string( pBuffer, nBufferSize ), " \t\n");

	
//   // we expect side length for the first line
	
//   if(vsTokens[0].size() > 1){
//     cerr << "[MicFileBase::ReadMicFileBase] Incorrect Mic File Format " << filename << endl;
//     exit(0);
//   }
	
//   fInitialSideLength = atof(vsTokens[0][0].c_str());
//   for( Size_Type i = 1; i < vsTokens.size(); i ++){
    
//     if(vsTokens[i].size() < 9 )
//     {
//       cerr << "[MicFileBase::ReadMicFile] Error: Incorrect number of columns, line: "
//            << i  << filename << endl;
//     }
//     else
//     {
//       DEBUG_ALERT( vsTokens[i].size() >= 10, "MicFileBase::Read:  Warning!  Old format MIC file with only 9 columns\n" );
//       VoxelType v;

//       v.nGeneration = atoi( vsTokens[i][4].c_str() );
//       v.fSideLength = fInitialSideLength / pow( 2,  v.nGeneration ) ;
//       v.nPhase = atoi( vsTokens[i][5].c_str() );

//       Float fX, fY, fZ;

//       fX = atof( vsTokens[i][0].c_str() );
//       fY = atof( vsTokens[i][1].c_str() );
//       fZ = atof( vsTokens[i][2].c_str() );

//       if ( atoi( vsTokens[i][3].c_str() ) == UP  ) 
//       {
//         // Winding order, counter clockwise
//         v.pVertex[0].Set( fX,                 fY, fZ );
//         v.pVertex[1].Set( fX + v.fSideLength, fY, fZ );
//         v.pVertex[2].Set( fX + v.fSideLength / (Float) 2.0 , fY + v.fSideLength / (Float) 2.0 * sqrt( (Float) 3.0 ) , fZ);
//         v.bVoxelPointsUp = true;
//       }
//       else
//       {
//         // Winding order, counter clockwise
//         v.pVertex[0].Set( fX,                 fY, fZ );
//         v.pVertex[1].Set( fX + v.fSideLength / (Float) 2.0 , fY - v.fSideLength / (Float) 2.0 * sqrt( (Float) 3.0 ) , fZ);
//         v.pVertex[2].Set( fX + v.fSideLength, fY, fZ );
//         v.bVoxelPointsUp = false;
//       }

//       SVector3 oOrientation;
//       oOrientation.Set(  DEGREE_TO_RADIAN( atof(vsTokens[i][6].c_str()) ),
//                          DEGREE_TO_RADIAN( atof(vsTokens[i][7].c_str()) ),
//                          DEGREE_TO_RADIAN( atof(vsTokens[i][8].c_str()) ) );
      
//       // note that orientation matrix assigned here - the choice of active vs. passive
//       // transform is made by the file format

//       v.oOrientMatrix.BuildActiveEulerMatrix( oOrientation.m_fX,
//                                               oOrientation.m_fY,
//                                               oOrientation.m_fZ );
      
//       if( vsTokens[i].size() > 9 )
//         v.fConfidence = atof( vsTokens[i][9].c_str() );
//       else
//         v.fConfidence = 0;
      
//       if( vsTokens[i].size() == 19 )   //  Defining our strain tensor to be ( I)
//       {
//         v.fCost                = atof( vsTokens[i][10].c_str() );
//         v.fPixelOverlapRatio   = atof( vsTokens[i][11].c_str() );
//         v.oDeformation.m[0][0] = atof( vsTokens[i][13].c_str() );
//         v.oDeformation.m[1][1] = atof( vsTokens[i][14].c_str() );
//         v.oDeformation.m[2][2] = atof( vsTokens[i][15].c_str() );
        
//         v.oDeformation.m[0][1] = v.oDeformation.m[1][0] = atof( vsTokens[i][16].c_str() );
//         v.oDeformation.m[1][2] = v.oDeformation.m[2][1] = atof( vsTokens[i][17].c_str() );
//         v.oDeformation.m[0][2] = v.oDeformation.m[2][0] = atof( vsTokens[i][18].c_str() );
//       }
//       else
//       {
//         v.oDeformation.SetIdentity();
//       }
      
//       oVoxelList.push_back(v);
//     }	

//   }

//   ClearBuffer();  // clear buffer after use
//   return true;
// }

//---------------------------------------------------------
//
//  Public: Write
//  
//  Action:  Write MicFileBase to a file.
//
//---------------------------------------------------------
// template<>
// bool MicFile<SVoxel>::Write(const string &filename) const
// {
//   std::ofstream osOutFile;
//   osOutFile.open( filename.c_str() );

//   if( !osOutFile.good() )
//     return false;

//   osOutFile << fInitialSideLength << std::endl;
  
//   // columns 1 -3: xyz location of left vertex of a triangle
//   // column 4: 1/2 for triangle type -- 1 for upward pointing or 2 for downward
//   // pointing
//   // column 5: generation number -- triangle side length = a0/(2**generation),
//   // generation = 0,1,2,...
//   // column 6: 0/1 for phase -- 0 mean it wasn't (or hasn't yet) fitted to an
//   // orientation, 1 means it was
//   // columns 7 -9: Eulers in degrees
//   // column 10: confidence parameter: fraction of simulated Bragg peaks that hit
//   // experimental peaks
  
//   for( Size_Type i = 0; i < oVoxelList.size(); i ++ )
//   {
//     Int nDir;
//     Int nLeftMostVertex;

//     //  Need to find the left most voxel
//     //  to follow MIC files specification.
//     nLeftMostVertex = GetLeftMostVertex( oVoxelList[i] );

//     if ( oVoxelList[i].bVoxelPointsUp ){
//       nDir = UP;
//     }else{
//       nDir = DOWN;
//     }
    
//     osOutFile << oVoxelList[i].pVertex[nLeftMostVertex].m_fX << " "
//               << oVoxelList[i].pVertex[nLeftMostVertex].m_fY  << " "
//               << oVoxelList[i].pVertex[nLeftMostVertex].m_fZ << " "

//               << nDir << " "
//               << oVoxelList[i].nGeneration << " "
//               << oVoxelList[i].nPhase << " "
//               << RadianToDegree( oVoxelList[i].oOrientMatrix.GetEulerAngles() ) << " "
//               << oVoxelList[i].fConfidence
      
//       //---------------------------------------------------------
//       // newly added - to be put into property map
//       //---------------------------------------------------------
//               << " " << oVoxelList[i].fCost
//               << " " << oVoxelList[i].fPixelOverlapRatio
//               << " " << oVoxelList[i].fFittingTime
      
//               << " " << oVoxelList[i].oDeformation.m[0][0]
//               << " " << oVoxelList[i].oDeformation.m[1][1]
//               << " " << oVoxelList[i].oDeformation.m[2][2]
//               << " " << oVoxelList[i].oDeformation.m[1][0]
//               << " " << oVoxelList[i].oDeformation.m[1][2]
//               << " " << oVoxelList[i].oDeformation.m[2][0]
      
//               << std::endl;;
//   }

//   osOutFile.close();

//   return true;
// }
