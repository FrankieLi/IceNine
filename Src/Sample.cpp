//============================================================================== 
// Copyright (c) 2014, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by S. F. Li (li31@llnl.gov)
// LLNL-CODE-657639
// All rights reserved.
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the Lawrence Livermore National Lab nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL LAB BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//============================================================================== 

//------------------------------------------------------------------------------------
//  Author:  S. F. Li (Frankie)
//  e-mail:  li31@llnl.gov; sfli@cmu.edu 
//------------------------------------------------------------------------------------
////////////////////////////////////////////////////////////
//
//  Filename:  Sample.cpp
//  Author:    Frankie Li
//
//
//
////////////////////////////////////////////////////////////

#include "Sample.h"

//----------------------------------------------------------------------
//
//
//  Public: CSample::~CSample
//
//----------------------------------------------------------------------
CSample::CSample(): v3SampleLoc(0, 0, 0),  v3SampleOrientation(0, 0, 0),
                    voCrystalStructList(), voRecipVectors()
                    
{
  oSampleToLabMatrix.SetIdentity();
}

//----------------------------------------------------------------------
//
//  Public:  CSample::SetLocation
//  It only sets location to specific values (without doing anything else)
//
//----------------------------------------------------------------------
void CSample::SetLocation(const SVector3 &loc)
{
  v3SampleLoc = loc;
  oSampleToLabMatrix.SetTranslation(-loc.m_fX, -loc.m_fY, -loc.m_fZ);
}

//----------------------------------------------------------------------
//
//  Public: CSample::SetOrientation
//  It only set orientations to specific values (without doing anything else)
//
//----------------------------------------------------------------------
void CSample::SetOrientation(Float phi, Float theta, Float psi)
{
  v3SampleOrientation.Set(phi, theta, psi);	
  oSampleToLabMatrix.SetPassiveEulerMatrix( phi, theta, psi );  // A passive rotation is used to go from the
                                                                // sample back to the lab. 
}

//----------------------------------------------------------------------
//
//   SetOrientation
//   Only sets the orientation part
//----------------------------------------------------------------------
void CSample::SetOrientation(const SMatrix3x3 &oOrientation)
{
  oSampleToLabMatrix.m[0][0] = oOrientation.m[0][0];
  oSampleToLabMatrix.m[1][0] = oOrientation.m[1][0];
  oSampleToLabMatrix.m[2][0] = oOrientation.m[2][0];

  oSampleToLabMatrix.m[0][1] = oOrientation.m[0][1];
  oSampleToLabMatrix.m[1][1] = oOrientation.m[1][1];
  oSampleToLabMatrix.m[2][1] = oOrientation.m[2][1];

  oSampleToLabMatrix.m[0][2] = oOrientation.m[0][2];
  oSampleToLabMatrix.m[1][2] = oOrientation.m[1][2];
  oSampleToLabMatrix.m[2][2] = oOrientation.m[2][2];
}

//----------------------------------------------------------------------
//  SetSampleSymmetry
//----------------------------------------------------------------------
void CSample::SetSampleSymmetry( ESymmetryT eSym )
{
  eSampleSymmetry = eSym;
}
  
//----------------------------------------------------------------------
//
// Public:  CSample::Rotate
//
// Note:  We're using active rotation to go rotate the sample, and therefore anything that goes from the
// lab frame to the sample frame uses an active rotation
//
//----------------------------------------------------------------------
void CSample::Rotate(Float phi, Float theta, Float psi)
{
  SMatrix4x4 tmpMatrix;
  tmpMatrix.BuildActiveEulerMatrix(phi, theta, psi);
  oSampleToLabMatrix = tmpMatrix * oSampleToLabMatrix;
}

//----------------------------------------------------------------------
//
// Public:  CSample::Rotate ( oAxis, fOmega )
//
//
//----------------------------------------------------------------------
void CSample::Rotate(const SVector3 & oAxis, Float fOmega )
{
  SMatrix4x4 oRotMatrix;
  oRotMatrix.BuildRotationAboutAxis( oAxis, fOmega );
  oSampleToLabMatrix = oRotMatrix * oSampleToLabMatrix ;      // A ACTIVE rotation is used to go from the
                                                              // sample back to the lab.  
}

//----------------------------------------------------------------------
//
//  Rotate Z  -- a specialized function for speed -- rotate about Z 
//  axis counter clockwise
//  
//  Profiler shows that the original function uses approximately 10% of the total time.
//  That's why it is hand optimized.
//----------------------------------------------------------------------
void CSample::RotateZ(Float fOmega )
{
  Float fCosOmega = cos( fOmega );
  Float fSinOmega = sin( fOmega );
  
  Float m00 = fCosOmega;
  Float m01 = -fSinOmega;
  Float m10 = fSinOmega;
  Float m11 = fCosOmega;
  
  Float mTmp[2][3];
  
  // A ACTIVE rotation is used to go from the
  // sample back to the lab.  

  mTmp[0][0] = oSampleToLabMatrix.m[0][0];
  mTmp[1][0] = oSampleToLabMatrix.m[1][0];
  mTmp[1][1] = oSampleToLabMatrix.m[1][1];
  mTmp[0][1] = oSampleToLabMatrix.m[0][1];
  mTmp[1][2] = oSampleToLabMatrix.m[1][2];
  mTmp[0][2] = oSampleToLabMatrix.m[0][2];
  
  oSampleToLabMatrix.m[0][0] = m00 * mTmp[0][0] + m01 * mTmp[1][0];
  oSampleToLabMatrix.m[0][1] = m00 * mTmp[0][1] + m01 * mTmp[1][1];
  oSampleToLabMatrix.m[0][2] = m00 * mTmp[0][2] + m01 * mTmp[1][2];
  
  oSampleToLabMatrix.m[1][0] = m10 * mTmp[0][0] + m11 * mTmp[1][0];
  oSampleToLabMatrix.m[1][1] = m10 * mTmp[0][1] + m11 * mTmp[1][1];
  oSampleToLabMatrix.m[1][2] = m10 * mTmp[0][2] + m11 * mTmp[1][2];

}
//----------------------------------------------------------------------
//
// Public:  CSample:Translate
// 
//
//----------------------------------------------------------------------
void CSample::Translate(const SVector3 &loc)
{
  SMatrix4x4 tmpMatrix;
  v3SampleLoc += loc;
  tmpMatrix.BuildTranslation( loc.m_fX, loc.m_fY, loc.m_fZ );    // Building active translate matrix
  oSampleToLabMatrix = oSampleToLabMatrix * tmpMatrix;   // A ACTIVE rotation is used to go from the
                                                         // sample back to the lab. 
}

//----------------------------------------------------------------------
//
// Public:  CSample::LoadSample
//
//  TODO:  Generalize this  - there are more ways to read a sample than Mic file
//----------------------------------------------------------------------
bool CSample::LoadSample( string filename, int GridType )
{
  pSampleMic = MicIOFactory::Create( GridType );
  return pSampleMic->Read( filename );
}

//----------------------------------------------------------------------
//
//  CSample::AddCrystalStructure
//
//
//----------------------------------------------------------------------
void CSample::AddCrystalStructure(const CUnitCell &c)
{
  voCrystalStructList.push_back(c);
}

//-----------------------------------------------------------------
//
//  ToLabFrame:  This is a coordinate transformation function
//  Output:      oRes is v, represented in the Lab coordinates
//
//-----------------------------------------------------------------
void CSample::ToLabFrame(SVector4 &oRes, const SVector4 &v) const
{
  oRes = v;
  oRes.Transform( oSampleToLabMatrix );
}

//-----------------------------------------------------------------
//  ToLabFrame
// 
//
//-----------------------------------------------------------------
SVector3 CSample::ToLabFrame( const SVector3 &v ) const
{
  SVector3 oRes = v;  
  oRes.Transform( oSampleToLabMatrix );
  return oRes;
}


//----------------------------------------------------------------------
//
// Accessors
//
//
//----------------------------------------------------------------------

//-----------------------------------------------------------------
//  GetLocation
//----------------------------------------------------------------- 
SVector3 CSample::GetLocation( ) const
{ 
  SVector3 oRet = v3SampleLoc;
  return oRet;
}

//-----------------------------------------------------------------
//  GetOrientation
//-----------------------------------------------------------------    
SVector3 CSample::GetOrientation( ) const
{
  SVector3 oRet = v3SampleOrientation;
  return oRet;
}
	
//-----------------------------------------------------------------
//  GetOrientation
//-----------------------------------------------------------------    
SMatrix3x3 CSample::GetOrientationMatrix () const
{
  SMatrix3x3 oRes;
  for( Int i = 0; i < 3; i ++ )
    for( Int j = 0; j < 3; j ++ )
      oRes.m[i][j] = oSampleToLabMatrix.m[i][j];
  
  return oRes;
}

//-----------------------------------------------------------------
//  GetSampleSymmetry
//
//  WARNING - not using nPhase.  Will NEED TO FIX THIS!
//
//-----------------------------------------------------------------    
LatticeSymmetry::CSymmetry *
CSample::GetSampleSymmetry( int nPhase ) const
{
  CSymmetry * pSym = NULL;
  switch( GetSampleSymmetry() )
  {
    case LatticeSymmetry::eCubic:
      pSym =  & ( LatticeSymmetry::CCubicSymmetry::Get() );
      break;
    case LatticeSymmetry::eHexagonal:
      pSym =  & ( LatticeSymmetry::CHexagonalSymmetry::Get() );
      break;
    case LatticeSymmetry::eTetragonal:
      pSym =  & ( LatticeSymmetry::CTetragonalSymmetry::Get() );
      break;
    default:
      cerr << "WARNING!!!!  Not using any symmetry!!" << std::endl;
      break;
  }
  return pSym;
}

//-----------------------------------------------------------------
//  Save
//-----------------------------------------------------------------    
Bool CSample::Save   ( CSerializer   & oSerialBuf ) const
{
  Bool bSuccess = true;
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( v3SampleLoc );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( v3SampleOrientation );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( oSampleToLabMatrix );

  bSuccess = bSuccess && pSampleMic->Save( oSerialBuf );
  
  bSuccess = bSuccess && oSerialBuf.InsertComplexVector( voCrystalStructList );
  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( voRecipVectors );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( eSampleSymmetry );
  return bSuccess;
}

//-----------------------------------------------------------------
//  Restore
//-----------------------------------------------------------------    
Bool CSample::Restore( CDeserializer & oSerialBuf)
{
  Bool bSuccess = true;
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & v3SampleLoc );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & v3SampleOrientation );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & oSampleToLabMatrix );

  bSuccess = bSuccess && pSampleMic->Restore( oSerialBuf );
  
  bSuccess = bSuccess && oSerialBuf.GetComplexVector( voCrystalStructList );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( voRecipVectors );
  
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( &eSampleSymmetry );
  return bSuccess;
}
