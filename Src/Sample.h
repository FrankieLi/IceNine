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
//  Filename:  Sample.h
//  Author:    Frankie Li
//
//
//  Purpose:   Abstraction for sample.  It should contain the discretized
//             subdivision of the sample (pixel or voxelization).  It also
//             contains the sample geometry and orientations.  
//
////////////////////////////////////////////////////////////

#ifndef _SAMPLE_H_
#define _SAMPLE_H_

//#include "MicMesh.h"
#include "Voxel.h"
#include "Types.h"
#include "3dMath.h"
#include "MicIO.h"
#include "MathDebug.h"
#include "CrystalStructure.h"
#include "Serializer.h"
#include "Symmetry.h"
#include <memory>

using namespace GeneralLib;

//-----------------------------------------------------------------------
//  
//  Class:  CSample
//  Description:  Contains the data for the orientations of the sample,
//                as well as the sample orientation, location, and 
//                size in the lab coordinates
//
//-----------------------------------------------------------------------
class CSample
{
public:
  typedef std::shared_ptr< MicIOBase > MicIOPtrT;
private:
  
  // Sample orientation and location in Lab coordinates
  SVector3           v3SampleLoc;
  SVector3           v3SampleOrientation;
  SMatrix4x4         oSampleToLabMatrix;
  
  // Sample material property
  vector<CUnitCell>  voCrystalStructList;
  vector<SVector3>   voRecipVectors;
  MicIOPtrT          pSampleMic;
  // CMic               oSampleMic; 

  LatticeSymmetry::ESymmetryT  eSampleSymmetry;
  
  //----------------------------------
  //  In the future, use a more adaptive data structure,
  //  such as the Quadtree.
  //----------------------------------
public:	

  //----------------------------------
  //  Consturctors
  //----------------------------------
  CSample();

  //----------------------------------
  // Sets the location of the sample (fills out only the translation part of the matrix)
  //----------------------------------
  void SetLocation(const SVector3 &loc);

  //----------------------------------
  // Set orientation as euler rotation
  //----------------------------------
  void SetOrientation(Float phi, Float theta, Float psi);
  void SetOrientation(const SMatrix3x3 &oOrientation);

  //----------------------------------
  //  SetSampleSymmetry
  //----------------------------------
  void SetSampleSymmetry( LatticeSymmetry::ESymmetryT eSym );
  
  //----------------------------------
  // Rotates the sample (multiply by a rotation matrix)
  //----------------------------------
  void Rotate(Float phi, Float theta, Float psi);

  //----------------------------------
  // Rotate( oAxis, omega ) - axis rotation according to right hand rule
  //----------------------------------
  void Rotate(const SVector3 & oAxis, Float fOmega );

  //----------------------------------
  // Rotate( oAxis, omega ) - axis rotation according to right hand rule
  //----------------------------------
  void RotateZ( Float fOmega );
  //----------------------------------
  // Translates the sample (muliply by a translation matrix)
  //----------------------------------
  void Translate(const SVector3 &loc);

  //----------------------------------
  // Read the sample from the filename
  //----------------------------------
  bool LoadSample( string filename, int GridType );

  //----------------------------------
  // ToLabFrame:  This is a coordinate transformation function
  //----------------------------------
  void ToLabFrame(SVector4 &oRes, const SVector4 &v) const;

  //----------------------------------
  // ToLabFrame  (The 3 Vector version, i.e., no translation)
  //----------------------------------
  SVector3 ToLabFrame( const SVector3 &v ) const;

  //----------------------------------
  // Add crystal structure to the sample
  //----------------------------------
  void AddCrystalStructure( const CUnitCell &c );

  //-----------------------------------------
  //  This is a hack.. need to separate MicIO from Sample.  Sample
  //  should only hold voxel iterators.  (maybe... not sure)  HACK
  //-----------------------------------------
  void InitializeMic( int GridType )
  {
    pSampleMic = MicIOFactory::Create( GridType );
  }
  
  //-----------------------------------------------------------------
  //  Accessors
  //-----------------------------------------------------------------    
  SVector3                  GetLocation            () const;
  SVector3                  GetOrientation         () const;
  SMatrix3x3                GetOrientationMatrix   () const;
  const vector<CUnitCell> & GetStructureList       () const { return voCrystalStructList; }
  LatticeSymmetry::ESymmetryT  GetSampleSymmetry    () const { return eSampleSymmetry; }

  //----------------------------------
  //  Not fully implemented yet
  //----------------------------------
  LatticeSymmetry::CSymmetry * GetSampleSymmetry    ( int nPhase ) const;
  
  const MicIOPtrT                GetMic        () const { return pSampleMic; }
  MicIOPtrT                      GetMic        ()       { return pSampleMic; }
//   const CMic &                GetMic        () const { return oSampleMic; }
//   CMic &                     GetMic        ()       { return oSampleMic; }


  Bool Save   ( CSerializer   & oSerialBuf ) const;
  Bool Restore( CDeserializer & oSerialBuf);
};


#endif
