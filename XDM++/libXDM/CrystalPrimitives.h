////////////////////////////////////////////////////////////////
//
//  File:    CrystalPrimitives.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Definition of primitives, such as reciprocal vectors
//  and atomic positions to be used in crystal structures.
//
/////////////////////////////////////////////////////////////////


#ifndef _CRYSTAL_PRIMITIVES_H_
#define _CRYSTAL_PRIMITIVES_H_

#include "Types.h"
#include "3dMath.h"
using namespace GeneralLib;


//----------------------------------
// translation vectors
//----------------------------------
class CAtom
{
public:
  Float fEffectiveZ;
  SVector3 v;
};

//----------------------------------
//
//  CRecpVector
//
//----------------------------------
class CRecpVector
{
public:
  SVector3 oDir;
  SVector3 v;
  Float fMag;
  Float d;  // distance of separation between planes
  Float fIntensity;
  Int h, k, l;
};


#endif
