////////////////////////////////////////////////////////////////
//
//  File:    CrystalStructure.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Crystal structure and x-ray diffraction related classes
//
//
/////////////////////////////////////////////////////////////////

#ifndef _CRYSTAL_STRUCTURE_H
#define _CRYSTAL_STRUCTURE_H

#include "Types.h"
#include "3dMath.h"
#include <string>
#include <fstream>
#include "Debug.h"
#include <vector>
#include <math.h>
#include "Error.h"
#include "Symmetry.h"
#include "CrystalPrimitives.h"
#include "Serializer.h"

using namespace GeneralLib;
using namespace LatticeSymmetry;
using std::ostream;
using std::vector;



//----------------------------------
//
//  class CUnitCell
//
//----------------------------------
class CUnitCell
{

private:
  
  //----------------------------------
  // CalculateIntensity
  //  oReciprocalVector is the reciprocal lattice vectork, G_hkl
  // 
  //----------------------------------
  Float CalculateIntensity( const SVector3 & oReciprocalLatticeVector ) const;

  //----------------------------------
  // Calculate primative translation vectors
  // Right now this function must be called before CalculateRecipVectors
  //----------------------------------
  void CalculatePrimitiveVectors();

  //----------------------------------
  // Calculate reciprocal vectors
  //----------------------------------
  void CalculateRecipVectors();

  //----------------------------------
  // Generate the list of reflection vector list, save to oReflectionVector
  //----------------------------------
  void GenerateReflectionVectorList( Int nHMax, Int nKMax, Int nLMax, Float fMaxQ );
  
  //----------------------------------
  // Generate the list of reflection vector list, save to oReflectionVector
  //----------------------------------
  vector<CRecpVector> FilterByMinIntensity( const vector<CRecpVector> &oRecpList,
                                            Float fMinIntensity );

  //----------------------------------
  //  Calculate the metric tensor based on the angles given
  //----------------------------------
  void CalculateRecpMetricTensor( );
  
  // To determine the reciprocal vector list
  Int nRecipListMaxH, nRecipListMaxK, nRecipListMaxL;
  Float fRecipListfMaxQ;
  Float fRecipMinIntensity;

  Bool bLimitSet;
  Bool bLimitChanged;
  
  // basis vector in cartesian coordinates
	
  //  The variables below follow naming convention for
  //  crystal structures.
  
  
  // number of atoms
  UInt nNumAtoms;

  // max intenisty
  Float fMaxIntensity;

  // angles between basis vectors
  Float fAlpha;
  Float fBeta;
  Float fGamma;
  
  // length between unit cells
  Float fLengthA;
  Float fLengthB;
  Float fLengthC;

  SMatrix3x3 oRecpLatMetricTensor;    // metric tensor of the reciprocal lattice
  



  //----------------------------------
  // TODO:  Make oTranslationVectors and oReciprocalVectors private.
  // Make them only accessible via access functions.  Make insertion
  // a function as well.
  //
  // Lattice Translation Vector
  // Positions of atoms in the unit cell described by these translation
  // vectors.  They are expressed in terms of the primitive vectors
  //
  //  i.e., AtomicPositionVector = oTranslactionVectors[i].m_fX * oPrimitiveVector[0]
  //                             + oTranslactionVectors[i].m_fY * oPrimitiveVector[1]
  //                             + oTranslactionVectors[i].m_fZ * oPrimitiveVector[2]
  //----------------------------------
  vector<CAtom> oTranslationVector;
  
  // recipricol vectors and location of
  // Recipricol Vectors
  // This is given by the length of the unit cell
  // and the angles between the three vectors
  vector<SVector3> oReciprocalVector;
  
  //----------------------------------
  //  Simply a list of three primitive vectors in cartesian coordinates
  //  TODO:  Think about a metric instead.
  //----------------------------------
  vector<SVector3> oPrimitiveVector;
  
  //----------------------------------
  // recipricol vectors and location of
  // Recipricol Vectors
  // This is given by the length of the unit cell
  // and the angles between the three vectors
  //----------------------------------
  vector<CRecpVector> oReflectionVectorList;
  vector<CRecpVector> oUniqueReflectionVectorList;

  
public:
  
  //----------------------------------
  //  C O N S T R U C T O R S
  //----------------------------------
  CUnitCell();
  CUnitCell(Float a, Float b, Float c, 
            Float lengthA, Float lengthB, Float lengthC);
  CUnitCell(const CUnitCell &c);


  //----------------------------------
  //  M U T A T O R S
  //----------------------------------
  
  //----------------------------------
  //  SetUnitCellLength
  //
  //  Set length of unit cell
  //----------------------------------
  void SetUnitCellLength( Float _fLengthA, Float _fLengthB, Float _fLengthC );

  //----------------------------------
  //  SetUnitCellLength
  //
  //  Set angles between unit cell according
  //  to typical crystallographic convention
  //----------------------------------
  void SetUnitCellBasisAngles( Float _fAlpha, Float _fBeta, Float _fGamma );

  //----------------------------------
  //  SetNumAtoms
  //----------------------------------
  void SetNumAtoms( Int _nNumAtoms );   // TO DEPRECATE

  //----------------------------------
  //  AddAtom
  //  Add a new atom by inputting a new positon of effective
  //  Z number and translation vector.  
  //----------------------------------
  void AddAtom( Float fEffectiveZ, const SVector3 & oPos );
  
  //----------------------------------
  //  This function takes care of producing both the primitive and the
  //  reciprocal vectors
  //----------------------------------
  void InitializeCoordinateSystem();
  
  //----------------------------------
  // SetReflectionVectorLimits
  //----------------------------------
  void SetReflectionVectorLimits( Int nHMax, Int nKMax, Int nLMax, Float fMaxQ, Float fMinIntensity );
  
  //----------------------------------
  //   A C C E S S O R S
  //----------------------------------

  //----------------------------------
  // calculate the list of reflection vectors *up to* max hkl
  // !! in the future, do this by intensity
  //----------------------------------
  const vector<CRecpVector> & GetReflectionVectorList(Int nHMax, Int nKMax, Int nLMax,
                                                      Float fMaxQ, Float fMinIntensityFraction );
  
  //----------------------------------
  // return the list of reflection vectors - the list of vectors is calculated only if
  // SetReflectionVectorLimits is recently called
  // 
  // Must first call SetReflectionVectorLimits();
  //----------------------------------
  const vector<CRecpVector> &  GetReflectionVectorList( ) const;
  void  SetUniqueReflectionVectorList( const CSymmetry &oSym );
  const vector<CRecpVector> &  GetUniqueReflectionVectorList( ) const;
  const vector<SVector3>    &  GetReciprocalVectorList() const { return oReciprocalVector; }
  const SMatrix3x3          &  GetRecipocalMetric() const { return oRecpLatMetricTensor ; }
  
  const SVector3 & GetPrimitiveVector( int n ) const
  {
    DEBUG_ASSERT( n >= 0 && n < 3, "n out of range for PrimitiveVector \n" );
    return oPrimitiveVector[n];
  }
  
  CRecpVector GetReciprocalVector( int h, int k, int l ) const;

  //----------------------------------
  //  fMinIntensity is the minimum fraction of the max - i.e., it can only between [0, 1]
  //----------------------------------
  void WriteReflectionVectorList(const string &filename, Float fMinIntensity ) const;


  //----------------------------------
  //  Save and Restore
  //----------------------------------
  Bool Save   ( CSerializer   & oSerialBuf ) const;
  Bool Restore( CDeserializer & oSerialBuf);

  //----------------------------------
  //  GetCompactSize
  //----------------------------------
  Size_Type GetCompactSize();
};

#endif
