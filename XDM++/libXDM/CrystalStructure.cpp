////////////////////////////////////////////////////////////////
//
//  File:    CrystalStructure.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//  Crystal structure and x-ray diffraction related classes
//
//
/////////////////////////////////////////////////////////////////


#include "CrystalStructure.h"

#include <iostream>

using namespace std;


//------------------------------------------------------------
//   Default constructor
//------------------------------------------------------------
CUnitCell::CUnitCell():
  bLimitSet(false), bLimitChanged(false),
  nNumAtoms(0),fMaxIntensity(0),
  fAlpha(-1), fBeta(-1), fGamma(-1), fLengthA(-1),  
  fLengthB(-1), fLengthC(-1), oTranslationVector(),
  oReciprocalVector()
{
}

//------------------------------------------------------------
//   Constructor
//------------------------------------------------------------
CUnitCell::CUnitCell(Float a, Float b, Float c, 
                     Float lengthA, Float lengthB, Float lengthC):
  bLimitSet(false), bLimitChanged(false),
  nNumAtoms(0),fMaxIntensity(0),
  fAlpha(a), fBeta(b), fGamma(c), fLengthA(lengthA), 
  fLengthB(lengthB), fLengthC(lengthC), oTranslationVector(),
  oReciprocalVector()
{
}

//------------------------------------------------------------
//   Copy Constructor
//------------------------------------------------------------
CUnitCell::CUnitCell(const CUnitCell &c)
{
   *this = c;
}

//------------------------------------------------------------
//  CalculateRecpMetricTensor
//------------------------------------------------------------
void CUnitCell::CalculateRecpMetricTensor()
{
  oRecpLatMetricTensor.m[0][0] = fLengthA * fLengthA;
  oRecpLatMetricTensor.m[1][1] = fLengthB * fLengthB;
  oRecpLatMetricTensor.m[2][2] = fLengthC * fLengthC;
  
  oRecpLatMetricTensor.m[1][0] = oRecpLatMetricTensor.m[0][1] = fLengthA * fLengthB * cos( fGamma );
  oRecpLatMetricTensor.m[2][0] = oRecpLatMetricTensor.m[0][2] = fLengthA * fLengthC * cos( fBeta  );
  oRecpLatMetricTensor.m[2][1] = oRecpLatMetricTensor.m[1][2] = fLengthB * fLengthC * cos( fAlpha  );  
}

//------------------------------------------------------------
//  SetUnitCellLength
//
//  Set length of unit cell
//------------------------------------------------------------
void CUnitCell::SetUnitCellLength( Float _fLengthA, Float _fLengthB, Float _fLengthC )
{
  fLengthA = _fLengthA;
  fLengthB = _fLengthB;
  fLengthC = _fLengthC;
}

//------------------------------------------------------------
//  SetUnitCellLength
//
//  Set angles between unit cell according
//  to typical crystallographic convention
//------------------------------------------------------------
void CUnitCell::SetUnitCellBasisAngles( Float _fAlpha, Float _fBeta, Float _fGamma )
{
  fAlpha = _fAlpha;
  fBeta  = _fBeta;
  fGamma = _fGamma;
}

//------------------------------------------------------------
//  SetNumAtoms
//------------------------------------------------------------
void CUnitCell::SetNumAtoms( Int _nNumAtoms )   // TO DEPRECATE
{
  nNumAtoms = _nNumAtoms;
}

//------------------------------------------------------------
//  AddAtom
//  Add a new atom by inputting a new positon of effective
//  Z number and translation vector.  
//------------------------------------------------------------
void CUnitCell::AddAtom( Float _fEffectiveZ, const SVector3 & _oPos )
{
  CAtom oAtom;
  oAtom.v = _oPos;
  oAtom.fEffectiveZ = _fEffectiveZ;
  oTranslationVector.push_back( oAtom );
}     

//------------------------------------------------------------
//  GetReciprocalVector
//
//
//------------------------------------------------------------
CRecpVector CUnitCell::GetReciprocalVector( int h, int k, int l ) const
{
  SVector3 oHDir, oKDir, oLDir;
  
  oHDir = oReciprocalVector[0] * h;
  oKDir = oReciprocalVector[1] * k;
  oLDir = oReciprocalVector[2] * l;

  CRecpVector oRefRecipVector;
  oRefRecipVector.v = oHDir + oKDir + oLDir;
  oRefRecipVector.fMag= oRefRecipVector.v.GetLength(); 
  
  
  oRefRecipVector.h = h;
  oRefRecipVector.k = k;
  oRefRecipVector.l = l;
  
  oRefRecipVector.oDir = oRefRecipVector.v; 
  oRefRecipVector.oDir.Normalize();
  
  SVector3 oHKL( h, k, l );
  oRefRecipVector.d = Dot( oHKL, oRecpLatMetricTensor * oHKL );
  oRefRecipVector.d = sqrt( oRefRecipVector.d );
  oRefRecipVector.fIntensity = CalculateIntensity( oRefRecipVector.v );

  return oRefRecipVector;
}

//------------------------------------------------------------
//
//  Private:  CalculatePrimitiveVectors
//
// Generate the list of reflection vector list, save to oReflectionVectorList
//------------------------------------------------------------
void CUnitCell::GenerateReflectionVectorList( Int nHMax, Int nKMax, Int nLMax, Float fMaxQ )
{
   
  oReflectionVectorList.clear();
  fMaxIntensity = 0;
  for(Int h = -nHMax; h <= nHMax; h ++)
  {
    for(Int k = -nKMax; k <= nKMax; k ++)
    {
      for(Int l = -nLMax; l <= nLMax; l++)
      {
        // one of them has to be nonzero.
        if( abs(k) > 0 || abs(h) > 0 || abs(l) > 0 )
        {
          CRecpVector oRefRecipVector = GetReciprocalVector( h, k, l );
          if ( oRefRecipVector.v.GetLength() < fMaxQ )
          {
            if( oRefRecipVector.fIntensity > fMaxIntensity )
              fMaxIntensity = oRefRecipVector.fIntensity;
            oReflectionVectorList.push_back( oRefRecipVector );
          }

//           std::cout << "---------------------------- " << std::endl;
//           SVector3 TestVec = oRecpLatMetricTensor * oHKL;
//           std::cout << h << " " << k << " " << l << " "
//                     << oRefRecipVector.v << std::endl
//                     << TestVec << std::endl;
//           std::cout << "---------------------------- " << std::endl;
        }
      }
    }
  }

  //  std::cout << "---------------------------- " << std::endl;
  //  std::cout << GetReciprocalVector( 1, 0, 0 ).v << std::endl;
  //  std::cout << GetReciprocalVector( 0, 1, 0 ).v << std::endl;
  //  std::cout << GetReciprocalVector( 0, 0, 1 ).v << std::endl;
  //  std::cout << "---------------------------- " << std::endl;
}

//------------------------------------------------------------
//
//  FilterByMinIntensity
//
//------------------------------------------------------------
vector<CRecpVector> CUnitCell::FilterByMinIntensity( const vector<CRecpVector> &oRecpList, Float fMinIntensity )
{
  vector<CRecpVector> oReducedList;
  
  for( Size_Type i = 0; i < oRecpList.size(); i ++ )
  {
    if ( oRecpList[i].fIntensity > fMinIntensity )
      oReducedList.push_back( oRecpList[i] );
  }

  return oReducedList;
}
//------------------------------------------------------------
//
//  Private:  CalculatePrimitiveVectors
//
// Calculate primative translation vectors
// Right now this function must be called before CalculateRecipVector
//
//------------------------------------------------------------
void CUnitCell::CalculatePrimitiveVectors()
{

  RUNTIME_ASSERT( oPrimitiveVector.size() == 0, 
                  "[CUnitCell::CalculatePrimitiveVectors]: Attempt to reinitialiez Primitive vectors already initialized once!\n" );

  oPrimitiveVector.resize( 3 );
  
  oPrimitiveVector[0].m_fX = 1;
  oPrimitiveVector[0].m_fY = 0;
  oPrimitiveVector[0].m_fZ = 0;
	
  oPrimitiveVector[1].m_fX = cos(fGamma);
  oPrimitiveVector[1].m_fY = sin(fGamma);
  oPrimitiveVector[1].m_fZ = 0;
  
  oPrimitiveVector[2].m_fX = cos(fBeta);
  oPrimitiveVector[2].m_fY = (cos(fAlpha) - cos(fBeta) * cos(fGamma))/sin(fGamma);
  oPrimitiveVector[2].m_fZ = sqrt(
                                  1. - 
                                    oPrimitiveVector[2].m_fX * oPrimitiveVector[2].m_fX
                                  - oPrimitiveVector[2].m_fY * oPrimitiveVector[2].m_fY
                                  );

  
  oPrimitiveVector[0] *= fLengthA;
  oPrimitiveVector[1] *= fLengthB;
  oPrimitiveVector[2] *= fLengthC;


  //  std::cout << "Primitive vectors ----------------------------- " << std::endl;
  //  std::cout << oPrimitiveVector[0] << " " << fLengthA << std::endl
  //            << oPrimitiveVector[1] << " " << fLengthB << std::endl
  //            << oPrimitiveVector[2] << " " << fLengthC << std::endl
  //            << std::endl; 
}

//------------------------------------------------------------
//
//  Private:  CalculateRecipVectors
//  
//------------------------------------------------------------
void CUnitCell::CalculateRecipVectors()
{
  RUNTIME_ASSERT( oReciprocalVector.size() == 0, 
                  "[CUnitCell::CalculateRecipVectors]: Attempt to reinitialize Reciprocal vectors already initialized once!\n" );

  oReciprocalVector.resize( 3 );
  
  SVector3 a1_x_a2 = Cross( oPrimitiveVector[1], oPrimitiveVector[2] );
  
  // 2 pi (a1 x a2) / (a0 dot (a1 x a2) )
  oReciprocalVector[0] = (Float) 2 * PI * a1_x_a2 / ( Dot( oPrimitiveVector[0], a1_x_a2 ) );

  SVector3 a2_x_a0 = Cross( oPrimitiveVector[2], oPrimitiveVector[0] ); 

  // 2 pi (a2 x a0) / (a1 dot (a2 x a0) )
  oReciprocalVector[1] = (Float) 2 * PI * a2_x_a0 / ( Dot( oPrimitiveVector[1], a2_x_a0 ) );

  SVector3 a0_x_a1 = Cross( oPrimitiveVector[0], oPrimitiveVector[1] ); 

  // 2 pi (a0 x a1) / (a2 dot (a0 x a1) )
  oReciprocalVector[2] = (Float) 2 * PI * a0_x_a1 / ( Dot( oPrimitiveVector[2], a0_x_a1 ) );


  //  std::cout << "Reciprocal vectors ----------------------------- " << std::endl;
  //  std::cout << oReciprocalVector[0] << std::endl
  //            << oReciprocalVector[1] << std::endl
  //            << oReciprocalVector[2] << std::endl;  
}

//------------------------------------------------------------
//
//  Public:  InitializeCoordinateSystem
//
//------------------------------------------------------------
void CUnitCell::InitializeCoordinateSystem()
{
  CalculatePrimitiveVectors();
  CalculateRecipVectors();
  CalculateRecpMetricTensor();
}

//------------------------------------------------------------
// SetReflectionVectorLimits
//------------------------------------------------------------
void CUnitCell::SetReflectionVectorLimits( Int nHMax, Int nKMax, Int nLMax, Float fMaxQ, Float fMinIntensityFraction )
{
  nRecipListMaxH = nHMax;
  nRecipListMaxK = nKMax;
  nRecipListMaxL = nLMax;
  fRecipListfMaxQ = fMaxQ;
  
  GenerateReflectionVectorList( nHMax, nKMax, nLMax, fMaxQ );  // need to initialize max intensity
  RUNTIME_ASSERT( oReciprocalVector.size() == 3, 
                  "[CUnitCell::GetReflectionVectorList]: Recipricol Vector not initialied\n");
  fRecipMinIntensity = fMaxIntensity * fMinIntensityFraction;
  oReflectionVectorList = FilterByMinIntensity( oReflectionVectorList, fRecipMinIntensity );
}


//------------------------------------------------------------
// return the list of reflection vectors - the list of vectors is calculated only if
// SetReflectionVectorLimits is recently called
// 
// Must first call SetReflectionVectorLimits();
//------------------------------------------------------------
const vector<CRecpVector> & CUnitCell::GetReflectionVectorList( ) const
{
  return oReflectionVectorList;
}

//------------------------------------------------------------
//
//  SetUniqueReflectionVectorList
//
//------------------------------------------------------------
void CUnitCell::SetUniqueReflectionVectorList( const CSymmetry &oSym )
{
  oUniqueReflectionVectorList.clear();
  GetUniqueRecpVectors( oUniqueReflectionVectorList, oSym, oReflectionVectorList );
}
//------------------------------------------------------------
//
//  GetUniqueReflectionVectorList
//
//------------------------------------------------------------
const vector<CRecpVector> & CUnitCell::GetUniqueReflectionVectorList( ) const
{
  return oUniqueReflectionVectorList;
}

//------------------------------------------------------------
//
//  Function:  GetReflectionVectorList
//
//
//  WARNING:  This function is NOT lazied.   It both sets the new limits and calculates the reflection
//  vectors
//
//
//  Currently, it returns reflection vector up to [nHMax, nKMax, nLMax] without regards to
//  symmetry.
//
//  TODO:  Add in symmetry considerations
//
//------------------------------------------------------------
const vector<CRecpVector> & CUnitCell::GetReflectionVectorList(Int nHMax, Int nKMax, Int nLMax, Float fMaxQ, Float fMinIntensityFraction )
{
  RUNTIME_ASSERT(oReciprocalVector.size() == 3, 
                 "[CUnitCell::GetReflectionVectorList]: Recipricol Vector not initialied\n");

  SetReflectionVectorLimits( nHMax, nKMax, nLMax, fMaxQ, fMinIntensityFraction );
  return oReflectionVectorList;
}

//------------------------------------------------------------
//  WriteReflectionVectorList
//  fMinIntensity is between [0, 1] (Fractional minimum)
//------------------------------------------------------------
void CUnitCell::WriteReflectionVectorList( const string &filename, Float fMinFraction ) const
{
  ofstream oOutStream( filename.c_str() );

  for ( Size_Type i = 0; i < oReflectionVectorList.size(); i ++){

    if ( oReflectionVectorList[i].fIntensity > fMinFraction * fMaxIntensity )
    {
      SVector3 Q = oReflectionVectorList[i].v;
    
      oOutStream << oReflectionVectorList[i].h << " " 
                 << oReflectionVectorList[i].k << " " 
                 << oReflectionVectorList[i].l << " " 
        
                 << oReflectionVectorList[i].v.m_fX << " "
                 << oReflectionVectorList[i].v.m_fY << " " << oReflectionVectorList[i].v.m_fZ << " "
                 << oReflectionVectorList[i].fIntensity << " " <<   Q.GetLength() << std::endl;
    }
  }
}

//------------------------------------------------------------
//  CalculateIntensity
//  oReciprocalVector is the reciprocal lattice vectork, G_hkl
//
//
//------------------------------------------------------------
Float CUnitCell::CalculateIntensity( const SVector3 & oReciprocalVector ) const
{

  Float fSImaginary = 0;
  Float fSReal = 0;
  
  for(Size_Type i = 0; i < oTranslationVector.size(); i ++)
  {
    SVector3 oAtomPosition = oTranslationVector[i].v.m_fX * oPrimitiveVector[0]
      + oTranslationVector[i].v.m_fY * oPrimitiveVector[1]
      + oTranslationVector[i].v.m_fZ * oPrimitiveVector[2];

    Float fRDotK = Dot( oReciprocalVector, oAtomPosition );
    
    // must break into imaginary and real parts;  Amplitude = A + iB, with the
    fSReal      += oTranslationVector[i].fEffectiveZ * cos( fRDotK );
    fSImaginary += oTranslationVector[i].fEffectiveZ * sin( fRDotK );
		
  }
  // return intensity, which is f^2 * (A^2 + B^2)
  return ( fSReal * fSReal + fSImaginary * fSImaginary );
}


//----------------------------------
//  Save
//----------------------------------
Bool CUnitCell::Save   ( CSerializer   & oSerialBuf ) const
{
  Bool bSuccess = oSerialBuf.InsertCompactObj( nRecipListMaxH );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nRecipListMaxK );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nRecipListMaxL );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fRecipListfMaxQ );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fRecipMinIntensity);

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( bLimitSet );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( bLimitChanged );
  
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( nNumAtoms );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fMaxIntensity );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fAlpha );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fBeta );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fGamma );

  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fLengthA );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fLengthB );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( fLengthC );
  bSuccess = bSuccess && oSerialBuf.InsertCompactObj( oRecpLatMetricTensor );
  
  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( oTranslationVector );
  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( oReciprocalVector );
  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( oPrimitiveVector );

  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( oReflectionVectorList );
  bSuccess = bSuccess && oSerialBuf.InsertCompactVector( oUniqueReflectionVectorList );

  return bSuccess;
}

//----------------------------------
//  Restore
//----------------------------------
Bool CUnitCell::Restore( CDeserializer & oSerialBuf)
{
  
  Bool bSuccess = oSerialBuf.GetCompactObj( & nRecipListMaxH );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nRecipListMaxK );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nRecipListMaxL );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fRecipListfMaxQ );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fRecipMinIntensity);

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & bLimitSet );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & bLimitChanged );
  
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & nNumAtoms );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fMaxIntensity );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fAlpha );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fBeta );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fGamma );

  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fLengthA );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fLengthB );
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & fLengthC );
  
  bSuccess = bSuccess && oSerialBuf.GetCompactObj( & oRecpLatMetricTensor );

  bSuccess = bSuccess && oSerialBuf.GetCompactVector( oTranslationVector );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( oReciprocalVector );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( oPrimitiveVector );

  bSuccess = bSuccess && oSerialBuf.GetCompactVector( oReflectionVectorList );
  bSuccess = bSuccess && oSerialBuf.GetCompactVector( oUniqueReflectionVectorList );
  
  return bSuccess;
}

//----------------------------------
//  GetCompactSize
//----------------------------------
Size_Type CUnitCell::GetCompactSize()
{
  Size_Type nCompactSize = 0;

  nCompactSize = sizeof( nRecipListMaxH ) + sizeof(nRecipListMaxK) + sizeof( nRecipListMaxL )
    + sizeof( fRecipListfMaxQ ) + sizeof( fRecipMinIntensity ) + sizeof( bLimitSet )
    + sizeof( bLimitChanged ) + sizeof( nNumAtoms ) + sizeof( fMaxIntensity )
    + sizeof( fAlpha ) + sizeof( fBeta ) + sizeof( fGamma )
	+ sizeof( fLengthA ) + sizeof( fLengthB ) + sizeof( fLengthC );

  
  nCompactSize += sizeof( CAtom ) * oTranslationVector.size();
  nCompactSize += sizeof( SVector3 ) * oReciprocalVector.size();
  nCompactSize += sizeof( SVector3 ) * oPrimitiveVector.size();
  nCompactSize += sizeof( CRecpVector ) * oReflectionVectorList.size();
  nCompactSize += sizeof( CRecpVector ) * oUniqueReflectionVectorList.size();

  return nCompactSize;
}
