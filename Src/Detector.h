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
//
//   Detector.h
//   Author:   Frankie Li  (sfli@cmu.edu)
//
//   Purpose:
//    Implementation of detector geometry and physical properties.
//    
//    Detector is a container that holds images.  It has physical properties
//    such as orientation, location, pixel resolutions, and so forth, and therefore
//    needs to be initialized before use.  Since each detector holds a number of images
//    (Think of them as cameras holding film), each of the picture from the film is taken
//    at the same detector setup.
//
//    Note that added capability is included in this class.  For example, the ability to
//    search quickly each of the images held in this detector are to have labels.
//
//
//   TODO:
//         Quadtree peak searching (?  don't know if this will be useful, though bounding boxes are nice)  
//                                 ( Will require insertion of bounding boxes. )
//
//
//
////////////////////////////////////////////////////////////


#ifndef DETECTOR_H_
#define DETECTOR_H_

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "Types.h"
#include "3dMath.h"
#include "InitFilesIO.h"
#include "Pixel.h"
#include "ConfigFile.h"
#include <vector>
#include <string>
#include "BBox.h"
#include "Debug.h"
#include <fstream>
#include "ImageData.h"
#include "DetectorFile.h"


using namespace GeneralLib;
using std::string;
using std::ofstream;

// Forward declaration
namespace InitFileIO
{
   class CDetectorInfo;
}


//--------------------------
//  Class CDetector
//
//  Detector image of the experiment.  Data is read in 
//  as collection of peaks, reduced from tif image
//  files.
//--------------------------
class CDetector
{
	
private:
 
  //--------------------------
  // 
  // Detector information
  //
  //
  //--------------------------
  SVector3 oDetectorRotationCenter;  // Rotation center of the detector plane.  i.e., the fixed point
                                     // under rotation of the detector given in a vector in the lab frame
                                     // In another words, this is where the detector intersects the beam.
   //--------------------------
  //
  //  Matrices used to transform between 
  //  different coordinates
  //--------------------------
  SMatrix3x3 oLabToDetectorMatrix;
  SMatrix3x3 oPrevOrientationMatrix;   // For state saving purpose
  
  //--------------------------
  // DETECTOR Coordinate Setup
  //--------------------------
  
  // define the plane of the detector with three points
  SVector3 oDetectorPt1, oDetectorPt2, oDetectorPt3;
  CPlane oLabFrameImagePlane;
	
  // Pixel oDetectorCenter;
  Float fBeamCenterJ;           // TODO:  Figure out if this affects precision limit
  Float fBeamCenterK;
  
  //--------------------------
  //  Defines the coordinate origin of the CCD
  //
  //  Note that all three of the vectors below are suseptible to detector tilts
  //
  //  Everythign is defined in the lab frame.  Note that oImageCoordinateOrigin does NOT
  //  include displacement of the detector.  This is because the displacement vector
  //  cannot be dependent on the detector tilt, unlike the origin of the image frame.
  //
  //
  //  All of these are defined in the lab frame
  //--------------------------
  SVector3 oImageBasisJ;  // the two directions form a basis for the image 
  SVector3 oImageBasisK;  // their direction indicates how rows and columns are counted
  SVector3 oImageCoordinateOrigin; // coordinate origin of the image frame (i.e., where (0, 0)
                                   // is located in the image.

  //--------------------------
  //
  //  A detector frame version of the vectors above
  //
  //--------------------------
  SVector3 oDetFrameImageBasisJ;
  SVector3 oDetFrameImageBasisK;
  SVector3 oDetFrameImageCoordinateOrigin;
  
  
  Int nNumRows;
  Int nNumCols;

  Float fPixelHeight;
  Float fPixelWidth;
  Float fPixelHalfHeight;
  Float fPixelHalfWidth;

  void CalculateImagePlane();

  void SetDetectorPlane( const SVector3 &v1,
                         const SVector3 &v2,
                         const SVector3 &v3 )
  {
    oDetectorPt1 = v1;
    oDetectorPt2 = v2;
    oDetectorPt3 = v3;
    CalculateImagePlane();
  }
  
 
  CDetector(){};   // This is not allowed
public:
  
 
  
  //--------------------------
  //
  //  CDetector
  //
  //  A word on detector Local coordinates
  //  In the local coordinate of the detector, the vertical component
  //  is Y, and the horizontal component is X
  //  
  //  fBeamCenterX fBeamCenterY are the x,y location in mm of the 
  //  of the beam center on the detector
  //
  //  fPixelXLength fPixelYLength is the physical size of the pixel
  //--------------------------
  CDetector(const BBox2D &range, 
            Float fBeamCenterX, Float fBeamCenterY,
            Float fPixelXLength, Float fPixelYLength,
            UInt nRows, UInt nCols,
            const SVector3 & oImageBasisJ,
            const SVector3 & oImageBasisK,
            const SVector3 & oImageCoordinateOrigin );

	
  void SetImageParameter(const BBox2D &range, 
                         Float fBeamCenterX, Float fBeamCenterY,
                         Float fPixelXLength, Float fPixelYLength,
                         UInt nCols, UInt nRows, 
                         const SVector3 & oImageBasisJ,
                         const SVector3 & oImageBasisK,
                         const SVector3 & oImageCoordinateOrigin );
	
  void SetLocation( const SVector3 &loc );
  void SetOrientation( Float phi, Float theta, Float psi );
  void SetOrientation( const SMatrix3x3 &oOrientation );

  //
  //  D E T E C T O R   T R A N S F O R M A T I O N S
  //
  void Rotate( Float phi, Float theta, Float psi );
  void Translate( const SVector3 &loc );
  
  //
  //  C O O R D I N A T E   T R A N S F O R M A T I O N S
  // 
  
  //--------------------------
  //
  //  Given a SVector3 oPos, the position of a point on the detector
  //  plane in the lab frame.  
  //  Output the row and col pixel value of the point.
  //
  //--------------------------
  void LabToPixel( Float & nRow, Float & nCol, const SVector3 &oPos ) const;
  void LabToDetectorCoordinate(Float &fJ, Float &fK,  const SVector3 &oPos ) const;
  SVector3 DetectorToLabCoordinate( Float fJ, Float fK ) const;
  
  //--------------------------
  //  PixelToLabCoordinate
  //  Col = X pixel, Row = Y pixel
  //--------------------------
  SVector3 PixelToLabCoordinate( const Float & fXPixel, const Float & fYPixel ) const;
  Float    ColPixelToImageJ( Float p) const;
  Float    RowPixelToImageK( Float p) const;
  Int      ToRowPixel(Float f)  const;
  Int      ToColPixel(Float f) const;
  
  
  // 
  //  A C C E S S O R S
  //

  SVector3 GetLocation( ) const;
  SVector3 GetRotationCenter( ) const;
  Float    GetBeamCenterJ () const;
  Float    GetBeamCenterK () const;
  Float    GetPixelWidth  () const;
  Float    GetPixelHeight () const;
  
  
  SMatrix3x3 GetOrientationMatrix( ) const;
  Float GetDetectorWidth( ) const;
  Float GetDetectorHeight( ) const;
  Int GetNumCols( ) const;
  Int GetNumRows( ) const;
  

  //-----------------
  //  GetDetectorCoordinateOrigin() -- return coordinate origin in lab frame
  //-----------------
  SVector3 GetDetectorCoordinateOrigin( ) const;
  void GetCoordinateBasis( SVector3 & oJBasis, SVector3 & oKBasis ) const;
  CPlane GetLabDetectorPlane();


  //-----------------
  //  Detector Frame version of some of the above accessors 
  //-----------------
  SVector3 GetDetFrameCoordinateOrigin() const;
  void GetDetFrameCoordBasis( SVector3 & oJBasis, SVector3 & oKBasis ) const;
  
  
  //
  //  U T I L I T I E S
  //

  //--------------------------
  //  Return the maximum extent direction
  //  in the number of pixels
  //  (This is a bounding box approximation)
  //
  //--------------------------
  Int GetPixelExtent( const SVoxel & oVoxel ) const;
  
  //--------------------------
  //  InRange -- return true of Col, Row is in the range of the detector
  //
  //  Note:  Col = x pixel, Row = y pixel 
  //--------------------------
  template< typename IndexType >
  inline bool InRange( IndexType Col, IndexType Row ) const
  {
    return ( Row < nNumRows && Row > 0 && Col < nNumCols && Col > 0);
  }
  
  //--------------------------
  //  TO CHANGE
  //
  //  AddDirectBeam
  // 
  //  fBeamHeight and fBeamWidth are in mm
  //--------------------------
  void AddDirectBeam( CImageData &oOutputImage, Float fBeamHeight, Float fBeamWidth );


  //--------------------------
  //  Intersect:  Returns true if interescts,
  //              false otherwise.  fT contains
  //              the value of the parameter t where oRay
  //              intersects the detector plane.
  //--------------------------
  bool Intersects( const CRay &oRay, Float &fT ) const;

  
};

//--------------------------
//
//  CXDMDetectorFactory
////  A simple factory to generate detectors.  Note that
//  this is specific to HEDM experiment.  (A new factor for
//  a new experiment.)
//
//  Purpose:  To keep the design general, often parameters used
//            in the simulation (internal representation) are different
//            from the actual experiment.  The reason range from performance
//            to convience.  Therefore, initialization from configuration
//            files may not be straightforward.  Therefore, it'd be wise
//            to have a centralized location where all initialization are
//            done.  Similarly, modification of some of the parameters may
//            also lead to significant changes in the overall structure.
//            One example would be a change of detector pixel size.  This
//            has a cascading effect in various geometric parameters of
//            the detector.  This is certainly not true in general.  Certain
//            things can be easily modified by using CDetector's mutators.
//            However, experiment specific changes is not known to the detector.
//            For example, convention of coordinate axis, direction of x-ray
//            beam, and so forth.  The rationale of not having this factory
//            inside CDetector is that we want different experiment to be able
//            to setup their corrosponding factories.
//--------------------------
class CXDMDetectorFactory
{
  
private:
  CXDMDetectorFactory();

  //--------------------------
  // GetCoordOrigin
  //--------------------------
  static SVector3 GetCoordOrigin( Float fBeamCenterJ, Float FBeamCenterK,
                                  Float fPixelWidth,  Float fPixelHeight,
                                  const SVector3 & oJUnitVector,
                                  const SVector3 & oKUnitVector );
public:


  //--------------------------
  //  SDetParameters
  //  Purpose:  Define a set of parameters
  //            that can be modified using this
  //            factory.
  //--------------------------
  struct SDetParameters
  {
    SMatrix3x3 oOrientation;
    SVector3   oPosition;
    Float      fBeamCenterJ;
    Float      fBeamCenterK;
    Float      fPixelWidth;    // J length
    Float      fPixelHeight;   // K length
    Int        nNumJPixels;
    Int        nNumKPixels;
  };
  
  //--------------------------
  //  MakeDetector
  //
  //  Purpose:  Construct a detector used in HEDM experiment from scratch
  //            
  //  Note that the orientations are Euler angles in radians.
  //--------------------------
  static CDetector MakeDetector( Int nNumJPixels, Int nNumKPixels,
                                 const SVector3 & oDetLabFrameLoc,
                                 Float fBeamCenterJ, Float fBeamCenterK,
                                 const SVector3 &oJUnitVector,
                                 const SVector3 &oKUnitVector,
                                 Float fPixelWidth, Float fPixelHeight,
                                 const SMatrix3x3 & oDetLabFrameOrientation );
  
  //--------------------------
  //  ModifyImageParameters
  //--------------------------
  static void ModifyImageParameters( CDetector & oDetector,
                                     const SVector3 & oDetLabFrameLoc,
                                     Float fBeamCenterJ, Float fBeamCenterK,
                                     Float fPixelWidth,  Float fPixelHeight,
                                     const SMatrix3x3 & oDetLabFrameOrientation );

  //--------------------------
  //  GetImageParameters
  //--------------------------
  static SDetParameters GetImageParameters( const CDetector & oDetector );

  //--------------------------
  //  GetImageInfo
  //  Purpose:  Return the specification of the images in a CDetectorInfo
  //            format.
  //--------------------------
  static InitFileIO::CDetectorInfo GetImageInfo( const CDetector & oDetector );
  
};

#endif
