////////////////////////////////////////////////////////////////
//
//  File:    GrainAnalysis.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//  
//
//  Grain analysis library
//
//
/////////////////////////////////////////////////////////////////
#ifndef GRAINANALYSIS_H_
#define GRAINANALYSIS_H_

#include "MicMesh.h"
#include "3dMath.h"
#include "Symmetry.h"
#include "Quaternion.h"

using namespace GeneralLib;


namespace GrainAnalysis
{
  //--------------------------------------------------------------------------------------------------------
  //   PointToPointOperation
  //   -- generalization of all point to point operations between two "MIC" files - or in fact, two spatially
  //      searchable 2D structure that implements Begin(), End(), and FindOverlap
  //
  //  -- change to use back_inserter later
  //--------------------------------------------------------------------------------------------------------
  template < class Inserter, class Symmetry, class Iterator, class Operator, class PropertyMap >
  void PointToPointOperation( Inserter oResultInserter, const Symmetry & oSym,
                              const CMicMesh & oMap1, CMicMesh oMap2,
                              Iterator pSamplePointBegin, Iterator pSamplePointEnd,
                              Operator FOp, PropertyMap PropMap )
  {
    typedef CMicMesh::ShapePtrList ShapePtrList;    // modify this ?
    for( Iterator pCur = pSamplePointBegin; pCur != pSamplePointEnd; ++pCur )
    {
      ShapePtrList oList1 = oMap1.FindOverlap( PropMap.GetVoxel( *pCur ) );
      ShapePtrList oList2 = oMap2.FindOverlap( PropMap.GetVoxel( *pCur ) );
      FOp( oResultInserter, oSym, oList1, oList2 );
    }
  }
  
  //-----------------------------------
  // Trait Classes
  //-----------------------------------
  namespace XDMTraits
  {
    template< class ShapeMesh >
    struct ShapeTraits
    {
      typedef ShapeMesh           Type;
      typedef typename ShapeMesh::ShapeT   Shape; 
      typedef typename ShapeMesh::ShapePtr ShapePtr;
    };
  }
  
  namespace Operators
  {
    template < class ShapePtr >
    struct FShapePtrToQuatMap
    {
      SQuaternion GetQuaternion( ShapePtr p ) const
      {
        SQuaternion q;
        q.Set( p->oOrientMatrix );
        return q;
      }
    };

    
    template < class ShapeTrait >
    struct DereferenceMap
    {
      typename ShapeTrait::Shape GetVoxel( typename ShapeTrait::ShapePtr p ) const
      { return *p; }
    };

    template < class ShapeTrait >
    struct IdentityMap
    {
      typename ShapeTrait::Shape GetVoxel( typename ShapeTrait::Shape p ) const
      { return p; }
    };

    
    //--------------------------------------------------------------------------------------------------------
    //  LocalAveragedMisorientation
    //
    //  Given list1, list2, calculate the average orientation of the two list, then calculate the misorientation
    //  between the two.  (We are currently not doing a weighted average.)
    //
    //  Parameters:  List1, List2 -- data structure containing pointers to Shapes, as implemented by MicMesh
    //
    //
    //--------------------------------------------------------------------------------------------------------
    template< class PropertyMap >
    struct LocalAveragedMisorientation
    {
      template< class Inserter, class Symmetry, class ShapePtrList >
      void operator()( Inserter oResultInserter, const Symmetry & oSym,
                       const ShapePtrList & oList1,
                       const ShapePtrList & oList2 ) const
      {
        SQuaternion qAve1 = LatticeSymmetry::Average( oSym, oList1.begin(), oList1.end(),
                                                      PropertyMap() );
        SQuaternion qAve2 = LatticeSymmetry::Average( oSym, oList2.begin(), oList2.end(),
                                                      PropertyMap() );
        oResultInserter = GetMisorientation( oSym, qAve1, qAve2 );
      }
    };
  } // namespace operators
    
}


#include "GrainAnalysis.tmpl.cpp"

#endif /*GRAINANALYSIS_H_*/
