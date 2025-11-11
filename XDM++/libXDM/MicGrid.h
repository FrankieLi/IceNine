/////////////////////////////////////////////////////////////////
//
//  File:     MicGrid.h
//  Authors:  Frankie Li
//  e-mail:   sfli@cmu.edu
//
//  Matrix representation of MicVolume.
//  Note:  This is a special case of general MicMesh, designed for
//         high performance at the cost of flexibility.  We demand that
//         the MIC file forms a grid instead of a general mesh.
//         This is NOT a multi-scale representation.
//
//
/////////////////////////////////////////////////////////////////
#ifndef _MIC_GRID_H_
#define _MIC_GRID_H_

#include "XDMVoxel.h"
#include "MicIO.h"
#include <boost/multi_array.hpp>
#include <tuple>
#include "MicMesh.h"
#include "IteratorAdapter.h"
#include <sstream>
#include <cmath>
namespace MicAnalysis
{

  namespace Details
  {
    struct null_deleter { void operator()( void const*) const{} };
    
    template<class ShapePtr>
    struct CopyCreate
    {
      template< class T>
      ShapePtr operator()( const T *s ) const
      {
        T* p = new T(*(const_cast<T*>(s)) );
        return ShapePtr( p );
      }
      /*
        template< class T>
        ShapePtr operator()( T *s )
        {
        T* p = new T( *s );
        return ShapePtr( p );
        }
      */
    };
    
    template<class ShapePtr>
    struct ReferenceCopyNullDelete
    {
      /*
        template< class T >
        ShapePtr operator()( const T *s ) const
        {
        ShapePtr p = ShapePtr( const_cast<T*>(s), null_deleter() ); 
        return p;
        }
      */
      template< class T >
      ShapePtr operator()( T *s )
      {
        ShapePtr p = ShapePtr( s, null_deleter() ); 
        return p;
      }
      

      //       template< class T >
      //       ShapePtr operator()( const T & s ) const
      //       {
      //         ShapePtr p = ShapePtr( & const_cast<T>(s), null_deleter() ); 
      //         return p;
      //       }
    };
  }
  //--------------------------
  //  CMicGrid
  //  Requirement: Shape class with centers defined.
  //
  //  Note:  A MicGrid only accomodates comenserating grids.
  //      
  //--------------------------
  class CMicGrid
  {
  public:

    typedef XDMUtility::CMicMesh::ShapeT     Shape;   // to be templated
    typedef Shape ShapeT; //  for compatibility
    
    
    typedef std::shared_ptr<Shape>         ShapePtr;
    typedef std::vector< ShapePtr >          ShapePtrList;
    typedef ShapePtrList::iterator           ShapePtrIter;
    typedef ShapePtrList::const_iterator     ConstShapePtrIter;

    typedef boost::multi_array<ShapePtr, 2>  ShapeGridT;
    typedef ShapeGridT::index                GridIndexT;
    typedef ShapeGridT::iterator             Iterator;
    typedef ShapeGridT::const_iterator       Const_Iterator;
    typedef boost::multi_array_types::index_range GridRangeT; 
    typedef PBRMath::BBox2D                  BBox;
    
    typedef MicFile<Shape>                   Mic;

  private:
    //------------------
    //  InsertChildren
    //------------------
    void InsertChildren( const ShapePtr & pVoxel );
    
    //------------------
    // PositionToIndices
    //------------------
    void PositionToIndices( Int & nI, Int & nJ, const SVector3 & oPos ) const;

    //------------------
    //  GetIndicies
    //------------------
    void GetIndices( Int & nI, Int & nJ, const Shape & oShape) const;
    
    
    SVector3   oOrigin;
    SVector3   oBasisI;
    SVector3   oBasisJ;
    ShapeGridT oGrid;
    Float      fMinSideLength;
    Float      fInitialSideLength;
    
    Int        nISize;
    Int        nJSize;
    Int        nMaxGeneration;
   
  public:
    
    // Default ctor
    CMicGrid() : oOrigin( 0, 0, 0 ),
                 oBasisI( 0, 0, 0 ), oBasisJ( 0, 0, 0 ),
                 oGrid( boost::extents[ 0 ][ 0 ] ),
                 fMinSideLength( 0 ), nMaxGeneration( 0 )  {}
    
    // Copy ctor
    CMicGrid( const CMicGrid & oRHS )
      : oOrigin( oRHS.oOrigin ),
        oBasisI( oRHS.oBasisI ),
        oBasisJ( oRHS.oBasisJ ),
        oGrid  ( oRHS.oGrid   ),
        fMinSideLength( oRHS.fMinSideLength ),
        fInitialSideLength( oRHS.fInitialSideLength ),
        nISize( oRHS.nISize ),
        nJSize( oRHS.nJSize ),
        nMaxGeneration( oRHS.nMaxGeneration ) { }
    
    // ctor
    CMicGrid( const SVector3 oOrigin_,
              const SVector3 oBasisI_,
              const SVector3 oBasisJ_,
              Float fMinSideLength_,
              Float fInitSideLength_,
              Int nISize_,
              Int nJSize_ ):
      oOrigin( oOrigin_ ),
      oBasisI( oBasisI_ ),
      oBasisJ( oBasisJ_ ),
      fMinSideLength( fMinSideLength_ ),
      fInitialSideLength( fInitSideLength_ ),
      nISize( nISize_ ),
      nJSize( nJSize_ )
    {
      oGrid.resize  ( boost::extents[ nISize ][ nJSize ] );
      nMaxGeneration = 0;
    }
    
    //---------------------------------------------
    //  Assignment operator
    //---------------------------------------------
    CMicGrid & operator=( const CMicGrid & oRHS )
    {
      oOrigin = oRHS.oOrigin;
      oBasisI = oRHS.oBasisI;
      oBasisJ = oRHS.oBasisJ;
      fMinSideLength = oRHS.fMinSideLength;
      fInitialSideLength = oRHS.fInitialSideLength;
      nISize = oRHS.nISize;
      nJSize = oRHS.nJSize;
      oGrid.resize( boost::extents[ nISize ][ nJSize ] );
      oGrid = oRHS.oGrid;

      return *this;
    }
    
    CMicGrid( const BBox2D & oBox, Float fMinSideLength_,
              Float fInitialSideLength_ ):
      oBasisI( 1, 0, 0 ), oBasisJ( 0, 1, 0 ),
      fMinSideLength( fMinSideLength_ ),
      fInitialSideLength( fInitialSideLength_ )
      
    {
      nISize = ceil( Float( 2 ) * ( oBox.pMax.x - oBox.pMin.x ) / fMinSideLength ) + 2;
      nJSize = ceil( ( oBox.pMax.y - oBox.pMin.y ) / ( fMinSideLength * sqrt(3) / 2  ) ) + 2;
      
      oOrigin.m_fX = oBox.pMin.x + fMinSideLength / Float( 4 );
      oOrigin.m_fY = oBox.pMin.y - fMinSideLength;
      oGrid.resize( boost::extents[ nISize ][ nJSize ] );
      nMaxGeneration  = round( log( fInitialSideLength / fMinSideLength ) / log( Float(2)) );
    }
    
    //-----------------------
    //  Insert
    //-----------------------
    void Insert( const CMic & oMic );

    //-----------------------
    //  Generalized Insert
    //-----------------------
    template< class IteratorT, class FConstructor >
    void InsertReplace( Float fInitialSideLength_, IteratorT pFirst, IteratorT pEnd,
                        FConstructor PtrConstructor )
    {
      oBasisI.Set( 1, 0, 0 );
      oBasisJ.Set( 0, 1, 0 );
      Int nMaxGen = 0;
      oOrigin.Set( std::numeric_limits<Float>::max(), std::numeric_limits<Float>::max(), 0 );
      Int NumVoxels = 0;
      for(IteratorT pCur = pFirst; pCur != pEnd; ++ pCur )
      {
        NumVoxels ++;
        nMaxGen = std::max( nMaxGen, pCur->nGeneration );
        for( int i = 0; i < 3; i ++ )
        {
          SVector3 CurVertex = pCur->Vertex(i);
          oOrigin.m_fX = std::min( oOrigin.m_fX, CurVertex.m_fX );
          oOrigin.m_fY = std::min( oOrigin.m_fY, CurVertex.m_fY );
        }
      }
      std::cout << "oOrigin pure " << oOrigin << " Num Voxels " << NumVoxels <<  std::endl;
      fInitialSideLength = fInitialSideLength_;
      fMinSideLength     = fInitialSideLength_ / pow( Float( 2 ), Float( nMaxGen ) ); 
      nMaxGeneration     = nMaxGen;
      nJSize = nISize = 2 * pow( Float(2), Float( nMaxGen ) ) + 1;
      nISize *= 2;
      
      oOrigin.m_fX -= fMinSideLength / Float( 4 );
      oOrigin.m_fY -= fMinSideLength / Float( 2 );

      //----------------------------------------------------
      // DEBUG
      /*      std::cout << "MicGrid Origin " << oOrigin << std::endl;
              std::cout << " nI, nJ " << nISize << " " << nJSize << std::endl;
              std::cout << fMinSideLength << " " << fInitialSideLength << std::endl;
              std::cout << nMaxGeneration << std::endl;
      */
      oGrid.resize  ( boost::extents[ nISize ][ nJSize ] );
      for( IteratorT pCur = pFirst; pCur != pEnd; ++ pCur )
      {
        ShapePtr pToInsert = PtrConstructor( &(*pCur) ); 
        Insert( pToInsert );
      }
    }

    //-----------------------
    //  Insert - Default version
    //
    //
    //-----------------------
    template< class IteratorT >
    void InsertReplace( Float fInitialSideLength, IteratorT pFirst, IteratorT pEnd )
    {
      typename Details::CopyCreate<ShapePtr> PtrConstructor;
      InsertReplace(fInitialSideLength, pFirst, pEnd, PtrConstructor );
    }
    
    //-----------------------
    //  InsertReference - only references
    //  are inserted, and therefore saves memory.  However,
    //  the references must NOT be freed in the time that this
    //  object, MicGrid exists.
    //-----------------------
    template< class IteratorT >
    void InsertReferenceReplace( Float fInitialSideLength, IteratorT pFirst, IteratorT pEnd )
    {
      typename Details::ReferenceCopyNullDelete<ShapePtr> PtrConstructor;
      InsertReplace( fInitialSideLength, pFirst, pEnd, PtrConstructor );
    }
    //-----------------------
    //  Insert
    //-----------------------
    ShapePtr Insert( const ShapePtr & pShape );
    ShapePtr Insert( const Shape & pShape );

    //-----------------------
    //  GetNeighbors
    //-----------------------
    vector<ShapePtr> GetNeighbors( const ShapePtr & pShape ) const;
    vector<ShapePtr> GetNeighbors( const Shape    & oShape ) const;

    //-----------------------------------
    //  FindOverlap(const &SVoxel)
    //
    //  Purpose:  Find all voxels within the voxel tree that overlaps the voxel, oCenter 
    //            All (Data) voxels in oRes must overlap oCenter.
    //
    //-----------------------------------
    ShapePtrList FindOverlap( const BBox & oCenter ) const;

    //-----------------------------------
    //  FindOverlap
    //  Purpose:  Find all overlapping voxel with oCenter, given a relative
    //            error limit.  (Fractional)   Note that this limits the resolution
    //            to fRelError * fSideLength
    //-----------------------------------
    ShapePtrList FindOverlap( const ShapeT & oCenter, Float fRelError = 0.001 ) const;

    void Find( vector<ShapePtr> & oRes, const BBox & oSearchBox ) const;
    
    //-----------------------
    //  A C C E S S O R S
    //-----------------------
    ShapePtr operator() ( Int nI, Int nJ )      const;
    ShapePtr operator() ( const ShapeT   & v ) const;
    ShapePtr operator() ( const ShapePtr & v ) const;

    // ---------------- Inconsistent with operatr( v ), use with caution!!!
    ShapePtr operator() ( const SVector3 & oPos ) const;    
    
    //------------------
    //  InRange ( check to see if the indices are in range )
    //------------------
    bool InRange( Int nI, Int nJ )     const
    {
      return ( nI >= 0 && nI < nISize
               && nJ >=0 && nJ < nJSize );
    }

    //------------------
    //  IsValid ( check validity of smart pointers )
    //------------------
    static bool IsValid( const ShapePtr & p )
    {
      return ( p != ShapePtr() );
    }
    
    Mic GetMic() const;
    
    Float MinSideLength()     const { return fMinSideLength; }
    Float InitialSideLength() const { return fInitialSideLength; }
    const SVector3 & Origin() const { return oOrigin; }
    const SVector3 & BasisI() const { return oBasisI; }
    const SVector3 & BasisJ() const { return oBasisJ; }

    Iterator Begin() { return oGrid.begin(); }
    Iterator End()   { return oGrid.end();   }
    
    Const_Iterator Begin() const { return oGrid.begin(); }
    Const_Iterator End()   const { return oGrid.end();   }

    Int Size1() const { return nISize; }
    Int Size2() const { return nJSize; }
    
    //-----------------------
    //  Advanced Options
    //-----------------------
    // Set basis vectors;
    void SetBasisI( const SVector3 & v ) { oBasisI = v; oBasisI.Normalize(); }
    void SetBasisJ( const SVector3 & v ) { oBasisJ = v; oBasisJ.Normalize(); }
    void SetMinSideLength( Float f )     { fMinSideLength     = f ; }
    void SetInitialSideLength( Float f )     { fInitialSideLength = f ; }

    
    void InitializeWithCopy     ( const Mic & MicInput )
    {
      RUNTIME_ASSERT( 0, "InitializeWithCopy:  This option is NOT implemented\n");
    }

    //----------------------------------------------------------------------
    //  SampleSideLength is defined differently for different types of voxel.
    //  That's why it's left as a "SideLength"
    //----------------------------------------------------------------------
    void InitializeWithReference( Mic & MicInput ) 
    {
      SetInitialSideLength  ( MicInput.GetInitialSideLength() );
      SetMinSideLength      ( MicInput.GetMinSideLength() );
      InsertReferenceReplace( MicInput.GetInitialSideLength(), MicInput.VoxelListBegin(), MicInput.VoxelListEnd() );
    }
    
  };

  //---------------------------------------------------------------------------------------
  //  Iterator Adapter
  //---------------------------------------------------------------------------------------
  //  class CMicGridIterator : public IteratorAdapter< CMicGrid::ShapePtr, 
  
  //---------------------------------------------------------------------------------------
  //  CGridShapeLocator -- override some of the insertion functions of ShapeLocator
  //---------------------------------------------------------------------------------------
  class CGridShapeLocator : public CGeneralShapeLocator< CMicGrid >
  {
  public:
    typedef CMicGrid                          LayerContainer;
    typedef std::shared_ptr<LayerContainer> LayerContainerPtr;
  private:
    Float fMinSideLength;
    Float fInitialSideLength;
   
  public:

    CGridShapeLocator()
    {
      RUNTIME_ASSERT( 0, "Defualt Constructor of CGridShapeLocator is NOT implemented\n");
    }
    
    //-----------------------------------
    //  CGeneralShapeLocator
    //-----------------------------------
    explicit CGridShapeLocator( const BBox2D & oPlanarBounds,
                                Float fMinSideLength_, Float fInitialSideLength_,
                                Float fError ):
      CGeneralShapeLocator< CMicGrid >( oPlanarBounds, fError, 0 ),
      fMinSideLength( fMinSideLength_ ), fInitialSideLength( fInitialSideLength_ )
    {
      FErrorEstimator.fError = fError;
    }

    
    //---------------------------------------------------------------------------------------
    //  Overriding Insertions
    //---------------------------------------------------------------------------------------
    ShapePtr Insert( const ShapeT & oShape )
    {
      UpdateErrorLimit( oShape.GetBoundingBox() );
      Float fZ = oShape.Vertex(0).m_fZ;
      ShapeTree3D::iterator pIter = oObjLocator.find( fZ );
      if( pIter == oObjLocator.end() )
      {
        ShapeTreePtr pNewTree = ShapeTreePtr( new ShapeTreeT( oPlanarBBox, fMinSideLength, fInitialSideLength ) );
        bool bSuccess;
        std::tie( pIter, bSuccess) = oObjLocator.insert( std::pair<Float, ShapeTreePtr>( fZ, pNewTree ) );
      }
      return pIter->second->Insert( oShape );
    }

    //---------------------------------------------------------------------------------------
    //  Overriding Insertions
    //---------------------------------------------------------------------------------------
    ShapePtr Insert( const ShapePtr & pShape )
    {
      UpdateErrorLimit( pShape->GetBoundingBox() );
      Float fZ = pShape->Vertex(0).m_fZ;
      ShapeTree3D::iterator pIter = oObjLocator.find( fZ );
      if( pIter == oObjLocator.end() )
      {
        ShapeTreePtr pNewTree = ShapeTreePtr( new ShapeTreeT( oPlanarBBox, fMinSideLength, fInitialSideLength ) );
        bool bSuccess;
        std::tie( pIter, bSuccess) = oObjLocator.insert( std::pair<Float, ShapeTreePtr>( fZ, pNewTree ) );
      }
      return pIter->second->Insert( pShape );
    }
    
    //---------------------------------------------------------------------------------------
    //  Inserting new layer
    //---------------------------------------------------------------------------------------
    void Insert(  ShapeTreePtr & pNewLayer, Float fZ )
    {
      oObjLocator.insert( std::pair<Float, ShapeTreePtr >( fZ, pNewLayer ) );
    }
  };
  

  //------------------------------------------------------------
  //
  //  SMicGridVolume  -- simple volume construction for MicGrid
  //
  //------------------------------------------------------------
  struct SMicGridVolume
  {
    vector<CMicGrid> oGridList;
    Float fMinZ;
    Int   nSize1, nSize2;
    Float fZStep;
    typedef CMicGrid::ShapePtr ShapePtr;
    ShapePtr operator() ( const SVector3 & oPos ) const
    {
      Int nGridIndex = floor( ( oPos.m_fZ - fMinZ ) / fZStep );

      //  clamp to [0, oGridList.size() ) 
      if( nGridIndex < 0 )
        nGridIndex = 0;
      if( nGridIndex >= Int( oGridList.size() ) )
        nGridIndex = oGridList.size() - 1;
      
      return oGridList[ nGridIndex ]( oPos );
    }
    

    bool InRange( const int pIndex[3] ) const
    {
      return  ( pIndex[2] >= 0 && pIndex[2] < Int( oGridList.size() )
                && pIndex[0] >= 0 && pIndex[0] < nSize1
                && pIndex[1] >= 0 && pIndex[1] < nSize2  );
    }
    
    Float Differentiate( Int oIndices[3],
                         const SVector3 & oDiffStepSize,
                         Int nSpatialInd, Int l, Int m ) const;
    
  }; 
  
  //--------------------------------------------------------------------------------------------------------
  //  To be moved to the deeper part of the library later
  //--------------------------------------------------------------------------------------------------------
  typedef boost::multi_array< SMatrix3x3, 3 > TensorField3;
  
  //--------------------------------------------------------------------------------------------------------
  //
  //  CalculateNyeTensorField
  //
  //--------------------------------------------------------------------------------------------------------
  TensorField3 CalculateNyeTensorField( const SMicGridVolume & oMicVolume, const SVector3 & oStepSize );
  
  
  namespace Details
  {
    //--------------------------------------------------------------------------------------------------------
    //  Calculate Nye Tensor using finite difference method.
    //
    //  alpha_ij = epislon_ilk D_k G_lj
    //
    //  Note that Elastic Strain is assumed to be negligible
    //
    //  WARNING!  Beware of noisy data.
    //
    //--------------------------------------------------------------------------------------------------------
    SMatrix3x3 CalculateNyeTensor( const SMicGridVolume & oMicVolume,
                                   Int oIndices[3], const SVector3 & oDiffStepSize );
    

    
  }
  
  
  

  //---------------------------------------------------------------------------------------
  //  Rectilinear MicGrid
  //
  //   Probably will need to refactor both this and MicGrid to inherit from the same
  //   interface in the near future.  Not everything will be kept.
  //---------------------------------------------------------------------------------------
  template< class ShapeT >
  class RectMicGrid
  {
  public:
    typedef std::shared_ptr<ShapeT>        ShapePtr;
    typedef std::vector< ShapePtr >          ShapePtrList;
    typedef typename ShapePtrList::iterator  ShapePtrIter;
    typedef typename boost::multi_array<ShapePtr, 2>  ShapeGridT;
    typedef typename ShapeGridT::index                GridIndexT;
    
    typedef typename XDMUtility::SerializingIterator< ShapeGridT >             Iterator;
    typedef typename XDMUtility::Const_SerializingIterator< const ShapeGridT > Const_Iterator;
    
    typedef typename boost::multi_array_types::index_range GridRangeT; 
    typedef PBRMath::BBox2D                  BBox;

    
    typedef MicFile<ShapeT>                  Mic;

  private:
    //------------------
    // PositionToIndices
    //------------------
    void PositionToIndices( Int & nI, Int & nJ, const SVector3 & oPos ) const
    {
      SVector3 D = oPos - oOrigin;
      //      std::cout << fVoxelSideLength << " | " << D << " " << oOrigin <<  std::endl;
      nI = floor( std::max( D.m_fX / fVoxelSideLength, static_cast<Float>( 0 ) ) );   // bounded below
      nJ = floor( std::max( D.m_fY / fVoxelSideLength, static_cast<Float>( 0 ) ) );
    }
    
    //------------------
    //  GetIndicies
    //------------------
    void GetIndices( Int & nI, Int & nJ, const ShapeT& oShape) const
    {
      PositionToIndices( nI, nJ, oShape.GetCenter() );
    }

    
    SVector3   oOrigin;
    ShapeGridT oGrid;
    Float      fVoxelSideLength;
    Float      fSampleSpaceSideLength;


  public:

    RectMicGrid( ) :
      oOrigin( 0, 0, 0 ), oGrid(),
      fVoxelSideLength( 0 ), fSampleSpaceSideLength( 0 ) {}
   
    RectMicGrid( const SVector3 & Origin_, const Float SampleSideLength,
                 const Float VoxelSideLength )
      : fVoxelSideLength( VoxelSideLength ),
        fSampleSpaceSideLength( SampleSideLength )
    {
      int VoxelPerSide = std::ceil( SampleSideLength / VoxelSideLength ); 
      oGrid.resize  ( boost::extents[ VoxelPerSide ][ VoxelPerSide ] );
    }
    
    // use default copy ctor
    // use default operator=
    // use default cop de-ctor


    //-----------------------
    //  Generalized Insert -- should really be called insert replace....
    //-----------------------
    template< class IteratorT, class FConstructor >
    void InsertReplace( Float fInitialSideLength_, 
                        IteratorT pFirst, IteratorT pEnd,
                        FConstructor PtrConstructor )
    {
      oOrigin.Set( std::numeric_limits<Float>::max(), std::numeric_limits<Float>::max(), 0 );
      Float fMinSideLength = std::numeric_limits<Float>::max();
      for(IteratorT pCur = pFirst; pCur != pEnd; ++ pCur )
      {
        oOrigin.m_fX = std::min( oOrigin.m_fX, pCur->GetCenter().m_fX );
        oOrigin.m_fY = std::min( oOrigin.m_fY, pCur->GetCenter().m_fY );
        fMinSideLength = std::min( pCur->fSideLength, fMinSideLength ); // there should only be one side length
      }
      
      fSampleSpaceSideLength   = fInitialSideLength_;
      oOrigin.m_fX -= fVoxelSideLength / Float( 2 );
      oOrigin.m_fY -= fVoxelSideLength / Float( 2 );
      
      fVoxelSideLength = fMinSideLength;
      int VoxelPerSide = std::ceil( fSampleSpaceSideLength / fVoxelSideLength ); 
      
      oGrid.resize  ( boost::extents[ VoxelPerSide ][ VoxelPerSide ] );
      for( IteratorT pCur = pFirst; pCur != pEnd; ++ pCur )
      {
        ShapePtr pToInsert = PtrConstructor( &(*pCur) ); 
        Insert( pToInsert );
      }
    }

    
    void Insert( const CMic & oMic )  // legacy support
    {
      Insert( oMic.GetInitialSideLength(), oMic.VoxelListBegin(), oMic.VoxelListEnd() );
    }
    //-----------------------
    //  Insert - Default version
    //-----------------------
    template< class IteratorT >
    void InsertReplace( Float SampleSpaceSideLength, IteratorT pFirst, IteratorT pEnd )
    {
      typename Details::CopyCreate<ShapePtr> PtrConstructor;
      InsertReplace( SampleSpaceSideLength, pFirst, pEnd, PtrConstructor );
    }
    
    //-----------------------
    //  InsertReference - only references
    //  are inserted, and therefore saves memory.  However,
    //  the references must NOT be freed in the time that this
    //  object, MicGrid exists.
    //-----------------------
    template< class IteratorT >
    void InsertReferenceReplace( Float fInitialSideLength, IteratorT pFirst, IteratorT pEnd )
    {
      typename Details::ReferenceCopyNullDelete<ShapePtr> PtrConstructor;
      InsertReplace(fInitialSideLength, pFirst, pEnd, PtrConstructor );
    }

    //-----------------------
    //  Insert
    //-----------------------
    ShapeT   Insert( const ShapePtr & pShape )
    {
      Int nI, nJ;
      GetIndices( nI, nJ, *pShape );
      oGrid[nI][nJ] = pShape;
      return *pShape;
    }
    ShapePtr Insert( const ShapeT & oShape )
    {
      ShapePtr pShape = ShapePtr( new ShapeT( oShape ) );
      Insert( pShape );
      return pShape;
    }
    
    //-----------------------
    //  GetNeighbors
    //-----------------------
    vector<ShapePtr> GetNeighbors( const ShapePtr & pShape ) const
    {
      return GetNeighbors( *pShape );
    }
    
    vector<ShapePtr> GetNeighbors( const ShapeT   & oShape ) const
    {
      vector<ShapePtr> Results;
      Int nI, nJ;
      GetIndices( nI, nJ, oShape );

      int i, j;
      i = nI - 1; j = nJ;
      if( InRange( i, j ) && IsValid( oGrid[i][j] ) )
        Results.push_back( oGrid[i][j] );
      
      i = nI + 1; j = nJ;
      if( InRange( i, j ) && IsValid( oGrid[i][j] ) )
        Results.push_back( oGrid[i][j] );

      i = nI; j = nJ - 1;
      if( InRange( i, j ) && IsValid( oGrid[i][j] ) )
        Results.push_back( oGrid[i][j] );

      i = nI; j = nJ + 1;
      if( InRange( i, j ) && IsValid( oGrid[i][j] ) )
        Results.push_back( oGrid[i][j] );
      
      
      return Results;
    }

    //-----------------------------------
    //  FindOverlap(const &SVoxel)
    //
    //  Purpose:  Find all voxels within the voxel tree that overlaps the voxel, oCenter 
    //            All (Data) voxels in oRes must overlap oCenter.
    //
    //-----------------------------------
    ShapePtrList FindOverlap( const BBox & oCenter ) const
    {
      ShapePtrList oRes;
      Find( oRes, oCenter );
      return oRes;
    }
    
    //-----------------------------------
    //  FindOverlap
    //  Purpose:  Find all overlapping voxel with oCenter, given a relative
    //            error limit.  (Fractional)   Note that this limits the resolution
    //            to fRelError * fSideLength
    //-----------------------------------
    ShapePtrList FindOverlap( const ShapeT & SearchCenter, Float fRelError = 0.001 ) const
    {
      BBox SearchBox;
      SVector3 CenterPos = SearchCenter.GetCenter();
      SearchBox.pMin.x = CenterPos.m_fX - fVoxelSideLength / static_cast< Float > ( 2 );
      SearchBox.pMin.y = CenterPos.m_fY - fVoxelSideLength / static_cast< Float > ( 2 );
      SearchBox.pMax.x = CenterPos.m_fX + fVoxelSideLength / static_cast< Float > ( 2 );
      SearchBox.pMax.y = CenterPos.m_fY + fVoxelSideLength / static_cast< Float > ( 2 );

      return FindOverlap( SearchBox );
    }
    
    void Find( vector<ShapePtr> & oRes, const BBox & oSearchBox ) const
    {
      SVector3 oMin( oSearchBox.pMin.x, oSearchBox.pMin.y, 0 );
      SVector3 oMax( oSearchBox.pMax.x, oSearchBox.pMax.y, 0 );
      
      Int nMinI, nMaxI, nMinJ, nMaxJ;
      
      PositionToIndices( nMinI, nMinJ, oMin );
      PositionToIndices( nMaxI, nMaxJ, oMax );
      
      for( Int i = nMinI; i <= nMaxI; i ++ )
        for( Int j = nMinJ; j <= nMaxJ; j ++ )
          if( InRange( i, j ) && IsValid( oGrid[i][j] ) )
            oRes.push_back( oGrid[i][j] );
      
    }
    
    //-----------------------
    //  A C C E S S O R S
    //-----------------------
    ShapePtr operator() ( Int nI, Int nJ )        const  { return oGrid[ nI ][ nJ ]; }

    ShapePtr operator() ( const ShapeT & v ) const
    {
      return operator()( v.GetCenter());
    }

    ShapePtr operator() ( const ShapePtr & v ) const
    {
      return operator()( v->GetCenter());
    }

    ShapePtr operator() ( const SVector3 & oPos ) const
    {
      Int nI, nJ;
      PositionToIndices( nI, nJ, oPos );
      return operator()( nI, nJ );
    }

    
    //------------------
    //  InRange ( check to see if the indices are in range )
    //------------------
    bool InRange( Int nI, Int nJ )     const
    {
      return ( nI >= 0 && nI < Size1()
               && nJ >=0 && nJ < Size2() );
    }

    //------------------
    //  IsValid ( check validity of smart pointers )
    //------------------
    static bool IsValid( const ShapePtr & p )
    {
      return ( p != ShapePtr() );
    }

    //------------------
    Mic GetMic() const
    {
      Mic oRes;
      oRes.Origin           = oOrigin;
      oRes.VoxelSideLength  = fVoxelSideLength;
      oRes.SampleSideLength = fSampleSpaceSideLength;
      
      for( Int i = 0; i < Size1(); i ++ )
        for( Int j = 0; j < Size2(); j ++ )
          if( IsValid( oGrid[i][j] ) )
            oRes.AddVoxel( *oGrid[i][j] );
      return oRes;
    }

    void Initialize( const SVector3 & Origin_, const Float SampleSideLength,
                     const Float VoxelSideLength )
    {
      oOrigin = Origin_;
      fVoxelSideLength       = VoxelSideLength;
      fSampleSpaceSideLength = SampleSideLength;
      int VoxelPerSide = std::ceil( SampleSideLength / VoxelSideLength ); 
      oGrid.resize  ( boost::extents[ VoxelPerSide ][ VoxelPerSide ] );
    }

    void InitializeAndFill( const SVector3 & Origin_, const Float SampleSideLength,
                            const Float VoxelSideLength )
    {
      fVoxelSideLength       = VoxelSideLength;
      fSampleSpaceSideLength = SampleSideLength;
      int VoxelPerSide = std::ceil( SampleSideLength / VoxelSideLength ); 
      oGrid.resize  ( boost::extents[ VoxelPerSide ][ VoxelPerSide ] );

      SVector3 OriginShift( 0.5 * fVoxelSideLength, 0.5 * fVoxelSideLength, 0 );
      for( Int i = 0; i < Size1(); i ++ )
      {
        for( Int j = 0; j < Size2(); j ++ )
        {
          ShapeT Voxel;
          SVector3 Center = Origin_ + SVector3( i * fVoxelSideLength, j * fVoxelSideLength, 0) + OriginShift;
          Voxel.SetCenter( Center, fVoxelSideLength );
          Insert( Voxel );
        }
      }
    }


    
    Float MinSideLength()     const { return fVoxelSideLength; }
    Float InitialSideLength() const { return fSampleSpaceSideLength; }
    
    const SVector3 & Origin() const { return oOrigin; }

    Const_Iterator begin() const { return Const_Iterator( oGrid, Size1(), Size2(), 0, 0 ); }
    Const_Iterator end  () const { return Const_Iterator( oGrid, Size1(), Size2(), Size1(), Size2() ); }

    Iterator begin() { return Iterator( oGrid, Size1(), Size2(), 0, 0 ); }
    Iterator end  () { return Iterator( oGrid, Size1(), Size2(), Size1(), Size2() ); }
    
    Int Size1() const { return oGrid.shape()[0]; }
    Int Size2() const { return oGrid.shape()[1]; }
    
    //-----------------------
    //  Advanced Options
    //-----------------------
    void SetMinSideLength( Float f )     {  std::cerr << " Warning:  SetMinSideLength called from RectMicGrid, which doesn't do anything" << std::endl; }
    void SetInitialSideLength( Float f ) {  std::cerr << " Warning:  SetInitialSideLength called from RectMicGrid, which doesn't do anything" << std::endl; }



    //----------------------------------------------------------------------
    //  SampleSideLength is defined differently for different types of voxel.
    //  That's why it's left as a "SideLength"
    //----------------------------------------------------------------------
    void InitializeWithReference( Mic & MicInput ) 
    {
      Initialize( MicInput.Origin, MicInput.SampleSideLength, MicInput.VoxelSideLength );
      typedef typename Mic::VoxelType_iterator Iter;
      for( Iter pCur = MicInput.VoxelListBegin(); pCur != MicInput.VoxelListEnd(); ++ pCur )
      {
        Details::ReferenceCopyNullDelete<ShapePtr> PtrConstructor;
        ShapePtr pToInsert = PtrConstructor( &(*pCur) ); 
        Insert( pToInsert );
      }
    }
    
  };
  
}

#endif

