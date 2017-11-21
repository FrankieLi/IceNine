//-----------------------------------------------------------
//
//  VTKPropertyMap
//  Author:  Frankie Li
//
//  e-mail:  sfli@cmu.edu
//
//-----------------------------------------------------------

#ifndef _VTK_PROPERTY_MAP_H
#define _VTK_PROPERTY_MAP_H

#include <map>
#include <CGAL/IO/File_medit.h>
#include <boost/shared_ptr.hpp>
#include "MicVolume.h"
#include <typeinfo>
#include <fstream>
namespace Pandora
{
  namespace Details
  {
    //-------------------------------------------------
    //  ConstructCellVoxelPtrMap
    //-------------------------------------------------
    template< class Map, class C3T3 >
    void ConstructCellVoxelPtrMap( Map & CellToRawDataMap,
                                   const MicAnalysis::CMicVolume & oVolume,
                                   const C3T3 & c3t3,
                                   const GeneralLib::SVector3 & Origin,
                                   const GeneralLib::SVector3 & ScaleVector );

    //-------------------------------------------------
    //  ConstructCellVoxelPtrMap
    //-------------------------------------------------
    template< class Map, class CellToVoxelPtrMap,
              class C3T3, class Symmetry, class Converter >
    void ConstructCellToGrainAverageMap( Map & CellToGrainAverageMap,
                                         const CellToVoxelPtrMap & VoxelPtrMap,
                                         const C3T3 & c3t3,
                                         const Symmetry & oSym,
                                         Converter Convert );
    
    //-------------------------------------------------
    //  ConstructVertexVoxelPtrMap
    //-------------------------------------------------
    template< class Map, class C3T3 >
    void ConstructVertexVoxelPtrMap( Map & VertexToRawDataMap,
                                     const MicAnalysis::CMicVolume & oVolume,
                                     const C3T3 & c3t3,
                                     const GeneralLib::SVector3 & Origin,
                                     const GeneralLib::SVector3 & ScaleVector );

    //-------------------------------------------------
    //  ConstructCellToMisorientationMaps
    //-------------------------------------------------
    template< class ThetaMap, class RodMap,
              class VHToVoxelMap, class C3T3, class Symmetry >
    void ConstructCellToMisorientationMaps( ThetaMap & CellToThetaMap,
                                            RodMap & CellToRodMap,
                                            const VHToVoxelMap & VHMap,
                                            const C3T3 & c3t3,
                                            const GeneralLib::SVector3 & ScaleVector,
                                            const Symmetry & oSym);
    
    template< class Map, class C3T3, class ShapePtrToFieldFn, class BBox >
    void ConstructVertexFieldMap( Map & VertexToRawDataMap,
                                  const MicAnalysis::CMicVolume & oVolume,
                                  const C3T3 & c3t3,
                                  const GeneralLib::SVector3 & Origin,
                                  const GeneralLib::SVector3 & ScaleVector,
                                  ShapePtrToFieldFn ShapePtrToField,
                                  const BBox & SmoothingVolume );
    
    //-------------------------------------------------
    //  FindClosest
    //-------------------------------------------------
    template< class ShapePtrVector >
    Int FindClosest( const GeneralLib::SVector3 & oPos,
                     const ShapePtrVector & NgbVoxels )
    {
      
      float fMinR = 10000;
      int   nMinR = 0;
      for( int i = 0; i < NgbVoxels.size(); i ++ )
      {
        float r = (NgbVoxels[i]->GetCenter() - oPos).GetLength() ;
        if(  r < fMinR )
        {
          fMinR = r;
          nMinR = i;
        }
      }
      return nMinR;
    }

    //-------------------------------------------------
    // QuatToRod
    //-------------------------------------------------
    GeneralLib::SVector3 QuatToRod( const GeneralLib::SQuaternion & q )
    {
      SVector3 rod( q.m_fX, q.m_fY, q.m_fZ );
      rod /= q.m_fW;
      return rod;
    }
    
    //-------------------------------------------------
    //  Clearly this is possibly the worse way to do things
    //-------------------------------------------------
    GeneralLib::SVector3 ToCubicColorFZ( const GeneralLib::SQuaternion & q )
    {
      typedef LatticeSymmetry::CCubicSymmetry Symmetry; 
      typedef GeneralLib::SQuaternion SQuaternion;
      typedef GeneralLib::SVector3    SVector3;

      float fColorEdge = (sqrt( 2.0 ) - 1) + 0.001;
      float CubicRange = static_cast< float >( 2 ) * fColorEdge;
      SVector3 Shift( fColorEdge, fColorEdge, fColorEdge );
      
      const std::vector<SQuaternion> & oSym = Symmetry::Get().GetQuatOperatorList();
      for( int i = 0; i < oSym.size(); i ++ )
      {
        SQuaternion qFZ =  oSym[i] * q  ;
        SVector3 rod = QuatToRod( qFZ );
        float x = std::fabs( rod.m_fX );
        float y = std::fabs( rod.m_fY );
        float z = std::fabs( rod.m_fZ );
        if( x <= fColorEdge && y <= fColorEdge && z <= fColorEdge
            && ( x+ y + z ) <= 1.0 )
        {
          rod += Shift;
          rod /= CubicRange;
          return SVector3( rod.m_fZ, rod.m_fY, rod.m_fX );
        }
      }
      return QuatToRod( q );
    }

    
    struct QuatToCubicColor
    {
      GeneralLib::SVector3 operator()( const GeneralLib::SQuaternion & q )
      { return ToCubicColorFZ( q );  }
    };
  }
  
  //-------------------------------------------------
  //  VtkPropMap
  //
  //  This is an erasure for VtkPropertyMap.  By having
  //  an erasure, we can remove type information without
  //  using void*, which is dangerous at best.  This particular
  //  erasure uses the idea of concept, and it originally
  //  came from:
  //  http://www.cplusplus.com/forum/articles/18756/
  //
  //-------------------------------------------------
  class VtkPropMap
  {
  private:

    struct VtkPropMapConcept
    {
      virtual ~VtkPropMapConcept() {}
      virtual std::string GetVtkRank()     const = 0;
      virtual std::string GetName()        const = 0;
      virtual std::string GetVtkDataType() const = 0;
      virtual int         GetNumTuples()   const = 0;
      virtual void OutputToVtk( std::ostream & os ) const = 0;
      virtual bool        DefaultLookupTable()      const = 0;
    };
    
    //------------------------------------------
    //  The definition of the model is solely for the
    //  purpose of being a copy-constructable.  I.e.,
    //  type information is not needed.
    //------------------------------------------
    template< typename T >
    struct VtkPropMapModel : VtkPropMapConcept
    {
    
      VtkPropMapModel( const T& t )
        : ConcreteMap( t )
      {}
      
      // Model Accessor - no one should have to read this
      virtual std::string GetVtkRank()              const { return ConcreteMap.GetVtkRank(); }
      virtual std::string GetName()                 const { return ConcreteMap.GetName(); }
      virtual std::string GetVtkDataType()          const { return ConcreteMap.GetVtkDataType(); }
      virtual int         GetNumTuples()            const { return ConcreteMap.GetNumTuples(); }
      virtual void OutputToVtk( std::ostream & os ) const { return ConcreteMap.OutputToVtk( os ); }
      virtual bool        DefaultLookupTable()      const { return ConcreteMap.DefaultLookupTable(); }
      virtual ~VtkPropMapModel() {}
    private:
      T ConcreteMap;
    };

    boost::shared_ptr<VtkPropMapConcept> MapModelHandle;
    
  public:
    
    template< typename T >
    VtkPropMap( const T& ConcreteMap ) :
      MapModelHandle( new VtkPropMapModel<T>( ConcreteMap ) )
    { }
    
    // Model Accessor - no one should have to read this
    
    std::string GetVtkRank()              const { return MapModelHandle->GetVtkRank(); }
    std::string GetName()                 const { return MapModelHandle->GetName(); }
    std::string GetVtkDataType()          const { return MapModelHandle->GetVtkDataType(); }
    int         GetNumTuples()            const { return MapModelHandle->GetNumTuples(); }
    void OutputToVtk( std::ostream & os ) const { return MapModelHandle->OutputToVtk( os ); }
    bool        DefaultLookupTable()      const { return MapModelHandle->DefaultLookupTable(); }
  };   // end VtkPropMap


  //-------------------------------------------------
  //  HandleFieldMap
  //
  //  --  A model of VtkPropMapConcept
  //  --  HandleIterator is a special type of iterator
  //      that is convertable between iterator and handle.
  //      It's a CGAL thing.
  //-------------------------------------------------
  template < class Handle, class Property,
             class HandleIterator,
             class FieldSelector,
             class PropertyMapT = std::map< Handle, Property >   >
  class HandleFieldMap
  {
  private:
    std::string Rank;
    std::string Name;
    std::string DataType;
    int         NumTuples;
    
    //   typedef typename std::map< Handle, Property > PropertyMap;
    typedef PropertyMapT PropertyMap;
  protected:

    const PropertyMap & PMap;
    FieldSelector  Selector;
    HandleIterator pBegin;
    HandleIterator pEnd;
    bool           bDefaultLookup;
  public:

    virtual ~HandleFieldMap() {}
    HandleFieldMap( const std::string & sRank,
                    const std::string & sName,
                    const std::string & sDataType,
                    int N,
                    const PropertyMap & PMap_,
                    FieldSelector  SFn_,
                    HandleIterator pBegin_, HandleIterator pEnd_,
                    bool bDefaultLookup_ = true )
      : Rank( sRank ), Name( sName ),
        DataType( sDataType ), NumTuples( N ),
        PMap( PMap_ ), Selector( SFn_ ),
        pBegin( pBegin_ ), pEnd( pEnd_ ), bDefaultLookup( bDefaultLookup_ )
    {   }
    
    std::string GetVtkRank()     const { return Rank; }
    std::string GetName()        const { return Name; }
    std::string GetVtkDataType() const { return DataType; }
    int         GetNumTuples()   const { return NumTuples; }
    bool  DefaultLookupTable()   const { return  bDefaultLookup;  }
    
    //-------------------------------------------------
    //  GetProperty
    //  if h is not in the property map, the behavior
    //  will be undefined
    //-------------------------------------------------
    Property GetProperty( const Handle & h ) const
    {
      typename PropertyMap::const_iterator pCur = PMap.find( h );
      if( pCur != PMap.end() )
        return pCur->second;
      else
        return Property();
    }
    
    //-------------------------------------------------
    // OutputProperty  -- Overriden.
    //-------------------------------------------------
    void OutputToVtk( std::ostream & os ) const
    {
      for( HandleIterator pCur = pBegin; pCur != pEnd; ++ pCur )
      {
        Handle h = pCur;
        os << Selector( GetProperty( h ) ) << std::endl;
      }
    }

  };
  
  
  //-------------------------------------------------
  //  By definition, ID must be >= 0
  //-------------------------------------------------
  template < class Handle >
  class HandleIDMap
  {
  private:
    std::map< Handle, int > IDMap;
    int CurrentID;
    
  public:
    
    HandleIDMap(): IDMap(), CurrentID(0) {}
    HandleIDMap( int nID ): IDMap(), CurrentID( nID ) {}
    
    template< class InputIterator >
    void insert( InputIterator pFirst, InputIterator pEnd )
    {
      for( InputIterator pCur = pFirst; pCur != pEnd; ++ pCur )
      {
        IDMap.insert( std::make_pair( *pCur, CurrentID ) );
        ++ CurrentID;
      }
    }
    
    int operator()( const Handle & vh ) const
    {
      typename std::map< Handle, int >::iterator pResult = IDMap.find( vh );
      if( pResult != IDMap.end() )
        return pResult->second;
      else
        return -1;
    }
  };


  //--------------------------------------------------------
  //  VoxelFieldSelector
  //
  //  -- Selectors (functors) for different voxel fields
  //--------------------------------------------------------
  namespace VoxelFieldSelector
  {

    template< class T, class U >
    struct FirstSelector
    {
      T operator()( const std::pair<T, U> & o ) const
      {
        return o.first;
      }
    };
    
    template< class T, class U >
    struct SecondSelector
    {
      U operator()( const std::pair<T, U> & o ) const
      {
        return o.second;
      }
      
    };
    
    struct IdentitySelector
    {
      template< typename T >
      T operator()( const T & o ) const
      {
        return o;
      }
    };
    
    struct ConfidenceSelector
    {
      template< typename ShapePtr >
      float operator()( const ShapePtr & p ) const
      {
        if( p != ShapePtr() )
          return p->fConfidence;
        else
          return -1;
      }
    };
    
    struct EulerAngleSelector
    {
      template< typename ShapePtr >
      GeneralLib::SVector3 operator()( const ShapePtr & p ) const
      {
        if( p != ShapePtr() )
          return p->oOrientMatrix.GetEulerAngles();
        else
          return GeneralLib::SVector3( -1, -1, -1 );
      }
    };

    struct RodriguzSelector
    {
      typedef GeneralLib::SVector3    SVector3;
      typedef GeneralLib::SQuaternion SQuaternion; 
      template< typename ShapePtr >
      SVector3 operator()( const ShapePtr & p ) const
      {
        if ( p == ShapePtr() )
          return SVector3( 0, 0, 0 );
        
        SQuaternion q;
        q.Set( p->oOrientMatrix );
        SVector3 rod( q.m_fX, q.m_fY, q.m_fZ );
        rod /= q.m_fW;
        return rod;
      }
    };
    
    struct RodColorSelector
    {
      RodriguzSelector RFMap;
      typedef GeneralLib::SVector3    SVector3;
      typedef GeneralLib::SQuaternion SQuaternion; 
      template< typename ShapePtr >
      SVector3 operator()( const ShapePtr & p ) const
      {
        SVector3 rod = RFMap( p );
        float fColorEdge = (sqrt( 2.0 ) - 1) + 0.001;

        rod += SVector3(  fColorEdge, fColorEdge, fColorEdge );
        rod /= ( 2.0 * fColorEdge );

        return rod;
      }
    };


    struct CubicFZRodColorSelector
    {
      typedef GeneralLib::SVector3    SVector3;
      typedef GeneralLib::SQuaternion SQuaternion;
      
      template< typename ShapePtr >      
      SVector3 operator()( const ShapePtr & p ) const
      {
        if( p == ShapePtr() )
          return SVector3( 0, 0, 0 );
        SQuaternion q;
        q.Set( p->oOrientMatrix );
        return Details::ToCubicColorFZ( q );
      }
    };

    
    //-------------------------------------------
    //  CellMisorientationSelector
    //------------------------------------------- 
    template< class VertexHandleToVoxelMap,
              class C3T3,
              class Symmetry >
    class CellMisorientationSelector
    {
    private:
      CellMisorientationSelector();
      const VertexHandleToVoxelMap & VHMap;
      const C3T3                   & c3t3;
      typedef typename VertexHandleToVoxelMap::mapped_type    ShapePtr;
      typedef typename VertexHandleToVoxelMap::const_iterator VHMapIter;
    public:
      CellMisorientationSelector( const VertexHandleToVoxelMap & VMap,
                                  const C3T3 & c3t3_ )
        : VHMap( VMap ), c3t3( c3t3_ ) {}
      
      typedef GeneralLib::SVector3                        SVector3;
      typedef GeneralLib::SQuaternion                     SQuaternion;
      typedef typename C3T3::Cell_handle                  Cell_handle;
      typedef typename C3T3::Vertex_handle                Vertex_handle;
      typedef typename C3T3::Triangulation::Vertex::Point Point3;
      
      //------------------------------------
      //  operator()
      //------------------------------------
      std::pair< SQuaternion, float > operator()( Cell_handle ch,
                                                  const SVector3 & ScaleVector ) const
      {
        //------------------------------------
        //  Calculating cell center
        SVector3 CellCenter( 0, 0, 0 );
        for( int i = 0; i < 4; i ++ )
        {
          Point3 p = ch->vertex(i)->point();
          CellCenter += SVector3(p.x(), p.y(), p.z() );
          
        }
        CellCenter /= static_cast<float>( 4 );   

        //------------------------------------
        //  Get all vertices 
        std::vector< Vertex_handle > VertexList;
        for( int i = 0; i < 4; i ++ )
          if( c3t3.in_dimension( ch->vertex(i) ) == 3 )
            VertexList.push_back( ch->vertex(i) );

        std::vector<SQuaternion>   CountedOrientList;
        std::vector<Vertex_handle> CountedVertexList;
        CountedOrientList.reserve( VertexList.size() );
        CountedVertexList.reserve( VertexList.size() );
        for( int i = 0; i < VertexList.size(); ++ i )
        {
          VHMapIter pFound = VHMap.find( VertexList[i] );
          if( pFound != VHMap.end() )
          {
            if( pFound->second != ShapePtr() )
            {
              SQuaternion q;
              q.Set( pFound->second->oOrientMatrix );
              CountedOrientList.push_back( q );
              CountedVertexList.push_back( VertexList[i] );
            }
          }
        }
        
        SQuaternion qAve = LatticeSymmetry::Average( Symmetry::Get(),
                                                     CountedOrientList );
        float Max_dTheta_dR = 0;
        SQuaternion Max_dQ;
        bool  bDegenerate = false;
        for( int i = 0; i < CountedVertexList.size(); ++ i )
        {
          float fMisorient =
            LatticeSymmetry::GetMisorientation( Symmetry::Get(),
                                                qAve, CountedOrientList[i] );
          SVector3 SamplePoint( CountedVertexList[i]->point().x(),
                                CountedVertexList[i]->point().y(),
                                CountedVertexList[i]->point().z() );
          SVector3 oDiff = SamplePoint - CellCenter;
          oDiff.m_fX *= ScaleVector.m_fX;
          oDiff.m_fY *= ScaleVector.m_fY;
          oDiff.m_fZ *= ScaleVector.m_fZ;
          
          float dR = oDiff.GetLength();
          if( dR > 0 )
          {
            float dTheta_dR = fMisorient / dR;
            if( dTheta_dR > Max_dTheta_dR )
            {
              Max_dTheta_dR = dTheta_dR;
              Max_dQ = qAve.Inverse() * CountedOrientList[i];
            }
          }
          else
          {
            bDegenerate = true;
          }
        }
        if( bDegenerate )
        {
          Max_dQ.Set( -2, -2, -2, -2 );
          return std::make_pair( Max_dQ, static_cast<float>( -5 ) );
        }
        return std::make_pair( Max_dQ,
                               Max_dTheta_dR );
      }
    };


    
    //---------------------------------------
    //  ShapePtrListToLocalDeformationEst
    //---------------------------------------
    template< class Symmetry >
    class ShapePtrListToOrientationGradient
    {
    public:
      
      typedef GeneralLib::SQuaternion                  SQuaternion;
      typedef typename std::pair< SQuaternion, float > ReturnT;

      //---------------------------------------
      //  Precondition - ShapePtr are all valid
      //---------------------------------------
      template< class ShapePtr, class VertexHandle >
      ReturnT operator()( const std::vector<ShapePtr> & ShapeList,
                          VertexHandle vh )
      {
        SQuaternion Max_dQ;
        Max_dQ.Set(-2, -2, -2, -2);
        float       Max_dTheta_dR = -5;

        if( ShapeList.size() == 0 )
          return ReturnT( Max_dQ, Max_dTheta_dR );
        
        vector<SQuaternion> OrientList;
        SVector3 Center( 0, 0, 0 );
        for( int i = 0; i < ShapeList.size(); i ++ )
        {
          SQuaternion q;
          q.Set( ShapeList[i]->oOrientMatrix );
          OrientList.push_back( q );
          Center += ShapeList[i]->GetCenter();
        }
        Center /= static_cast<float>( ShapeList.size() );
        SQuaternion qAve = LatticeSymmetry::Average( Symmetry::Get(),
                                                     OrientList );
        
        for( int i = 0; i < ShapeList.size(); i ++ )
        {
          SVector3 oDiff = ShapeList[i]->GetCenter() - Center;
          float dR = oDiff.GetLength();
          if( dR > 0 )
          {
            float fMisorient =
              LatticeSymmetry::GetMisorientation( Symmetry::Get(),
                                                  qAve, OrientList[i] );   
            float dTheta_dR = fMisorient / dR;
            if( dTheta_dR > Max_dTheta_dR )
            {
              //      std::cout << "Accept " << std::endl;
              Max_dTheta_dR = dTheta_dR;
              Max_dQ = ReduceToFundamentalZone( Symmetry::Get(),
                                                qAve.Inverse() * OrientList[i] );
            }
          }
        }
        return ReturnT( Max_dQ, Max_dTheta_dR );
      }
      
    };
    
    //---------------------------------------
    //  ShapePtrListToLocalDeformationEst
    //---------------------------------------
    template< class Symmetry >
    class ShapePtrListToLocalAverageMisorientation
    {
    public:
      
      typedef GeneralLib::SQuaternion                  SQuaternion;
      typedef GeneralLib::SVector3                     SVector3;
      typedef typename std::pair< SQuaternion, float > ReturnT;

      //---------------------------------------
      //  Precondition - ShapePtr are all valid
      //---------------------------------------
      template< class ShapePtr, class VertexHandle >
      ReturnT operator()( const std::vector<ShapePtr> & ShapeList,
                          VertexHandle vh )
      {
        SQuaternion Max_dQ;
        Max_dQ.Set(-2, -2, -2, -2);
        float       Ave_dTheta = -5;

        if( ShapeList.size() < 2 )
          return ReturnT( Max_dQ, Ave_dTheta );

        SVector3 VertexCenter( vh->point().x(),
                               vh->point().y(),
                               vh->point().z() );
        
        int CenterIndex = Details::FindClosest( VertexCenter, ShapeList );
        
        int nCounted = 0;
        float fTotalMisorient = 0;
        const float Threshold = DEGREE_TO_RADIAN( 15 );
        for( int i = 0; i < ShapeList.size(); i ++ )
        {
          if( i != CenterIndex )
          {
            float fMisorient =
              LatticeSymmetry::GetMisorientation( Symmetry::Get(),
                                                  ShapeList[i]->oOrientMatrix,
                                                  ShapeList[ CenterIndex ]->oOrientMatrix );   
            if( fMisorient < Threshold )
            {
              fTotalMisorient += fMisorient;
              nCounted ++;
            }
          }
        }
        if( nCounted > 0 )
          Ave_dTheta = fTotalMisorient / static_cast<float>( nCounted );
        else
          Ave_dTheta = 0;
        
        return ReturnT( Max_dQ, Ave_dTheta );
      }
      
    };
  }
}

//--------------------------------------------------------
//
//
//--------------------------------------------------------
namespace Pandora
{
  namespace Details
  {
    //------------------------------------------
    //  ConstructCellOrientMap
    //------------------------------------------
    template< class Map, class C3T3 >
    void ConstructCellVoxelPtrMap( Map & CellToRawDataMap,
                                   const MicAnalysis::CMicVolume & oVolume,
                                   const C3T3 & c3t3,
                                   const GeneralLib::SVector3 & Origin,
                                   const GeneralLib::SVector3 & ScaleVector )
    {
      typedef typename C3T3::Cell_iterator Cell_iter;
      typedef typename C3T3::Cell_handle   Cell_handle;
      typedef typename C3T3::Triangulation::Vertex::Point Point3;
      typedef typename GeneralLib::SVector3 SVector3;
      typedef typename MicAnalysis::CMicVolume::ShapePtr ShapePtr;
      
      SVector3 oDiag( oVolume.MinSideLength(),
                      oVolume.MinSideLength(),
                      oVolume.MinSideLength() );

      for( Cell_iter pCur = c3t3.cells_begin(); pCur != c3t3.cells_end(); ++ pCur )
      {
        SVector3 pCenter( 0, 0, 0 );
        for( int i = 0; i < 4; i ++ )
        {
          Point3 p = pCur->vertex(i)->point();
          pCenter += SVector3(p.x(), p.y(), p.z() );
          
        }
        pCenter /= static_cast<float>( 4 );

        GeneralLib::SVector3 CorrectedPosition( pCenter.m_fX * ScaleVector.m_fX,
                                                pCenter.m_fY * ScaleVector.m_fY,
                                                pCenter.m_fZ * ScaleVector.m_fZ);
        CorrectedPosition += Origin;
        MicAnalysis::CMicVolume::BBox3D oSearchBox;
        oSearchBox.m_oBoxMin = CorrectedPosition - oDiag;
        oSearchBox.m_oBoxMax = CorrectedPosition + oDiag;
        std::vector< ShapePtr > oRes = oVolume.Find( oSearchBox );
        Cell_handle ch = pCur;
        if( oRes.size() > 0 )
        {
          int nClosest = FindClosest( CorrectedPosition, oRes );
          CellToRawDataMap.insert( std::make_pair( ch, oRes[ nClosest ] ) ); 
        }
        else
        {
          CellToRawDataMap.insert( std::make_pair( ch, ShapePtr() ) ); 
        }
      }
    }


    //-------------------------------------------------
    //  ConstructCellVoxelPtrMap
    //-------------------------------------------------
    template< class Map, class CellToVoxelPtrMap,
              class C3T3, class Symmetry, class ConverterT >
    void ConstructCellToGrainAverageMap( Map & CellToGrainAverageMap,
                                         const CellToVoxelPtrMap & VoxelPtrMap,
                                         const C3T3 & c3t3,
                                         const Symmetry & oSym, ConverterT Convert )
      
    {
      typedef typename C3T3::Subdomain_index  Subdomain_index;
      typedef typename C3T3::Cell_iterator    Cell_iterator;
      typedef typename C3T3::Cell_handle                  Cell_handle;
      typedef typename CellToVoxelPtrMap::mapped_type    ShapePtr;
              
      typedef GeneralLib::SQuaternion         SQuaternion;
      typedef typename std::vector<SQuaternion>   QuatList;
      typedef typename std::map< Subdomain_index, QuatList > GrainID2OrientationMap;
      
      GrainID2OrientationMap IDToOrientListMap;
      CGAL::Default_cell_index_pmap<C3T3>     CellToIDMap( c3t3  );

      
      for( Cell_iterator ch = c3t3.cells_begin(); ch != c3t3.cells_end(); ++ ch )
      {
        typename CellToVoxelPtrMap::const_iterator pFound = VoxelPtrMap.find( ch );
        if( pFound != VoxelPtrMap.end() )
        {
          SQuaternion q;
          if( pFound->second != ShapePtr() )
          {
            Subdomain_index id = CellToIDMap.subdomain_index( ch );
            q.Set( pFound->second->oOrientMatrix );
            if( IDToOrientListMap.find( id ) != IDToOrientListMap.end() )
            {
              IDToOrientListMap[id].push_back(q);
            }
            else
            {
              QuatList qList;
              qList.push_back( q );
              IDToOrientListMap[id] = qList;
            }
          }
        }
      }
      typename std::map< Subdomain_index, SQuaternion> ID2AveOrientMap;
      typedef typename GrainID2OrientationMap::iterator Iter;
      for( Iter pCur = IDToOrientListMap.begin(); pCur != IDToOrientListMap.end();
           ++ pCur )
      {
        SQuaternion q  = LatticeSymmetry::Average( oSym, pCur->second );
        ID2AveOrientMap.insert( std::make_pair(  pCur->first, q ) );
      }
      
      for( Cell_iterator pCur = c3t3.cells_begin(); pCur != c3t3.cells_end(); ++ pCur )
      {
        Cell_handle ch = pCur;
        SQuaternion qAve = ID2AveOrientMap[ CellToIDMap.subdomain_index( ch ) ];
        CellToGrainAverageMap.insert( std::make_pair( ch, Convert( qAve ) ) );
      }
    }

    
    //-------------------------------------------------
    //  ConstructVertexVoxelPtrMap
    //-------------------------------------------------
    template< class Map, class C3T3 >
    void ConstructVertexVoxelPtrMap( Map & VertexToRawDataMap,
                                     const MicAnalysis::CMicVolume & oVolume,
                                     const C3T3 & c3t3,
                                     const GeneralLib::SVector3 & Origin,
                                     const GeneralLib::SVector3 & ScaleVector )
    {
      typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
      typedef typename C3T3::Vertex_handle   Vertex_handle;
      typedef typename C3T3::Triangulation::Vertex Vertex;
      typedef typename C3T3::Triangulation::Vertex::Point Point3;
      typedef typename GeneralLib::SVector3 SVector3;
      typedef typename MicAnalysis::CMicVolume::ShapePtr ShapePtr;
      
      SVector3 oDiag( oVolume.MinSideLength(),
                      oVolume.MinSideLength(),
                      oVolume.MinSideLength() );
      
      for( Finite_vertices_iterator pCur = c3t3.triangulation().finite_vertices_begin();
           pCur != c3t3.triangulation().finite_vertices_end(); ++ pCur )
      {
        Point3 p = pCur->point();
        
        GeneralLib::SVector3 CorrectedPosition( p.x() * ScaleVector.m_fX,
                                                p.y() * ScaleVector.m_fY,
                                                p.z() * ScaleVector.m_fZ);
        CorrectedPosition += Origin;
        
        MicAnalysis::CMicVolume::BBox3D oSearchBox;
        oSearchBox.m_oBoxMin = CorrectedPosition - oDiag;
        oSearchBox.m_oBoxMax = CorrectedPosition + oDiag;
        std::vector< ShapePtr > oRes = oVolume.Find( oSearchBox );
        
        Vertex_handle vh = pCur;
        if( oRes.size() > 0 )
        {
          int nClosest = FindClosest( CorrectedPosition, oRes );
          VertexToRawDataMap.insert( std::make_pair( vh, oRes[ nClosest ] ) ); 
        }
        else
        {
          VertexToRawDataMap.insert( std::make_pair( vh, ShapePtr() ) ); 
        }
      }
    }

    //-------------------------------------------------
    //  ConstructCellToMisorientationMaps
    //-------------------------------------------------
    template< class ThetaMap, class RodMap,
              class VHToVoxelMap, class C3T3, class Symmetry >
    void ConstructCellToMisorientationMaps( ThetaMap & CellToThetaMap,
                                            RodMap & CellToRodMap,
                                            const VHToVoxelMap & VHMap,
                                            const C3T3 & c3t3,
                                            const GeneralLib::SVector3 & ScaleVector,
                                            const Symmetry & oSym )
    {
      typedef typename C3T3::Cell_iterator Cell_iterator;
      typedef typename VoxelFieldSelector
        ::CellMisorientationSelector< VHToVoxelMap, C3T3, Symmetry > C2MisT;
      typedef GeneralLib::SQuaternion SQuaternion;
      
      C2MisT MisorientationFn( VHMap, c3t3 );
      for( Cell_iterator pCur = c3t3.cells_begin(); pCur != c3t3.cells_end(); ++ pCur )
      {
        std::pair< SQuaternion, float >  oRes = MisorientationFn( pCur, ScaleVector );
        CellToThetaMap[ pCur ] = RADIAN_TO_DEGREE( oRes.second );
        SVector3 rod( oRes.first.m_fX, oRes.first.m_fY, oRes.first.m_fZ );
        rod /= oRes.first.m_fW;
        CellToRodMap[ pCur ] = rod;
      }
    }

    //-------------------------------------------------
    //  ConstructVertexFieldMap
    //
    //  Requirement:  Map::mapped_type must be the same
    //                as ShapePtrToFieldFn's return type
    //
    //-------------------------------------------------
    template< class Map, class C3T3, class ShapePtrToFieldFn, class BBox >
    void ConstructVertexFieldMap( Map & VertexToRawDataMap,
                                  const MicAnalysis::CMicVolume & oVolume,
                                  const C3T3 & c3t3,
                                  const GeneralLib::SVector3 & Origin,
                                  const GeneralLib::SVector3 & ScaleVector,
                                  ShapePtrToFieldFn ShapePtrToField,
                                  const BBox & SmoothingVolume )
    {
      typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
      typedef typename C3T3::Vertex_handle   Vertex_handle;
      typedef typename C3T3::Triangulation::Vertex Vertex;
      typedef typename C3T3::Triangulation::Vertex::Point Point3;
      typedef typename GeneralLib::SVector3 SVector3;
      typedef typename MicAnalysis::CMicVolume::ShapePtr ShapePtr;
      
      
      for( Finite_vertices_iterator pCur = c3t3.triangulation().finite_vertices_begin();
           pCur != c3t3.triangulation().finite_vertices_end(); ++ pCur )
      {
        Point3 p = pCur->point();
        GeneralLib::SVector3 CorrectedPosition( p.x() * ScaleVector.m_fX,
                                                p.y() * ScaleVector.m_fY,
                                                p.z() * ScaleVector.m_fZ);
        CorrectedPosition += Origin;
        MicAnalysis::CMicVolume::BBox3D SearchBox = SmoothingVolume;
        SearchBox.m_oBoxMin += CorrectedPosition;
        SearchBox.m_oBoxMax += CorrectedPosition;
        
        std::vector< ShapePtr > oRes = oVolume.Find( SearchBox );
        Vertex_handle vh = pCur;
        VertexToRawDataMap.insert( std::make_pair( vh, ShapePtrToField( oRes, vh ) ) ); 
      }
    }

    
  }  // end namespace Details

}


#endif
