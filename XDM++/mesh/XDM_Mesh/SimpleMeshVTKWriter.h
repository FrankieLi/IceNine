//-----------------------------------------------------------
//
//  SimpleMeshVTKWriter
//  Author:  Frankie Li
//
//  e-mail:  sfli@cmu.edu
//
//  Note:  Some parts of this file is modeled after File_medit.h
//  from CGAL/IO (as included)
//
//-----------------------------------------------------------

#ifndef _SIMPLE_MESH_VTK_WRITER_H
#define _SIMPLE_MESH_VTK_WRITER_H

#include <iostream>
#include <fstream>
#include <map>
#include <CGAL/IO/File_medit.h>
namespace Pandora
{

  namespace Details
  {
    template< class Complex_3_In_Triangulation_3, class VertexToIDMap >
    void WriteVerticesToPolydataVtk( std::ostream & os,
                                     const Complex_3_In_Triangulation_3 & c3t3,
                                     VertexToIDMap & VertexIDPMap );
    
    template< class Complex_3_In_Triangulation_3, class VertexToIDMap >
    void WriteFacetsToPolydataVtk( std::ostream & os,
                                   const Complex_3_In_Triangulation_3 & c3t3,
                                   const VertexToIDMap & VertexIDPMap );

    template< class Complex_3_In_Triangulation_3, class VertexToIDMap >
    void WriteCellsToUnstructuredVtk( std::ostream & os,
                                      const Complex_3_In_Triangulation_3 & c3t3,
                                      const VertexToIDMap & VertexIDPMap );
  }
  
  //-----------------------------------------------------------
  //  WriteMeshToPolydataVtk
  //
  //  Purpose:  Given a complex_3_in_triangulation_3 object and
  //            an iterator of property maps, write to a Polydata
  //            format Vtk file.
  //-----------------------------------------------------------
  template < class Complex_3_In_Triangulation_3, class PropertyMapIterator >
  void WriteMeshToPolydataVtk( std::ostream & os,
                               const Complex_3_In_Triangulation_3 & c3t3,
                               PropertyMapIterator PMapFirst,
                               PropertyMapIterator PMapLast );
  
  //-----------------------------------------------------------
  //  WriteMeshToUnstructuredVtk
  //
  //  Purpose:  Given a complex_3_in_triangulation_3 object and
  //            an iterator of property maps, write to an unstructured
  //            format Vtk file.
  //-----------------------------------------------------------
  template < class Complex_3_In_Triangulation_3, class PropertyMapIterator >
  void WriteMeshToUnstructuredVtk( std::ostream & os,
                                   const Complex_3_In_Triangulation_3 & c3t3,
                                   PropertyMapIterator PointDataPMapFirst,
                                   PropertyMapIterator CointDataPMapLast,
                                   PropertyMapIterator CellDataPMapFirst,
                                   PropertyMapIterator CellDataPMapLast,
                                   PropertyMapIterator FieldDataPMapFirst,
                                   PropertyMapIterator FieldDataPMapLast );

  
}//  end namespace Pandora



//--------------------------------------------------------------------
//
//  I M P L E M E N T A T I O N
//
//--------------------------------------------------------------------
namespace Pandora
{

  namespace Details
  {
    //-----------------------------------------------------------
    // WriteVerticesToPolydataVtk
    //
    //-----------------------------------------------------------
    template< class Complex_3_In_Triangulation_3, class VertexToIDMap >
    void WriteVerticesToPolydataVtk( std::ostream & os,
                                     const Complex_3_In_Triangulation_3 & c3t3,
                                     VertexToIDMap & VertexIDPMap )
    {
      typedef typename Complex_3_In_Triangulation_3::Triangulation Tr;
      typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Point Point_3;

      const Tr & tr = c3t3.triangulation();
      int nCurID = 0;
      for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
           vit != tr.finite_vertices_end();
           ++vit)
      {
        Vertex_handle v = vit;
        VertexIDPMap[v] = nCurID;
        ++ nCurID;
        Point_3 p = vit->point();
        os << CGAL::to_double(p.x()) << " "
           << CGAL::to_double(p.y()) << " "
           << CGAL::to_double(p.z()) << std::endl;
      }
    }
    
    //-----------------------------------------------------------
    // WriteFacetsToPolydataVtk
    //
    //-----------------------------------------------------------
    template< class Complex_3_In_Triangulation_3, class VertexToIDMap >
    void WriteFacetsToPolydataVtk( std::ostream & os,
                                   const Complex_3_In_Triangulation_3 & c3t3,
                                   const VertexToIDMap & VertexIDPMap )
    {
      typedef typename Complex_3_In_Triangulation_3::Triangulation Tr;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Complex_3_In_Triangulation_3::Facet_iterator Facet_iterator;
      
      for( Facet_iterator fit = c3t3.facets_begin();
           fit != c3t3.facets_end();
           ++fit)
      {
        os << 3 << " ";
        for (int i=0; i<4; i++)
        {
          if (i != fit->second)   // facet is represeted by a cell and an index opposite the facet
                                  // Therefore, fit->second (opp index) must not be included
          {
            const Vertex_handle& vh = (*fit).first->vertex(i);
            typename VertexToIDMap::const_iterator pFound = VertexIDPMap.find( vh );
            if( pFound != VertexIDPMap.end() )
              os << pFound->second << " ";
            else
              os << -1 << " ";  // error condition
          }
        }
        os << std::endl;
      }
    }

    
    //-----------------------------------------------------------
    // WriteCellsToUnstructuredData
    //
    //-----------------------------------------------------------
    template< class Complex_3_In_Triangulation_3, class VertexToIDMap >
    void WriteCellsToUnstructuredVtk( std::ostream & os,
                                      const Complex_3_In_Triangulation_3 & c3t3,
                                      const VertexToIDMap & VertexIDPMap )
    {
      typedef typename Complex_3_In_Triangulation_3::Triangulation Tr;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Complex_3_In_Triangulation_3::Cell_iterator Cell_iterator;
      
      for( Cell_iterator cit = c3t3.cells_begin() ;
           cit != c3t3.cells_end() ;
           ++cit )
      {
        os << "4 ";
        for (int i=0; i<4; i++)
        {
          const Vertex_handle& vh = cit->vertex(i);
          typename VertexToIDMap::const_iterator pFound = VertexIDPMap.find( vh );
          if( pFound != VertexIDPMap.end() )
            os <<  pFound->second << " ";
          else
            os << -1 << " ";  // error condition
          
        }      
        os << std::endl;
      }
    }
    
  }   //-----------------------------------------------------------

  //-----------------------------------------------------------
  //  WriteMeshToPolydataVtk
  //
  //  Purpose:  Given a complex_3_in_triangulation_3 object and
  //            an iterator of property maps, write to a Polydata
  //            format Vtk file.
  //-----------------------------------------------------------
  template < class Complex_3_In_Triangulation_3, class PropertyMapIterator >
  void WriteMeshToPolydataVtk( std::ostream & os,
                               const Complex_3_In_Triangulation_3 & c3t3,
                               PropertyMapIterator PMapFirst,
                               PropertyMapIterator PMapLast )
  {
    // --- header information
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << "Polydata file from Pandora" << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET POLYDATA" << std::endl;
    // --- end header information

    typedef typename Complex_3_In_Triangulation_3::Facet_iterator Facet_iterator;
    typedef typename Complex_3_In_Triangulation_3::Facet Facet;
    typedef typename Complex_3_In_Triangulation_3::Surface_index  Surface_index;
    typedef typename Complex_3_In_Triangulation_3::Triangulation::Vertex_handle Vertex_handle;
    CGAL::Default_facet_index_pmap<Complex_3_In_Triangulation_3>   FacetToIDMap( c3t3  );
    
    typedef typename Complex_3_In_Triangulation_3::Triangulation Tr;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename std::map<Vertex_handle, int> VertexPMap;
    VertexPMap VertexToIDMap;

    // section header information
    os << "POINTS " << c3t3.triangulation().number_of_vertices() << " float" << std::endl;
    Details::WriteVerticesToPolydataVtk( os, c3t3, VertexToIDMap );
    os << "POLYGONS " << c3t3.number_of_facets() << " " << c3t3.number_of_facets()  * 4 << std::endl;
    Details::WriteFacetsToPolydataVtk( os, c3t3, VertexToIDMap );
    os << "CELL_DATA " << c3t3.number_of_facets() << std::endl;
    os << "SCALARS SurfaceID int 1 " << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    for( Facet_iterator fit = c3t3.facets_begin();
         fit != c3t3.facets_end();
         ++fit)
    {
      os <<  FacetToIDMap.surface_index( *fit ) << std::endl;
    }
      
    for( PropertyMapIterator pCur = PMapFirst; pCur != PMapLast; ++ pCur )
    {
      os << pCur->GetVtkRank() << " " << pCur->GetName() << " " << pCur->GetVtkDataType() << std::endl;
      if( pCur->DefaultLookupTable() )
        os << "LOOKUP_TABLE default" << std::endl;
      pCur->OutputToVtk( os );
    }
  }

   
  //-----------------------------------------------------------
  //  WriteMeshToUnstructuredVtk
  //
  //  Purpose:  Given a complex_3_in_triangulation_3 object and
  //            an iterator of property maps, write to an unstructured
  //            format Vtk file.
  //
  //  Note:     Currently this only outputs volumetric data.
  //
  //
  //-----------------------------------------------------------
  template < class Complex_3_In_Triangulation_3,
             class PropertyMapIterator >
  void WriteMeshToUnstructuredVtk( std::ostream & os,
                                   const Complex_3_In_Triangulation_3 & c3t3,
                                   PropertyMapIterator PointDataPMapFirst,
                                   PropertyMapIterator PointDataPMapLast,
                                   PropertyMapIterator CellDataPMapFirst,
                                   PropertyMapIterator CellDataPMapLast,
                                   PropertyMapIterator FieldDataPMapFirst,
                                   PropertyMapIterator FieldDataPMapLast )
  {
    typedef typename Complex_3_In_Triangulation_3::Cell_iterator Cell_iterator;
    typedef typename Complex_3_In_Triangulation_3::Cell_handle   Cell_handle;
    typedef typename Complex_3_In_Triangulation_3::Subdomain_index  Subdomain_index;
    typedef typename Complex_3_In_Triangulation_3::Triangulation::Vertex_handle Vertex_handle;
    // CGAL::Default_cell_index_pmap<Complex_3_In_Triangulation_3>   CellToIDMap( c3t3  );
    
    typedef typename Complex_3_In_Triangulation_3::Triangulation Tr;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename std::map<Vertex_handle, int> VertexPMap;
    VertexPMap VertexToIDMap;

    
    // --- header information
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << "Unstructured Data converged from Pandora" << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET UNSTRUCTURED_GRID" << std::endl;
    // --- end header information
    
    // section header information
    os << "POINTS " << c3t3.triangulation().number_of_vertices() << " float" << std::endl;
    Details::WriteVerticesToPolydataVtk( os, c3t3, VertexToIDMap );  // Polydata and Unstructure share the same point format

    os << "CELLS " << c3t3.number_of_cells()
       << " " << c3t3.number_of_cells() * 5 << std::endl;
    Details::WriteCellsToUnstructuredVtk( os, c3t3, VertexToIDMap );
       
    // Write cells
    os << "CELL_TYPES " << c3t3.number_of_cells() << std::endl;
    for( int i = 0; i < c3t3.number_of_cells(); ++ i )
      os << 10 << std::endl;
    
    if( PointDataPMapFirst != PointDataPMapLast )
    {
      //  Grain ID  -- should make into property map
      //--------------------------------------------------
      os << "POINT_DATA " << c3t3.triangulation().number_of_vertices() << std::endl;
      //--------------------------------------------------
      for( PropertyMapIterator pCur = PointDataPMapFirst;
           pCur != PointDataPMapLast; ++ pCur )
      {
        os << pCur->GetVtkRank() << " " << pCur->GetName()
           << " " << pCur->GetVtkDataType() << std::endl;
        if( pCur->DefaultLookupTable() )
          os << "LOOKUP_TABLE default" << std::endl;
        pCur->OutputToVtk( os );
      }
    }
    //  Grain ID  -- should make into property map
    //--------------------------------------------------
    os << "CELL_DATA " << c3t3.number_of_cells() << std::endl;
    os << "SCALARS GrainID int 1" << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;

    for( Cell_iterator ch = c3t3.cells_begin(); ch != c3t3.cells_end(); ++ ch )
    {
      os << c3t3.subdomain_index( ch ) << std::endl;
    }
    std::cout << "Before writing cell " << std::endl;
    //      os << CellToIDMap.subdomain_index( ch ) << std::endl;
    //--------------------------------------------------
    if( CellDataPMapFirst != CellDataPMapLast )
    { 
      for( PropertyMapIterator pCur = CellDataPMapFirst; pCur != CellDataPMapLast; ++ pCur )
      {
        os << pCur->GetVtkRank() << " " << pCur->GetName()
           << " " << pCur->GetVtkDataType() << std::endl;
        if( pCur->DefaultLookupTable() )
          os << "LOOKUP_TABLE default" << std::endl;
        pCur->OutputToVtk( os );
        std::cout << "done with a cell data block " << std::endl;
      }
    }
    std::cout << "Before field data " << std::endl;
    if( FieldDataPMapFirst != FieldDataPMapLast )
    {
      os << "FIELD Grain_Field_Data " << (FieldDataPMapLast - FieldDataPMapFirst) << std::endl;
      for( PropertyMapIterator pCur = FieldDataPMapFirst;
           pCur != FieldDataPMapLast; ++ pCur )
      {
        std::cout << pCur->GetName() << std::endl;
        os << pCur->GetName() << " " << pCur->GetNumTuples()
           << " " << c3t3.number_of_cells()
           << " " << pCur->GetVtkDataType() << std::endl;
        
        pCur->OutputToVtk( os );
      }
    }
    std::cout << "Done writing " << std::endl;
  }
  
  
}  // end implementation of pandora


#endif
