//-----------------------------------------------------------
//
//  SimpleC3T3Writer.h
//  Author:  Frankie Li
//
//  e-mail:   sfli@cmu.edu
//  
//  Purpose:  Simple Complex_3_in_triangulation_3 writer/reader
//
//  Note:  Some parts of this file is modeled after File_medit.h
//  from CGAL/IO (as included)
//
//-----------------------------------------------------------


#ifndef SIMPLE_C3T3_WRITER_
#define SIMPLE_C3T3_WRITER_

#include <iostream>
#include <map>

namespace Pandora
{

  namespace Details
  {
    
  }  // Details

  namespace C3T3Writer
  {
  //--------------------------------------
  //  WriteC3T3
  //--------------------------------------
  template< class Complex_3_In_Triangulation_3 >
  void WriteC3T3( std::ostream & os,
                  const Complex_3_In_Triangulation_3 & c3t3 );
  
  //--------------------------------------
  //  ReadC3T3
  //--------------------------------------
  template< class Complex_3_In_Triangulation_3 >
  void ReadC3T3( std::ostream & os,
                 Complex_3_In_Triangulation_3 & c3t3 );
  }
}


namespace Pandora
{

  namespace Details
  {
    
  }  // Details

  namespace C3T3Writer
  {
    //--------------------------------------
    //  WriteC3T3
    //--------------------------------------
    template< class Complex_3_In_Triangulation_3 >
    void WriteC3T3( std::ostream & os,
                    const Complex_3_In_Triangulation_3 & c3t3 )
    {
      typedef typename Complex_3_In_Triangulation_3::Triangulation Tr;
      typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
      typedef typename Tr::Vertex_iterator Vertex_iterator;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Point Point_3;
      
      std::map< Vertex_handle, int > VertexIDPMap;
      
      const Tr & tr = c3t3.triangulation();
      int nCurID = 0;

      os << tr.number_of_vertices() << std::endl;
      for( Vertex_iterator vit = tr.vertices_begin();
           vit != tr.vertices_end();
           ++vit)
      {
        Vertex_handle v = vit;
        os << *v << " " << nCurID << std::endl;
        VertexIDPMap[v] = nCurID;
        ++ nCurID;
      }
      
      for( Cell_iterator pCur = c3t3_copy.cells_begin();
           pCur != c3t3_copy.cells_end(); ++pCur )
      {
        for( int i = 0; i < 3; i ++ )
          os << VertexIDPMap[ pCur->vertex( i ) ] << " "; 
        os << pCur->subdomain_index << std::endl;
      }
      
   

      
    }
    
    //--------------------------------------
    //  ReadC3T3
    //--------------------------------------
    template< class Complex_3_In_Triangulation_3 >
    void ReadC3T3( std::ostream & os,
                   Complex_3_In_Triangulation_3 & c3t3 )
    {
    }
  }
}


#endif
