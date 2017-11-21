#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

#include "Voxel.h"
#include "Quaternion.h"
#include <fstream>

//--------------------------------------------------------------------------------
//  Delauney Triagulation Setup
//
//
//--------------------------------------------------------------------------------
template < class GT, class Vb = CGAL::Triangulation_vertex_base_3<GT> >
class My_vertex_base
  : public Vb
{
public:
  typedef typename Vb::Vertex_handle  Vertex_handle;
  typedef typename Vb::Cell_handle    Cell_handle;
  typedef typename Vb::Point          Point;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_base<GT, Vb2>                        Other;
  };

  My_vertex_base() {}

  My_vertex_base(const Point& p)
    : Vb(p) {}

  My_vertex_base(const Point& p, Cell_handle c)
    : Vb(p, c) {}

  SVoxel Voxel;

};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_data_structure_3<My_vertex_base<K> >    Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;

typedef Delaunay::Vertex_handle    Vertex_handle;
typedef Delaunay::Point            Point;
//--------------------------------------------------------------------------------




//--------------------------------------------------------------------------------
//
//  QuadNode
//--------------------------------------------------------------------------------
struct QuadNode
{
  Point p;
  Point COM;
  std::vector<Point> Vertices;
  std::vector<int> IDList;
};

//--------------------------------------------------------------------------------
// QuadNodeIdentification
//
//  Quad nodes are identified by finding "cells," or tetrahedrons in the Delauney
//  triangulation.  Since the DT is the dual of the Voronoi diagram, then the tets of
//  a DT corresponds to a vertex in the Voronoi.  All vertices with 4 voronoi cells
//  corresponds to a quad node.  In the DT, the 4 voronoi cells show up as 4 vertices
//  of the tetrahedron with specified grain ID.  If all four IDs are different, then
//  this quad node will be a junction between 4 materials.
//
//--------------------------------------------------------------------------------
std::vector<QuadNode> QuadNodeIdentification( const Delaunay & T )
{
  typedef Delaunay::Vertex_iterator Vertex_iterator;
  typedef Delaunay::Cell_iterator Cell_iterator;

  std::vector<QuadNode> QuadNodeList;
  
  for( Cell_iterator CellIter = T.cells_begin(); CellIter != T.cells_end();
       ++ CellIter )
  {
    if( !T.is_infinite( CellIter ) )
    {
      std::set<int> IDSet;
      bool bInterior = true;;
      for( int i = 0; i < 4; i ++ )
      {
        IDSet.insert( CellIter->vertex(i)->Voxel.nGrainID );
        if( CellIter->vertex(i)->Voxel.nGrainID <= 0 )
          bInterior = false;
      }
      if( IDSet.size() == 4 && bInterior )
      {
        QuadNode QN;
        QN.p  = T.dual( CellIter );

        GeneralLib::SVector3 v(0, 0, 0);
        QN.COM = Point(0, 0, 0);
        for( int i = 0; i < 4; i ++ )
        {
          QN.Vertices.push_back( CellIter->vertex(i)->point() );
          QN.IDList.push_back( CellIter->vertex(i)->Voxel.nGrainID );
          v.m_fX +=  CellIter->vertex(i)->point().x();
          v.m_fY +=  CellIter->vertex(i)->point().y();
          v.m_fZ +=  CellIter->vertex(i)->point().z();
        }
        v /= 4.;
        QN.COM = Point( v.m_fX, v.m_fY, v.m_fZ );
        QuadNodeList.push_back( QN );
      }
    }
  }

  return QuadNodeList;
}


//--------------------------------------------------------------------------------
//  ProcessDelaunayTriangulation
//
//
//--------------------------------------------------------------------------------
void ProcessDelaunayTriagulation( const std::string & FilenamePrefix,
                                  const std::string & FilenamePostfix,
                                  int Start, int Stop,
                                  Delaunay & T )
{
  std::ofstream DebugOS("DebugInput.csv");
    std::stringstream Filename;
    //    Filename << FilenamePrefix << i << FilenamePostfix;
    Filename << FilenamePrefix << FilenamePostfix;

    std::cout << "Reading: " << Filename.str() << std::endl;
    std::ifstream is;
    is.open( Filename.str().c_str() );

    if( !is.good() )
    {
      std::cout << "Error:  file " << Filename.str() << " cannot be opend " << std::endl;
      exit(0);
    }

    //------------------------------ DEBUG
    std::string line;
    std::set<Point> unique_points;
    int NumPoints = 0;  
    while( std::getline( is, line ) )
    {
      if( line.size() > 2 )
      {
        std::stringstream ss;
        ss << line;
        int i, j, k, nPhase, nGrainID;
        float fConfidence;
        Point p;
        GeneralLib::SQuaternion q;
        
        ss >> std::skipws >> i
           >> std::skipws >> j
           >> std::skipws >> k
          
           >> std::skipws >> p
          
          
           >> std::skipws >> q.m_fW
           >> std::skipws >> q.m_fX
           >> std::skipws >> q.m_fY
           >> std::skipws >> q.m_fZ
           >> std::skipws >> nPhase
           >> std::skipws >> fConfidence
           >> std::skipws >> nGrainID;
        unique_points.insert( p );
        NumPoints ++;
      }
    }
    if( unique_points.size() != NumPoints )
    {
      std::cout << "Error : there are duplicate points " << std::endl;
      exit( 0 );
    }
    std::cout << "Total number of points " << NumPoints << std::endl;
    is.close();
    is.open( Filename.str().c_str() );
    while( std::getline( is, line ) )
    {
      if( line.size() > 2 )  // sanity check 
      {
        std::stringstream ss;
        ss << line;
        int i, j, k, nPhase, nGrainID;
        float fConfidence;
        Point p;
        GeneralLib::SQuaternion q;
        
        ss >> std::skipws >> i
           >> std::skipws >> j
           >> std::skipws >> k
          
           >> std::skipws >> p
          
          
           >> std::skipws >> q.m_fW
           >> std::skipws >> q.m_fX
           >> std::skipws >> q.m_fY
           >> std::skipws >> q.m_fZ
           >> std::skipws >> nPhase
           >> std::skipws >> fConfidence
           >> std::skipws >> nGrainID;

        Vertex_handle vh      = T.insert( p );
        vh->Voxel.nPhase      = nPhase;
        vh->Voxel.fConfidence = fConfidence;
        vh->Voxel.nGrainID    = nGrainID;

        DebugOS << p << " " << nPhase << " " << fConfidence << " " << nGrainID << " " << std::endl;
      }
    }
    
    // }
  DebugOS.close();
}


int main()
{
  Delaunay T;


  if( T.is_valid() )
  {
    std::cout << "Triangulation is valid delaunay " << std::endl;
  }


  typedef Delaunay::Vertex_iterator Vertex_iterator;
  typedef Delaunay::Cell_iterator Cell_iterator;
  

  std::string InputPrefix;
  std::string InputPostfix;
  int Start;
  int Stop;
  
  
  std::cout << "Input prefix   " << std::endl;
  std::cin >> std::skipws >> InputPrefix;
  
  std::cout << "Input postfix: " << std::endl;
  std::cin >> std::skipws >> InputPostfix;
  
  std::cout << "Start " << std::endl;
  std::cin >> std::skipws >> Start;
  
  std::cout << "Stop " << std::endl;
  std::cin >> std::skipws >> Stop;

  Delaunay  DT;
  ProcessDelaunayTriagulation( InputPrefix,
                               InputPostfix,
                               Start, Stop,
                               DT );
  
  std::vector<QuadNode> QNList = QuadNodeIdentification( DT );
  

  std::cout << "Number of quad nodes: " << QNList.size() << std::endl;


  std::ofstream TestOS("CrapQN.txt");
  for( int i = 0; i < QNList.size(); i ++ )
  {
    TestOS << QNList[i].p << " "
           << QNList[i].IDList.size() << " ";
    for( int j = 0; j < QNList[i].IDList.size(); j ++ )
      TestOS << QNList[i].IDList[j] << " ";

    
    TestOS << QNList[i].COM << " ";
    for( int j = 0; j < QNList[i].IDList.size(); j ++ )
      TestOS << QNList[i].Vertices[j] << " ";
    TestOS << std::endl;
  }
      
      
  
  return 0;
  
}
