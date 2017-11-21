#include <iostream>
#include <cmath>
#include <fstream>


int main()
{
  int nXC = 25;
  int nYC = 25;
  int nZC = 25;
  float r = 20;
  
  std::ofstream os("TestSphere.dx");
  for( int k = 0; k < 50; k ++ )
    for( int j = 0; j < 50; j ++ )
      for( int i = 0; i < 50; i ++ )
        
      {
        float CurrentR2 = ( i - nXC ) * ( i - nXC )
          + ( j - nYC ) * ( j - nYC )
          + ( k - nZC ) * ( k - nZC );
        
        if(  CurrentR2  < r * r )
        {
          if( i < nXC && j <= nYC )
            os << 5 << " ";
          if( i < nXC && j > nYC )
            os << 4 << " ";
          if( i >= nXC && j > nYC )
            os << 2 << " ";
          if( i >= nXC && j <= nYC )
            os << 1 << " ";
          
        }
        else
          os << -1 << " ";
        
      }
  os.close();
}
