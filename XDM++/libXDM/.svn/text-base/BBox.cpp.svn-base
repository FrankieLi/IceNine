#include "BBox.h"

namespace PBRMath
{

  
BBox2D Union(const BBox2D &b, const Point &p) {
	BBox2D ret = b;
	ret.pMin.x = min(b.pMin.x, p.x);
	ret.pMin.y = min(b.pMin.y, p.y);
	ret.pMax.x = max(b.pMax.x, p.x);
	ret.pMax.y = max(b.pMax.y, p.y);
	return ret;
}

BBox2D Union(const BBox2D &b, const BBox2D &b2) {
	BBox2D ret;
	ret.pMin.x = min(b.pMin.x, b2.pMin.x);
	ret.pMin.y = min(b.pMin.y, b2.pMin.y);
	ret.pMax.x = max(b.pMax.x, b2.pMax.x);
	ret.pMax.y = max(b.pMax.y, b2.pMax.y);
	return ret;
}

//------------------------------------------------
//  operator<<
//------------------------------------------------
std::ostream & operator<< ( std::ostream & os, const BBox2D &m )
{
  os << m.pMin.x << " "
     << m.pMin.y << " "
     << m.pMax.x << " "
     << m.pMax.y; 
  return os;
}


}
