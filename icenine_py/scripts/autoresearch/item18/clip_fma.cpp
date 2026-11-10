// Review item 18 (D-20261109-N4): does the C++ build fuse the clip's multiply-add?
// The clip's intersection is `p0 + slope * (m - p0)` in single precision
// (XDM++/libXDM/SutherlandHodgman.h, BoundaryHor::Intersect and BoundaryVert::Intersect).
// Compile to assembly with the project's Release flags and count the instructions, from the
// repository root:
//   c++ -std=c++11 -O3 -DNDEBUG -Wno-everything -IXDM++/libXDM -S -o clip_fma.s \
//       icenine_py/scripts/autoresearch/item18/clip_fma.cpp
//   grep -E "fmadd|fmul|fadd" clip_fma.s
// Result of 2026-11-09 (Apple clang 21.0.0, arm64): one `fmadd` per function, no separate
// multiply and add. With -ffp-contract=off: one `fmul` and one `fadd` per function.
#include <functional>
#include <vector>
#define private public
#include "SutherlandHodgman.h"
#undef private
typedef SutherlandHodgman SH;
__attribute__((noinline)) PointF vert(const SH::BoundaryVert<std::less<REAL> >& b, const PointF& p0, const PointF& p1) { return b.Intersect(p0, p1); }
__attribute__((noinline)) PointF hor(const SH::BoundaryHor<std::less<REAL> >& b, const PointF& p0, const PointF& p1) { return b.Intersect(p0, p1); }
