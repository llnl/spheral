//------------------------------------------------------------------------------
// Compute the second moment about the give position for a polytope
//
// Note, these methods currently assume the polytopes are convex.
//------------------------------------------------------------------------------

#include "SmoothingScale/polySecondMoment.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {

using FastMath::pow2;

//------------------------------------------------------------------------------
// 1D -- nothing to do
//------------------------------------------------------------------------------
Dim<1>::SymTensor
polySecondMoment(const Dim<1>::FacetedVolume& poly,
                 const Dim<1>::Vector& center) {
  return Dim<1>::SymTensor(1);
}

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
Dim<2>::SymTensor
polySecondMoment(const Dim<2>::FacetedVolume& poly,
                 const Dim<2>::Vector& center) {
  Dim<2>::SymTensor result;
  const auto& facets = poly.facets();
  for (const auto& f: facets) {
    const auto v1 = f.point1() - center;
    const auto v2 = f.point2() - center;
    const auto thpt = std::abs(v1.x()*v2.y() - v2.x()*v1.y())/12.0;
    result[0] += (v1.x()*v1.x() + v1.x()*v2.x() + v2.x()*v2.x())*thpt;
    result[1] += (v1.x()*v1.y() + v2.x()*v2.y() + 0.5*(v2.x()*v1.y() + v1.x()*v2.y()))*thpt;
    result[2] += (v1.y()*v1.y() + v1.y()*v2.y() + v2.y()*v2.y())*thpt;
  }
  return result;
}

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
Dim<3>::SymTensor
polySecondMoment(const Dim<3>::FacetedVolume& poly,
                 const Dim<3>::Vector& center) {
  using Scalar = Dim<3>::Scalar;
  using Vector = Dim<3>::Vector;
  using SymTensor = Dim<3>::SymTensor;
  SymTensor result;
  std::vector<std::array<Vector, 3>> tris;
  Vector v1, v2, v3;
  Scalar thpt, x1, x2, x3, y1, y2, y3, z1, z2, z3;
  const auto& facets = poly.facets();
  for (const auto& f: facets) {
    f.decompose(tris);
    for (const auto& tri: tris) {
      v1 = tri[0] - center;
      v2 = tri[1] - center;
      v3 = tri[2] - center;
      x1 = v1.x(); y1 = v1.y(); z1 = v1.z();
      x2 = v2.x(); y2 = v2.y(); z2 = v2.z();
      x3 = v3.x(); y3 = v3.y(); z3 = v3.z();
      thpt = std::abs(x3*y2*z1 - x2*y3*z1 - x3*y1*z2 + x1*y3*z2 + x2*y1*z3 - x1*y2*z3);
      result[0] += thpt * pow2(x1 + x2) + x3*(x1 + x2 + x3);                                              // xx
      result[1] += thpt * 0.5*(x1*(2.0*y1 + y2 + y3) + x2*(y1 + 2.0*y2 + y3) + x3*(y1 + y2 + 2.0*y3));    // xy
      result[2] += thpt * 0.5*(x1*(2.0*z1 + z2 + z3) + x2*(z1 + 2.0*z2 + z3) + x3*(z1 + z2 + 2.0*z3));    // xz
      result[3] += thpt * pow2(y1 + y2) + y3*(y1 + y2 + y3);                                              // yy
      result[4] += thpt * 0.5*(y1*(2.0*z1 + z2 + z3) + y2*(z1 + 2.0*z2 + z3) + y3*(z1 + z2 + 2.0*z3));    // yz
      result[5] += thpt * pow2(z1 + z2) + z3*(z1 + z2 + z3);                                              // zz
    }
  }
  result /= 60.0;
  return result;
}

}
