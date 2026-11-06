#ifndef __Spheral_HExtent__
#define __Spheral_HExtent__

#include "Utilities/DBC.hh"
#include <cmath>

namespace Spheral {

// Axis-aligned bounding extent of the compact support described by H.
template<typename Dimension>
inline typename Dimension::Vector
HExtent(const typename Dimension::SymTensor& H, const double kernelExtent) {
  using Vector = typename Dimension::Vector;
  const auto Hdet = H.Determinant();
  CHECK(Hdet > 0.0);
  if constexpr (Dimension::nDim == 1) {
    return Vector(kernelExtent/H.xx());
  } else if constexpr (Dimension::nDim == 2) {
    const auto M = H.square();
    return kernelExtent/Hdet*Vector(std::sqrt(M.yy()), std::sqrt(M.xx()));
  } else {
    const auto M = H.square();
    return kernelExtent/Hdet*Vector(std::sqrt(M.yy()*M.zz() - M.yz()*M.zy()),
                                    std::sqrt(M.xx()*M.zz() - M.xz()*M.zx()),
                                    std::sqrt(M.xx()*M.yy() - M.xy()*M.yx()));
  }
}

}

#endif
