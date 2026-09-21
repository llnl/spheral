#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
NeighborBase(const NeighborSearchType searchType,
             const typename Dimension::Scalar kernelExtent):
  mSearchType(searchType),
  mKernelExtent(kernelExtent) {
}

//------------------------------------------------------------------------------
// Calculate the maximum radial extent of a given smoothing tensor.
//------------------------------------------------------------------------------
// Make the general case the SPH approximation: assume the argument is a 
// scalar representing the inverse smoothing scale.
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
typename Dimension::Vector
NeighborBase<Dimension>::
HExtent(const typename Dimension::Scalar& H, 
        const typename Dimension::Scalar& kernelExtent) {
  REQUIRE(H > 0.0);
  return Vector(kernelExtent/H);
}

//------------------------------------------------------------------------------
// Specializations for symmetric tensor smoothing transformation (ASPH).
// The ASPH tensor has units of inverse length.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
Dim<1>::Vector
NeighborBase< Dim<1> >::
HExtent(const Dim<1>::SymTensor& H,
        const Dim<1>::Scalar& kernelExtent) {
  REQUIRE(H.Determinant() > 0.0);
  return Vector(kernelExtent/H.xx());
}

template<>
SPHERAL_HOST_DEVICE
inline
Dim<2>::Vector
NeighborBase< Dim<2> >::
HExtent(const Dim<2>::SymTensor& H,
        const Dim<2>::Scalar& kernelExtent) {
  const auto Hdet = H.Determinant();
  const auto M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Vector(std::sqrt(M.yy()), std::sqrt(M.xx()));
}
  
template<>
SPHERAL_HOST_DEVICE
inline
Dim<3>::Vector
NeighborBase< Dim<3> >::
HExtent(const Dim<3>::SymTensor& H,
        const Dim<3>::Scalar& kernelExtent) {
  const auto Hdet = H.Determinant();
  const auto M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Vector(std::sqrt(M.yy()*M.zz() - M.yz()*M.zy()),
                                  std::sqrt(M.xx()*M.zz() - M.xz()*M.zx()),
                                  std::sqrt(M.xx()*M.yy() - M.xy()*M.yx()));
}

}
