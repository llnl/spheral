//---------------------------------Spheral++----------------------------------//
// GeomVector_fwd -- forward declaration for Geometric Vector Class.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomVector_fwd_hh__
#define __Spheral_GeomVector_fwd_hh__

#include "config.hh"

namespace Spheral {
  template<int nDim> class GeomVector;

  // Forward declare the global functions.
  template<int nDim> SPHERAL_HOST_DEVICE GeomVector<nDim> operator*(const double val, const GeomVector<nDim>& vec);
  template<int nDim> SPHERAL_HOST_DEVICE GeomVector<nDim> elementWiseMin(const GeomVector<nDim>& lhs,
                                                                         const GeomVector<nDim>& rhs);
  template<int nDim> SPHERAL_HOST_DEVICE GeomVector<nDim> elementWiseMax(const GeomVector<nDim>& lhs,
                                                                         const GeomVector<nDim>& rhs);

  template<> SPHERAL_HOST_DEVICE GeomVector<1> elementWiseMin(const GeomVector<1>& lhs,
                                                              const GeomVector<1>& rhs);
  template<> SPHERAL_HOST_DEVICE GeomVector<2> elementWiseMin(const GeomVector<2>& lhs,
                                                              const GeomVector<2>& rhs);
  template<> SPHERAL_HOST_DEVICE GeomVector<3> elementWiseMin(const GeomVector<3>& lhs,
                                                              const GeomVector<3>& rhs);

  template<> SPHERAL_HOST_DEVICE GeomVector<1> elementWiseMax(const GeomVector<1>& lhs,
                                                              const GeomVector<1>& rhs);
  template<> SPHERAL_HOST_DEVICE GeomVector<2> elementWiseMax(const GeomVector<2>& lhs,
                                                              const GeomVector<2>& rhs);
  template<> SPHERAL_HOST_DEVICE GeomVector<3> elementWiseMax(const GeomVector<3>& lhs,
                                                              const GeomVector<3>& rhs);

  template<int nDim> std::istream& operator>>(std::istream& is, GeomVector<nDim>& vec);
  template<int nDim> std::ostream& operator<<(std::ostream& os, const GeomVector<nDim>& vec);
}

#endif
