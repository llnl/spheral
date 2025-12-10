#include <math.h>
#include "Utilities/DBC.hh"
#include "Geometry/GeomVector.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return as a reference descendent class.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
Descendant&
Kernel<Dimension, Descendant>::asDescendant() const {
  return static_cast<Descendant&>(const_cast<Kernel<Dimension, Descendant>&>(*this));
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance or position.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::operator()(const double& etaij, 
                                          const typename Dimension::Scalar& Hdet) const {
  REQUIRE(etaij >= 0.0);
  return asDescendant().kernelValue(etaij, Hdet);
}

template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::operator()(const typename Dimension::Vector& etaj,
                                          const typename Dimension::Vector& etai,
                                          const typename Dimension::Scalar& Hdet) const {
  return asDescendant().kernelValue((etai - etaj).magnitude(), Hdet);
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::grad(const double& etaij,
                                    const typename Dimension::Scalar& Hdet) const {
  REQUIRE(etaij >= 0.0);
  return asDescendant().gradValue(etaij, Hdet);
}

template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::grad(const typename Dimension::Vector& etaj,
                                    const typename Dimension::Vector& etai,
                                    const typename Dimension::Scalar& Hdet) const {
  return asDescendant().gradValue((etai - etaj).magnitude(), Hdet);
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::grad2(const double& etaij,
                                     const typename Dimension::Scalar& Hdet) const {
  REQUIRE(etaij >= 0.0);
  return asDescendant().grad2Value(etaij, Hdet);
}

template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::grad2(const typename Dimension::Vector& etaj,
                                     const typename Dimension::Vector& etai,
                                     const typename Dimension::Scalar& Hdet) const {
  return asDescendant().grad2Value((etai - etaj).magnitude(), Hdet);
}

//------------------------------------------------------------------------------
// All kernels must redefine this method to give the W value for a given
// distance and H determinant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::kernelValue(double etaij, const double Hdet) const {
  return asDescendant().kernelValue(etaij, Hdet);
}

//------------------------------------------------------------------------------
// All kernels must redefine this method to give the gradient value for a given
// distance and H determinant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::gradValue(double etaij, const double Hdet) const {
  return asDescendant().gradValue(etaij, Hdet);
}

//------------------------------------------------------------------------------
// All kernels must redefine this method to give the second derivative value for
// a given distance and H determinant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::grad2Value(double etaij, const double Hdet) const {
  return asDescendant().grad2Value(etaij, Hdet);
}

//------------------------------------------------------------------------------
// Return the volume normalization.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::volumeNormalization() const {
  return mVolumeNormalization;
}

//------------------------------------------------------------------------------
// Set the volume normalization.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
void
Kernel<Dimension, Descendant>::setVolumeNormalization(double volumeNormalization) {
  REQUIRE(volumeNormalization > 0.0);
  mVolumeNormalization = volumeNormalization;
}

//------------------------------------------------------------------------------
// Return the kernel extent.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::kernelExtent() const {
  return mKernelExtent;
}

//------------------------------------------------------------------------------
// Set the kernel extent.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
void
Kernel<Dimension, Descendant>::setKernelExtent(double extent) {
  REQUIRE(extent > 0.0);
  mKernelExtent = extent;
}

//------------------------------------------------------------------------------
// Return the inflection point.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
double
Kernel<Dimension, Descendant>::inflectionPoint() const {
  return mInflectionPoint;
}

//------------------------------------------------------------------------------
// Set the inflection point.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE inline
void
Kernel<Dimension, Descendant>::setInflectionPoint(double x) {
  mInflectionPoint = x;
}

}
