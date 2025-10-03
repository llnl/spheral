#include "Geometry/Dimension.hh"

#include <math.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// SuperGaussianKernel<Dimension>::SuperGaussianKernel():
//   Kernel<Dimension, SuperGaussianKernel<Dimension> >() {
// }

template<>
inline
SuperGaussianKernel< Dim<1> >::SuperGaussianKernel():
  Kernel<Dim<1>, SuperGaussianKernel< Dim<1> > >() {
  setVolumeNormalization(1.0/std::sqrt(M_PI));
  setKernelExtent(3.0);
}

template<>
inline
SuperGaussianKernel< Dim<2> >::SuperGaussianKernel():
  Kernel<Dim<2>, SuperGaussianKernel< Dim<2> > >() {
  setVolumeNormalization(1.0/M_PI);
  setKernelExtent(3.0);
}

template<>
inline
SuperGaussianKernel< Dim<3> >::SuperGaussianKernel():
  Kernel<Dim<3>, SuperGaussianKernel< Dim<3> > >() {
  setVolumeNormalization(1.0/std::pow(M_PI, 1.5));
  setKernelExtent(3.0);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SuperGaussianKernel<Dimension>::~SuperGaussianKernel() {
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SuperGaussianKernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);
  double etaij2 = etaij*etaij;
  return this->volumeNormalization()*Hdet*(mKW - etaij2)*exp(-etaij2);
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SuperGaussianKernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);
  double etaij2 = etaij*etaij;
  return -2.0*this->volumeNormalization()*Hdet*
    etaij*(mKW + 1.0 - etaij2)*exp(-etaij2);
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SuperGaussianKernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);
  double etaij2 = etaij*etaij;
  return 2.0*this->volumeNormalization()*Hdet*
    ((2.0*mKW + 5.0)*etaij2 - 4.0*etaij2*etaij2 - mKW - 1.0)*
    exp(-etaij2);
}

#if defined(SPHERAL1D)
template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;
#endif

#if defined(SPHERAL2D)
template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;
#endif

#if defined(SPHERAL3D)
template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;
template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;
#endif

}
