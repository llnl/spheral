#include "Geometry/Dimension.hh"

#include <math.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// PiGaussianKernel<Dimension>::PiGaussianKernel():
//   Kernel<Dimension, PiGaussianKernel<Dimension> >(),
//     mK(0.0),
//     mKV(0.0) {
// }

template<>
inline
PiGaussianKernel< Dim<1> >::PiGaussianKernel():
  Kernel<Dim<1>, PiGaussianKernel< Dim<1> > >(),
    mK(0.0),
    mKV(0.0) {
  setVolumeNormalization(2.0/3.6256099082);
  setKernelExtent(2.0);
}

template<>
inline
PiGaussianKernel< Dim<2> >::PiGaussianKernel():
  Kernel<Dim<2>, PiGaussianKernel< Dim<2> > >(),
    mK(0.0),
    mKV(0.0) {
  setVolumeNormalization(2.0/std::pow(M_PI, 1.5));
  setKernelExtent(2.0);
}

template<>
inline
PiGaussianKernel< Dim<3> >::PiGaussianKernel():
  Kernel<Dim<3>, PiGaussianKernel< Dim<3> > >(),
    mK(0.0),
    mKV(0.0) {
  setVolumeNormalization(1.0/(M_PI*1.2254167024));
  setKernelExtent(2.0);
}

//------------------------------------------------------------------------------
// Construct with the given value for K.
//------------------------------------------------------------------------------
template<>
inline
PiGaussianKernel< Dim<1> >::PiGaussianKernel(double K):
  Kernel<Dim<1>, PiGaussianKernel< Dim<1> > >(),
  mK(K) {
  setVolumeNormalization(2.0/3.6256099082);
  setKernelExtent(2.0);
  mKV = std::pow(K, 1.0/4.0);
}

template<>
inline
PiGaussianKernel< Dim<2> >::PiGaussianKernel(double K):
  Kernel<Dim<2>, PiGaussianKernel< Dim<2> > >(),
  mK(K) {
  setVolumeNormalization(2.0/std::pow(M_PI, 1.5));
  setKernelExtent(2.0);
  mKV = std::pow(K, 2.0/4.0);
}

template<>
inline
PiGaussianKernel< Dim<3> >::PiGaussianKernel(double K):
  Kernel<Dim<3>, PiGaussianKernel< Dim<3> > >(),
  mK(K) {
  setVolumeNormalization(1.0/(M_PI*1.2254167024));
  setKernelExtent(2.0);
  mKV = std::pow(K, 3.0/4.0);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PiGaussianKernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);
  CHECK(mK > 0.0);
  CHECK(mKV > 0.0);
  return this->volumeNormalization()*mKV*Hdet*exp(-mK*pow(etaij, 4));
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PiGaussianKernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);
  CHECK(mK > 0.0);
  CHECK(mKV > 0.0);
  return -4.0*mK*pow(etaij, 3)*kernelValue(etaij, Hdet);
}

//------------------------------------------------------------------------------
// Return the second derivative for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PiGaussianKernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);
  CHECK(mK > 0.0);
  CHECK(mKV > 0.0);
  return 4.0*mK*etaij*etaij*(4.0*mK*pow(etaij, 4) - 3.0)*
    kernelValue(etaij, Hdet);
}

//------------------------------------------------------------------------------
// Return the exponential constant K.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PiGaussianKernel<Dimension>::getK() const {
  return mK;
}

//------------------------------------------------------------------------------
// Set the exponential constant K.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
PiGaussianKernel<Dimension>::setK(double K) {
  CHECK(K > 0.0);
  mK = K;
  mKV = pow(K, double(Dimension::nDim)/4.0);
}

}
