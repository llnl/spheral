#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// W4SplineKernel<Dimension>::W4SplineKernel():
//   Kernel<Dimension, W4SplineKernel<Dimension> >() {
// }

template<>
inline
W4SplineKernel< Dim<1> >::W4SplineKernel():
  Kernel<Dim<1>, W4SplineKernel< Dim<1> > >() {
  setVolumeNormalization(1.0);
  setKernelExtent(2.0);
}

template<>
inline
W4SplineKernel< Dim<2> >::W4SplineKernel():
  Kernel<Dim<2>, W4SplineKernel< Dim<2> > >() {
  setVolumeNormalization(30.0/7.0*M_PI);
  setKernelExtent(2.0);
}

template<>
inline
W4SplineKernel< Dim<3> >::W4SplineKernel():
  Kernel<Dim<3>, W4SplineKernel< Dim<3> > >() {
  setVolumeNormalization(5.0/6.0*M_PI);
  setKernelExtent(2.0);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
W4SplineKernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  const double& eta = etaij;
  const double eta2 = eta*eta;
  const double eta3 = eta*eta2;
  if (eta < 1.0) {
    return this->volumeNormalization()*Hdet*(1.0 - 2.5*eta2 + 1.5*eta3);
  } else if (etaij < 2.0) {
    return this->volumeNormalization()*Hdet*(2.0 - 4.0*etaij +
                                       2.5*eta2 - 0.5*eta3);
//     return this->volumeNormalization()*Hdet*0.5*pow(2.0 - etaij, 2)*(1.0 - etaij);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
W4SplineKernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  const double& eta = etaij;
  const double eta2 = eta*eta;
  if (eta < 1.0) {
    return  -this->volumeNormalization()*Hdet*(5.0*eta - 4.5*eta2);
  } else if (etaij < 2.0) {
    return -this->volumeNormalization()*Hdet*(4.0 - 5.0*etaij + 1.5*eta2);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
W4SplineKernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  if (etaij < 1.0) {
    return  -this->volumeNormalization()*Hdet*(5.0 - 9.0*etaij);
  } else if (etaij < 2.0) {
    return this->volumeNormalization()*Hdet*(5.0 - 3.0*etaij);
  } else {
    return 0.0;
  }
}

}
