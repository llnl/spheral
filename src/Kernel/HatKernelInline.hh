#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
HatKernel< Dim<1> >::HatKernel(double eta0, double W0):
  Kernel<Dim<1>, HatKernel< Dim<1> > >(),
  mEta0(eta0),
  mW0(W0),
  mSlope(0.0) {
  REQUIRE(eta0 > 0.0);
  REQUIRE(W0 > 0.0);
  mSlope = W0/eta0;
  setVolumeNormalization(1.0/(W0*eta0));
  setKernelExtent(eta0);
  setInflectionPoint(0.0);
}

template<>
inline
HatKernel< Dim<2> >::HatKernel(double eta0, double W0):
  Kernel<Dim<2>, HatKernel< Dim<2> > >(),
  mEta0(eta0),
  mW0(W0),
  mSlope(0.0) {
  REQUIRE(eta0 > 0.0);
  REQUIRE(W0 > 0.0);
  mSlope = W0/eta0;
  setVolumeNormalization(3.0/(M_PI*W0*eta0*eta0));
  setKernelExtent(eta0);
  setInflectionPoint(0.0);
}

template<>
inline
HatKernel< Dim<3> >::HatKernel(double eta0, double W0):
  Kernel<Dim<3>, HatKernel< Dim<3> > >(),
  mEta0(eta0),
  mW0(W0),
  mSlope(0.0) {
  REQUIRE(eta0 > 0.0);
  REQUIRE(W0 > 0.0);
  mSlope = W0/eta0;
  setVolumeNormalization(3.0/(M_PI*W0*eta0*eta0*eta0));
  setKernelExtent(eta0);
  setInflectionPoint(0.0);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
HatKernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  if (etaij < mEta0) {
    return this->volumeNormalization()*Hdet*(mW0 - mSlope*etaij);
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
HatKernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  if (etaij < mEta0) {
    return  -this->volumeNormalization()*Hdet*mSlope;
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
HatKernel<Dimension>::grad2Value(double /*etaij*/, const double /*Hdet*/) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Const access to the intercept values.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double 
HatKernel<Dimension>::eta0() const {
  return mEta0;
}

template<typename Dimension>
inline
double 
HatKernel<Dimension>::W0() const {
  return mW0;
}

}
