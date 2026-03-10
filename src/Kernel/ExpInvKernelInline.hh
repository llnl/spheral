#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE inline
ExpInvKernel< Dim<1> >::ExpInvKernel():
  Kernel<Dim<1>, ExpInvKernel< Dim<1> > >() {
  setVolumeNormalization(12.570510753090891498);
  // setVolumeNormalization(588.13226316830821361);
  setKernelExtent(3.0);
  setInflectionPoint(0.0);
}

template<>
SPHERAL_HOST_DEVICE inline
ExpInvKernel< Dim<2> >::ExpInvKernel():
  Kernel<Dim<2>, ExpInvKernel< Dim<2> > >() {
  setVolumeNormalization(14.453960393997373757*M_PI);
  //  setVolumeNormalization(28.256642009101042845*M_PI);
  setKernelExtent(3.0);
  setInflectionPoint(0.0);
}

template<>
SPHERAL_HOST_DEVICE inline
ExpInvKernel< Dim<3> >::ExpInvKernel():
  Kernel<Dim<3>, ExpInvKernel< Dim<3> > >() {
  setVolumeNormalization(17.824612984468913623*M_PI);
  // setVolumeNormalization(20.151784602855109085*M_PI);
  setKernelExtent(3.0);
  setInflectionPoint(0.0);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
double
ExpInvKernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  const double A = this->volumeNormalization();
  const double x = std::abs(etaij) + 0.5;
  return A*Hdet*exp(1.0/x);
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
double
ExpInvKernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  const double A = this->volumeNormalization();
  const double x = std::abs(etaij) + 0.5;
  return -A*Hdet/(x*x)*exp(1.0/x);
}

//------------------------------------------------------------------------------
// Return the second derivative for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
double
ExpInvKernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  const double A = this->volumeNormalization();
  const double x = std::abs(etaij) + 0.5;
  return A*Hdet*(1.0/(x*x*x*x) + 1.0/(x*x*x))*exp(1.0/x);
}

}
