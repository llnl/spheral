#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
QuarticSplineKernel< Dim<1> >::QuarticSplineKernel():
  Kernel<Dim<1>, QuarticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(5.0/8.0);
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
QuarticSplineKernel< Dim<2> >::QuarticSplineKernel():
  Kernel<Dim<2>, QuarticSplineKernel< Dim<2> > >() {
  setVolumeNormalization(5.0/(4.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
QuarticSplineKernel< Dim<3> >::QuarticSplineKernel():
  Kernel<Dim<3>, QuarticSplineKernel< Dim<3> > >() {
  setVolumeNormalization(105.0/(128.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
QuarticSplineKernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  if (etaij < 2.0) {
    return this->volumeNormalization()*Hdet*(1.0 -
                                             1.5*etaij*etaij +
                                             etaij*etaij*etaij -
                                             3.0/16.0*FastMath::pow4(etaij));
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
QuarticSplineKernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  if (etaij < 2.0) {
    return  -this->volumeNormalization()*Hdet*(-3.0*etaij +
                                               3.0*etaij*etaij - 
                                               0.75*FastMath::pow3(etaij));
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
QuarticSplineKernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  if (etaij < 2.0) {
    return this->volumeNormalization()*Hdet*(-3.0 + 6.0*etaij - 
                                             2.25*etaij*etaij);
  } else {
    return 0.0;
  }
}

}
