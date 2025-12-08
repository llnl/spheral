#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
WendlandC2Kernel< Dim<1> >::WendlandC2Kernel():
  Kernel<Dim<1>, WendlandC2Kernel< Dim<1> > >() {
  setVolumeNormalization(5.0/4.0);
  setKernelExtent(1.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
WendlandC2Kernel< Dim<2> >::WendlandC2Kernel():
  Kernel<Dim<2>, WendlandC2Kernel< Dim<2> > >() {
  setVolumeNormalization(7.0/(M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(1.0/4.0);
}

template<>
inline
WendlandC2Kernel< Dim<3> >::WendlandC2Kernel():
  Kernel<Dim<3>, WendlandC2Kernel< Dim<3> > >() {
  setVolumeNormalization(21.0/(2.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(1.0/4.0);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC2Kernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);

    if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(pow(1.0-etaij,3)*(1.0+3.0*etaij))*(etaij < 1.0);
    else
      return this->volumeNormalization()*Hdet*(pow(1.0-etaij,4)*(1.0+4.0*etaij))*(etaij < 1.0);

}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC2Kernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);

    if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(-12.0*pow(1.0-etaij,2)*etaij)*(etaij < 1.0);
    else
      return  this->volumeNormalization()*Hdet*(20.0*pow(etaij-1.0,3)*etaij)*(etaij < 1.0);


}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC2Kernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);


    const double eta2 = etaij*etaij;
    if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(-12.0*(3.0*eta2-4.0*etaij+1.0))*(etaij < 1.0);
    else
      return this->volumeNormalization()*Hdet*(20.0*pow(etaij-1.0,2)*
                                               (4.0*etaij-1.0))*(etaij < 1.0);

}

}
