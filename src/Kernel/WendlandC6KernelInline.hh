#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
WendlandC6Kernel< Dim<1> >::WendlandC6Kernel():
  Kernel<Dim<1>, WendlandC6Kernel< Dim<1> > >() {
  setVolumeNormalization(55.0/32.0);
  setKernelExtent(1.0);
  setInflectionPoint(0.239983);
}

template<>
inline
WendlandC6Kernel< Dim<2> >::WendlandC6Kernel():
  Kernel<Dim<2>, WendlandC6Kernel< Dim<2> > >() {
  setVolumeNormalization(78.0/(7.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.215286);
}

template<>
inline
WendlandC6Kernel< Dim<3> >::WendlandC6Kernel():
  Kernel<Dim<3>, WendlandC6Kernel< Dim<3> > >() {
  setVolumeNormalization(1365.0/(64.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.215286);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC6Kernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);

    double eta2 = etaij*etaij;
    if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(pow(1.0-etaij,7)*
                                             (1.0+7.0*etaij+19.0*eta2+21.0*etaij*eta2))*(etaij < 1.0);
    else
      return this->volumeNormalization()*Hdet*(pow(1.0-etaij,8)*
                                             (1.0+8.0*etaij+25.0*eta2+32.0*etaij*eta2))*(etaij < 1.0);

}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC6Kernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);

    double eta2 = etaij*etaij;
    if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(-6.0*pow(1.0-etaij,6)*etaij*
                                             (3.0+18.0*etaij+35.0*eta2))*(etaij < 1.0);
    else
      return  this->volumeNormalization()*Hdet*(22.0*pow(etaij-1.0,7)*
                                              etaij*(16.0*eta2+7.0*etaij+1.0))*(etaij < 1.0);


}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC6Kernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);


    const double eta2 = etaij*etaij;
    if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(-18.0*pow(etaij-1.0,5)*
                                             (105.0*eta2*etaij+13.0*eta2-5.0*etaij-1.0))*(etaij < 1.0);
    else
      return this->volumeNormalization()*Hdet*(22.0*pow(etaij-1.0,6)*
                                               (160.0*eta2*etaij+15.0*eta2-6.0*etaij-1.0))*(etaij < 1.0);

}

}
