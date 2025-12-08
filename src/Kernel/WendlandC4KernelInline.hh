#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
WendlandC4Kernel< Dim<1> >::WendlandC4Kernel():
  Kernel<Dim<1>, WendlandC4Kernel< Dim<1> > >() {
  setVolumeNormalization(3.0/2.0);
  setKernelExtent(1.0);
  setInflectionPoint(0.275978);
}

template<>
inline
WendlandC4Kernel< Dim<2> >::WendlandC4Kernel():
  Kernel<Dim<2>, WendlandC4Kernel< Dim<2> > >() {
  setVolumeNormalization(9.0/(M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.235571);
}

template<>
inline
WendlandC4Kernel< Dim<3> >::WendlandC4Kernel():
  Kernel<Dim<3>, WendlandC4Kernel< Dim<3> > >() {
  setVolumeNormalization(495.0/(32.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.235571);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC4Kernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);

  double eta2 = etaij*etaij;
  if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(pow(1.0-etaij,5)*
                                             (1.0+5.0*etaij+8.0*eta2))*(etaij < 1.0);
  else
      return this->volumeNormalization()*Hdet*(pow(1.0-etaij,6)*
                                             (1.0+6.0*etaij+(35.0/3.0)*eta2))*(etaij < 1.0);
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC4Kernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(-14.0*pow(1.0-etaij,4)*etaij*
                                             (1.0+4.0*etaij))*(etaij < 1.0);
  else
      return  this->volumeNormalization()*Hdet*((56.0/3.0)*pow(etaij-1.0,5)*
                                              etaij*(5.0*etaij+1.0))*(etaij < 1.0);


}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC4Kernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
    const double eta2 = etaij*etaij;
    if(Dimension::nDim == 1)
      return this->volumeNormalization()*Hdet*(-14.0*pow(etaij-1.0,3)*
                                             (24.0*eta2-3.0*etaij-1.0))*(etaij < 1.0);
    else
      return this->volumeNormalization()*Hdet*((56.0/3.0)*pow(etaij-1.0,4)*
                                               (35.0*eta2-4.0*etaij-1.0))*(etaij < 1.0);



}

}
