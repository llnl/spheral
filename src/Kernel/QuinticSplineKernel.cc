#include "QuinticSplineKernel.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
QuinticSplineKernel<Dimension>::~QuinticSplineKernel() {
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
double
QuinticSplineKernel<Dimension>::kernelValue(double eta, const double Hdet) const {
  REQUIRE(eta >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (eta < 1.0/3.0) {
    return this->volumeNormalization()*Hdet*(    FastMath::pow5(1.0 - eta) - 
                                              6.0*FastMath::pow5(2.0/3.0 - eta) +
                                             15.0*FastMath::pow5(1.0/3.0 - eta));
  } else if (eta < 2.0/3.0) {
    return this->volumeNormalization()*Hdet*(    FastMath::pow5(1.0 - eta) - 
                                              6.0*FastMath::pow5(2.0/3.0 - eta));
  } else if (eta < 1.0) {
    return this->volumeNormalization()*Hdet*(    FastMath::pow5(1.0 - eta));
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
double
QuinticSplineKernel<Dimension>::gradValue(double eta, const double Hdet) const {
  REQUIRE(eta >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (eta < 1.0/3.0) {
    return this->volumeNormalization()*Hdet*(-5.0*FastMath::pow4(1.0 - eta) +
                                             30.0*FastMath::pow4(2.0/3.0 - eta) - 
                                             75.0*FastMath::pow4(1.0/3.0 - eta));
  } else if (eta < 2.0/3.0) {
    return this->volumeNormalization()*Hdet*(-5.0*FastMath::pow4(1.0 - eta) +
                                             30.0*FastMath::pow4(2.0/3.0 - eta));
  } else if (eta < 1.0) {
    return this->volumeNormalization()*Hdet*(-5.0*FastMath::pow4(1.0 - eta));
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
double
QuinticSplineKernel<Dimension>::grad2Value(double eta, const double Hdet) const {
  REQUIRE(eta >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (eta < 1.0/3.0) {
    return this->volumeNormalization()*Hdet*( 20.0*FastMath::pow3(1.0 - eta) -
                                             120.0*FastMath::pow3(2.0/3.0 - eta) + 
                                             300.0*FastMath::pow3(1.0/3.0 - eta));
  } else if (eta < 2.0/3.0) {
    return this->volumeNormalization()*Hdet*( 20.0*FastMath::pow3(1.0 - eta) -
                                             120.0*FastMath::pow3(2.0/3.0 - eta));
  } else if (eta < 1.0) {
    return this->volumeNormalization()*Hdet*( 20.0*FastMath::pow3(1.0 - eta));
  } else {
    return 0.0;
  }
}

}
