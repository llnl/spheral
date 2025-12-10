#include "Geometry/Dimension.hh"
#include "VolumeIntegrationFunctions.hh"
#include "Utilities/DBC.hh"

#include <math.h>
#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given order (allowed values = 1, 3, 5, 7, & 9).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NSincPolynomialKernel<Dimension>::NSincPolynomialKernel(const int order):
  Kernel<Dimension, NSincPolynomialKernel<Dimension> >(),
  mOrder(order),
  mAij((order + 1)/2) {

  // We only have equations for order = 1, 3, 5, 7, or 9.
  REQUIRE(order == 1 || order == 3 || order == 5 || order == 7 || order == 9);
  if (order != 1 && order != 3 && order != 5 && order != 7 && order !=9) {
    std::cerr << "NSincPolynomialKernel ERROR: only support order = [1, 3, 5, 7, 9], "
              << order << " specified." << std::endl;
  }

  this->setKernelExtent((order + 1)/2);
  this->setInflectionPoint(0.0);
  
  // Set the polynomial coefficients.
  setPolynomialCoefficients(order, mAij);

  // Numerically evaluate the volume normalization.
  this->setVolumeNormalization(1.0);
  double volumeIntegral = simpsonsVolumeIntegral<Dimension, NSincPolynomialKernel<Dimension> >(*this, 0.0, this->kernelExtent(), 10000);
  CHECK(volumeIntegral > 0.0);
  this->setVolumeNormalization(1.0/volumeIntegral);
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NSincPolynomialKernel<Dimension>::kernelValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  const int numPolynomials = (mOrder + 1)/2;
  CHECK(numPolynomials > 0 && numPolynomials < (int)mAij.size());

  CHECK(this->kernelExtent() > 0.0);
  const int iPolynomial = int(etaij);

  if (iPolynomial >= numPolynomials) {
    return 0.0;
  } else {
    double result = 0.0;
    for (int i = 0; i <= mOrder; ++i)
      result += mAij[iPolynomial][i]*pow(etaij, i);
    result *= this->volumeNormalization()*Hdet;
    return result;
  }
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NSincPolynomialKernel<Dimension>::gradValue(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  const int numPolynomials = (mOrder + 1)/2;
  CHECK(numPolynomials > 0 && numPolynomials < (int)mAij.size());

  CHECK(this->kernelExtent() > 0.0);
  const int iPolynomial = int(etaij);

  if (iPolynomial >= numPolynomials) {
    return 0.0;
  } else {
    double result = 0.0;
    for (int i = 1; i <= mOrder; ++i) 
      result += mAij[iPolynomial][i]*i*pow(etaij, i - 1);
    result *= this->volumeNormalization()*Hdet;
    return result;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NSincPolynomialKernel<Dimension>::grad2Value(double etaij, const double Hdet) const {
  CHECK(etaij >= 0.0);
  CHECK(Hdet >= 0.0);

  const int numPolynomials = (mOrder + 1)/2;
  CHECK(numPolynomials > 0 && numPolynomials < (int)mAij.size());

  CHECK(this->kernelExtent() > 0.0);
  const int iPolynomial = int(etaij);

  if (iPolynomial >= numPolynomials) {
    return 0.0;
  } else {
    double result = 0.0;
    for (int i = 2; i <= mOrder; ++i) 
      result += mAij[iPolynomial][i]*i*(i - 1)*pow(etaij, i - 2);
    result *= this->volumeNormalization()*Hdet;
    return result;
  }
}

}
