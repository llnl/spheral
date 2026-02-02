//---------------------------------Spheral++----------------------------------//
// monotonicallyLimitedGradient
//
// The methods for applying monotonic limiting to the gradient operations.
//
// Created by JMO, Wed Feb 13 09:32:33 PST 2008
//----------------------------------------------------------------------------//

#include "Geometry/MathTraits.hh"
#include "Geometry/innerProduct.hh"
#include "Utilities/SpheralFunctions.hh"

#include <vector>
#include <map>

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the necessary limit such that the projected value falls within
// the specified delta.
// Note that for vectors and tensors, we just use component wise limiting.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
computeMonotonicLimiting(const typename Dimension::Vector& nhat,
                         const typename Dimension::Scalar& dFproj,
                         const typename Dimension::Scalar& dF,
                         const double fuzz = 1.0e-15) {
  REQUIRE(fuzz > 0.0);
  REQUIRE(fuzzyEqual(nhat.magnitude2(), 1.0, 1.0e-10));
  return std::min(1.0, std::max(0.0, std::max(dF/(std::abs(dFproj) + fuzz)*sgn(dFproj), fuzz/(dFproj*dFproj + fuzz))));
}

template<typename Dimension>
inline
double
computeMonotonicLimiting(const typename Dimension::Vector& nhat,
                         const typename Dimension::Vector& dFproj,
                         const typename Dimension::Vector& dF) {
  return computeMonotonicLimiting<Dimension>(nhat, dFproj.dot(nhat), dF.dot(nhat));
}

template<typename Dimension>
inline
double
computeMonotonicLimiting(const typename Dimension::Vector& nhat,
                         const typename Dimension::Tensor& dFproj,
                         const typename Dimension::Tensor& dF) {
  return computeMonotonicLimiting<Dimension>(nhat, dFproj.dot(nhat), dF.dot(nhat));
}


//------------------------------------------------------------------------------
// Build the tensor limiter based on point-wise projected limits.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::SymTensor
buildLimiter(const std::vector<std::pair<double, typename Dimension::Vector> >& phiScales) {

  using std::vector;
  using std::pair;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename vector<pair<double, Vector> >::const_iterator PhiItr;

  const double fuzz = 1.0e-20;

  // Pre-conditions.
  REQUIRE(phiScales.size() > 0);
  BEGIN_CONTRACT_SCOPE
  {
    for (size_t k = 0; k != phiScales.size(); ++k) {
      REQUIRE(fuzzyEqual(phiScales[k].second.magnitude2(), 1.0, 1.0e-12));
    }
  }
  END_CONTRACT_SCOPE

  // Iterate over the input scalar phis and build the second moment tensor.
  SymTensor result;
  for (PhiItr itr = phiScales.begin(); itr != phiScales.end(); ++itr) {
    const double phi0 = itr->first;
    const Vector& vec0 = itr->second;
    result += 1.0/(phi0 + 1.0e-8) * vec0.selfdyad();
  }
  result /= phiScales.size();
  const double thpt = std::max(fuzz, result.maxAbsElement());
  CHECK(thpt > 0.0);
  result /= thpt;
  result = (fuzzyEqual(result.Determinant(), 0.0) ? 
            SymTensor::zero() :
            thpt * (result.Inverse()));

  // Scale the thing to make sure all our limits are met.
  double f = 1.0;
  for (PhiItr itr = phiScales.begin(); itr != phiScales.end(); ++itr) {
    const double phi0 = itr->first;
    const Vector& vec0 = itr->second;
    const double phi1 = (result.dot(vec0)).magnitude();
    f = std::min(f, std::min(1.0, phi0/(phi1 + fuzz)));
  }
  result *= f;

  // Post-condtions.
  BEGIN_CONTRACT_SCOPE
  {
    for (PhiItr itr = phiScales.begin(); itr != phiScales.end(); ++itr) {
      const double phi0 = itr->first;
      const Vector& vec0 = itr->second;
      ENSURE(fuzzyLessThanOrEqual(result.dot(vec0).magnitude(), phi0, 1.0e-12));
    }
  }
  END_CONTRACT_SCOPE

  return result;
}

//------------------------------------------------------------------------------
// Return the monotonically limited gradient for the given node, based on a 
// pre-computed gradient and connectivity.
// Scalar version.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
typename MathTraits<Dimension, DataType>::GradientType
scalarLimitedGradient(const DataType& fieldi,
                      const typename MathTraits<Dimension, DataType>::GradientType& gradi,
                      const typename Dimension::Vector& ri,
                      const std::vector<typename Dimension::Vector>& rNeighbors,
                      const std::vector<DataType>& fieldNeighbors) {

  using std::vector;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  const unsigned numNeighbors = rNeighbors.size();
  REQUIRE2(fieldNeighbors.size() == numNeighbors, "Bad sizing:  " << fieldNeighbors.size() << " " << numNeighbors);

  // Build up the set of per neighbors constraints (phi).
  double phimin = 1.0;
  unsigned k;
  Vector rji, drHat;
  DataType dF, dFproj;
  for (k = 0; k != numNeighbors; ++k) {
    const Vector& rj = rNeighbors[k];
    const DataType& fieldj = fieldNeighbors[k];
    rji = rj - ri;
    drHat = rji.unitVector();
    dF = fieldj - fieldi;
    dFproj = Geometry::innerProduct<Dimension>(gradi, rji);
    phimin = std::min(phimin, computeMonotonicLimiting<Dimension>(drHat, dFproj, dF));
  }

  // Return the limited gradient, and we're done.
  return phimin*gradi;
}

//------------------------------------------------------------------------------
// Return the monotonically limited gradient for the given node, based on a 
// pre-computed gradient and connectivity.
// Tensor version.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
typename MathTraits<Dimension, DataType>::GradientType
tensorLimitedGradient(const DataType& fieldi,
                      const typename MathTraits<Dimension, DataType>::GradientType& gradi,
                      const typename Dimension::Vector& ri,
                      const std::vector<typename Dimension::Vector>& rNeighbors,
                      const std::vector<DataType>& fieldNeighbors) {

  using std::vector;
  using std::pair;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  const unsigned numNeighbors = rNeighbors.size();
  REQUIRE2(fieldNeighbors.size() == numNeighbors, "Bad sizing:  " << fieldNeighbors.size() << " " << numNeighbors);

  // Build up the set of per neighbors constraints (phi).
  vector<pair<double, Vector> > phiScales;
  unsigned k;
  Vector rji, drHat;
  DataType dF, dFproj;
  for (k = 0; k != numNeighbors; ++k) {
    const Vector& rj = rNeighbors[k];
    const DataType& fieldj = fieldNeighbors[k];
    rji = rj - ri;
    drHat = rji.unitVector();
    dF = fieldj - fieldi;
    dFproj = Geometry::innerProduct<Dimension>(gradi, rji);
    phiScales.push_back(pair<double, Vector>(computeMonotonicLimiting<Dimension>(drHat, dFproj, dF),
                                             drHat));
  }

  // Build the tensor limiter.
  const SymTensor phi = buildLimiter<Dimension>(phiScales);

  // Return the limited gradient, and we're done.
  return phi*gradi;
}

}
