#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the kernel and gradient values for a set of normalized distances.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TableKernel<Dimension>::kernelAndGradValues(const std::vector<Scalar>& etaijs,
                                            const std::vector<Scalar>& Hdets,
                                            std::vector<Scalar>& kernelValues,
                                            std::vector<Scalar>& gradValues) const {
  // Preconditions.
  const auto n = etaijs.size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(Hdets.size() == n);
    for (auto i = 0u; i < n; ++i) {
      REQUIRE(etaijs[i] >= 0.0);
      REQUIRE(Hdets[i] >= 0.0);
    }
  }
  END_CONTRACT_SCOPE

  // Prepare the results.
  kernelValues.resize(n);
  gradValues.resize(n);

  // Fill those suckers in.
  for (auto i = 0u; i < n; ++i) {
    const auto i0 = this->mInterp.lowerBound(etaijs[i]);
    kernelValues[i] = Hdets[i]*this->mInterp(etaijs[i], i0);
    gradValues[i] = Hdets[i]*this->mGradInterp(etaijs[i], i0);
  }
}

}

