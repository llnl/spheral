//---------------------------------Spheral++----------------------------------//
// IncrementASPHHtensor
//
// Specialized version of UpdatePolicyBase for time integrating the H tensor.
//
// Created by JMO, Mon Oct  7 13:31:02 PDT 2024
//----------------------------------------------------------------------------//
#include "IncrementASPHHtensor.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
IncrementASPHHtensor<Dimension>::
IncrementASPHHtensor(const bool fixShape,
                     const bool radialOnly,
                     std::shared_ptr<RadialFunctorType> radialFunctorPtr):
  UpdatePolicyBase<Dimension>(),
  mFixShape(fixShape),
  mRadialOnly(radialOnly),
  mRadialFunctorPtr(radialFunctorPtr) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IncrementASPHHtensor<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::H and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Get the state we're updating.
  auto       H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto DHDt = derivs.fields(prefix() + HydroFieldNames::H, SymTensor::zero());
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());  // Only needed if we're using radial scaling
  const auto numFields = H.numFields();
  CHECK(DHDt.numFields() == numFields);
  CHECK(pos.numFields() == numFields);

  // Walk the NodeLists
  for (auto k = 0u; k < numFields; ++k) {
    const auto& nodeList = H[k]->nodeList();
    const auto  hminInv = 1.0/nodeList.hmin();
    const auto  hmaxInv = 1.0/nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto n = nodeList.numInternalNodes();

    // Walk the nodes and update H (with limiting)
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
    
      // Check for special update rules
      if (mFixShape) {

        // Fix the shape (only volume scaling allowed)
        auto fi = Dimension::rootnu((H(k,i) + multiplier*DHDt(k,i)).Determinant()/H(k,i).Determinant());
        H(k,i) *= fi;

      } else if (mRadialOnly) {

        // Force only the radial component of H to be scaled
        const auto nhat = mRadialFunctorPtr->radialUnitVector(k, i, pos(k,i));
        const auto T = rotationMatrix(nhat);
        H(k,i).rotationalTransform(T);          // Should have one eigenvector aligned with the x' axis in this frame
        auto DHDti = DHDt(k,i);
        DHDti.rotationalTransform(T);
        H(k,i)[0] += multiplier * DHDti[0];
        H(k,i).rotationalTransform(T.Transpose());

      } else {

        H(k,i) += multiplier * DHDt(k,i);

      }

      // Apply limiting
      const auto hev = H(k,i).eigenVectors();
      const auto hminEffInv = min(hminInv, max(hmaxInv, hev.eigenValues.minElement())/hminratio);
      H(k,i) = constructSymTensorWithBoundedDiagonal(hev.eigenValues, hmaxInv, hminEffInv);
      H(k,i).rotationalTransform(hev.eigenVectors);
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
IncrementASPHHtensor<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  const auto rhsPtr = dynamic_cast<const IncrementASPHHtensor<Dimension>*>(&rhs);
  if (rhsPtr == nullptr) return false;

  // Ok, now do we agree on min & max?
  return (fixShape() == rhsPtr->fixShape() and
          radialOnly() == rhsPtr->radialOnly());
}

}
