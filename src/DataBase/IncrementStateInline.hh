//---------------------------------Spheral++----------------------------------//
// IncrementState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//

#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(std::initializer_list<std::string> depends,
               const bool wildCardDerivs):
  FieldUpdatePolicy<Dimension, Value>(depends),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const bool wildCardDerivs):
  FieldUpdatePolicy<Dimension, Value>({}),
  mWildCardDerivs(wildCardDerivs) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
void
IncrementState<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  // Get the state we're updating.
  auto& f = state.field(key, Value());

  // Find all the available matching derivative Field keys.
  const auto incrementKey = prefix() + fieldKey;
  const auto allkeys = derivs.keys();
  KeyType dfKey, dfNodeListKey;
  auto numDeltaFields = 0u;
  CONTRACT_VAR(numDeltaFields);
  for (const auto& key: allkeys) {
    StateBase<Dimension>::splitFieldKey(key, dfKey, dfNodeListKey);
    if (dfNodeListKey == nodeListKey and
        dfKey.compare(0, incrementKey.size(), incrementKey) == 0) {
      ++numDeltaFields;

      // This delta field matches the base of increment key, so apply it.
      const auto& df = derivs.field(key, Value());
      const auto  n = f.numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        f(i) += multiplier*(df(i));
      }
    }
  }

  // If we're not allowing wildcard update, there should have only be one match.
  VERIFY2(mWildCardDerivs or numDeltaFields == 1,
          "IncrementState ERROR: unable to find unique match for derivative field key " << incrementKey);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
bool
IncrementState<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const IncrementState<Dimension, Value>*>(&rhs);
  return rhsPtr != nullptr;
}

//------------------------------------------------------------------------------
// Serialize the derivative data.  The most common usage here will be for
// descendant classes to override this and call this base method with a modified
// key.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
void
IncrementState<Dimension, Value>::
serializeDerivatives(std::vector<double>& buf,
                     const KeyType& key,
                     const StateDerivatives<Dimension>& derivs) const {

  // Are we allowing multiple derivative contributions?  If not this is simpler.
  if (not mWildCardDerivs) {
    const auto incrementKey = prefix() + key;
    const auto& f = derivs.template field<Value>(incrementKey);
    CHECK(f.fixedSizeDataType());
    const std::vector<char> rawbuf = packFieldValues(f);
    const auto ndvals = f.numInternalElements() * f.numValsInDataType();
    const auto nraw = rawbuf.size();
    CHECK(nraw == ndvals*sizeof(double));
    const auto istart = buf.size();
    buf.resize(istart + ndvals);
    CHECK((buf.size() - istart)*sizeof(double) == nraw);
    std::memcpy(&buf[istart], &rawbuf[0], nraw);

  } else {
    // We potentially have multiple derivative fields accumulating,
    // so we need a temporary buffer to accumulate the result.
    KeyType fKey, nodeListKey, dfKey, dfNodeListKey;
    StateBase<Dimension>::splitFieldKey(key, fKey, nodeListKey);
    const auto incrementKey = prefix() + fKey;

    std::vector<Value> vals;
    const auto allkeys = derivs.keys();
    for (const auto& dkey: allkeys) {
      StateBase<Dimension>::splitFieldKey(dkey, dfKey, dfNodeListKey);
      if (dfNodeListKey == nodeListKey and
          dfKey.compare(0, incrementKey.size(), incrementKey) == 0) {
        // This delta field matches the base of increment key, so apply it.
        const auto& df = derivs.template field<Value>(dkey);
        const auto n = df.numInternalElements();
        vals.resize(n, DataTypeTraits<Value>::zero());
        for (auto i = 0u; i < n; ++i) vals[i] += df[i];
      }
    }

    const auto n = vals.size();
    if (n > 0u) {
      const auto ndvals = n * DataTypeTraits<Value>::numElements(vals[0]);
      const auto nraw = ndvals * sizeof(double);
      std::vector<char> rawbuf;
      for (auto i = 0u; i < n; ++i) packElement(vals[i], rawbuf);
      CHECK(rawbuf.size() == nraw);
      const auto istart = buf.size();
      buf.resize(istart + ndvals);
      CHECK((buf.size() - istart)*sizeof(double) == nraw);
      std::memcpy(&buf[istart], &rawbuf[0], nraw);
    }
  }      
}

}

