//---------------------------------Spheral++----------------------------------//
// IncrementPairFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#include "DEM/IncrementPairFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"

#include <vector>
using std::vector;
using std::string;


namespace Spheral {

//------------------------------------------------------------------------------
// Constructors mostly wrappers to the base
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
IncrementPairFieldList<Dimension, Value>::
IncrementPairFieldList(std::initializer_list<std::string> depends):
  UpdatePolicyBase<Dimension>(depends) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
IncrementPairFieldList<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Get the state we're updating.
  auto f = state.fields(fieldKey, Value());
  const auto numNodeLists = f.size();

  // Find all the available matching derivative FieldList keys.
  const auto incrementKey = prefix() + fieldKey;
  // cerr << "IncrementPairFieldList: [" << fieldKey << "] [" << incrementKey << "] : " << endl;
  const auto allkeys = derivs.fullFieldKeys();
  vector<string> incrementKeys;
  for (const auto& key: allkeys) {
    // if (std::regex_search(key, std::regex("^" + incrementKey))) {
    if (key.compare(0, incrementKey.size(), incrementKey) == 0) {
      incrementKeys.push_back(key);
    }
  }
  CHECK(not incrementKeys.empty());
  // Update by each of our derivative fields.
  for (const auto& key: incrementKeys) {
    const auto df = derivs.fields(key, Value());
    CHECK(df.size() == f.size());
    for (auto k = 0u; k != numNodeLists; ++k) {
      const auto numNodes = f[k]->numInternalElements();
      for (auto i = 0u; i != numNodes; ++i) {
        const auto numContacts = df(k,i).size();
        if (f(k,i).size()!=numContacts)f(k,i).resize(numContacts);
        for (auto j = 0u; j != numContacts; ++j) {
          f(k, i)[j] += multiplier*(df(k, i)[j]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
IncrementPairFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const IncrementPairFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != nullptr;
}


}

