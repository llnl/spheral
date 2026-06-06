//---------------------------------Spheral++----------------------------------//
// State -- Accumulate and cart the state for a set of physics packages around.
//
// Created by JMO, Fri Aug 27 10:56:40 2004
//----------------------------------------------------------------------------//
#include "State.hh"
#include "StateBase.hh"
#include "StateDerivatives.hh"
#include "DataBase.hh"
#include "UpdatePolicyBase.hh"
#include "Physics/Physics.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"

#include <string>
using std::vector;
using std::map;
using std::set;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Helper to compare a set of state to a key.
// We assume that wildcards are only applied at the end of a Key!
//------------------------------------------------------------------------------
// template<typename Key, typename Iterator>
// bool keyPresent(const Key& key,
//                 const Key& wildcard,
//                 const Iterator setBegin,
//                 const Iterator setEnd) {
//   const size_t pos = key.find(wildcard);
//   if (pos == string::npos) {
//     Key fieldKey, nodeListKey;
//     Iterator itr = setBegin;
//     while (itr != setEnd) {
//       StateBase<Dim<1> >::splitFieldKey(*itr, fieldKey, nodeListKey);
//       if (fieldKey == key) return true;
//       ++itr;
//     }
//     return false;
//     //return (find(setBegin, setEnd, key) != setEnd);
//   } else {
//     const Key testKey = key.substr(0, pos);
//     Iterator itr = setBegin;
//     while (itr != setEnd) {
//       if (testKey == itr->substr(0, pos)) return true;
//       ++itr;
//     }
//     return false;
//   }
// }

// //------------------------------------------------------------------------------
// // Helper to compare a set of state to a key.
// // We assume that wildcards are only applied at the end of a Key!
// //------------------------------------------------------------------------------
// bool compareKeysByField(const std::string& a,
//                         const std::string& b) {
//   std::string af, bf, an, bn;
//   StateBase<Dim<1> >::splitFieldKey(a, af, an);
//   StateBase<Dim<1> >::splitFieldKey(b, bf, bn);
//   return af < bf;
// }

}

//------------------------------------------------------------------------------
// Construct with the state for the given set of Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>::
State(DataBase<Dimension>& dataBase,
      typename State<Dimension>::PackageList& physicsPackages):
  StateBase<Dimension>(),
  mPolicyMap(),
  mTimeAdvanceOnly(false) {
  // Iterate over the physics packages, and have them register their state.
  for (auto pkg: physicsPackages) pkg->registerState(dataBase, *this);
  auto cmp = dataBase.connectivityMapPtr();
  if (cmp) this->enrollConnectivityMap(cmp);
}

//------------------------------------------------------------------------------
// Construct with the state for the given set of Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>::
State(DataBase<Dimension>& dataBase,
      typename State<Dimension>::PackageIterator physicsPackageBegin,
      typename State<Dimension>::PackageIterator physicsPackageEnd):
  StateBase<Dimension>(),
  mPolicyMap(),
  mTimeAdvanceOnly(false) {
  // Iterate over the physics packages, and have them register their state.
  for (auto pkg: range(physicsPackageBegin, physicsPackageEnd)) pkg->registerState(dataBase, *this);
  auto cmp = dataBase.connectivityMapPtr();
  if (cmp) this->enrollConnectivityMap(cmp);
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension>
bool
State<Dimension>::
operator==(const StateBase<Dimension>& rhs) const {
  return StateBase<Dimension>::operator==(rhs);
}

//------------------------------------------------------------------------------
// The set of keys for all registered policies.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<typename State<Dimension>::KeyType>
State<Dimension>::
policyKeys() const {
  vector<KeyType> result;
  for (const auto& itr: mPolicyMap) result.push_back(itr.first);
  ENSURE(result.size() == mPolicyMap.size());
  return result;
}

//------------------------------------------------------------------------------
// Return the policy for the given key.
//------------------------------------------------------------------------------
template<typename Dimension>
typename State<Dimension>::PolicyPointer
State<Dimension>::
policy(const typename State<Dimension>::KeyType& key) const {
  KeyType fieldKey, nodeKey;
  this->splitFieldKey(key, fieldKey, nodeKey);
  const auto outerItr = mPolicyMap.find(fieldKey);
  if (outerItr == mPolicyMap.end()) return PolicyPointer();
  // VERIFY2(outerItr != mPolicyMap.end(),
  //         "State ERROR: attempted to retrieve non-existent policy for key " << key);
  const auto& key2policies = outerItr->second;
  const auto innerItr = key2policies.find(key);
  if (innerItr == key2policies.end()) return PolicyPointer();
  // VERIFY2(innerItr != policies.end(),
  //         "State ERROR: attempted to retrieve non-existent policy for key " << key);
  return innerItr->second;
}

//------------------------------------------------------------------------------
// Return all the policies for the given field key.
//------------------------------------------------------------------------------
template<typename Dimension>
std::map<typename State<Dimension>::KeyType, typename State<Dimension>::PolicyPointer>
State<Dimension>::
policies(const typename State<Dimension>::KeyType& fieldKey) const {
  const auto outerItr = mPolicyMap.find(fieldKey);
  if (outerItr == mPolicyMap.end()) return std::map<KeyType, PolicyPointer>();
  // VERIFY2(outerItr != mPolicyMap.end(),
  //         "State ERROR: attempted to retrieve non-existent policy for key " << key);
  return outerItr->second;
}

//------------------------------------------------------------------------------
// Remove the policy associated with the given key.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
removePolicy(const typename State<Dimension>::KeyType& key) {
  KeyType fieldKey, nodeKey;
  this->splitFieldKey(key, fieldKey, nodeKey);
  typename PolicyMapType::iterator outerItr = mPolicyMap.find(fieldKey);
  VERIFY2(outerItr != mPolicyMap.end(),
          "State ERROR: attempted to remove non-existent policy for field key " << fieldKey);
  std::map<KeyType, PolicyPointer>& policies = outerItr->second;
  typename std::map<KeyType, PolicyPointer>::iterator innerItr = policies.find(key);
  if (innerItr == policies.end()) {
    cerr << "State ERROR: attempted to remove non-existent policy for inner key " << key << endl
         << "Known keys are: " << endl;
    for (auto itr = policies.begin(); itr != policies.end(); ++itr) cerr << " --> " << itr->first << endl;
    VERIFY(innerItr != policies.end());
  }
  policies.erase(innerItr);
  if (policies.size() == 0) mPolicyMap.erase(outerItr);
}

//------------------------------------------------------------------------------
// Remove the policy associated with a Field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
removePolicy(FieldBase<Dimension>& field) {
  this->removePolicy(StateBase<Dimension>::key(field));
}

//------------------------------------------------------------------------------
// Remove the policy associated with a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
removePolicy(FieldListBase<Dimension>& fieldList,
             const bool clonePerField) {
  if (clonePerField) {
    for (auto fieldPtrItr = fieldList.begin_base();
         fieldPtrItr < fieldList.end_base();
         ++fieldPtrItr) this->removePolicy(**fieldPtrItr);
  } else {
    this->removePolicy(StateBase<Dimension>::key(fieldList));
  }
}

//------------------------------------------------------------------------------
// Update the state with the given derivatives object, according to the per
// state field policies.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
update(StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt,
       const bool dependentOnly) {

  // cerr << "################################################################################" << endl;

  // Prepare lists of the keys to be completed.
  vector<KeyType> fieldsToBeCompleted;
  map<KeyType, set<KeyType>> stateToBeCompleted;
  for (auto& fieldkey_policies: mPolicyMap) {
    auto& fieldkey = fieldkey_policies.first;               // Just Field names
    for (auto& [fullkey, policyPtr]: fieldkey_policies.second) {
      // if ((not dependentOnly) or policyPtr->dependent()) {  // Fully encoded Field + NodeList names
        stateToBeCompleted[fieldkey].insert(fullkey);
      // }
    }
  }
  for (const auto& key_fullkey: stateToBeCompleted) fieldsToBeCompleted.push_back(key_fullkey.first);
  CHECK(fieldsToBeCompleted.size() == stateToBeCompleted.size());

  // Keep a copy of all the Field keys we are evolving
  const std::set<KeyType> allFieldKeys(fieldsToBeCompleted.begin(), fieldsToBeCompleted.end());

  // Iterate until all state has been updated.
  while (not stateToBeCompleted.empty()) {

    // Walk the remaining state to be completed.
    vector<KeyType> stateToRemove;
    for (auto& [fieldKey, remainingKeys]: stateToBeCompleted) { // (fieldKey, set<full_fieldKeys>)

      // Walk the remaining individual keys for this fieldKey.
      for (auto& key: remainingKeys) {
        const PolicyPointer policyPtr = mPolicyMap[fieldKey][key];

        // Check if any of the dependencies for this state have not been completed yet.
        bool allDependenciesMet = true;
        {
          const auto& dependencies = policyPtr->dependencies();
          auto itr = dependencies.begin();
          while (itr < dependencies.end() and
                 find(fieldsToBeCompleted.begin(), fieldsToBeCompleted.end(), *itr) == fieldsToBeCompleted.end()) itr++;
          allDependenciesMet = (itr == dependencies.end());
        }
        if (allDependenciesMet) {

          // We also require that any FieldList policies fire before NodeList specific 
          // versions of the same Field names.
          KeyType fieldKey, nodeListKey;
          this->splitFieldKey(key, fieldKey, nodeListKey);
          bool fire = (nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
          if (not fire) {
            const KeyType wildKey = this->buildFieldKey(fieldKey, UpdatePolicyBase<Dimension>::wildcard());
            fire = (find(remainingKeys.begin(), remainingKeys.end(), wildKey) == remainingKeys.end());
          }

          if (fire) {
            // cerr <<" --> Update " << key << endl;
            if (dependentOnly and policyPtr->independent()) {
              if (mTimeAdvanceOnly) {
                policyPtr->updateAsIncrement(key, *this, derivs, 0.0, t, 0.0);
              } else {
                policyPtr->update(key, *this, derivs, 0.0, t, 0.0);
              }
            } else {
              if (mTimeAdvanceOnly) {
                policyPtr->updateAsIncrement(key, *this, derivs, multiplier, t, dt);
              } else {
                policyPtr->update(key, *this, derivs, multiplier, t, dt);
              }
            }
            // List this field as completed.
            stateToRemove.push_back(key);
          }
        }
      }
    }

    // Check that the some state has been updated on this iteration.  If not, then *somebody*
    // has specified a circular dependency tree!
    if (stateToRemove.empty()) {
      std::stringstream message;
      message << "State::update ERROR: someone has specified a circular state dependency.\n"
              << "Remaining State:\n";
      for (const auto& itr: stateToBeCompleted) message << "   " << itr.first << "\n";
      message << "State dependencies:\n";
      for (const auto& itr: mPolicyMap) {
        const KeyType fieldKey = itr.first;
        const map<KeyType, PolicyPointer>& keysAndPolicies = itr.second;
        for (const auto& pitr: keysAndPolicies) {
          const KeyType key = pitr.first;
          const PolicyPointer policyPtr = pitr.second;
          message << key << " : ";
          for (const auto& dep: policyPtr->dependencies()) message << dep << "  ++  ";
          message << "\n";
        }
      }
      VERIFY2(not stateToRemove.empty(), message.str());
    }
    
    // Remove the completed state.
    for (const auto& key: stateToRemove) {
      KeyType fieldKey, nodeListKey;
      this->splitFieldKey(key, fieldKey, nodeListKey);
      CHECK(stateToBeCompleted.find(fieldKey) != stateToBeCompleted.end());
      stateToBeCompleted[fieldKey].erase(key);
      if (stateToBeCompleted[fieldKey].empty()) stateToBeCompleted.erase(fieldKey);
    }
    fieldsToBeCompleted = vector<KeyType>();
    for (const auto& itr: stateToBeCompleted) fieldsToBeCompleted.push_back(itr.first);
    CHECK(fieldsToBeCompleted.size() == stateToBeCompleted.size());
  }
}

//------------------------------------------------------------------------------
// Serialize independent data to a buffer for use in implicit time integration
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
serializeIndependentData(std::vector<double>& buf) const {
  buf.clear();

  // Look for any state with update policies that indicate an independent variable
  for (const auto& stuff: mPolicyMap) {
    const auto& keys2policies = stuff.second;
    for (const auto& key2policy: keys2policies) {
      const auto& key = key2policy.first;
      const auto& policyPtr = key2policy.second;
      if (policyPtr->independent()) {
        // cerr << "Serializing " << key << endl;
        policyPtr->serializeData(buf, key, *this);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Deserialize independent data from a buffer; for use in implicit time
// integration
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
deserializeIndependentData(const std::vector<double>& buf) const {
  
  // Look for any state with update policies that indicate an independent variable
  size_t offset = 0u;
  for (const auto& stuff: mPolicyMap) {
    const auto& keys2policies = stuff.second;
    for (const auto& key2policy: keys2policies) {
      const auto& key = key2policy.first;
      const auto& policyPtr = key2policy.second;
      if (policyPtr->independent()) {
        // cerr << "Deserializing " << key << endl;
        offset = policyPtr->deserializeData(buf, key, *this, offset);
        CHECK(offset <= buf.size());
      }
    }
  }
  CHECK(offset == buf.size());
}

//------------------------------------------------------------------------------
// Serialize the derivatives associated with independent data to a buffer for
// use in implicit time integration
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
serializeDerivatives(std::vector<double>& buf,
                     const StateDerivatives<Dimension>& derivs) const {
  buf.clear();

  // Look for any state with update policies that indicate an independent variable
  for (const auto& stuff: mPolicyMap) {
    const auto& keys2policies = stuff.second;
    for (const auto& key2policy: keys2policies) {
      const auto& key = key2policy.first;
      const auto& policyPtr = key2policy.second;
      if (policyPtr->independent()) {
        // cerr << "Serializing DERIV " << key << endl;
        policyPtr->serializeDerivatives(buf, key, derivs);
      }
    }
  }
}

}

