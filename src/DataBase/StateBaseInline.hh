#include "boost/algorithm/string.hpp"
#include "DataBase/UpdatePolicyBase.hh"
#include "RK/RKCorrectionParams.hh"
#include "RK/ReproducingKernel.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/range.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Enroll an arbitrary type
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename T>
inline
void
StateBase<Dimension>::
enroll(const KeyType& key, T& thing) {
  // std::cerr << "StateBase::enroll " << key << std::endl;
  mStorage[key] = std::ref(thing);
}

template<typename Dimension>
template<typename T>
inline
void
StateBase<Dimension>::
enrollCopy(const KeyType& key, T& thing) {
  mCache.emplace_back(thing);
  mStorage[key] = std::ref(std::any_cast<T&>(mCache.back()));
}

//------------------------------------------------------------------------------
// Return the Field for the given key.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
Field<Dimension, Value>&
StateBase<Dimension>::
field(const KeyType& key) const {
  FieldBase<Dimension>& fb = this->template get<FieldBase<Dimension>>(key);
  auto* fptr = dynamic_cast<Field<Dimension, Value>*>(&fb);
  VERIFY2(fptr != nullptr,
          "StateBase::field ERROR: field type incorrect for key " << key);
  return *fptr;
}

template<typename Dimension>
template<typename Value>
inline
Field<Dimension, Value>&
StateBase<Dimension>::
field(const typename StateBase<Dimension>::KeyType& key, 
      const Value&) const {
  return this->template field<Value>(key);
}

//------------------------------------------------------------------------------
// Return all the Fields of the given Value element type.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
std::vector<Field<Dimension, Value>*>
StateBase<Dimension>::
allFields() const {
  std::vector<Field<Dimension, Value>*> result;
  KeyType fieldName, nodeListName;
  for (auto [key, aref]: mStorage) {
    if (aref.type() == typeid(std::reference_wrapper<FieldBase<Dimension>>)) {
      auto fb = std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(aref);
      auto* fptr = dynamic_cast<Field<Dimension, Value>*>(&fb.get());
      if (fptr != nullptr) result.push_back(fptr);
    }
  }
  return result;
}

template<typename Dimension>
template<typename Value>
inline
std::vector<Field<Dimension, Value>*>
StateBase<Dimension>::
allFields(const Value&) const {
  return this->template allFields<Value>();
}

//------------------------------------------------------------------------------
// Return a FieldList containing all registered fields of the given name.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
FieldList<Dimension, Value>
StateBase<Dimension>::
fields(const std::string& name,
       bool allowNone) const {
  FieldList<Dimension, Value> result;
  KeyType fieldName, nodeListName;
  for (auto [key, aref]: mStorage) {
    splitFieldKey(key, fieldName, nodeListName);
    if (fieldName == name) {
      CHECK(nodeListName != "");
      if (aref.type() == typeid(std::reference_wrapper<FieldBase<Dimension>>)) {
        auto fb = std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(aref);
        auto* fptr = dynamic_cast<Field<Dimension, Value>*>(&fb.get());
        if (fptr != nullptr) result.appendField(*fptr);
      }
    }
  }
  CHECK2(allowNone or result.size() > 0u, "Error: found no fields for key " << name);
  return result;
}

template<typename Dimension>
template<typename Value>
inline
FieldList<Dimension, Value>
StateBase<Dimension>::
fields(const std::string& name, const Value& dummy, bool allowNone) const {
  return this->template fields<Value>(name, allowNone);
}

//------------------------------------------------------------------------------
// Extract an arbitrary type
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
Value&
StateBase<Dimension>::
get(const typename StateBase<Dimension>::KeyType& key) const {
  Value* ptr = this->template getPtr<Value>(key);
  VERIFY2(ptr != nullptr, "StateBase ERROR: failed to return type for key " << key);
  return *ptr;
}

// Same thing passing a dummy argument to help with template type
template<typename Dimension>
template<typename Value>
inline
Value&
StateBase<Dimension>::
get(const typename StateBase<Dimension>::KeyType& key,
    const Value&) const {
  return this->get<Value>(key);
}

//------------------------------------------------------------------------------
// Extract an arbitrary type as a pointer
// Does not throw, but rather returns nullptr if failing
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
Value*
StateBase<Dimension>::
getPtr(const typename StateBase<Dimension>::KeyType& key) const {
  auto itr = mStorage.find(key);
  if (itr == mStorage.end() or
      itr->second.type() != typeid(std::reference_wrapper<Value>)) return nullptr;
  return &(std::any_cast<std::reference_wrapper<Value>>(itr->second).get());
}

// Same thing passing a dummy argument to help with template type
template<typename Dimension>
template<typename Value>
inline
Value*
StateBase<Dimension>::
getPtr(const typename StateBase<Dimension>::KeyType& key,
       const Value&) const {
  return this->getPtr<Value>(key);
}

//------------------------------------------------------------------------------
// Assign the Fields matching the given name of this State object to be equal to
// the values in another.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
void
StateBase<Dimension>::
assignFields(const StateBase<Dimension>& rhs, const std::string name) {
  auto lhsfields = this->fields(name, Value());
  auto rhsfields = rhs.fields(name, Value());
  lhsfields.assignFields(rhsfields);
}

//------------------------------------------------------------------------------
// Construct the lookup key for the given field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename StateBase<Dimension>::KeyType
StateBase<Dimension>::
key(const FieldBase<Dimension>& field) {
  return buildFieldKey(field.name(), field.nodeListPtr()->name());
}

//------------------------------------------------------------------------------
// Construct the lookup key for the given FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename StateBase<Dimension>::KeyType
StateBase<Dimension>::
key(const FieldListBase<Dimension>& fieldList) {
  REQUIRE(fieldList.begin_base() != fieldList.end_base());
  return buildFieldKey((*fieldList.begin_base())->name(), UpdatePolicyBase<Dimension>::wildcard());
}

//------------------------------------------------------------------------------
// Internal methods to encode the convention for combining Field and NodeList
// names into a single unique key.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename StateBase<Dimension>::KeyType
StateBase<Dimension>::
buildFieldKey(const std::string& fieldName,
              const std::string& nodeListName) {
  return fieldName + "|" + nodeListName;
}

template<typename Dimension>
inline
void
StateBase<Dimension>::
splitFieldKey(const KeyType& key,
              KeyType& fieldName,
              KeyType& nodeListName) {
  std::vector<std::string> keys;
  boost::split(keys, key, boost::is_any_of("|"));
  if (keys.size() > 1) {
    fieldName = keys[0];
    nodeListName = keys[1];
  } else if (keys.size() == 1) {
    fieldName = keys[0];
    nodeListName = "";
  } else {
    fieldName = "";
    nodeListName = "";
  }
}

}
