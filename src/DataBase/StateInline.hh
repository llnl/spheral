namespace Spheral {

//------------------------------------------------------------------------------
// Functors in a detail namespace to help with partial specialization
//------------------------------------------------------------------------------
namespace Detail {
    
template<typename Dimension, typename T>
struct EnrollAny {
  void operator()(State<Dimension>& state,
                  const typename State<Dimension>::KeyType& key,
                  T& thing) {
    dynamic_cast<StateBase<Dimension>*>(&state)->enroll(key, thing);
  }
};

template<typename Dimension, typename T>
struct EnrollAny<Dimension, std::shared_ptr<T>> {
  void operator()(State<Dimension>& state,
                  const typename State<Dimension>::KeyType& key,
                  std::shared_ptr<T>& thing) {
    auto UPP = std::dynamic_pointer_cast<UpdatePolicyBase<Dimension>>(thing);
    if (UPP) {
      state.enroll(key, UPP);
    } else {
      dynamic_cast<StateBase<Dimension>*>(&state)->enroll(key, thing);
    }
  }
};

}

//------------------------------------------------------------------------------
// Enroll the given field and policy.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(FieldBase<Dimension>& field,
       typename State<Dimension>::PolicyPointer policy) {
  this->enroll(field);
  this->enroll(this->key(field), policy);
}

//------------------------------------------------------------------------------
// Enroll the given FieldList and policy.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(FieldListBase<Dimension>& fieldList,
       typename State<Dimension>::PolicyPointer policy) {
  if (policy->clonePerField()) {
    // std::cerr << "Registering FieldList " << this->key(fieldList) << " with cloning policy" << std::endl;
    for (auto fptr: range(fieldList.begin_base(), fieldList.end_base())) {
      this->enroll(*fptr, policy);
    }
  } else {
    // std::cerr << "Registering FieldList " << this->key(fieldList) << " with SINGLE policy" << std::endl;
    // this->enroll(this->key(fieldList), fieldList);
    this->enroll(fieldList);  // enrolls each field without a policy
    this->enroll(this->key(fieldList), policy);
  }
}

//------------------------------------------------------------------------------
// Enroll the given policy.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(const typename State<Dimension>::KeyType& key,
       typename State<Dimension>::PolicyPointer policy) {
  KeyType fieldKey, nodeKey;
  this->splitFieldKey(key, fieldKey, nodeKey);
  mPolicyMap[fieldKey][key] = policy;
}

//------------------------------------------------------------------------------
// Return the policy for the specified field.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
typename State<Dimension>::PolicyPointer
State<Dimension>::
policy(const Field<Dimension, Value>& field) const {
  const KeyType key = StateBase<Dimension>::key(field);
  return this->policy(key);
}
  
//------------------------------------------------------------------------------
// Enroll an arbitrary type
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename T>
inline
void
State<Dimension>::
enroll(const KeyType& key, T& thing) {
  Detail::EnrollAny<Dimension, T>()(*this, key, thing);
}

}
