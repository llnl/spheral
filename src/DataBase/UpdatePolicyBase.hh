//---------------------------------Spheral++----------------------------------//
// UpdatePolicyBase -- Base/interface class for the policies on how state 
// variables are to be updated.
//
// Note that we currenly only support specifying [0,6] dependencies.  There's 
// no problem adding more constructors to expand the allowed number of 
// dependencies, this is just where I've arbitrarily drawn the line for now.
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_UpdatePolicyBase_hh__
#define __Spheral_UpdatePolicyBase_hh__

#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <initializer_list>

// Since std::make_shared has problems with our use of initializer_lists, we
// include a convenience function here for making policies.
#include "make_policy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class UpdatePolicyBase {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef std::string KeyType;

  // Constructors, destructor.
  UpdatePolicyBase(std::initializer_list<std::string> depends = {});
  virtual ~UpdatePolicyBase() {};
  
  // The methods you overload to define behavior mapping new state field 
  // from derivative methods.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) = 0;

  // An alternate method to be called when you want to specify that the derivative information
  // should be assumed to not necessarily be properly time-centered, and therefore you should 
  // only use time advancement ideas, no "replace" or more sophisticated approaches.
  // Default to just calling the generic method.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) {
    this->update(key, state, derivs, multiplier, t, dt);
  }

  // Require descendents to define an equivalence operator.
  virtual bool operator==(const UpdatePolicyBase& rhs) const = 0;
  bool operator!=(const UpdatePolicyBase& rhs) const;

  // Should this policy be cloned per Field when registering for a FieldList?
  virtual bool clonePerField() const { return false; }

  // Return the set of field names that this state depends upon (if any).
  const std::vector<std::string>& dependencies() const;

  // Allow the addition of new dependencies.
  void addDependency(const std::string& depend);

  // The wildcard string for comparing dependency keys.
  static const std::string wildcard() { return "*"; }

  //............................................................................
  // Methods to support implicit time advancement
  // Test is this policy is for independent or dependent state, where independent means
  // this should be treated as implicitly advanced state to be solved for
  virtual bool independent() const;
  virtual bool dependent() const { return not this->independent(); }

  // Serialize the state we're updating to a std::vector<double> -- needed for packing State data in implicit time solve
  virtual void serializeData(std::vector<double>& buf,
                             const KeyType& key,
                             const State<Dimension>& state) const;
  virtual size_t deserializeData(const std::vector<double>& buf,
                                 const KeyType& key,
                                 const State<Dimension>& state,
                                 const size_t offset) const;
  virtual void serializeDerivatives(std::vector<double>& buf,
                                    const KeyType& key,
                                    const StateDerivatives<Dimension>& derivs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<std::string> mDependencies;

  UpdatePolicyBase(const UpdatePolicyBase& rhs);
  UpdatePolicyBase& operator=(const UpdatePolicyBase& rhs);
};

}

#include "UpdatePolicyBaseInline.hh"

#endif
