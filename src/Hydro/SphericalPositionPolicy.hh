//---------------------------------Spheral++----------------------------------//
// SphericalPositionPolicy
//
// Specializes the IncrementFieldListPolicy for use updating the position in
// spherical coordinates.  This position does not allow points to pass through
// the origin.
//
// Created by JMO, Tue Mar 29 11:12:31 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalPositionPolicy_hh__
#define __Spheral_SphericalPositionPolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SphericalPositionPolicy: public UpdatePolicyBase<Dim<1>> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Dimension = Dim<1>;
  using Vector = Dimension::Vector;
  using UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  SphericalPositionPolicy();
  virtual ~SphericalPositionPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Should this policy be cloned per Field when registering for a FieldList?
  virtual bool clonePerField() const { return true; }

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  SphericalPositionPolicy(const SphericalPositionPolicy& rhs);
  SphericalPositionPolicy& operator=(const SphericalPositionPolicy& rhs);
};

}

#endif
