//---------------------------------Spheral++----------------------------------//
// MFVIncrementVelocityPolicy -- specialized policy for hydros that allow for mass
//                      flux between nodes. The momentum time derivative
//                      is used to update the velocity. The "hydro acceleration"
//                      is also added in to be compatible w/ phys packages
//                      that apply a pure acceleration. 
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
// TODO : HydroAcceleration needs to be added in
//----------------------------------------------------------------------------//

#ifndef __Spheral_MFVIncrementVelocityPolicy_hh__
#define __Spheral_MFVIncrementVelocityPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

template<typename Dimension>
class MFVIncrementVelocityPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Vector> {
public:

  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename FieldUpdatePolicy<Dimension, Vector>::KeyType;

  // Constructors, destructor.
  MFVIncrementVelocityPolicy(std::initializer_list<std::string> depends={});
   ~MFVIncrementVelocityPolicy() = default;
  
  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "delta "; }
  
  // Forbidden methods
  MFVIncrementVelocityPolicy(const MFVIncrementVelocityPolicy& rhs) = delete;
  MFVIncrementVelocityPolicy& operator=(const MFVIncrementVelocityPolicy& rhs) = delete;
};

}

#endif
