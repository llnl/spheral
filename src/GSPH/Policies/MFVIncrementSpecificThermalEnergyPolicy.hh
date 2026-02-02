//---------------------------------Spheral++----------------------------------//
// MFVIncrementSpecificThermalEnergyPolicy -- This is a specialized increment
//            policy for the specific thermal energy for schemes that allow
//            for flux between nodes. The specific thermal energy is updated
//            based on the time derivative of thermal energy. The mass and 
//            time derivative are needed to got from thermal to specific 
//            thermal. 
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
// TODO: the edge case handing for m->0 needs to be improved to robustly
//       handle void when full Eulerian.
//----------------------------------------------------------------------------//

#ifndef __Spheral_MFVIncrementSpecificThermalEnergyPolicy_hh__
#define __Spheral_MFVIncrementSpecificThermalEnergyPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

template<typename Dimension>
class MFVIncrementSpecificThermalEnergyPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:

  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  MFVIncrementSpecificThermalEnergyPolicy(std::initializer_list<std::string> depends={});
  ~MFVIncrementSpecificThermalEnergyPolicy() = default;
  
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
  MFVIncrementSpecificThermalEnergyPolicy(const MFVIncrementSpecificThermalEnergyPolicy& rhs) = delete;
  MFVIncrementSpecificThermalEnergyPolicy& operator=(const MFVIncrementSpecificThermalEnergyPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const std::string mStateKey;
  const std::string mDerivativeKey;
};

}

#endif
