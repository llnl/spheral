//---------------------------------Spheral++----------------------------------//
// InverseEquivalentStressDeviatorPolicy -- updates the inverse equivalent stress
//      deviator for the FSISPH hydro package. This allows yield applied in
//      a pairwise manner during the eval derivs loop.
//
// J.M. Pearl 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_InverseEquivalentStressDeviatorPolicy_hh__
#define __Spheral_InverseEquivalentStressDeviatorPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class InverseEquivalentStressDeviatorPolicy: public FieldUpdatePolicy<Dimension,typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension,Scalar>::KeyType;

  // Constructors, destructor.
  InverseEquivalentStressDeviatorPolicy();
  virtual ~InverseEquivalentStressDeviatorPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  InverseEquivalentStressDeviatorPolicy(const InverseEquivalentStressDeviatorPolicy& rhs);
  InverseEquivalentStressDeviatorPolicy& operator=(const InverseEquivalentStressDeviatorPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class InverseEquivalentStressDeviatorPolicy;
}

#endif
