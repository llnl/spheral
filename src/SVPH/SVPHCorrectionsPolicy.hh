//---------------------------------Spheral++----------------------------------//
// SVPHCorrectionsPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the SVPHCorrections in the state.
//
// Created by JMO, Fri Aug  9 15:24:04 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SVPHCorrectionsPolicy_hh__
#define __Spheral_SVPHCorrectionsPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/DataBase.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class SVPHCorrectionsPolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  SVPHCorrectionsPolicy(const DataBase<Dimension>& dataBase,
                        const TableKernel<Dimension>& kernel);
  virtual ~SVPHCorrectionsPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  SVPHCorrectionsPolicy() = delete;
  SVPHCorrectionsPolicy(const SVPHCorrectionsPolicy& rhs) = delete;
  SVPHCorrectionsPolicy& operator=(const SVPHCorrectionsPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dimension>& mDataBase;
  const TableKernel<Dimension>& mKernel;
};

}

#endif
