//---------------------------------Spheral++----------------------------------//
// StateDerivatives -- Accumulate and cart around the derivatives/changes to
// the state for a set of physics packages.
//
// Created by JMO, Fri Aug 27 15:50:37 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StateDerivatives_hh__
#define __Spheral_StateDerivatives_hh__

#include "StateBase.hh"
#include "Field/Field.hh"
#include "Field/NodeIteratorBase.hh"

#include <vector>
#include <map>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class DataBase;
template<typename Dimension> class Physics;

template<typename Dimension>
class StateDerivatives: public StateBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PackageList = std::vector<Physics<Dimension>*>;
  using PackageIterator = typename PackageList::iterator;

  using KeyType = typename StateBase<Dimension>::KeyType;

  // Constructors, destructor
  StateDerivatives(DataBase<Dimension>& dataBase, PackageList& physicsPackage);
  StateDerivatives(DataBase<Dimension>& dataBase,
                   PackageIterator physicsPackageBegin,
                   PackageIterator physicsPackageEnd);
  StateDerivatives() = default;
  StateDerivatives(const StateDerivatives& rhs) = default;
  StateDerivatives& operator=(const StateDerivatives& rhs) = default;
  virtual ~StateDerivatives() = default;

  // Test if two StateDerivatives have equivalent fields.
  virtual bool operator==(const StateBase<Dimension>& rhs) const override;

  // Force all derivative FieldLists to zero.
  void Zero();

private:
  //--------------------------- Private Interface ---------------------------//
  using StateBase<Dimension>::mStorage;
};

}

#endif

