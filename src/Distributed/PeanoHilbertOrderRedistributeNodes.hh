//---------------------------------Spheral++----------------------------------//
// PeanoHilbertOrderRedistributeNodes
// Attempt to redistribute nodes such that they are laid out in memory
// in a PeanoHilbert ordering.  Note that this involves renumbering the nodes of 
// each NodeList, not just redistributing them between processors.
//
// Springel V. (2005), MNRAS
//
// Created by JMO, Tue Apr  8 14:25:18 PDT 2008
//----------------------------------------------------------------------------//
#ifndef PeanoHilbertOrderRedistributeNodes_HH
#define PeanoHilbertOrderRedistributeNodes_HH

#include "SpaceFillingCurveRedistributeNodes.hh"

namespace Spheral {
  template<typename Dimension> class DataBase;
}

namespace Spheral {

template<typename Dimension>
class PeanoHilbertOrderRedistributeNodes: public SpaceFillingCurveRedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using Key = typename KeyTraits::Key;

  // Constructors
  PeanoHilbertOrderRedistributeNodes(const double minNodesPerDomainFraction = 0.5,
                                     const double maxNodesPerDomainFraction = 1.5,
                                     const bool workBalance = true,
                                     const bool localReorderOnly = false);
  virtual ~PeanoHilbertOrderRedistributeNodes() = default;

  // Hash the positions.
  virtual
  FieldList<Dimension, Key> 
  computeHashedIndices(const DataBase<Dimension>& dataBase) const override;

  // No copy or assignment operations.
  PeanoHilbertOrderRedistributeNodes(const PeanoHilbertOrderRedistributeNodes& nodes) = delete;
  PeanoHilbertOrderRedistributeNodes& operator=(const PeanoHilbertOrderRedistributeNodes& rhs) = delete;
};

}

#endif
