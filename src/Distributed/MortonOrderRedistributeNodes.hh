//---------------------------------Spheral++----------------------------------//
// MortonOrderRedistributeNodes
//
// Attempt to redistribute nodes such that they are laid out in memory
// in a Morton ordering.  Note that this involves renumbering the nodes of 
// each NodeList, not just redistributing them between processors.
//
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Tue Mar 25 14:19:18 PDT 2008
//----------------------------------------------------------------------------//
#ifndef MortonOrderRedistributeNodes_HH
#define MortonOrderRedistributeNodes_HH

#include "SpaceFillingCurveRedistributeNodes.hh"

namespace Spheral {
  template<typename Dimension> class DataBase;
}

namespace Spheral {

template<typename Dimension>
class MortonOrderRedistributeNodes: public SpaceFillingCurveRedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using Key = KeyTraits::Key;

  // Constructors
  MortonOrderRedistributeNodes(const double minNodesPerDomainFraction = 0.5,
                               const double maxNodesPerDomainFraction = 1.5,
                               const bool workBalance = true,
                               const bool localReorderOnly = false);
  virtual ~MortonOrderRedistributeNodes() = default;

  // Hash the positions.
  virtual
  FieldList<Dimension, Key> 
  computeHashedIndices(const DataBase<Dimension>& dataBase) const override;

  // No copy or assignment operations.
  MortonOrderRedistributeNodes(const MortonOrderRedistributeNodes& nodes) = delete;
  MortonOrderRedistributeNodes& operator=(const MortonOrderRedistributeNodes& rhs) = delete;

};

}

#endif
