//---------------------------------Spheral++----------------------------------//
// TreeDistributedBoundary -- Implementation of the Distributed Boundary
// condition for use with TreeNeighbor based NodeLists.
//
// Created by JMO, Mon Aug 27 21:57:51 PDT 2001
//----------------------------------------------------------------------------//
#ifndef __Spheral_TreeDistributedBoundary__
#define __Spheral_TreeDistributedBoundary__

#include "DistributedBoundary.hh"
#include "Neighbor/TreeNeighbor.hh"

#include <vector>

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class NodeList;

template<typename Dimension>
class TreeDistributedBoundary: public DistributedBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using CellKey = typename TreeNeighbor<Dimension>::CellKey;
  using NodeID = typename TreeNeighbor<Dimension>::NodeID;
  using DomainBoundaryNodes = typename DistributedBoundary<Dimension>::DomainBoundaryNodes;

  // This method returns the singleton instance of the object.
  static TreeDistributedBoundary& instance();

  // Destructor.
  virtual ~TreeDistributedBoundary() = default;

  //**********************************************************************
  // Set the ghost nodes based on the NodeLists in the given DataBase.
  virtual void setAllGhostNodes(DataBase<Dimension>& dataBase) override;
  //**********************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Private methods.
  const TreeNeighbor<Dimension>* getTreeNeighborPtr(const DataBase<Dimension>& dataBase) const;
  void broadcastTree(const std::vector<std::vector<CellKey>>& localTree);
  void buildSendNodes(const DataBase<Dimension>& dataBase,
                      const std::vector<std::vector<CellKey>>& localTree);

  // Disabled methods.
  TreeDistributedBoundary();
  TreeDistributedBoundary(const TreeDistributedBoundary&);
  TreeDistributedBoundary& operator=(const TreeDistributedBoundary&);
};

}

#endif
