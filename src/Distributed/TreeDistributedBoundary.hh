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
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename TreeNeighbor<Dimension>::CellKey CellKey;
  typedef typename DistributedBoundary<Dimension>::DomainBoundaryNodes DomainBoundaryNodes;

  // This method returns the singleton instance of the object.
  static TreeDistributedBoundary& instance();

  // Destructor.
  virtual ~TreeDistributedBoundary();

  //**********************************************************************
  // Set the ghost nodes based on the NodeLists in the given DataBase.
  virtual void setAllGhostNodes(DataBase<Dimension>& dataBase) override;
  //**********************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Disabled methods.
  TreeDistributedBoundary();
  TreeDistributedBoundary(const TreeDistributedBoundary&);
  TreeDistributedBoundary& operator=(const TreeDistributedBoundary&);

  // Private methods.
  const TreeNeighbor<Dimension>* getTreeNeighborPtr(const NodeList<Dimension>* nodeListPtr) const;
  std::vector<std::vector<CellKey>> flattenTrees(const DataBase<Dimension>& dataBase) const;
  void broadcastTree(const std::vector<std::vector<CellKey>>& localTree);
  void buildSendNodes(const DataBase<Dimension>& dataBase,
                      const std::vector<std::vector<CellKey>>& localTree);
};

}

#endif
