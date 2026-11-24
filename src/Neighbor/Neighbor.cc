//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
// Split into Value/View pattern Tue Sep 22 13:55:08 PDT 2026
//----------------------------------------------------------------------------//
#include "Neighbor/Neighbor.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>::
Neighbor(DataBase<Dimension>& dataBase):
  mDataBase(dataBase),
  mNodeExtent(dataBase.newGlobalFieldList(Vector::zero(), "node extents")) {
}

//------------------------------------------------------------------------------
// Update the neighbor info for just the nodes of a NodeList
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
updateNodes(const NodeList<Dimension>& nodes) {
  const auto k = mDataBase.nodeListIndex(nodes);
  CHECK(k < mDataBase.numNodeLists());
  const auto n = nodes.numNodes();
  std::vector<NodeID> ids(n);
  for (auto i = 0u; i < n; ++i) ids[i] = std::make_pair(k, i);
  this->updateNodes(ids);
}

}
