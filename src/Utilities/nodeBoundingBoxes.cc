//---------------------------------Spheral++----------------------------------//
// nodeBoundingBoxes
//
// Compute minimum volume bounding boxes for nodes and their extent.
//
// Created by JMO, Sun Jan 24 16:13:16 PST 2010
//----------------------------------------------------------------------------//
#include "nodeBoundingBoxes.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/DataBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// The bounding boxes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, std::pair<typename Dimension::Vector, typename Dimension::Vector> >
nodeBoundingBoxes(const NodeList<Dimension>& nodes,
                  const typename Dimension::Scalar kernelExtent) {
  using Vector = typename Dimension::Vector;
  using Box = std::pair<Vector, Vector>;

  Field<Dimension, Box> result("NodeList bounding boxes", nodes);
  const auto& positions = nodes.positions();
  const auto& Hfield = nodes.Hfield();
  for (auto i = 0u; i < nodes.numNodes(); ++i) {
    result(i) = boundingBox<Dimension>(positions(i), Hfield(i), kernelExtent);
  }
  return result;
}


//------------------------------------------------------------------------------
// The bounding boxes for all nodes in a DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, std::pair<typename Dimension::Vector, typename Dimension::Vector> >
nodeBoundingBoxes(const DataBase<Dimension>& dataBase) {
  using Vector = typename Dimension::Vector;
  using Box = std::pair<Vector, Vector>;

  FieldList<Dimension, Box> result = dataBase.newGlobalFieldList(Box(), "Bounding boxes");
  const auto positions = dataBase.globalPosition();
  const auto Hfield = dataBase.globalHfield();
  const auto kernelExtent = dataBase.kernelExtent();
  for (auto [nodeListi, nptr]: enumerate(dataBase.nodeListPtrs())) {
    const auto numNodes = nptr->numNodes();
    for (auto i = 0u; i < numNodes; ++i) {
      result(nodeListi, i) = boundingBox<Dimension>(positions(nodeListi, i), Hfield(nodeListi, i), kernelExtent);
    }
  }
  return result;
}

}

