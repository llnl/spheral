//---------------------------------Spheral++----------------------------------//
// globalBoundingVolumes
//
// Compute minimum bounding volumes (convex hull) for the nodes and their 
// sampling extents in the given DataBase.
// The two values computed are:
// 1.  nodeVolume -- the box containing the node positions.
// 2.  sampleVolume -- the box containing the node extents.
//
// Created by JMO, Sun Jan 31 19:53:36 PST 2010
//----------------------------------------------------------------------------//
#include "DataBase/DataBase.hh"
#include "Distributed/allReduce.hh"
#include "Geometry/Dimension.hh"

#include <vector>
#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Add the sampling position box bounds to a vector.
//------------------------------------------------------------------------------
inline
void
appendSamplingPositions(const Dim<1>::Vector& position,
                        const Dim<1>::Vector& extent,
                        vector<Dim<1>::Vector>& result) {
  typedef Dim<1>::Vector Vector;
  result.push_back(Vector(position.x() - extent.x()));
  result.push_back(Vector(position.x() + extent.x()));
}

inline
void
appendSamplingPositions(const Dim<2>::Vector& position,
                        const Dim<2>::Vector& extent,
                        vector<Dim<2>::Vector>& result) {
  typedef Dim<2>::Vector Vector;
  const double xi = position.x();
  const double yi = position.y();
  const double xexti = extent.x();
  const double yexti = extent.y();
  result.push_back(Vector(xi - xexti, yi - yexti));
  result.push_back(Vector(xi + xexti, yi - yexti));
  result.push_back(Vector(xi - xexti, yi + yexti));
  result.push_back(Vector(xi + xexti, yi + yexti));
}

inline
void
appendSamplingPositions(const Dim<3>::Vector& position,
                        const Dim<3>::Vector& extent,
                        vector<Dim<3>::Vector>& result) {
  typedef Dim<3>::Vector Vector;
  const double xi = position.x();
  const double yi = position.y();
  const double zi = position.z();
  const double xexti = extent.x();
  const double yexti = extent.y();
  const double zexti = extent.z();
  result.push_back(Vector(xi - xexti, yi - yexti, zi - zexti));
  result.push_back(Vector(xi + xexti, yi - yexti, zi - zexti));
  result.push_back(Vector(xi - xexti, yi + yexti, zi - zexti));
  result.push_back(Vector(xi + xexti, yi + yexti, zi - zexti));
  result.push_back(Vector(xi - xexti, yi - yexti, zi + zexti));
  result.push_back(Vector(xi + xexti, yi - yexti, zi + zexti));
  result.push_back(Vector(xi - xexti, yi + yexti, zi + zexti));
  result.push_back(Vector(xi + xexti, yi + yexti, zi + zexti));
}

//------------------------------------------------------------------------------
// Compute the minimum volume box containing all the points in a Field of
// positions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
globalBoundingBox(const Field<Dimension, typename Dimension::Vector>& positions,
                  typename Dimension::Vector& xmin,
                  typename Dimension::Vector& xmax,
                  const bool ghost) {
  typedef typename Dimension::Vector Vector;

  // Find our local bounds.
  xmin = DBL_MAX;
  xmax = -DBL_MAX;
  const unsigned n = ghost ? positions.numElements() : positions.numInternalElements();
  for (unsigned i = 0; i != n; ++i) {
    const Vector& xi = positions(i);
    xmin = elementWiseMin(xmin, xi);
    xmax = elementWiseMax(xmax, xi);
  }

  // Now find the global bounds across all processors.
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    xmin(i) = allReduce(xmin(i), SPHERAL_OP_MIN);
    xmax(i) = allReduce(xmax(i), SPHERAL_OP_MAX);
  }
}

//------------------------------------------------------------------------------
// Compute the minimum volume box containing all the points in a FieldList of
// positions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
globalBoundingBox(const FieldList<Dimension, typename Dimension::Vector>& positions,
                  typename Dimension::Vector& xmin,
                  typename Dimension::Vector& xmax,
                  const bool ghost) {
  typedef typename Dimension::Vector Vector;

  // Find our local bounds.
  xmin = DBL_MAX;
  xmax = -DBL_MAX;
  for (unsigned nodeList = 0; nodeList != positions.numFields(); ++nodeList) {
    const unsigned n = ghost ? positions[nodeList]->numElements() : positions[nodeList]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Vector& xi = positions(nodeList, i);
      xmin = elementWiseMin(xmin, xi);
      xmax = elementWiseMax(xmax, xi);
    }
  }

  // Now find the global bounds across all processors.
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    xmin(i) = allReduce(xmin(i), SPHERAL_OP_MIN);
    xmax(i) = allReduce(xmax(i), SPHERAL_OP_MAX);
  }
}

//------------------------------------------------------------------------------
// Compute the minimum volume FacetedVolume containing all the nodes in the 
// DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
globalBoundingVolumes(const DataBase<Dimension>& dataBase,
                      typename Dimension::ConvexHull& nodeVolume,
                      typename Dimension::ConvexHull& sampleVolume) {
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::ConvexHull ConvexHull;

  // Extract the node positions as a single vector.
  const size_t numNodes = dataBase.numNodes();
  const size_t numSamples = (1U << Dimension::nDim) * numNodes;
  const FieldList<Dimension, Vector> positions = dataBase.globalPosition();
  const FieldList<Dimension, Vector> extents = dataBase.globalNodeExtent();
  vector<Vector> nodePositions, samplePositions;
  nodePositions.reserve(numNodes);
  samplePositions.reserve(numSamples);
  for (size_t fieldi = 0u; fieldi < positions.numFields(); ++fieldi) {
    for (size_t i = 0u; i < positions[fieldi]->numElements(); ++i) {
      nodePositions.push_back(positions(fieldi, i));
      appendSamplingPositions(positions(fieldi, i), extents(fieldi, i), samplePositions);
    }
  }
  CHECK(nodePositions.size() == numNodes);
  CHECK(samplePositions.size() == numSamples);

  // Build the node bounding volume.
  nodeVolume = ConvexHull(nodePositions);

  // Throw away sampling positions that are inside the node hull.
  vector<size_t> kill;
  for (size_t i = 0; i != numSamples; ++i) {
    if (nodeVolume.contains(samplePositions[i])) kill.push_back(i);
  }
  removeElements(samplePositions, kill);

  // Build the sampling bounding volume.
  sampleVolume = ConvexHull(samplePositions);
}

}

