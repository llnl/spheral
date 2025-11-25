//---------------------------------Spheral++----------------------------------//
// LineZone -- 1-D zone class.
//
// Created by JMO, Thu Oct 14 21:21:02 PDT 2010
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Mesh::Zone(...)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<1> >::Zone::
Zone(const Mesh<Dim<1> >& mesh,
     const unsigned ID,
     const vector<int>& faceIDs):
  mMeshPtr(&mesh),
  mID(ID),
  mNodeIDs(),
  mEdgeIDs(),
  mFaceIDs(faceIDs) {

  // Pre-conditions.
  REQUIRE(faceIDs.size() == 2);
  REQUIRE(faceIDs[0] <  0);
  REQUIRE(faceIDs[1] >= 0);

  // Iterate over the input faces and build the unique sets of nodes and edges.
  int i;
  for (vector<int>::const_iterator faceIDitr = faceIDs.begin();
       faceIDitr != faceIDs.end();
       ++faceIDitr) {
    i = *faceIDitr;
    if (i < 0) i = ~i;
    const Face& face = mMeshPtr->face(i);
    CHECK(face.numNodes() == 1);
    CHECK(face.numEdges() == 1);

    // Add this face's edges and nodes to the zones sets.
    mNodeIDs.push_back(face.nodeIDs()[0]);
    mEdgeIDs.push_back(face.edgeIDs()[0]);
  }
  CHECK(mNodeIDs.size() == 2);
  CHECK(mEdgeIDs.size() == 2);

  // Order the faces, edges, & nodes as (left, right).
  if (mMeshPtr->node(mNodeIDs[0]).position().x() > mMeshPtr->node(mNodeIDs[1]).position().x()) std::swap(mNodeIDs[0], mNodeIDs[1]);
  if (mMeshPtr->edge(mEdgeIDs[0]).position().x() > mMeshPtr->edge(mEdgeIDs[1]).position().x()) std::swap(mEdgeIDs[0], mEdgeIDs[1]);

  // Post-conditions.
  ENSURE(mNodeIDs.size() == 2);
  ENSURE(mEdgeIDs.size() == 2);
  ENSURE(mFaceIDs.size() == 2);
  ENSURE(mMeshPtr->node(mNodeIDs[0]).position().x() < mMeshPtr->node(mNodeIDs[1]).position().x());
  ENSURE(mMeshPtr->edge(mEdgeIDs[0]).position().x() < mMeshPtr->edge(mEdgeIDs[1]).position().x());
  ENSURE(mMeshPtr->face(mFaceIDs[0]).position().x() < mMeshPtr->face(~mFaceIDs[1]).position().x());
}

//------------------------------------------------------------------------------
// LineZone::position
//------------------------------------------------------------------------------
template<>
Mesh<Dim<1> >::Vector
Mesh<Dim<1> >::Zone::
position() const {
  return 0.5*(mMeshPtr->mNodePositions[mNodeIDs[0]] +
              mMeshPtr->mNodePositions[mNodeIDs[1]]);
}

//------------------------------------------------------------------------------
// LineZone::volume
//------------------------------------------------------------------------------
template<>
double
Mesh<Dim<1> >::Zone::
volume() const {
  return (mMeshPtr->mNodePositions[mNodeIDs[1]].x() -
          mMeshPtr->mNodePositions[mNodeIDs[0]].x());
}

}
