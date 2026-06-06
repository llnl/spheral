//---------------------------------Spheral++----------------------------------//
// PolyhedralFace -- 3-D face class.
//
// Created by JMO, Wed Jan  5 21:09:11 PST 2011
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Mesh::Face(...)
// Note we assume that the edges passed in are already sorted counter-clockwise
// from the desired orientation.
//------------------------------------------------------------------------------
template<>
Mesh<Dim<3> >::Face::
Face(const Mesh<Dim<3> >& mesh,
     const unsigned ID,
     const int zone1ID,
     const int zone2ID,
     const vector<unsigned>& edgeIDs):
  mMeshPtr(&mesh),
  mID(ID),
  mZone1ID(zone1ID),
  mZone2ID(zone2ID),
  mNodeIDs(),
  mEdgeIDs(edgeIDs) {
  REQUIRE(not (mZone1ID ==  (int)UNSETID and mZone2ID ==  (int)UNSETID));
  REQUIRE(not (mZone1ID == ~(int)UNSETID and mZone2ID ==  (int)UNSETID));
  REQUIRE(not (mZone1ID ==  (int)UNSETID and mZone2ID == ~(int)UNSETID));
  REQUIRE(not (mZone1ID == ~(int)UNSETID and mZone2ID == ~(int)UNSETID));
  REQUIRE(mEdgeIDs.size() > 2);
  REQUIRE(*std::max_element(mEdgeIDs.begin(), mEdgeIDs.end()) < mMeshPtr->mEdges.size());

  // Build the node IDs.  As stated above, this relies on the edges already being 
  // sorted in the desired counter-clockwise order when passed in.
  const unsigned n1 = mMeshPtr->mEdges[edgeIDs[0]].node1ID();
  const unsigned n2 = mMeshPtr->mEdges[edgeIDs[0]].node2ID();
  const unsigned n3 = mMeshPtr->mEdges[edgeIDs[1]].node1ID();
  const unsigned n4 = mMeshPtr->mEdges[edgeIDs[1]].node2ID();
  CHECK2((n1 == n3 ? 1 : 0) +
         (n1 == n4 ? 1 : 0) +
         (n2 == n3 ? 1 : 0) +
         (n2 == n4 ? 1 : 0) == 1, n1 << " " << n2 << " " << n3 << " " << n4);
  mNodeIDs.push_back(n1 == n3 ? n1 :
                     n1 == n4 ? n1 : 
                                n2);
  for (vector<unsigned>::const_iterator itr = edgeIDs.begin() + 1;
       itr != edgeIDs.end();
       ++itr) {
    const unsigned n1 = mMeshPtr->mEdges[*itr].node1ID();
    const unsigned n2 = mMeshPtr->mEdges[*itr].node2ID();
    CHECK((n1 == mNodeIDs.back() ? 1 : 0) +
          (n2 == mNodeIDs.back() ? 1 : 0) == 1);
    mNodeIDs.push_back(n1 == mNodeIDs.back() ? n2 : n1);
  }
  CHECK(mNodeIDs.size() == mEdgeIDs.size());
}

//------------------------------------------------------------------------------
// PolyhedralFace::position
//------------------------------------------------------------------------------
template<>
Dim<3>::Vector
Mesh<Dim<3> >::Face::
position() const {
  REQUIRE(mNodeIDs.size() > 2);
  Vector result;
  for (vector<unsigned>::const_iterator itr = mNodeIDs.begin();
       itr != mNodeIDs.end();
       ++itr) result += mMeshPtr->mNodePositions[*itr];
  result /= mNodeIDs.size();
  return result;
}

//------------------------------------------------------------------------------
// PolyhedralFace::area
// This logic depends on a face being convex!
//------------------------------------------------------------------------------
template<>
double
Mesh<Dim<3> >::Face::
area() const {
  REQUIRE(mNodeIDs.size() > 2);
  double result = 0.0;
  const Vector origin = position();
  unsigned i, j;
  for (i = 0; i != mNodeIDs.size(); ++i) {
    j = (i + 1) % mNodeIDs.size();
    result += ((mMeshPtr->mNodePositions[mNodeIDs[i]] - origin).cross
               (mMeshPtr->mNodePositions[mNodeIDs[j]] - origin).magnitude());
  }
  ENSURE(result > 0.0);
  return 0.5*result;
}

//------------------------------------------------------------------------------
// PolyhedralFace::unitNormal
//------------------------------------------------------------------------------
template<>
Dim<3>::Vector
Mesh<Dim<3> >::Face::
unitNormal() const {
  REQUIRE(mNodeIDs.size() > 2);
  Vector result;
  const Vector origin = position();
  unsigned i, j;
  for (i = 0; i != mNodeIDs.size(); ++i) {
    j = (i + 1) % mNodeIDs.size();
    result += ((mMeshPtr->mNodePositions[mNodeIDs[i]] - origin).cross
               (mMeshPtr->mNodePositions[mNodeIDs[j]] - origin));
  }
  return result.unitVector();
}

}
