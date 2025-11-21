//---------------------------------Spheral++----------------------------------//
// PolygonalFace -- 2-D face class.
//
// Created by JMO, Thu Oct 14 21:21:02 PDT 2010
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Mesh::Face(...)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<2> >::Face::
Face(const Mesh<Dim<2> >& mesh,
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
//   REQUIRE(mZone1ID < mMeshPtr->mZones.size() or mZone1ID == UNSETID);
//   REQUIRE(mZone2ID < mMeshPtr->mZones.size() or mZone2ID == UNSETID);
  REQUIRE(not (mZone1ID ==  (int)UNSETID and mZone2ID ==  (int)UNSETID));
  REQUIRE(not (mZone1ID == ~(int)UNSETID and mZone2ID ==  (int)UNSETID));
  REQUIRE(not (mZone1ID ==  (int)UNSETID and mZone2ID == ~(int)UNSETID));
  REQUIRE(not (mZone1ID == ~(int)UNSETID and mZone2ID == ~(int)UNSETID));
  REQUIRE(mEdgeIDs.size() == 1);
  REQUIRE(mEdgeIDs[0] < mMeshPtr->mEdges.size());
  mNodeIDs.push_back(mMeshPtr->mEdges[mEdgeIDs[0]].node1ID());
  mNodeIDs.push_back(mMeshPtr->mEdges[mEdgeIDs[0]].node2ID());
  REQUIRE(mNodeIDs.size() == 2);
  REQUIRE(mNodeIDs[0] < mMeshPtr->mNodePositions.size() and 
          mNodeIDs[1] < mMeshPtr->mNodePositions.size());
}

//------------------------------------------------------------------------------
// PolygonalFace::position
//------------------------------------------------------------------------------
template<>
Dim<2>::Vector
Mesh<Dim<2> >::Face::
position() const {
  REQUIRE(mEdgeIDs.size() == 1);
  REQUIRE(mEdgeIDs[0] == mID);
  REQUIRE(mID < mMeshPtr->mEdges.size());
  return mMeshPtr->mEdges[mEdgeIDs[0]].position();
}

//------------------------------------------------------------------------------
// PolygonalFace::area
//------------------------------------------------------------------------------
template<>
double
Mesh<Dim<2> >::Face::
area() const {
  REQUIRE(mEdgeIDs.size() == 1);
  return mMeshPtr->mEdges[mEdgeIDs[0]].length();
}

//------------------------------------------------------------------------------
// PolygonalFace::unitNormal
//------------------------------------------------------------------------------
template<>
Dim<2>::Vector
Mesh<Dim<2> >::Face::
unitNormal() const {
  REQUIRE(mNodeIDs.size() == 2);
  const Vector dx = (mMeshPtr->mNodePositions[mNodeIDs[1]] - mMeshPtr->mNodePositions[mNodeIDs[0]]).unitVector();
  return Vector(dx.y(), -(dx.x()));
}

}
