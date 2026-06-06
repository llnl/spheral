//---------------------------------Spheral++----------------------------------//
// PolygonalZone -- 2-D zone class.
//
// Created by JMO, Tue Nov 16 15:30:39 PST 2010
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/DBC.hh"

#include <vector>
#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Comparator to help sorting the faces counter-clockwise about a position.
//------------------------------------------------------------------------------
template<typename Element, typename Vector>
struct CounterClockwiseCompareElements {
  CounterClockwiseCompareElements(const vector<Element>& elements,
                                  const unsigned originID):
    mElements(elements),
    mCentroid(elements[originID].position()) {}
  bool operator()(const unsigned id1, const unsigned id2) {
    REQUIRE(id1 < mElements.size());
    REQUIRE(id2 < mElements.size());
    return ((mElements[id1].position() - mCentroid).cross(mElements[id2].position() - mCentroid).z() > 0.0);
//     return isgn0((mElements[id1].position() - mCentroid).cross(mElements[id2].position() - mCentroid).z());
  }
  const vector<Element>& mElements;
  Vector mCentroid;
};

//------------------------------------------------------------------------------
// Mesh::Zone(...)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<2> >::Zone::
Zone(const Mesh<Dim<2> >& mesh,
     const unsigned ID,
     const vector<int>& faceIDs):
  mMeshPtr(&mesh),
  mID(ID),
  mNodeIDs(),
  mEdgeIDs(),
  mFaceIDs(faceIDs) {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(mFaceIDs.size() > 2);
    for (const int i: mFaceIDs) {
      int j = (i < 0 ? ~i : i);
      CONTRACT_VAR(j);
      REQUIRE(j < (int)mMeshPtr->mFaces.size());
      REQUIRE(mMeshPtr->mFaces[j].mEdgeIDs.size() == 1);
      REQUIRE((int)mMeshPtr->mFaces[j].mEdgeIDs[0] == j);
    }
  }
  END_CONTRACT_SCOPE
  
  // Copy the face IDs as the edge IDs (they are degenerate after all!).
  for (const int faceID: mFaceIDs) mEdgeIDs.push_back(faceID < 0 ? ~faceID : faceID);
  CHECK(mEdgeIDs.size() == mFaceIDs.size());

  // We need the nodes sorted counter-clockwise around the zone.  The faces
  // are already sorted, so we can leverage that info.  Note we are using the
  // fact that faces and edges are degenerate here.
  int i;
  unsigned n1, n2;
  for (const int faceID: mFaceIDs) {
    i = (faceID < 0 ? ~faceID : faceID);
    n1 = mMeshPtr->mEdges[i].node1ID();
    n2 = mMeshPtr->mEdges[i].node2ID();
    mNodeIDs.push_back(faceID < 0 ? n2 : n1);
  }
  CHECK(mNodeIDs.size() == mEdgeIDs.size());

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    ENSURE2(mNodeIDs.size() == mFaceIDs.size(), mID << " " << mNodeIDs.size() << " " << mFaceIDs.size());
    ENSURE2(mEdgeIDs.size() == mFaceIDs.size(), mID << " " << mEdgeIDs.size() << " " << mFaceIDs.size());
    ENSURE2(mFaceIDs.size() > 2, mID << " " << mFaceIDs.size());
    for (i = 0; i != (int)mFaceIDs.size(); ++i) {
      ENSURE2(mNodeIDs[i] < mMeshPtr->mNodes.size(), mID << " " << mNodeIDs[i] << " " <<  mMeshPtr->mNodes.size());
      ENSURE2(mEdgeIDs[i] < mMeshPtr->mEdges.size(), mID << " " << mEdgeIDs[i] << " " <<  mMeshPtr->mEdges.size());
      ENSURE2(Mesh<Dim<2> >::positiveID(mFaceIDs[i]) < mMeshPtr->mFaces.size(),
              mID << " " << mFaceIDs[i] << " " << mMeshPtr->mFaces.size());
    }

    // Make sure the elements are unique!
    vector<unsigned> nodeIDs(mNodeIDs);
    sort(nodeIDs.begin(), nodeIDs.end());
    if (!(unique(nodeIDs.begin(), nodeIDs.end()) == nodeIDs.end())) {
      cerr << "Blago!  : " << mID << " : Faces : ";
      copy(mFaceIDs.begin(), mFaceIDs.end(), std::ostream_iterator<int>(cerr, " "));
      cerr << endl << " Face nodes :";
      for (vector<int>::const_iterator itr = mFaceIDs.begin(); itr != mFaceIDs.end(); ++itr) {
        const unsigned i = Mesh<Dim<2> >::positiveID(*itr);
        cerr << " (" << mMeshPtr->mFaces[i].mNodeIDs[0] << " " << mMeshPtr->mFaces[i].mNodeIDs[1] << ")";
      }
      cerr << endl;
    }
    ENSURE(unique(nodeIDs.begin(), nodeIDs.end()) == nodeIDs.end());
    vector<unsigned> edgeIDs(mEdgeIDs);
    sort(edgeIDs.begin(), edgeIDs.end());
    ENSURE(unique(edgeIDs.begin(), edgeIDs.end()) == edgeIDs.end());
    vector<unsigned> faceIDs;
    for (const int i: mFaceIDs) faceIDs.push_back(Mesh<Dim<2> >::positiveID(i));
    sort(faceIDs.begin(), faceIDs.end());
    ENSURE(unique(faceIDs.begin(), faceIDs.end()) == faceIDs.end());

    // Make sure elements are listed counter-clockwise.
    CounterClockwiseCompareElements<Node, Vector> nodeComparator(mMeshPtr->mNodes, mNodeIDs[0]);
    CounterClockwiseCompareElements<Edge, Vector> edgeComparator(mMeshPtr->mEdges, mEdgeIDs[0]);
    CounterClockwiseCompareElements<Face, Vector> faceComparator(mMeshPtr->mFaces, faceIDs[0]);
    for (unsigned i = 0; i < faceIDs.size() - 1; ++i) {
      ENSURE2(nodeComparator(mNodeIDs[i], mNodeIDs[i + 1]),
              nodeComparator(mNodeIDs[i], mNodeIDs[i + 1]) << " "
              << mMeshPtr->mNodes[mNodeIDs[0]].position() << " "
              << mMeshPtr->mNodes[mNodeIDs[i]].position() << " "
              << mMeshPtr->mNodes[mNodeIDs[i + 1]].position());
      ENSURE2(edgeComparator(mEdgeIDs[i], mEdgeIDs[i + 1]),
              edgeComparator(mEdgeIDs[i], mEdgeIDs[i + 1]) << " "
              << mMeshPtr->mEdges[mEdgeIDs[0]].position() << " "
              << mMeshPtr->mEdges[mEdgeIDs[i]].position() << " "
              << mMeshPtr->mEdges[mEdgeIDs[i + 1]].position());
      int id1 = Mesh<Dim<2> >::positiveID(mFaceIDs[i]),
          id2 = Mesh<Dim<2> >::positiveID(mFaceIDs[i + 1]);
      CONTRACT_VAR(id1);
      CONTRACT_VAR(id2);
      ENSURE2(faceComparator(id1, id2),
              faceComparator(id1, id2) << " "
              << mMeshPtr->mFaces[Mesh<Dim<2> >::positiveID(mFaceIDs[0])].position() << " "
              << mMeshPtr->mFaces[id1].position() << " "
              << mMeshPtr->mFaces[id2].position());
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// PolygonalZone::position
//------------------------------------------------------------------------------
template<>
Dim<2>::Vector
Mesh<Dim<2> >::Zone::
position() const {
  unsigned i, j;
  double weight, weightSum = 0.0;
  Vector result;
  for (i = 0; i != mNodeIDs.size(); ++i) {
    j = (i + 1) % mNodeIDs.size();
    weight = (mMeshPtr->mNodePositions[mNodeIDs[j]] - mMeshPtr->mNodePositions[mNodeIDs[i]]).magnitude();
    weightSum += weight;
    result += weight*(mMeshPtr->mNodePositions[mNodeIDs[j]] + mMeshPtr->mNodePositions[mNodeIDs[i]]);
  }
  CHECK(weightSum > 0.0);
  result /= 2.0*weightSum;
  return result;
}

//------------------------------------------------------------------------------
// PolygonalZone::volume
//------------------------------------------------------------------------------
template<>
double
Mesh<Dim<2> >::Zone::
volume() const {
  double result = 0.0;
  const Vector centroid = this->position();
  unsigned i, j;
  for (i = 0; i != mNodeIDs.size(); ++i) {
    j = (i + 1) % mNodeIDs.size();
    const Vector x1 = mMeshPtr->mNodePositions[mNodeIDs[i]];
    const Vector x2 = mMeshPtr->mNodePositions[mNodeIDs[j]];
    const double dv = (x1 - centroid).cross(x2 - centroid).z();
    // CHECK2(fuzzyGreaterThanOrEqual(dv, 0.0, 1.0e-8),
    //        "Negative triangle!  " << dv << " " << centroid << " " << this->convexHull().contains(centroid) << " " << this->convexHull().distance(centroid));
    result += dv;
  }
  result *= 0.5;
  ENSURE2(result > 0.0, "Nonpositive volume!  " << result);
  return result;
}

}
