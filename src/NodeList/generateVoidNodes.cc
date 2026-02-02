//------------------------------------------------------------------------------
// generateVoidNodes
//
// This algorithm tries to analyze how continuous a node distribution is, and 
// if it determines there is an edge to the distribution creates new void nodes
// outside that surface.
// We assume here that the caller has already created all the boundary ghost 
// nodes.
//------------------------------------------------------------------------------
#include "generateVoidNodes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/BSplineKernel.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/pointDistances.hh"
#include "Mesh/Mesh.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Overloaded method to find the closest point on a face for each dimension.
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
closestPointOnFace(const Mesh<Dim<1> >& /*mesh*/,
                   const Mesh<Dim<1> >::Face& face,
                   const Dim<1>::Vector /*ri*/) {
  return face.position();
}

inline
Dim<2>::Vector
closestPointOnFace(const Mesh<Dim<2> >& mesh,
                   const Mesh<Dim<2> >::Face& face,
                   const Dim<2>::Vector ri) {
  const vector<unsigned>& nodeIDs = face.nodeIDs();
  REQUIRE(nodeIDs.size() == 2);
  return closestPointOnSegment(ri, 
                               mesh.node(nodeIDs[0]).position(),
                               mesh.node(nodeIDs[1]).position());
}

inline
Dim<3>::Vector
closestPointOnFace(const Mesh<Dim<3> >& mesh,
                   const Mesh<Dim<3> >::Face& face,
                   const Dim<3>::Vector ri) {
  return closestPointOnPlane(ri, 
                             mesh.node(face.nodeIDs()[0]).position(),
                             face.unitNormal());
}

//------------------------------------------------------------------------------
// The actual method.
//------------------------------------------------------------------------------
template<typename Dimension>
void generateVoidNodes(const vector<typename Dimension::Vector>& generators,
                       const vector<typename Dimension::SymTensor>& Hs,
                       const Mesh<Dimension>& mesh,
                       const typename Dimension::Vector& /*xmin*/,
                       const typename Dimension::Vector& /*xmax*/,
                       const unsigned numInternal,
                       const double nPerh,
                       const double threshold,
                       NodeList<Dimension>& voidNodes) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ConvexHull ConvexHull;

  // Pre-conditions.
  VERIFY(generators.size() >= numInternal);
  VERIFY(generators.size() == Hs.size());
  VERIFY(mesh.numZones() >= numInternal);
  VERIFY(mesh.numZones() <= generators.size());
  VERIFY(voidNodes.numNodes() == 0);
  VERIFY(nPerh > 0.0);
  VERIFY(threshold >= 0.0);

  unsigned i, j, k, ii, ivoid = 0;
  Vector ri, rj, rij, etai, etaj, fmhat;
  SymTensor Hi, Hj;
  Field<Dimension, Vector>& vpos = voidNodes.positions();
  Field<Dimension, SymTensor>& vH = voidNodes.Hfield();
  for (i = 0; i != numInternal; ++i) {
    ri = generators[i];
    Hi = Hs[i];

    // Generate the convex hull of this zones and it's neighbors.
    const vector<unsigned>& nodeIDs = mesh.zone(i).nodeIDs();
    vector<Vector> neighborNodes;
    for (k = 0; k != nodeIDs.size(); ++k) {
      const vector<unsigned>& zoneIDs = mesh.node(nodeIDs[k]).zoneIDs();
      for (j = 0; j != zoneIDs.size(); ++j) {
        if (zoneIDs[j] != Mesh<Dimension>::UNSETID) {
          CHECK2(zoneIDs[j] < mesh.numZones(), zoneIDs[j] << " " << mesh.numZones());
          const vector<unsigned>& otherNodeIDs = mesh.zone(zoneIDs[j]).nodeIDs();
          for (ii = 0; ii != otherNodeIDs.size(); ++ii) {
            neighborNodes.push_back(mesh.node(otherNodeIDs[ii]).position());
          }
        }
      }
    }
    const ConvexHull neighborHull(neighborNodes);

    // Look for any faces of this node's cell across which we 
    // should be generating a void point.
    const vector<int>& faceIDs = mesh.zone(i).faceIDs();
    for (const int faceID: faceIDs) {
      const typename Mesh<Dimension>::Face& face = mesh.face(faceID);
      j = Mesh<Dimension>::positiveID(face.oppositeZoneID(i));
      CHECK(j == Mesh<Dimension>::UNSETID or j < generators.size());
      if (j != Mesh<Dimension>::UNSETID) {
        rj = generators[j];
        Hj = Hs[j];
        rij = ri - rj;
        etai = Hi*rij;
        etaj = Hj*rij;
      }
      if (j == Mesh<Dimension>::UNSETID or 
          (etai.magnitude() > threshold and etaj.magnitude() > threshold)) {

        // Yep, generate a new void point.
        fmhat = (face.position() - ri).unitVector(); //  (closestPointOnFace(mesh, face, ri) - ri).unitVector();
        rj = ri + Hi.Inverse()/nPerh * fmhat;
        if (neighborHull.convexContains(rj)) {
// (mesh.zone(i).convexHull().convexContains(rj) or
//             (j != Mesh<Dimension>::UNSETID and mesh.zone(j).convexHull().convexContains(rj))) {
          voidNodes.numInternalNodes(ivoid + 1);
          vpos[ivoid] = rj;
          vH[ivoid] = Hi;
          ++ivoid;
        }
      }
    }
  }

  // That's it.
}

}
