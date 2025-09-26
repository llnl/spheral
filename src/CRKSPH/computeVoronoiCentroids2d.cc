//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on convex hulls.
//------------------------------------------------------------------------------
#include "polytope/polytope.hh"

#include "computeVoronoiCentroids.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
FieldList<Dim<2>, Dim<2>::Vector>
computeVoronoiCentroids(const FieldList<Dim<2>, Dim<2>::Vector>& position) {

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();

  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Copy the input positions to polytope's convention.
  // For convenience we gather all ghost node coordinates at the end of the
  // list.
  vector<double> coords;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = position[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Vector& ri = position(nodeListi, i);
      std::copy(ri.begin(), ri.end(), std::back_inserter(coords));
    }
  }
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned firstGhostNode = position[nodeListi]->nodeListPtr()->firstGhostNode();
    const unsigned n = position[nodeListi]->nodeListPtr()->numGhostNodes();
    for (unsigned i = 0; i != n; ++i) {
      const Vector& ri = position(nodeListi, firstGhostNode + i);
      std::copy(ri.begin(), ri.end(), std::back_inserter(coords));
    }
  }
  CHECK(coords.size() == 2*numGens);

  // Do the polytope tessellation.
  polytope::Tessellation<2, double> tessellation;
  {
#ifdef USE_MPI
    polytope::DistributedTessellator<2, double> tessellator
#if defined USE_TRIANGLE
      (new polytope::TriangleTessellator<double>(),
#else
      (new polytope::BoostTessellator<double>(),
#endif
       true,     // Manage memory for serial tessellator
       true);    // Build parallel connectivity
#else
#if defined USE_TRIANGLE
    polytope::TriangleTessellator<double> tessellator;
#else
    polytope::BoostTessellator<double> tessellator;
#endif
#endif
    tessellator.tessellate(coords, tessellation);
  }
  CHECK(tessellation.cells.size() == numGens);

  // Prepare the result.
  FieldList<Dim<2>, Vector> result = position;
  result.copyFields();

  // Copy the centroids from the tessellation to the result.
  const vector<double>& nodes = tessellation.nodes;
  const vector<vector<unsigned> >& faces = tessellation.faces;
  unsigned icell = 0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = result[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i, ++icell) {
      const Vector& ri = position(nodeListi, i);
      Vector& centroidi = result(nodeListi, i);
      centroidi.Zero();
      const vector<int>& faceIDs = tessellation.cells[icell];
      const unsigned nfaces = faceIDs.size();
      CHECK(nfaces >= 3);
      Scalar areasum = 0.0;
      for (unsigned j = 0; j != nfaces; ++j) {
        int iface = faceIDs[j] < 0 ? ~faceIDs[j] : faceIDs[j];
        CHECK(faces[iface].size() == 2);
        const Vector a = Vector(nodes[2*faces[iface][0]], nodes[2*faces[iface][0] + 1]);
        const Vector b = Vector(nodes[2*faces[iface][1]], nodes[2*faces[iface][1] + 1]);
        const Scalar area = abs(((a - ri).cross(b - ri)).z());
        areasum += area;
        centroidi += area*(ri + a + b);
      }
      centroidi /= 3*areasum;
    }
  }

  // That's it.
  return result;
}

}

