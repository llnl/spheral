//------------------------------------------------------------------------------
// Find the set of nodes that see through a pair of planes.
//------------------------------------------------------------------------------
#include "findNodesTouchingThroughPlanes.hh"
#include "Distributed/allReduce.hh"

namespace Spheral {
  
using std::vector;

template<typename Dimension>
std::vector<size_t>
findNodesTouchingThroughPlanes(const NodeList<Dimension>& nodeList,
                               const GeomPlane<Dimension>& enterPlane,
                               const GeomPlane<Dimension>& exitPlane,
                               const double hmultiplier) {
  vector<size_t> result;

  // Get the Neighbor object associated with the node list.
  auto& neighbor = nodeList.neighbor();

  // Here we switch between using the Neighbor magic or just doing an O(N)
  // search for anyone who overlaps the exit plane.
  if (false) {

    // Begin by identifying the set of master and neighbor nodes, where master
    // nodes see through the enter plane, and neighbors see through the exit plane.
    vector<int> masterList, coarseNeighbors, refineNeighbors;
    neighbor.setMasterList(enterPlane, exitPlane, masterList, coarseNeighbors);

    // Set the list of control nodes.
    // std::copy(neighbor.masterBegin(), neighbor.masterEnd(), std::back_inserter(controlNodes));
    // std::copy(neighbor.coarseNeighborBegin(), neighbor.coarseNeighborEnd(), std::back_inserter(controlNodes));
    const auto& pos = nodeList.positions();
    for (auto itr = coarseNeighbors.begin(); itr < coarseNeighbors.end(); ++itr) {
      if (exitPlane.signedDistance(pos(*itr)) >= 0.0) result.push_back(*itr);
    }

  } else {

    // Find the maximum smoothing scale of any point touching either plane.
    const auto  kernelExtent = hmultiplier * neighbor.kernelExtent();
    const auto  n = nodeList.numNodes();
    const auto& pos = nodeList.positions();
    const auto& H = nodeList.Hfield();
    auto hmax = 0.0;
    for (auto i = 0u; i != n; ++i) {
      const auto& ri = pos(i);
      const auto& Hi = H(i);
      const auto  hmaxi = 1.0/Hi.eigenValues().minElement();
      if (hmaxi > hmax and std::min(exitPlane.minimumDistance(ri), enterPlane.minimumDistance(ri)) < kernelExtent*hmaxi) hmax = hmaxi;
    }
    hmax = allReduce(hmax, SPHERAL_OP_MAX);

    // Now find all points within this range of the exit plane.
    if (hmax > 0.0) {
      for (auto i = 0u; i != n; ++i) {
        const auto& ri = pos(i);
        const auto  disti = exitPlane.signedDistance(ri)/hmax;
        // const GeomPlane<Dimension> exitPlanePrime(Hi*(exitPlane.point() - ri),
        //                                           (Hi*exitPlane.normal()).unitVector());
        // const Scalar disti = exitPlanePrime.signedDistance(Vector::zero());
        if (disti >= 0.0 and disti <= kernelExtent) result.push_back(i);
        // cerr << " --> " << i << " " << ri << " " << enterPlanePrime.minimumDistance(Vector::zero()) << " " << exitPlanePrime.minimumDistance(Vector::zero()) << endl;
      }
    }

  }

  return result;
}

}
