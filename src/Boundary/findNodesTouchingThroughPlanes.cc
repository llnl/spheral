//------------------------------------------------------------------------------
// Find the set of nodes that see through a pair of planes.
//------------------------------------------------------------------------------
#include "findNodesTouchingThroughPlanes.hh"
#include "Distributed/allReduce.hh"

namespace Spheral {
  
template<typename Dimension>
std::vector<size_t>
findNodesTouchingThroughPlanes(const NodeList<Dimension>& nodeList,
                               const GeomPlane<Dimension>& enterPlane,
                               const GeomPlane<Dimension>& exitPlane,
                               const double kernelExtent,
                               const double hmultiplier) {
  std::vector<size_t> result;

  // Find the maximum smoothing scale of any point touching either plane.
  const auto  etaMax = hmultiplier * kernelExtent;
  const auto  n = nodeList.numNodes();
  const auto& pos = nodeList.positions();
  const auto& H = nodeList.Hfield();
  auto hmax = 0.0;
  for (auto i = 0u; i < n; ++i) {
    const auto& ri = pos(i);
    const auto& Hi = H(i);
    const auto  hmaxi = 1.0/Hi.eigenValues().minElement();
    if (hmaxi > hmax and std::min(exitPlane.minimumDistance(ri), enterPlane.minimumDistance(ri)) < etaMax*hmaxi) hmax = hmaxi;
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
      if (disti >= 0.0 and disti <= etaMax) result.push_back(i);
      // cerr << " --> " << i << " " << ri << " " << enterPlanePrime.minimumDistance(Vector::zero()) << " " << exitPlanePrime.minimumDistance(Vector::zero()) << endl;
    }
  }

  return result;
}

}
