//---------------------------------Spheral++----------------------------------//
// ConnectivityMap_RAJA -- RAJA pair-construction backend for ConnectivityMap.
//
// ConnectivityMap retains the established CPU implementation and all common
// pair finalization.  This subclass overrides only the optional accelerated
// pair-construction hook.
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_ConnectivityMap_RAJA_hh_
#define _Spheral_NeighborSpace_ConnectivityMap_RAJA_hh_

#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

template<typename Dimension>
class ConnectivityMap_RAJA: public ConnectivityMap<Dimension> {
public:
  ConnectivityMap_RAJA() = default;
  ~ConnectivityMap_RAJA() override = default;

protected:
  bool tryBuildNodePairs(const double kernelExtent,
                         const bool ghostConnectivity,
                         std::vector<NodePairIdxType>& nodePairs) override;
};

}

#endif
