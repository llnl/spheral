#ifndef Spheral_NodePairListView_hh
#define Spheral_NodePairListView_hh

#include "Neighbor/NodePairIdxType.hh"
#include "Utilities/GPUUtils.hh"

#include <vector>
#include <unordered_map>

namespace Spheral {

class NodePairListView : public chai::CHAICopyable {
#ifdef SPHERAL_UNIFIED_MEMORY
  using ContainerType = SPHERAL_SPAN_TYPE<NodePairIdxType>;
#else
  using ContainerType = typename chai::ManagedArray<NodePairIdxType>;
#endif

public:
  SPHERAL_HOST_DEVICE NodePairListView() = default;
  SPHERAL_HOST_DEVICE virtual ~NodePairListView() = default;
  SPHERAL_HOST NodePairListView(ContainerType const &d) : mData(d) {}

  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) { return mData[i]; }

  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) const { return mData[i]; }

  SPHERAL_HOST_DEVICE
  size_t size() const { return mData.size(); }
  SPHERAL_HOST_DEVICE
  const NodePairIdxType* data() const { return mData.data(); }

  SPHERAL_HOST
  void move(chai::ExecutionSpace space) { GPUUtils::move(mData, space); }

  SPHERAL_HOST
  void touch(chai::ExecutionSpace space) { GPUUtils::touch(mData, space); }

protected:
  ContainerType mData;
};

}

#endif // Spheral_NodePairListView_hh
