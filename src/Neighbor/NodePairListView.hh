#ifndef Spheral_NodePairListView_hh
#define Spheral_NodePairListView_hh

#include "Neighbor/NodePairIdxType.hh"
#include "config.hh"
#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"
#include "Utilities/CHAI_MA_wrapper.hh"

#include <vector>
#include <unordered_map>

namespace Spheral {

class NodePairListView : public chai::CHAICopyable {
  using MAContainer = typename chai::ManagedArray<NodePairIdxType>;

public:
  SPHERAL_HOST_DEVICE NodePairListView() = default;
  SPHERAL_HOST_DEVICE virtual ~NodePairListView() = default;
  SPHERAL_HOST NodePairListView(MAContainer const &d) : mData(d) {}

  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) { return mData[i]; }

  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) const { return mData[i]; }

  SPHERAL_HOST_DEVICE
  size_t size() const { return mData.size(); }
  SPHERAL_HOST_DEVICE
  const NodePairIdxType* data() const { return mData.data(); }

  void move(chai::ExecutionSpace space) { mData.move(space); }

  SPHERAL_HOST
  void touch(chai::ExecutionSpace space) { mData.registerTouch(space); }

protected:
  MAContainer mData;
};

}

#endif // Spheral_NodePairListView_hh
