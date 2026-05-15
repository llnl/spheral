#ifndef Spheral_ConnectivityMapFlatView_hh
#define Spheral_ConnectivityMapFlatView_hh

#include "Utilities/GPUUtils.hh"

namespace Spheral {

class ConnectivityMapFlatView : public chai::CHAICopyable {
#ifdef SPHERAL_UNIFIED_MEMORY
  using SpanType = SPHERAL_SPAN_TYPE<int>;
#else
  using SpanType = typename chai::ManagedArray<int>;
#endif

public:
  SPHERAL_HOST_DEVICE ConnectivityMapFlatView() = default;
  SPHERAL_HOST_DEVICE virtual ~ConnectivityMapFlatView() = default;
  SPHERAL_HOST ConnectivityMapFlatView(const size_t numNodeLists,
                                       SpanType const& offsets,
                                       SpanType const& neighbors):
    mNumNodeLists(numNodeLists),
    mOffsets(offsets),
    mNeighbors(neighbors) {}

  SPHERAL_HOST_DEVICE size_t numNodeLists() const { return mNumNodeLists; }
  SPHERAL_HOST_DEVICE size_t numEntries() const { return (mOffsets.size() > 0u ? mOffsets.size() - 1u : 0u); }
  SPHERAL_HOST_DEVICE size_t numNodes() const { return (mNumNodeLists > 0u ? this->numEntries()/mNumNodeLists : 0u); }
  SPHERAL_HOST_DEVICE bool empty() const { return mOffsets.size() <= 1u; }

  SPHERAL_HOST_DEVICE size_t entryIndex(const size_t globalNodeIndex,
                                        const size_t neighborNodeList) const {
    return globalNodeIndex*mNumNodeLists + neighborNodeList;
  }

  SPHERAL_HOST_DEVICE size_t beginOffset(const size_t entryIndex) const { return mOffsets[entryIndex]; }
  SPHERAL_HOST_DEVICE size_t endOffset(const size_t entryIndex) const { return mOffsets[entryIndex + 1u]; }

  SPHERAL_HOST_DEVICE size_t size(const size_t entryIndex) const {
    return this->endOffset(entryIndex) - this->beginOffset(entryIndex);
  }

  SPHERAL_HOST_DEVICE size_t size(const size_t globalNodeIndex,
                                  const size_t neighborNodeList) const {
    return this->size(this->entryIndex(globalNodeIndex, neighborNodeList));
  }

  SPHERAL_HOST_DEVICE int operator()(const size_t entryIndex,
                                     const size_t localIndex) const {
    return mNeighbors[mOffsets[entryIndex] + localIndex];
  }

  SPHERAL_HOST_DEVICE int operator()(const size_t globalNodeIndex,
                                     const size_t neighborNodeList,
                                     const size_t localIndex) const {
    return (*this)(this->entryIndex(globalNodeIndex, neighborNodeList), localIndex);
  }

  SPHERAL_HOST_DEVICE const int* offsetsData() const { return mOffsets.data(); }
  SPHERAL_HOST_DEVICE const int* neighborsData() const { return mNeighbors.data(); }

  SPHERAL_HOST void move(chai::ExecutionSpace space) {
    GPUUtils::move(mOffsets, space);
    GPUUtils::move(mNeighbors, space);
  }

  SPHERAL_HOST void touch(chai::ExecutionSpace space) {
    GPUUtils::touch(mOffsets, space);
    GPUUtils::touch(mNeighbors, space);
  }

private:
  size_t mNumNodeLists = 0u;
  SpanType mOffsets;
  SpanType mNeighbors;
};

}

#endif // Spheral_ConnectivityMapFlatView_hh
