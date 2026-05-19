#ifndef Spheral_NestedGridNeighborView_hh
#define Spheral_NestedGridNeighborView_hh

#include "Utilities/GPUUtils.hh"

#include <stdint.h>

namespace Spheral {

template<typename Dimension>
class NestedGridNeighborView : public chai::CHAICopyable {
public:
  using CellKey = uint64_t;
#ifdef SPHERAL_UNIFIED_MEMORY
  using IntSpanType = SPHERAL_SPAN_TYPE<int>;
  using KeySpanType = SPHERAL_SPAN_TYPE<CellKey>;
#else
  using IntSpanType = typename chai::ManagedArray<int>;
  using KeySpanType = typename chai::ManagedArray<CellKey>;
#endif

  SPHERAL_HOST_DEVICE NestedGridNeighborView() = default;
  SPHERAL_HOST_DEVICE virtual ~NestedGridNeighborView() = default;
  SPHERAL_HOST NestedGridNeighborView(const size_t numLevels,
                                      const IntSpanType& levelOffsets,
                                      const KeySpanType& cellKeys,
                                      const IntSpanType& memberOffsets,
                                      const IntSpanType& members):
    mNumLevels(numLevels),
    mLevelOffsets(levelOffsets),
    mCellKeys(cellKeys),
    mMemberOffsets(memberOffsets),
    mMembers(members) {}

  SPHERAL_HOST_DEVICE size_t numLevels() const { return mNumLevels; }
  SPHERAL_HOST_DEVICE bool empty() const { return mCellKeys.size() == 0u; }

  SPHERAL_HOST_DEVICE size_t levelBegin(const size_t level) const { return mLevelOffsets[level]; }
  SPHERAL_HOST_DEVICE size_t levelEnd(const size_t level) const { return mLevelOffsets[level + 1u]; }
  SPHERAL_HOST_DEVICE size_t levelSize(const size_t level) const { return this->levelEnd(level) - this->levelBegin(level); }

  SPHERAL_HOST_DEVICE CellKey cellKey(const size_t cellIndex) const { return mCellKeys[cellIndex]; }

  SPHERAL_HOST_DEVICE size_t memberBegin(const size_t cellIndex) const { return mMemberOffsets[cellIndex]; }
  SPHERAL_HOST_DEVICE size_t memberEnd(const size_t cellIndex) const { return mMemberOffsets[cellIndex + 1u]; }
  SPHERAL_HOST_DEVICE size_t memberSize(const size_t cellIndex) const { return this->memberEnd(cellIndex) - this->memberBegin(cellIndex); }
  SPHERAL_HOST_DEVICE int member(const size_t cellIndex, const size_t memberIndex) const {
    return mMembers[mMemberOffsets[cellIndex] + memberIndex];
  }

  SPHERAL_HOST void move(chai::ExecutionSpace space) {
    GPUUtils::move(mLevelOffsets, space);
    GPUUtils::move(mCellKeys, space);
    GPUUtils::move(mMemberOffsets, space);
    GPUUtils::move(mMembers, space);
  }

  SPHERAL_HOST void touch(chai::ExecutionSpace space) {
    GPUUtils::touch(mLevelOffsets, space);
    GPUUtils::touch(mCellKeys, space);
    GPUUtils::touch(mMemberOffsets, space);
    GPUUtils::touch(mMembers, space);
  }

  SPHERAL_HOST_DEVICE void shallowCopy(NestedGridNeighborView const& other) const {
#ifdef SPHERAL_UNIFIED_MEMORY
    const_cast<size_t&>(mNumLevels) = other.mNumLevels;
    const_cast<IntSpanType&>(mLevelOffsets) = other.mLevelOffsets;
    const_cast<KeySpanType&>(mCellKeys) = other.mCellKeys;
    const_cast<IntSpanType&>(mMemberOffsets) = other.mMemberOffsets;
    const_cast<IntSpanType&>(mMembers) = other.mMembers;
#else
    const_cast<size_t&>(mNumLevels) = other.mNumLevels;
    mLevelOffsets.shallowCopy(other.mLevelOffsets);
    mCellKeys.shallowCopy(other.mCellKeys);
    mMemberOffsets.shallowCopy(other.mMemberOffsets);
    mMembers.shallowCopy(other.mMembers);
#endif
  }

private:
  size_t mNumLevels = 0u;
  IntSpanType mLevelOffsets;
  KeySpanType mCellKeys;
  IntSpanType mMemberOffsets;
  IntSpanType mMembers;
};

}

#endif // Spheral_NestedGridNeighborView_hh
