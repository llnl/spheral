#ifndef Spheral_TreeNeighborView_hh
#define Spheral_TreeNeighborView_hh

#include "Utilities/GPUUtils.hh"

#include <stdint.h>

namespace Spheral {

template<typename Dimension>
class TreeNeighborView : public chai::CHAICopyable {
public:
  using LevelKey = uint32_t;
  using CellKey = uint64_t;
#ifdef SPHERAL_UNIFIED_MEMORY
  using IntSpanType = SPHERAL_SPAN_TYPE<int>;
  using KeySpanType = SPHERAL_SPAN_TYPE<CellKey>;
#else
  using IntSpanType = typename chai::ManagedArray<int>;
  using KeySpanType = typename chai::ManagedArray<CellKey>;
#endif

  SPHERAL_HOST_DEVICE TreeNeighborView() = default;
  SPHERAL_HOST_DEVICE virtual ~TreeNeighborView() = default;
  SPHERAL_HOST TreeNeighborView(const size_t numLevels,
                                const IntSpanType& levelOffsets,
                                const KeySpanType& cellKeys,
                                const IntSpanType& daughterOffsets,
                                const IntSpanType& daughterIndices,
                                const IntSpanType& memberOffsets,
                                const IntSpanType& members):
    mNumLevels(numLevels),
    mLevelOffsets(levelOffsets),
    mCellKeys(cellKeys),
    mDaughterOffsets(daughterOffsets),
    mDaughterIndices(daughterIndices),
    mMemberOffsets(memberOffsets),
    mMembers(members) {}

  SPHERAL_HOST_DEVICE size_t numLevels() const { return mNumLevels; }
  SPHERAL_HOST_DEVICE size_t numCells() const { return mCellKeys.size(); }
  SPHERAL_HOST_DEVICE bool empty() const { return mCellKeys.size() == 0u; }

  SPHERAL_HOST_DEVICE size_t levelBegin(const size_t level) const { return mLevelOffsets[level]; }
  SPHERAL_HOST_DEVICE size_t levelEnd(const size_t level) const { return mLevelOffsets[level + 1u]; }
  SPHERAL_HOST_DEVICE size_t levelSize(const size_t level) const { return this->levelEnd(level) - this->levelBegin(level); }

  SPHERAL_HOST_DEVICE CellKey cellKey(const size_t cellIndex) const { return mCellKeys[cellIndex]; }

  SPHERAL_HOST_DEVICE size_t daughterBegin(const size_t cellIndex) const { return mDaughterOffsets[cellIndex]; }
  SPHERAL_HOST_DEVICE size_t daughterEnd(const size_t cellIndex) const { return mDaughterOffsets[cellIndex + 1u]; }
  SPHERAL_HOST_DEVICE size_t daughterSize(const size_t cellIndex) const { return this->daughterEnd(cellIndex) - this->daughterBegin(cellIndex); }
  SPHERAL_HOST_DEVICE int daughterIndex(const size_t cellIndex, const size_t daughterIndex) const {
    return mDaughterIndices[mDaughterOffsets[cellIndex] + daughterIndex];
  }

  SPHERAL_HOST_DEVICE size_t memberBegin(const size_t cellIndex) const { return mMemberOffsets[cellIndex]; }
  SPHERAL_HOST_DEVICE size_t memberEnd(const size_t cellIndex) const { return mMemberOffsets[cellIndex + 1u]; }
  SPHERAL_HOST_DEVICE size_t memberSize(const size_t cellIndex) const { return this->memberEnd(cellIndex) - this->memberBegin(cellIndex); }
  SPHERAL_HOST_DEVICE int member(const size_t cellIndex, const size_t memberIndex) const {
    return mMembers[mMemberOffsets[cellIndex] + memberIndex];
  }

  SPHERAL_HOST void move(chai::ExecutionSpace space) {
    GPUUtils::move(mLevelOffsets, space);
    GPUUtils::move(mCellKeys, space);
    GPUUtils::move(mDaughterOffsets, space);
    GPUUtils::move(mDaughterIndices, space);
    GPUUtils::move(mMemberOffsets, space);
    GPUUtils::move(mMembers, space);
  }

  SPHERAL_HOST void touch(chai::ExecutionSpace space) {
    GPUUtils::touch(mLevelOffsets, space);
    GPUUtils::touch(mCellKeys, space);
    GPUUtils::touch(mDaughterOffsets, space);
    GPUUtils::touch(mDaughterIndices, space);
    GPUUtils::touch(mMemberOffsets, space);
    GPUUtils::touch(mMembers, space);
  }

  SPHERAL_HOST_DEVICE void shallowCopy(TreeNeighborView const& other) const {
#ifdef SPHERAL_UNIFIED_MEMORY
    const_cast<size_t&>(mNumLevels) = other.mNumLevels;
    const_cast<IntSpanType&>(mLevelOffsets) = other.mLevelOffsets;
    const_cast<KeySpanType&>(mCellKeys) = other.mCellKeys;
    const_cast<IntSpanType&>(mDaughterOffsets) = other.mDaughterOffsets;
    const_cast<IntSpanType&>(mDaughterIndices) = other.mDaughterIndices;
    const_cast<IntSpanType&>(mMemberOffsets) = other.mMemberOffsets;
    const_cast<IntSpanType&>(mMembers) = other.mMembers;
#else
    const_cast<size_t&>(mNumLevels) = other.mNumLevels;
    mLevelOffsets.shallowCopy(other.mLevelOffsets);
    mCellKeys.shallowCopy(other.mCellKeys);
    mDaughterOffsets.shallowCopy(other.mDaughterOffsets);
    mDaughterIndices.shallowCopy(other.mDaughterIndices);
    mMemberOffsets.shallowCopy(other.mMemberOffsets);
    mMembers.shallowCopy(other.mMembers);
#endif
  }

private:
  size_t mNumLevels = 0u;
  IntSpanType mLevelOffsets;
  KeySpanType mCellKeys;
  IntSpanType mDaughterOffsets;
  IntSpanType mDaughterIndices;
  IntSpanType mMemberOffsets;
  IntSpanType mMembers;
};

}

#endif // Spheral_TreeNeighborView_hh
