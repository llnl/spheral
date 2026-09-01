//---------------------------------Spheral++----------------------------------//
// TreeNeighborView -- non-owning, device-capturable view of a flattened
// TreeNeighbor tree.
//
// TreeNeighbor remains responsible for owning and populating the arrays.
//----------------------------------------------------------------------------//
#ifndef __Spheral_TreeNeighborView_hh__
#define __Spheral_TreeNeighborView_hh__

#include "Threading/GPUUtils.hh"

#include <cstdint>
#include <algorithm>
#include <cmath>
#include <vector>

namespace Spheral {

template<typename Dimension> class TreeNeighbor;

// One flattened tree cell.  Daughter indices address the global cell array;
// member offsets address the flattened node-ID array.
struct TreeNeighborCellRecord {
  using CellKey = std::uint64_t;
  using Index = std::uint32_t;

  CellKey key = 0u;
  Index memberOffset = 0u;
  Index memberCount = 0u;
  Index daughterOffset = 0u;
  Index daughterCount = 0u;
};

template<typename Dimension>
class TreeNeighborView final: public chai::CHAICopyable {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using CellRecord = TreeNeighborCellRecord;
  using CellKey = typename CellRecord::CellKey;
  using Index = typename CellRecord::Index;
  using LevelKey = std::uint32_t;
  using NodeID = std::int32_t;
#ifdef SPHERAL_UNIFIED_MEMORY
  template<typename Value> using SpanType = SPHERAL_SPAN_TYPE<Value>;
#else
  template<typename Value> using SpanType = chai::ManagedArray<Value>;
#endif
  using CellSpan = SpanType<CellRecord>;
  using MemberSpan = SpanType<NodeID>;
  using DaughterSpan = SpanType<Index>;

  SPHERAL_HOST_DEVICE TreeNeighborView() = default;
  SPHERAL_HOST_DEVICE TreeNeighborView(const TreeNeighborView&) = default;
  SPHERAL_HOST_DEVICE TreeNeighborView(TreeNeighborView&&) = default;
  SPHERAL_HOST_DEVICE ~TreeNeighborView() = default;
  SPHERAL_HOST_DEVICE TreeNeighborView& operator=(const TreeNeighborView&) = default;
  SPHERAL_HOST_DEVICE TreeNeighborView& operator=(TreeNeighborView&&) = default;

  SPHERAL_HOST_DEVICE const CellRecord& cell(const Index i) const { return mCells[i]; }
  SPHERAL_HOST_DEVICE NodeID member(const Index i) const { return mMembers[i]; }
  SPHERAL_HOST_DEVICE Index daughter(const Index i) const { return mDaughters[i]; }

  SPHERAL_HOST_DEVICE Index numCells() const { return mCells.size(); }
  SPHERAL_HOST_DEVICE Index numMembers() const { return mMembers.size(); }
  SPHERAL_HOST_DEVICE Index numDaughters() const { return mDaughters.size(); }
  SPHERAL_HOST_DEVICE LevelKey numLevels() const { return mNumLevels; }
  SPHERAL_HOST_DEVICE bool empty() const { return mCells.size() == 0u; }

  // Walk the candidate members for one master query.  The fixed-size stack
  // keeps the traversal device-safe while retaining the host tree's
  // level-by-level range semantics.  The callback is invoked on the device;
  // candidates never need to be materialized on the host.
  template<typename Callback>
  SPHERAL_HOST_DEVICE
  void forEachCandidate(const Vector& position,
                        const typename Dimension::SymTensor& H,
                        Callback callback) const {
    const auto h = 1.0/H.eigenValues().minElement();
    const auto masterLevel = this->gridLevel(h);
    CellKey ixMaster, iyMaster, izMaster;
    this->buildCellIndices(masterLevel, position,
                           ixMaster, iyMaster, izMaster);

    if (this->empty() or this->numLevels() == 0u) return;

    // A tree cell has at most 2^nDim daughters.  A depth of 64 is a generous
    // bound for the current 21-bit tree and leaves room for malformed input
    // to be rejected by the host projection builder.
    constexpr Index maxStack = Index((1u << Dimension::nDim) * 64u);
    Index stackCells[maxStack];
    LevelKey stackLevels[maxStack];
    Index stackSize = 0u;

    constexpr Index root = 0u;
    const auto& rootCell = this->cell(root);
    CHECK(rootCell.daughterCount <= maxStack);
    for (Index offset = 0u; offset < rootCell.daughterCount; ++offset) {
      stackCells[stackSize] = this->daughter(rootCell.daughterOffset + offset);
      stackLevels[stackSize] = 1u;
      ++stackSize;
    }

    while (stackSize > 0u) {
      --stackSize;
      const auto cellIndex = stackCells[stackSize];
      const auto level = stackLevels[stackSize];
      CHECK(level < this->numLevels());

      const CellKey delta = (level <= masterLevel ?
                             CellKey(1u) : (CellKey(1u) << (level - masterLevel)));
      const CellKey ix = this->shiftKeyLevel(ixMaster, masterLevel, level);
      const CellKey iy = this->shiftKeyLevel(iyMaster, masterLevel, level);
      const CellKey iz = this->shiftKeyLevel(izMaster, masterLevel, level);
      const CellKey ixMin = (ix > delta ? ix - delta : CellKey(0u));
      const CellKey iyMin = (iy > delta ? iy - delta : CellKey(0u));
      const CellKey izMin = (iz > delta ? iz - delta : CellKey(0u));
      const CellKey ixMax = ((max1dKey - ix) > delta ?
                             ix + 2u*delta - 1u : max1dKey);
      const CellKey iyMax = ((max1dKey - iy) > delta ?
                             iy + 2u*delta - 1u : max1dKey);
      const CellKey izMax = ((max1dKey - iz) > delta ?
                             iz + 2u*delta - 1u : max1dKey);
      const auto& cellRecord = this->cell(cellIndex);

      if (this->keyInRange(cellRecord.key, ixMin, iyMin, izMin,
                           ixMax, iyMax, izMax)) {
        for (Index offset = 0u; offset < cellRecord.memberCount; ++offset) {
          callback(this->member(cellRecord.memberOffset + offset));
        }
        if (level + 1u < this->numLevels()) {
          CHECK(cellRecord.daughterCount <= maxStack - stackSize);
          for (Index offset = 0u; offset < cellRecord.daughterCount; ++offset) {
            stackCells[stackSize] = this->daughter(cellRecord.daughterOffset + offset);
            stackLevels[stackSize] = level + 1u;
            ++stackSize;
          }
        }
      }
    }
  }

  SPHERAL_HOST void move(chai::ExecutionSpace space) {
    GPUUtils::move(mCells, space);
    GPUUtils::move(mMembers, space);
    GPUUtils::move(mDaughters, space);
  }

private:
  friend class TreeNeighbor<Dimension>;

  SPHERAL_HOST
  void initialize(std::vector<CellRecord>& cells,
                  std::vector<NodeID>& members,
                  std::vector<Index>& daughters,
                  const LevelKey numLevels,
                  const Vector& xmin,
                  const Scalar boxLength,
                  const Scalar gridLevelConst0) {
    GPUUtils::initMAView(mCells, cells);
    GPUUtils::initMAView(mMembers, members);
    GPUUtils::initMAView(mDaughters, daughters);
    mNumLevels = numLevels;
    mXmin = xmin;
    mBoxLength = boxLength;
    mGridLevelConst0 = gridLevelConst0;
  }

  SPHERAL_HOST
  void release() {
    GPUUtils::freeMAView(mCells);
    GPUUtils::freeMAView(mMembers);
    GPUUtils::freeMAView(mDaughters);
  }

  SPHERAL_HOST_DEVICE
  LevelKey gridLevel(const Scalar h) const {
    const auto raw = int(mGridLevelConst0 - std::log(h)/std::log(2.0));
    return LevelKey(raw < 0 ? 0 : (raw >= int(num1dbits) ? num1dbits - 1u : raw));
  }

  SPHERAL_HOST_DEVICE
  void buildCellIndices(const LevelKey level,
                        const Vector& position,
                        CellKey& ix,
                        CellKey& iy,
                        CellKey& iz) const {
    const CellKey ncell = CellKey(1u) << level;
    const CellKey maxcell = ncell - 1u;
    ix = std::min(maxcell, CellKey(std::max(0.0, std::min(1.0, (position.x() - mXmin.x())/mBoxLength)) * ncell));
    iy = std::min(maxcell, CellKey(std::max(0.0, std::min(1.0, (position.y() - mXmin.y())/mBoxLength)) * ncell));
    iz = std::min(maxcell, CellKey(std::max(0.0, std::min(1.0, (position.z() - mXmin.z())/mBoxLength)) * ncell));
  }

  SPHERAL_HOST_DEVICE
  CellKey shiftKeyLevel(const CellKey ix,
                        const LevelKey level0,
                        const LevelKey level1) const {
    return (level1 <= level0 ? ix >> (level0 - level1) : ix << (level1 - level0));
  }

  SPHERAL_HOST_DEVICE
  bool keyInRange(const CellKey key,
                  const CellKey ixMin,
                  const CellKey iyMin,
                  const CellKey izMin,
                  const CellKey ixMax,
                  const CellKey iyMax,
                  const CellKey izMax) const {
    const auto ix = key & xkeymask;
    const auto iy = (key & ykeymask) >> num1dbits;
    const auto iz = (key & zkeymask) >> 2u*num1dbits;
    return (ix >= ixMin and ix <= ixMax and
            iy >= iyMin and iy <= iyMax and
            iz >= izMin and iz <= izMax);
  }

  static constexpr unsigned num1dbits = 21u;
  static constexpr CellKey max1dKey = CellKey(1u) << num1dbits;
  static constexpr CellKey xkeymask = max1dKey - 1u;
  static constexpr CellKey ykeymask = xkeymask << num1dbits;
  static constexpr CellKey zkeymask = ykeymask << num1dbits;

  CellSpan mCells;
  MemberSpan mMembers;
  DaughterSpan mDaughters;
  LevelKey mNumLevels = 0u;
  Vector mXmin;
  Scalar mBoxLength = 0.0;
  Scalar mGridLevelConst0 = 0.0;
};

} // namespace Spheral

#endif
