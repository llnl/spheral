//---------------------------------Spheral++----------------------------------//
// TreeNeighbor
//
// An SPH neighbor finder based on an oct-tree like structure with specialized
// cell membership criteria.
//
// Based on the algorithm originally described in 
// Owen, Villumsen, Shapiro, & Martel 1998, ApJS, 116, 155
//
// Created by J. Michael Owen, Fri Mar 23 15:50:35 PDT 2012
//----------------------------------------------------------------------------//
#include "TreeNeighbor.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/packElement.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/FastMath.hh"
#include "Boundary/mapPositionThroughPlanes.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <array>

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the vertex coordinates for a cell.
//------------------------------------------------------------------------------
// 1D
inline
std::array<Dim<1>::Vector, 2>
findCellVertices(const Dim<1>::Vector& xmin,
                 const double& boxLength,
                 const uint32_t& ilevel,
                 const uint64_t& ix,
                 const uint64_t& /*iy*/,
                 const uint64_t& /*iz*/) {
  using Vector = Dim<1>::Vector;
  const auto cellSize = boxLength/(1U << ilevel);
  return std::array<Vector, 2> result({xmin + Vector(ix       *cellSize),
                                       xmin + Vector((ix + 1U)*cellSize)});
}

// 2D
inline
std::array<Dim<2>::Vector, 4>
findCellVertices(const Dim<2>::Vector& xmin,
                 const double& boxLength,
                 const uint32_t& ilevel,
                 const uint64_t& ix,
                 const uint64_t& iy,
                 const uint64_t& /*iz*/) {
  using Vector = Dim<2>::Vector;
  const auto cellSize = boxLength/(1U << ilevel);
  return std::array<Vector, 4>({xmin + Vector(ix       *cellSize, iy       *cellSize),
                                xmin + Vector((ix + 1U)*cellSize, iy       *cellSize),
                                xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize),
                                xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize)});
}

// 3D
inline
std::array<Dim<3>::Vector, 8>
findCellVertices(const Dim<3>::Vector& xmin,
                 const double& boxLength,
                 const uint32_t& ilevel,
                 const uint64_t& ix,
                 const uint64_t& iy,
                 const uint64_t& iz) {
  using Vector = Dim<3>::Vector;
  const auto cellSize = boxLength/(1U << ilevel);
  return std::array<Vector, 8>({xmin + Vector(ix       *cellSize, iy       *cellSize, iz       *cellSize),
                                xmin + Vector((ix + 1U)*cellSize, iy       *cellSize, iz       *cellSize),
                                xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize, iz       *cellSize),
                                xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize, iz       *cellSize),
                                xmin + Vector(ix       *cellSize, iy       *cellSize, (iz + 1U)*cellSize),
                                xmin + Vector((ix + 1U)*cellSize, iy       *cellSize, (iz + 1U)*cellSize),
                                xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize, (iz + 1U)*cellSize),
                                xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize, (iz + 1U)*cellSize)});
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeNeighbor<Dimension>::
TreeNeighbor(DataBase<Dimension>& dataBase,
             const NeighborSearchType searchType,
             const double kernelExtent,
             const Vector& xmin,
             const Vector& xmax):
  Neighbor<Dimension>(dataBase, searchType, kernelExtent),
  mBoxLength((xmax - xmin).maxElement()),
  mGridLevelConst0(0.0),
  mXmin(xmin),
  mXmax(xmax),
  mTree() {
  this->reinitialize(xmin, xmax, (xmax - xmin).maxElement()/4.0);
}

//------------------------------------------------------------------------------
// Update the internal data for all nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeNeighbor<Dimension>::
updateNodes() {

  // Clear our internal data.
  mTree.clear();
  mTree.reserve(num1dbits); // This is necessary to avoid memory errors!

  // Grab the state.
  const auto numNodeLists = mDataBase.numNodeLists();
  const auto positions = mDataBase.globalPositions();
  const auto H = mDataBase.globalHfield();
  CHECK(positions.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);

  // Walk all the points and assign them to Cells
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = positions[k]->numElements();
    for (auto i = 0u; i < n; ++i) {
      this->addNodeToTree(positions(k,i), H(k,i), k, i, mTree);
    }
  }

  // Sort each Cells info
  for (auto klevel = 0u; klevel < mTree.size(); ++klevel) {
    for (auto& [key, cellt]: mTree[klevel]) {
      std::sort(cellt.daughters.begin(), cellt.daughters.end());
      std::sort(cellt.members.begin(), cellt.members.end());
    }
  }

  // Force the node extents to be calculated.
  this->setNodeExtents();

  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Update the internal data for a subset of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeNeighbor<Dimension>::
updateNodes(const std::vector<NodeID>&) {
  this->updateNodes();  // Punt and just rebuild everything.
}

//------------------------------------------------------------------------------
// reinitialize
// For TreeNeighbor we want to ensure we have a grid level a close fit to the
// target kernelExtent*htarget.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
reinitialize() {
  const auto etaMax = this->kernelExtent();
  mBoxLength = (mXmax - mXmin).maxElement();
  mGridLevelConst0 = log(mBoxLength/etaMax)/log(2.0);
  mTree.clear();
}

template<typename Dimension>
void
TreeNeighbor<Dimension>::
reinitialize(const typename Dimension::Vector& xmin,
             const typename Dimension::Vector& xmax,
             const Scalar /*htarget*/) {
  const auto etaMax = this->kernelExtent();
  mXmin = xmin;
  mXmax = xmax;
  mBoxLength = (mXmax - mXmin).maxElement();
  mGridLevelConst0 = log(mBoxLength/etaMax)/log(2.0);
  mTree.clear();

  // // Now optimize the box-size so we have a level near the target size.
  // const auto lvl = this->gridLevel(htarget) + 1U;
  // mBoxLength = (1U << lvl) * etaMax*htarget;
  // mGridLevelConst0 = log(mBoxLength/etaMax)/log(2.0);
  // cerr << lvl << " " << mBoxLength << endl;
}

//------------------------------------------------------------------------------
// Internal consistency checks.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TreeNeighbor<Dimension>::
valid() const {
  bool result = true;

  // Make sure each node is listed only once.
  unordered_map<NodeID, size_t> nodeCount;
  for (const auto& level: mTree) {
    for (const auto& [key, cell]: level) {
      for (const auto nid: cell.members) {
        auto itr = nodeCount.find(nid);
        if (itr == nodeCount.end()) {
          nodeCount[nid] = 1u;
        } else {
          ++(itr->second);
        }
      }
    }
  }

  for (const auto [i, count]: nodeCount) {
    if (count != 1u) {
      SpheralWarning << "TreeNeighor::valid failing test of nodes uniquely assigned to cell: " << i << " " << count << std::endl;
      result = false;
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// dumpTree
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeNeighbor<Dimension>::
dumpTree(const bool globalTree) const {
  return this->dumpTree(mTree, globalTree);
}

//------------------------------------------------------------------------------
// dumpTree
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeNeighbor<Dimension>::
dumpTree(const Tree& tree,
         const bool globalTree) const {
  std::stringstream ss;
  CellKey ix, iy, iz;
  unsigned nlevels = tree.size();
  if (globalTree) nlevels = allReduce(nlevels, SPHERAL_OP_MAX);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (auto ilevel = 0u; ilevel < nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(tree[ilevel].size());
    if (ilevel < tree.size()) {
      for (const auto& [key, cell]: tree[ilevel]) {
        cells.push_back(cell);
        this->serialize(cell, localBuffer);
      }
    }
#ifdef SPHERAL_ENABLE_MPI
    const auto numProcs = Process::getTotalNumberOfProcesses();
    const auto rank = Process::getRank();
    if (globalTree) {
      for (auto sendProc = 0u; sendProc < numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0u) {
          std::vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
          if (rank != sendProc) {
            std::vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    std::sort(cells.begin(), cells.end());
    cells.erase(std::unique(cells.begin(), cells.end()), cells.end());

    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : boxsize = " << mBoxLength/(1 << ilevel) << " : numCells = " << cells.size() << "\n";
    for (const auto& cell: cells) {
      extractCellIndices(cell.key, ix, iy, iz);
      ss << "    Cell key=" << cell.key << " : (ix,iy,iz)=(" << ix << " " << iy << " " << iz << ")\n"
         << "         daughters = ( ";
      for (const auto daughter: cell.daughters) ss << *ditr << " ";
      ss << ")\n"
         << "         nodes = [";
      for (const auto& nid: cell.members) ss << "(" << cell.members[k].first << " " cell.members.second << ") ";
      ss <<" ]\n";
    }
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// dumpTreeStatistics
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeNeighbor<Dimension>::
dumpTreeStatistics(const bool globalTree) const {
  return this->dumpTreeStatistics(mTree, globalTree);
}

//------------------------------------------------------------------------------
// dumpTreeStatistics
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeNeighbor<Dimension>::
dumpTreeStatistics(const Tree& tree,
                   const bool globalTree) const {
  std::stringstream ss;
  auto nlevels = tree.size();
  if (globalTree) nlevels = allReduce(nlevels, SPHERAL_OP_MAX);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (auto ilevel = 0u; ilevel < nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(tree[ilevel].size());
    if (ilevel < tree.size()) {
      for (const auto& [key, cell]: tree[ilevel]) {
        cells.push_back(cell);
        this->serialize(cell, localBuffer);
      }
    }
#ifdef SPHERAL_ENABLE_MPI
    const auto numProcs = Process::getTotalNumberOfProcesses();
    const auto rank = Process::getRank();
    if (globalTree) {
      for (auto sendProc = 0u; sendProc < numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
          if (rank != sendProc) {
            vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    std::sort(cells.begin(), cells.end());
    cells.erase(std::unique(cells.begin(), cells.end()), cells.end());

    // Count up how much of everything we have.
    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : boxsize = " << mBoxLength/(1 << ilevel) << " : numCells = " << cells.size() << "\n";
    size_t nparticles = 0u;
    for (const auto& cell: cells) nparticles += cell.members.size();
    ss << "         : nparts = " << nparticles << "\n";
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// Serialize to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
serialize(std::vector<char>& buffer) const {
  packElement(mBoxLength, buffer);
  packElement(mGridLevelConst0, buffer);
  packElement(mXmin, buffer);
  packElement(mXmax, buffer);
  const unsigned nlevels = mTree.size();
  packElement(nlevels, buffer);
  for (auto ilevel = 0u; ilevel < nlevels; ++ilevel) {
    const auto ncells = mTree[ilevel].size();
    packElement(ncells, buffer);
    for (const auto& [key, cell]: mTree[ilevel]) {
      packElement(key, buffer);
      serialize(cell, buffer);
    }
  }
}

//------------------------------------------------------------------------------
// Serialize a cell to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
serialize(const TreeNeighbor<Dimension>::Cell& cell,
          std::vector<char>& buffer) const {
  packElement(cell.key, buffer);
  packElement(cell.daughters, buffer);
  packElement(cell.members, buffer);
}

//------------------------------------------------------------------------------
// Deserialize a tree from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
deserialize(vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) {
  unpackElement(mBoxLength, bufItr, endItr);
  unpackElement(mGridLevelConst0, bufItr, endItr);
  unpackElement(mXmin, bufItr, endItr);
  unpackElement(mXmax, bufItr, endItr);
  unsigned nlevels;
  CellKey key;
  Cell cell;
  unpackElement(nlevels, bufItr, endItr);
  mTree.resize(nlevels);
  for (auto ilevel = 0u; ilevel < nlevels; ++ilevel) {
    unsigned ncells;
    unpackElement(ncells, bufItr, endItr);
    for (auto i = 0u; i < ncells; ++i) {
      cell.daughters.clear();
      cell.members.clear();
      unpackElement(key, bufItr, endItr);
      deserialize(cell, bufItr, endItr);
      mTree[ilevel][key] = cell;
    }
  }
}

//------------------------------------------------------------------------------
// Deserialize a cell from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
deserialize(typename TreeNeighbor<Dimension>::Cell& cell,
            vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) const {
  unpackElement(cell.key, bufItr, endItr);
  unpackElement(cell.daughters, bufItr, endItr);
  unpackElement(cell.members, bufItr, endItr);
}

//------------------------------------------------------------------------------
// Return the set of occupied cells.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<vector<typename TreeNeighbor<Dimension>::CellKey>>
TreeNeighbor<Dimension>::
occupiedCells() const {
  const auto numLevels = mTree.size();
  vector<vector<CellKey>> result(numLevels);
  for (auto ilevel = 0u; ilevel < numLevels; ++ilevel) {
    for (const auto& keyCell: mTree[ilevel]) {
      if (not keyCell.second.members.empty()) result[ilevel].push_back(keyCell.first);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Do any existing points occupy cells overlapping the given position and
// extent?
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TreeNeighbor<Dimension>::
occupied(const Vector& xi,
         const SymTensor& Hi) const {
  if (mTree.size() == 0) return false;
  const auto gridLeveli = std::min(gridLevel(Hi), unsigned(mTree.size()) - 1u);
  CellKey keyi, ix, iy, iz;
  for (auto ilevel = 0u; ilevel < gridLeveli; ++ilevel) {
    buildCellKey(ilevel, xi, keyi, ix, iy, iz);
    const auto keyCellItr = mTree[ilevel].find(keyi);
    if (keyCellItr != mTree[ilevel].end() and
        (not keyCellItr->second.members.empty())) return true;
  }
  buildCellKey(gridLeveli, xi, keyi, ix, iy, iz);
  return mTree[gridLeveli].find(keyi) != mTree[gridLeveli].end();
}

template<typename Dimension>
bool
TreeNeighbor<Dimension>::
occupied(const Vector& xi,
         const double hi) const {
  return this->occupied(xi, SymTensor::one() * 1.0/hi);
}

//------------------------------------------------------------------------------
// Add a daughter to a cell if not present.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
addDaughter(Cell& cell,
            const CellKey daughterKey) const {
  if (std::find(cell.daughters.begin(), cell.daughters.end(), daughterKey) == cell.daughters.end())
    cell.daughters.push_back(daughterKey);
  ENSURE(cell.daughters.size() <= (1U << Dimension::nDim));
}

//------------------------------------------------------------------------------
// Add a node to the internal Tree structure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
addNodeToTree(const Vector& xi,
              const SymTensor& Hi,
              const size_t nodeListID,
              const size_t i,
              Tree& tree) const {

  // Determine the level for this point.
  const LevelKey homeLevel = this->gridLevel(Hi);

  LevelKey ilevel = 0u;
  CellKey key, parentKey, ix, iy, iz;
  typename TreeLevel::iterator itr;

  // First walk the tree for all levels above our native level.
  while (ilevel <= homeLevel) {

    // Do we need to add another level to the tree?
    if (ilevel == tree.size()) tree.push_back(TreeLevel());

    // Create the key for the cell containing this particle on this level.
    buildCellKey(ilevel, xi, key, ix, iy, iz);
    itr = tree[ilevel].find(key);

    // Is this a new cell?
    if (itr == tree[ilevel].end()) {
      tree[ilevel][key] = Cell(key);
      itr = tree[ilevel].find(key);
    }

    // Link this cell as a daughter of its parent.
    if (ilevel > 0) {
      CHECK(tree[ilevel - 1].find(parentKey) != tree[ilevel - 1].end());
      addDaughter(tree[ilevel - 1][parentKey], key);
    }

    // Prepare for the next level.
    parentKey = key;
    ++ilevel;
  }

  // Add the node to it's final cell.
  itr->second.members.push_back(std::make_pair(nodeListID, i));
}

//------------------------------------------------------------------------------
// Flatten the Tree to a FlatTree for use in the View class
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
flattenTreeAndUpdateView() {

  // Clear the old data
  mFlatTree.clear();
  mFlatTreeDaughters.clear();
  mFlatTreeMembers.clear();
  mFlatTreeLevalOffsets.clear();

  // Walk the tree
  for (const auto& level: mTree) {
    mFlatTreeLevelOffsets.push_back(mFlatTreeLevelOffsets.size());

    // Walk the cells on this level
    for (const auto& [key, cell]: level) {
      const auto numDaughters = cell.daughters.size();
      const auto numMembers = cell.members.size();
      mFlatTree.emplace_back(FlatCell(key,
                                      mFlatTreeDaughters.size(),
                                      cell.daughters.size(),
                                      mFlatTreeMembers.size(),
                                      cell.members.size()));
      std::ranges::copy(cell.daughters, std::back_inserter(mFlatTreeDaughters));
      std::ranges::copy(cell.members, std::back_inserter(mFlatTreeMembers));
    }
  }
  mFlatTreeLevelOffsets.push_back(mFlatTreeLevelOffsets.size());
  CHECK(mFlatTreeLevelOffsets.size() == mTree.size() + 1u);

  // Update the TreeNeighborView data to reference the new FlatTree
#ifdef SPHERAL_UNIFIED_MEMORY
  ViewType::mFlatTreeView = mFlatTree;
  ViewType::mFlatTreeDaughtersIndicesView = mFlatTreeDaughtersIndices;
  ViewType::mFlatTreeMembersView = mFlatTreeMembers;
#else
  GPUUtils::initMAView(ViewType::mFlatTreeView, mFlatTree);
  GPUUtils::initMAView(ViewType::mFlatTreeDaughtersIndicesView, mFlatTreeDaughtersIndices);
  GPUUtils::initMAView(ViewType::mFlatTreeMembersView, mFlatTreeMembers);
  ViewType::mFlatTreeView.registerTouch(chai::CPU);
  ViewType::mFlatTreeDaughtersIndicesView.registerTouch(chai::CPU);
  ViewType::mFlatTreeMembersView.registerTouch(chai::CPU);
#endif
}

}
