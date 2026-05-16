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
#ifndef __Spheral_TreeNeighbor_hh__
#define __Spheral_TreeNeighbor_hh__

#include "Neighbor.hh"
#include "TreeNeighborView.hh"

#include <stdint.h>
#include <unordered_map>
#include <vector>

namespace Spheral {

template<typename Dimension>
class TreeNeighbor: public Neighbor<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef std::vector<int>::iterator iterator;
  typedef std::vector<int>::const_iterator const_iterator;

  // Data types we use to build the internal tree structure.
  typedef uint32_t LevelKey;
  typedef uint64_t CellKey;

  // Constructors and destructors
  TreeNeighbor(NodeList<Dimension>& nodeList, 
               const NeighborSearchType searchType,
               const double kernelExtent,
               const Vector& xmin,
               const Vector& xmax);
  virtual ~TreeNeighbor();

  // ********** Descendent Neighbor types must provide these methods. **********
  // Set or refine the neighbor lists for the given position and smoothing 
  // scale.
  virtual void setMasterList(int nodeID,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const override;
  virtual void setMasterList(const Vector& position,
                             const Scalar& H,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const override;
  virtual void setMasterList(const Vector& position,
                             const SymTensor& H,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const override;
  virtual size_t countMasterList(const Vector& position,
                                 const SymTensor& H,
                                 const bool ghostConnectivity = false) const override;
  virtual void fillMasterList(const Vector& position,
                              const SymTensor& H,
                              int* result,
                              const bool ghostConnectivity = false) const override;
  virtual size_t countCoarseNeighbors(const Vector& position,
                                      const SymTensor& H,
                                      const bool ghostConnectivity = false) const override;
  virtual void fillCoarseNeighbors(const Vector& position,
                                   const SymTensor& H,
                                   int* result,
                                   const bool ghostConnectivity = false) const override;
  virtual NeighborGroupDescriptor groupDescriptor(const Vector& position,
                                                  const SymTensor& H,
                                                  const uint64_t fallbackToken) const override;

  virtual void setRefineNeighborList(const Vector& position,
                                     const Scalar& H,
                                     const std::vector<int>& coarseNeighbors,
                                     std::vector<int>& refineNeighbors) const override;
  virtual void setRefineNeighborList(const Vector& position,
                                     const SymTensor& H,
                                     const std::vector<int>& coarseNeighbors,
                                     std::vector<int>& refineNeighbors) const override;

  // Set Neighbors for the given position.
  virtual void setMasterList(const Vector& position,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const override;
  virtual void setRefineNeighborList(const Vector& position,
                                     const std::vector<int>& coarseNeighbors,
                                     std::vector<int>& refineNeighbors) const override;

  // Set the neighbor lists based on proximity to planes.
  virtual void setMasterList(const GeomPlane<Dimension>& enterPlane,
                             const GeomPlane<Dimension>& exitPlane,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors) const override;

  // Force the update of internal data for the NodeList.
  virtual void updateNodes() override;
  virtual void updateNodes(const std::vector<int>& nodeIDs) override;
  //****************************************************************************

  // Reinitialize based on the desired target h and box dimensions.
  virtual void reinitialize() override;
  virtual void reinitialize(const Vector& xmin, const Vector& xmax, const Scalar htarget) override;

  // Checks for internal validity.
  virtual bool valid() const override;

  // Compute the grid level appropriate for the given smoothing scale.
  unsigned gridLevel(const double& h) const;        // units of length
  unsigned gridLevel(const SymTensor& H) const;     // units of 1/length

  // Return a dump of the tree structure as a string.
  std::string dumpTree(const bool globalTree) const;

  // Return a string describing the overall statistics of the tree.
  std::string dumpTreeStatistics(const bool globalTree) const;

  // Access internal attributes.
  const Vector& xmin() const;
  const Vector& xmax() const;
  double boxLength() const;

  // Cell size on the given level
  double cellSize(const LevelKey levelID) const;
  TreeNeighborView<Dimension> view() const {
    return TreeNeighborView<Dimension>(mTree.size(),
                                       mTreeLevelOffsetsSpan,
                                       mTreeCellKeysSpan,
                                       mTreeDaughterOffsetsSpan,
                                       mTreeDaughterIndicesSpan,
                                       mTreeMemberOffsetsSpan,
                                       mTreeMembersSpan);
  }

  // Methods to serialize/deserialize the state of TreeNeighbor.
  void serialize(std::vector<char>& buffer) const;
  void deserialize(std::vector<char>::const_iterator& itr,
                   const std::vector<char>::const_iterator& end);

  // Return the set of occupied cells on each level.
  std::vector<std::vector<CellKey>> occupiedCells() const;

  // Find the nearest cell center on the given grid level.
  Vector nearestCellCenter(const Vector& xi, const double hi) const;
  Vector nearestCellCenter(const Vector& xi, const SymTensor& Hi) const;

  // Test if the given position and extent is occupied, in the sense that
  // any tree cell is occupied at any coarser level down to it's native
  // level.
  bool occupied(const Vector& xi, const double hi) const;
  bool occupied(const Vector& xi, const SymTensor& Hi) const;

  // For our parallel algorithm it is useful to be able to set the master/coarse
  // information based on the given (level, cell).
  size_t countTreeMasterList(const Vector& position,
                             const SymTensor& H,
                             const bool ghostConnectivity) const;
  void fillTreeMasterList(const Vector& position,
                          const SymTensor& H,
                          int* result,
                          const bool ghostConnectivity) const;
  size_t countTreeCoarseNeighbors(const Vector& position,
                                  const SymTensor& H) const;
  void fillTreeCoarseNeighbors(const Vector& position,
                               const SymTensor& H,
                               int* result) const;
  size_t countTreeCoarseNeighbors(const LevelKey levelID,
                                  const CellKey cellID) const;
  void fillTreeCoarseNeighbors(const LevelKey levelID,
                               const CellKey cellID,
                               int* result) const;
  void setTreeMasterList(const Vector& position,
                         const SymTensor& H,
                         std::vector<int>& masterList,
                         std::vector<int>& coarseNeighbors,
                         const bool ghostConnectivity) const;
  void setTreeMasterList(const LevelKey levelID,
                         const CellKey cellID,
                         std::vector<int>& masterList,
                         std::vector<int>& coarseNeighbors,
                         const bool ghostConnectivity) const;

  // Expose hidden Neighbor methods
  using Neighbor<Dimension>::setMasterList;
  using Neighbor<Dimension>::setRefineNeighborList;

private:
  //--------------------------- Private Interface ---------------------------//
  static const unsigned num1dbits;                   // The number of bits we quantize 1D coordinates to.  We have to fit three of these in 64 bits.
  static const CellKey max1dKey;                     // The maximum number of cells this corresponds to in a direction.
  static const CellKey xkeymask, ykeymask, zkeymask; // Bit masks we can use to extract the coordinate specific indices from a cell key.

  //----------------------------------------------------------------------------
  // Cell holds the properties of cells in the tree.
  //----------------------------------------------------------------------------
  struct Cell {
    CellKey key;                     // Key for this cell.
    std::vector<CellKey> daughters;  // Keys of any daughter cells on level+1.
    std::vector<Cell*> daughterPtrs; // Pointers to the daughter cells.
    std::vector<int> members;        // Indices of nodes that are members of the cell.

    // Convenience constructors for OctTreeGravity::addNodeToTree.
    Cell(): key(0), daughters(), daughterPtrs(), members() {}
    Cell(const CellKey& keyi):
      key(keyi), daughters(), daughterPtrs(), members() {}

    // Throw in comparison operators for help sorting.
    bool operator==(const Cell& rhs) const { return key == rhs.key; }
    bool operator<(const Cell& rhs) const { return key < rhs.key; }
  };

  // Define the types we use to build the tree.
  typedef std::unordered_map<CellKey, Cell> TreeLevel;
  typedef std::vector<TreeLevel> Tree;

  // Default constructor -- disabled.
  TreeNeighbor();

  // Copy constructor -- disabled.
  TreeNeighbor(const TreeNeighbor&);

  // Assignment operator -- disabled.
  TreeNeighbor& operator=(const TreeNeighbor&);

  // Build a cell key based on a level and position.
  void buildCellKey(const LevelKey ilevel,
                    const Vector& xi,
                    CellKey& key,
                    CellKey& ix,
                    CellKey& iy,
                    CellKey& iz) const;

  // Extract the individual coordinate indices from a cell key.
  void extractCellIndices(const CellKey& key,
                          CellKey& ix,
                          CellKey& iy,
                          CellKey& iz) const;

  // Add a cell key to the daughters of a cell.
  void addDaughter(Cell& cell, const CellKey daughterKey) const;

  // Add a node to the internal tree.
  void addNodeToTree(const Vector& xi,
                     const SymTensor& Hi,
                     const unsigned i,
                     Tree& tree) const;

  // Construct all the daughterPtrs in a tree.
  void constructDaughterPtrs(Tree& tree) const;
  void rebuildView();

  // Actual method for setting the master list.
  void setTreeMasterList(const Vector& position,
                         const double& h,
                         std::vector<int>& masterList,
                         std::vector<int>& coarseNeighbors,
                         const bool ghostConnectivity) const;

  // Actual method for setting the refine list.
  void setTreeRefineNeighborList(const Vector& position,
                                 const SymTensor& H,
                                 const std::vector<int>& coarseNeighbors,
                                 std::vector<int>& refineNeighbors) const;

  // Return all the members in range of the specified cell.
  std::vector<int> findTreeNeighbors(const LevelKey& masterLevel,
                                     const CellKey& ix_master,
                                     const CellKey& iy_master,
                                     const CellKey& iz_master) const;

  // Shift a 1D key from one level to another.
  CellKey shiftKeyLevel(const CellKey& ix, const LevelKey& level0, const LevelKey& level1) const;

  // Determine if a full cell key is in the set of 1D key ranges.
  bool keyInRange(const CellKey& key, 
                  const CellKey& ix_min, const CellKey& iy_min, const CellKey& iz_min,
                  const CellKey& ix_max, const CellKey& iy_max, const CellKey& iz_max) const;

  // Find the minimum distance from a cell/level to a plane.
  double distanceToCell(const LevelKey& ilevel, const CellKey& key, const GeomPlane<Dimension>& plane) const;

  // Map a cell key/level through an entrance/exit plane pair, returning the set of cells
  // that overlap the mapped cell.
  std::vector<CellKey> mapKey(const LevelKey& ilevel, const CellKey& key, 
                              const GeomPlane<Dimension>& enterPlane,
                              const GeomPlane<Dimension>& exitPlane) const;

  // Seralize/deserialize cells.
  void serialize(const Cell& cell, std::vector<char>& buffer) const;
  void deserialize(Cell& cell, 
                   std::vector<char>::const_iterator& bufItr, 
                   const std::vector<char>::const_iterator& endItr) const;

  // Tree dump methods that work on the given Tree.
  std::string dumpTree(const Tree& tree,
                       const bool globalTree) const;
  std::string dumpTreeStatistics(const Tree& tree,
                                 const bool globalTree) const;

  // Private data.
  double mBoxLength, mGridLevelConst0;
  Vector mXmin, mXmax;
  Tree mTree;

  std::vector<int> mTreeLevelOffsets;
  std::vector<CellKey> mTreeCellKeys;
  std::vector<int> mTreeDaughterOffsets;
  std::vector<int> mTreeDaughterIndices;
  std::vector<int> mTreeMemberOffsets;
  std::vector<int> mTreeMembers;
#ifdef SPHERAL_UNIFIED_MEMORY
  using FlatTreeSpan = SPHERAL_SPAN_TYPE<int>;
  using FlatTreeCellKeySpan = SPHERAL_SPAN_TYPE<CellKey>;
#else
  using FlatTreeSpan = chai::ManagedArray<int>;
  using FlatTreeCellKeySpan = chai::ManagedArray<CellKey>;
#endif
  FlatTreeSpan mTreeLevelOffsetsSpan;
  FlatTreeCellKeySpan mTreeCellKeysSpan;
  FlatTreeSpan mTreeDaughterOffsetsSpan;
  FlatTreeSpan mTreeDaughterIndicesSpan;
  FlatTreeSpan mTreeMemberOffsetsSpan;
  FlatTreeSpan mTreeMembersSpan;
};

}

#endif
