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
#ifndef __Spheral_TreeNeighbor__
#define __Spheral_TreeNeighbor__

#include "Neighbor/NeighborSearchType.hh"

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <span>
#include <string>

namespace Spheral {

// Forward declarations
template<typename Dimension> class Neighbor;
template<typename Dimension> class TreeNeighborView;
template<typename Dimension> class GeomPlane;
template<typename Dimension> class DataBase;

template<typename Dimension>
class TreeNeighbor:
    public Neighbor<Dimension>,
    public TreeNeighborView<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using NodeID = typename Neighbor<Dimension>::NodeID;

  using ViewType = TreeNeighborView<Dimension>;

  // Data types we use to build the internal tree structure.
  using LevelKey = typename ViewType::LevelKey;
  using CellKey = typename ViewType::CellKey;

  // Constructors and destructors
  TreeNeighbor(DataBase<Dimension>& dataBase, 
               const NeighborSearchType searchType,
               const double kernelExtent,
               const Vector& xmin,
               const Vector& xmax);
  virtual ~TreeNeighbor() = default;

  // Disabled methods
  TreeNeighbor() = delete;
  TreeNeighbor(const TreeNeighbor&) = delete;
  TreeNeighbor(TreeNeighbor&&) = delete;
  TreeNeighbor& operator=(const TreeNeighbor&) = delete;

  //.............................................................................
  // The required Neighbor interface
  // Force the update of internal data for the NodeList.
  virtual void updateNodes() override;
  virtual void updateNodes(const std::span<NodeID>& nodeIDs) override;

  // Reinitialize based on the desired target h and box dimensions.
  virtual void reinitialize() override;
  virtual void reinitialize(const Vector& xmin,
                            const Vector& xmax,
                            const Scalar htarget) override;
  //.............................................................................

  // Checks for internal validity.
  bool valid() const;

  // Return a dump of the tree structure as a string.
  std::string dumpTree(const bool globalTree) const;

  // Return a string describing the overall statistics of the tree.
  std::string dumpTreeStatistics(const bool globalTree) const;

  // Methods to serialize/deserialize the state of TreeNeighbor.
  void serialize(std::vector<char>& buffer) const;
  void deserialize(std::vector<char>::const_iterator& itr,
                   const std::vector<char>::const_iterator& end);

  // Return the set of occupied cells on each level.
  std::vector<std::vector<CellKey>> occupiedCells() const;

private:
  //--------------------------- Private Interface ---------------------------//
  //----------------------------------------------------------------------------
  // Cell holds the properties of cells in the tree
  // This is the dynamic class used on the GPU for tree building
  //----------------------------------------------------------------------------
  struct Cell {
    CellKey key;                     // Key for this cell
    std::vector<CellKey> daughters;  // Keys of any daughter cells on level+1
    std::vector<NodeID> members;     // Indices (nodeListID, nodeID) of nodes that are members of the cell

    // Constructors
    Cell(): key(0), daughters(), members() {}
    Cell(const CellKey& keyi): key(keyi), daughters(), members() {}
    Cell(const Cell&) = default;
    Cell(const Cell&&) = default;
    Cell& operator=(const Cell&) = default;

    // Throw in comparison operators for help sorting.
    bool operator==(const Cell& rhs) const { return key == rhs.key; }
    bool operator<(const Cell& rhs) const { return key < rhs.key; }
  };

  // Define the types we use to build the tree.
  using TreeLevel = std::unordered_map<CellKey, Cell>;
  using Tree = std::vector<TreeLevel>;

  // Types for flattening a fully constructed Tree for use in the TreeNeighborView
  using FlatCell = typename ViewType::FlatCell;
  using FlatTree = std::vector<FlatCell>;
  using FlatTreeDaughters = std::vector<size_t>;
  using FlatTreeMembers = std::vector<NodeID>;

  // Add a cell key to the daughters of a cell.
  void addDaughter(Cell& cell, const CellKey daughterKey) const;

  // Add a node to the internal tree.
  void addNodeToTree(const Vector& xi,
                     const SymTensor& Hi,
                     const size_t nodeListID,
                     const size_t i,
                     Tree& tree) const;

  // Flatten a fully constructed Tree to FlatTree
  void flattenTreeAndUpdateView();

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
  Tree mTree;
  FlatTree mFlatTree;
  FlatTreeDaughters mFlatTreeDaughters;
  FlatTreeMembers mFlatTreeMembers;
  std::vector<size_t> mFlatTreeLevelOffsets;

  // Stuff we get from base classes
  using Neighbor<Dimension>::mDataBase;
  using Neighbor<Dimension>::mNodeExtent;

  using TreeNeighborView<Dimension>::num1dbits;
  using TreeNeighborView<Dimension>::max1dKey;
  using TreeNeighborView<Dimension>::xkeymask;
  using TreeNeighborView<Dimension>::ykeymask;
  using TreeNeighborView<Dimension>::zkeymask;
  using TreeNeighborView<Dimension>::log2inv;
  using TreeNeighborView<Dimension>::mSearchType;
  using TreeNeighborView<Dimension>::mKernelExtent;
  using TreeNeighborView<Dimension>::mBoxLength;
  using TreeNeighborView<Dimension>::mGridLevelConst0;
  using TreeNeighborView<Dimension>::mXmin;
  using TreeNeighborView<Dimension>::mXmax;
};

}

#include "Neighbor.hh"
#include "Neighbor/TreeNeighborView.hh"

#endif
