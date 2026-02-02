//---------------------------------Spheral++----------------------------------//
// Tree -- An implementation of the hierarchical binary tree for Spheral.
//         This instance subdivides space equally on recursively finer levels,
//         so it is a quad-tree in 2D and an oct-tree in 3D.  (This works in 1D
//         as well by subdividing in two's -- is that a binary tree?)
//
// Note: Tree currently hashes coordinates into 64 bit representations, which
// implies a minimum possible resolution as well as specifies that Tree's can 
// only be 64 levels deep.
//
// Created by JMO, 2012-03-12
//----------------------------------------------------------------------------//
#ifndef __Spheral_Tree__
#define __Spheral_Tree__

#include <stdint.h>
#include <bitset>
#include <unordered_map>

#include "NullCellValue.hh"
#include "UniqueNodeLeaf.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Tree
//------------------------------------------------------------------------------
template<typename Dimension, 
         typename CellValue = EmptyCell,
         typename LeafPolicy = UniqueNodeLeaf>
class Tree {
public:
  //---------------------------- Public Interface ----------------------------//
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef uint64_t CoordHash;
  typedef unsigned LevelKey;
  typedef std::bitset<Dimension::nDim*64> CellKey;
  typedef std::pair<LevelKey, CellKey> Key;

  //----------------------------------------------------------------------------
  // Cell holds the properties of cells in the tree.
  //----------------------------------------------------------------------------
  struct Cell {
    CellKey key;                     // Key for this cell.
    std::vector<CellKey> daughters;  // Keys of any daughter cells on level+1
    std::vector<Cell*> daughterPtrs; // Pointers to the daughter cells.
    std::vector<Vector> positions;   // Positions of the input points that terminate in this cell.
    CellValue value;                 // The value associated with this cell.

    // Constructors.
    Cell(): 
      key(0), daughters(), daughterPtrs(), positions(), value() {}

    Cell(const CellKey& keyi, const Vector& xi, const CellValue& valuei):
      key(keyi), daughters(), daughterPtrs(), positions(1, xi), value(valuei) {}

    Cell(const CellKey& keyi, const CellKey& daughter, const Vector& xi, const CellValue& valuei):
      key(keyi), daughters(1, daughter), daughterPtrs(), positions(1, xi), value(valuei) {}

    // Throw in comparison operators for help sorting.
    bool operator==(const Cell& rhs) const { return key == rhs.key; }
    bool operator<(const Cell& rhs) const { return key < rhs.key; }
  };

  //----------------------------------------------------------------------------
  // Construct from an iterator pair over positions.
  //----------------------------------------------------------------------------
  template<typename Iterator,
           typename CellValueFactory = NullCellFactory>
  Tree(const Iterator positionBegin,
       const Iterator positionEnd,
       CellValueFactory factory = CellValueFactory());

  //----------------------------------------------------------------------------
  // Construct from a FieldList of positions.
  //----------------------------------------------------------------------------
  template<typename CellValueFactory = NullCellFactory>
  Tree(const FieldList<Dimension, Vector>& positions,
       CellValueFactory factory = CellValueFactory());

  //----------------------------------------------------------------------------
  // The number of 1d bits (also levels) the tree can support.
  //----------------------------------------------------------------------------
  static unsigned num1dbits() { return 64U; }

  //--------------------------- Private Interface ----------------------------//
  // Define the types we use to build the tree.
  typedef std::unordered_map<CellKey, Cell> TreeLevel;
  typedef std::vector<TreeLevel> Tree;

  // Private data.
  Vector mXmin, mXmax;
  Tree mTree;

  // Build a cell key for a position on a level.
  void buildCellKey(const LevelKey ilevel, const Vector& xi, 
                    CellKey& key, CoordHash& ix, CoordHash& ix, CoordHash& ix) const;

  // Extract the individual indices corresponding to a cell key.
  void extractCellIndices(const CellKey& key, CoordHash& ix, CoordHash& iy, CoordHash& iz) const;

  // Add a daughter key to a cell if not present.
  void addDaughter(Cell& cell, const CellKey& daughterKey);
};

}

#include "TreeInline.hh"

#endif


