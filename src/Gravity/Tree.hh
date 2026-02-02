//---------------------------------Spheral++----------------------------------//
// Tree -- a implementation of a quad/oct (2D/3D) tree structure.
//
// Extracted from our original implementation in TreeGravity.
// This class carries along a few data members from that heritage (like the mass
// and velocity) which we retain for convenience. Those attributes can be ignored
// for purely geometrical applications.
//
// Created by JMO, Tue Oct  4 10:17:41 PDT 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_Tree__
#define __Spheral_Tree__

#include "Geometry/Dimension.hh"

#include <stdint.h>
#include <unordered_map>
#include <unordered_set>

namespace Spheral {

template<typename Dimension>
class Tree {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef uint64_t CellKey;

  //----------------------------------------------------------------------------
  // Cell holds the properties of cells in the tree.
  //----------------------------------------------------------------------------
  struct Cell {
    double M, Mglobal;               // total mass (and global sum)
    Vector xcm;                      // center of mass
    Vector vcm;                      // velocity of center of mass
    double rcm2cc2;                  // square of the distance between center of mass and geometric center
    CellKey key;                     // Key for this cell.
    std::vector<CellKey> daughters;  // Keys of any daughter cells on level+1
    std::vector<Cell*> daughterPtrs; // Pointers to the daughter cells.
    std::vector<double> masses;      // Masses of the nodes that terminate in this cell.
    std::vector<Vector> positions;   // Positions of the nodes that terminate in this cell.
    std::vector<Vector> velocities;  // Velocities of the nodes that terminate in this cell.

    // Convenience constructors for Tree::addNodeToTree.
    Cell(): M(0.0), Mglobal(0.0), xcm(), vcm(), rcm2cc2(0.0), key(0), daughters(), daughterPtrs(), masses(), positions(), velocities() {}
    Cell(const double mi, const Vector& xi, const Vector& vi, const CellKey& keyi):
      M(mi), Mglobal(mi), xcm(xi), vcm(vi), rcm2cc2(0.0), key(keyi), daughters(), daughterPtrs(), masses(1, mi), positions(1, xi), velocities(1, vi) {}
    Cell(const double mi, const Vector& xi, const Vector& vi, const CellKey& keyi, const CellKey& daughter):
      M(mi), Mglobal(mi), xcm(xi), vcm(vi), rcm2cc2(0.0), key(keyi), daughters(1, daughter), daughterPtrs(), masses(), positions(), velocities() {}

    // Throw in comparison operators for help sorting.
    bool operator==(const Cell& rhs) const { return key == rhs.key; }
    bool operator<(const Cell& rhs) const { return key < rhs.key; }
  };

  // Data types we use to build the internal tree structure.
  typedef uint32_t LevelKey;
  typedef std::pair<size_t, size_t> NodeID;
  typedef std::unordered_map<NodeID, std::vector<std::unordered_set<CellKey> > > CompletedCellSet;
  typedef std::unordered_map<CellKey, Cell> TreeLevel;

  static unsigned num1dbits;                   // The number of bits we quantize 1D coordinates to.  We have to fit three of these in 64 bits.
  static CellKey max1dKey;                     // The maximum number of cells this corresponds to in a direction.
  static CellKey xkeymask, ykeymask, zkeymask; // Bit masks we can use to extract the coordinate specific indices from a cell key.

  //! Constructor.
  Tree(const Vector& xmin, const Vector& xmax);

  //! Destructor.
  ~Tree();

  //! Coordinate boundaries of the tree
  const Vector& xmin() const;
  const Vector& xmax() const;

  //! Center of a tree cell
  Vector cellCenter(const LevelKey& level, const CellKey& key) const;

  //! Lower bound coordinate in cell
  Vector lowerBound(const LevelKey& level, const CellKey& key) const;

  //! Upper bound coordinate in cell
  Vector upperBound(const LevelKey& level, const CellKey& key) const;

  //! Cell size on the given level
  double cellSize(const LevelKey& level) const;

  //! Number of levels in the Tree
  size_t numLevels() const;
  void numLevels(const size_t n);

  //! Access a tree level by indexing
  TreeLevel& operator[](const size_t i);
  const TreeLevel& operator[](const size_t i) const;

  //! Build a cell key based on a level and position.
  void buildCellKey(const LevelKey ilevel,
                    const Vector& xi,
                    CellKey& key,
                    CellKey& ix,
                    CellKey& iy,
                    CellKey& iz) const;

  //! Extract the individual coordinate indices from a cell key.
  void extractCellIndices(const CellKey& key,
                          CellKey& ix,
                          CellKey& iy,
                          CellKey& iz) const;

  //! Add a cell key to the daughters of a cell.
  void addDaughter(Cell& cell, const CellKey daughterKey) const;

  //! Add a node to the internal tree.
  void addNodeToTree(const double mi,      // Including mass and velocity
                     const Vector& xi,
                     const Vector& vi);
  void addNodeToTree(const Vector& xi);    // Coordinate only

  //! Return a dump of the tree structure as a string.
  std::string dumpTree(const bool globalTree) const;

  //! Return a string describing the overall statistics of the tree.
  std::string dumpTreeStatistics(const bool globalTree) const;

  //! Methods to help serializing/deserializing Trees to buffers of char.
  void serialize(std::vector<char>& buffer) const;
  void serialize(const Cell& cell, std::vector<char>& buffer) const;

  //! Unpack a tree from a buffer.
  void deserialize(std::vector<char>::const_iterator& bufItr, const std::vector<char>::const_iterator& endItr);
  void deserialize(Cell& cell, std::vector<char>::const_iterator& bufItr, const std::vector<char>::const_iterator& endItr) const;

private:
  // Private data.
  double mBoxLength;
  Vector mXmin, mXmax;
  std::vector<TreeLevel> mLevels;

  // Construct all the daughterPtrs in a tree.
  void constructDaughterPtrs();

  // No default constructor, copying, or assignment
  Tree();
  Tree(const Tree& rhs);
  Tree& operator=(const Tree& rhs);
};

//------------------------------------------------------------------------------
// Define our static members.
//------------------------------------------------------------------------------
template<typename Dimension> unsigned Tree<Dimension>::num1dbits = 21U;
template<typename Dimension> uint64_t Tree<Dimension>::max1dKey = 1U << Tree<Dimension>::num1dbits;
template<typename Dimension> uint64_t Tree<Dimension>::xkeymask = (1U << Tree<Dimension>::num1dbits) - 1U;
template<typename Dimension> uint64_t Tree<Dimension>::ykeymask = Tree<Dimension>::xkeymask << Tree<Dimension>::num1dbits;
template<typename Dimension> uint64_t Tree<Dimension>::zkeymask = Tree<Dimension>::ykeymask << Tree<Dimension>::num1dbits;

}

#include "TreeInline.hh"

#endif
