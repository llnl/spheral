//---------------------------------Spheral++----------------------------------//
// NestedGridDistributedBoundary -- Implementation of the Distributed Boundary
// condition for use with NestedGridNeighbor based NodeLists.
//
// Created by JMO, Mon Aug 27 21:57:51 PDT 2001
//----------------------------------------------------------------------------//

#ifndef NestedGridDistributedBoundary_HH
#define NestedGridDistributedBoundary_HH

#include "DistributedBoundary.hh"

#include <vector>
#include <set>
#include <map>

namespace Spheral {
  template<typename Dimension> class DataBase;
  template<typename Dimension> class NestedGridNeighbor;
  template<typename Dimension> class GridCellIndex;
  template<typename Dimension> class NodeList;
}

namespace Spheral {

template<typename Dimension>
class NestedGridDistributedBoundary: public DistributedBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename DistributedBoundary<Dimension>::DomainBoundaryNodes DomainBoundaryNodes;

  // This method returns the singleton instance of the object.
  static NestedGridDistributedBoundary& instance();

  // Destructor.
  virtual ~NestedGridDistributedBoundary();

  //**********************************************************************
  // Apply the boundary condition to the given Field.
  // Set the ghost nodes based on the NodeLists in the given DataBase.
  virtual void setAllGhostNodes(DataBase<Dimension>& dataBase) override;
  virtual void reset(const DataBase<Dimension>& dataBase) override;
  //**********************************************************************

  // The last set of occupied grid cells (only the master process knows the full
  // set).
  const std::vector< std::vector< std::vector< GridCellIndex<Dimension> > > >& occupiedGridCells() const;

  // Determine the max number of occupied grid levels for all NodeLists in a
  // DataBase.
  int maxNumGridLevels(const DataBase<Dimension>& dataBase) const;

  // Determine the radius of influence in gridcells being used.
  int setGridCellInfluenceRadius(DataBase<Dimension>& dataBase,
                                 const int newGridCellInfluenceRadius) const;

  // Build the set of occupied grid cells for all NodeLists on this process.
  void flattenOccupiedGridCells(const DataBase<Dimension>& dataBase,
                                std::vector< std::vector< GridCellIndex<Dimension> > >& gridCells) const;

  // Pack the occupied grid cell set into a C style array syntax for messaging
  // with MPI.
  void packGridCellIndices(const std::vector< std::vector< GridCellIndex<Dimension> > >& gridCellSet,
                           std::vector<int>& packedGridCellIndices) const;

  // Unpack the occupied grid cell set from C style array syntax to the more
  // sensible set of occupied grid cells.
  void unpackGridCellIndices(const std::vector<int>& packedGridCellIndices,
                             const std::vector<int>& gridCellDimension,
                             std::vector< std::vector< GridCellIndex<Dimension> > >& gridCellSet) const;

  // Access the flag determing if we're applying box culling.
  bool boxCulling() const;
  void boxCulling(bool x);

  // The grid cell influence radius to use for communicating.
  int gridCellInfluenceRadius() const;
  void gridCellInfluenceRadius(int x);

private:
  //--------------------------- Private Interface ---------------------------//
  // List of the occupied grid cells, in the order [process][gridLevel][gridCell].
  std::vector< std::vector< std::vector< GridCellIndex<Dimension> > > > mOccupiedGridCells;

  // The grid cell influence radius we use for building ghost nodes.
  int mGridCellInfluenceRadius;

  // Flag indicating whether we're using the box culling algorithm for reducing the
  // number of ghost nodes.
  bool mBoxCulling;

  // Disabled methods.
  NestedGridDistributedBoundary();
  NestedGridDistributedBoundary(const NestedGridDistributedBoundary&);
  NestedGridDistributedBoundary& operator=(const NestedGridDistributedBoundary&);

  // Private methods.
  void distributeOccupiedGridCells();
  void buildSendNodes(const DataBase<Dimension>& dataBase);
};

}

#endif
