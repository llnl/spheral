//---------------------------------Spheral++----------------------------------//
// VoronoiRedistributeNodes
//
// This algorithm uses the Voronoi tessellation to decide how to domain 
// decompose our points.  The idea is to relax a set of generator points into
// the SPH node distribution -- the generators are attracted to the SPH points
// repelled by one and other.  These generator points then become the seeds to
// draw the Voronoi tessellation about, each cell of which then represents a 
// computational domain.
//
// Created by JMO, Fri Jan 15 09:56:56 PST 2010
//----------------------------------------------------------------------------//
#ifndef VoronoiRedistributeNodes_HH
#define VoronoiRedistributeNodes_HH

#include "RedistributeNodes.hh"
#include "Utilities/KeyTraits.hh"

#include <vector>
#include <map>

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class NodeList;
template<typename Dimension> class Boundary;

template<typename Dimension>
class VoronoiRedistributeNodes: public RedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef KeyTraits::Key Key;

  // Constructors
  VoronoiRedistributeNodes(const double dummy,
                           const bool workBalance,
                           const bool balanceGenerators,
                           const double tolerance,
                           const unsigned maxIterations);

  // Destructor
  virtual ~VoronoiRedistributeNodes();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  // This is the method required of all descendent classes.
  virtual void redistributeNodes(DataBase<Dimension>& dataBase,
                                 std::vector<Boundary<Dimension>*> boundaries =
                                 std::vector<Boundary<Dimension>*>()) override;

  // Given the set of DomainNodes with their domain assignments, compute the 
  // (work weighted) domain centroids.
  void computeCentroids(const std::vector<DomainNode<Dimension> >& nodes,
                        std::vector<Vector>& generators) const;

  // Assign the nodes to the given generator positions, simultaneously computing the total
  // generator work load.
  void assignNodesToGenerators(const std::vector<Vector>& generators,
                               const std::vector<int>& generatorFlags,
                               std::vector<double>& generatorWork,
                               std::vector<DomainNode<Dimension> >& nodes,
                               double& minWork,
                               double& maxWork,
                               unsigned& minNodes,
                               unsigned& maxNodes) const;

  // Look for any generator that has too much work, and unassign it's most 
  // distant nodes.
  void cullGeneratorNodesByWork(const std::vector<Vector>& generators,
                                const std::vector<double>& generatorWork,
                                const double targetWork,
                                std::vector<int>& generatorFlags,
                                std::vector<DomainNode<Dimension> >& nodes) const;

  // Find the adjacent generators in the Voronoi diagram.
  std::vector<size_t> findNeighborGenerators(const size_t igen,
                                             const std::vector<Vector>& generators) const;

  // Flag for whether we should compute the work per node or strictly balance by
  // node count.
  bool workBalance() const;
  void workBalance(bool val);

  // Should we try to work balance between generators?
  // node count.
  bool balanceGenerators() const;
  void balanceGenerators(bool val);

  // The tolerance we're using to check for convergence of the generators.
  double tolerance() const;
  void tolerance(double val);

  // The maximum number of iterations to try and converge the generator positions.
  unsigned maxIterations() const;
  void maxIterations(unsigned val);

private:
  //--------------------------- Private Interface ---------------------------//
  bool mWorkBalance, mBalanceGenerators;
  double mTolerance;
  unsigned mMaxIterations;

  // No copy or assignment operations.
  VoronoiRedistributeNodes(const VoronoiRedistributeNodes& nodes);
  VoronoiRedistributeNodes& operator=(const VoronoiRedistributeNodes& rhs);
};

}

#endif
