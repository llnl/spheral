//---------------------------------Spheral++----------------------------------//
// RedistributeNodes -- An abstract base class for methods that repartition
// the Spheral++ NodeLists among domains.
//
// Created by JMO, Tue Feb  4 14:23:11 PST 2003
//----------------------------------------------------------------------------//

#ifndef RedistributeNodes_HH
#define RedistributeNodes_HH

#include <string>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class NodeList;
template<typename Dimension> class Boundary;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> struct DomainNode;

template<typename Dimension>
class RedistributeNodes {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors
  RedistributeNodes();

  // Destructor
  virtual ~RedistributeNodes();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  // This is the method required of all descendent classes.
  virtual void redistributeNodes(DataBase<Dimension>& dataBase,
                                 std::vector<Boundary<Dimension>*> boundaries = std::vector<Boundary<Dimension>*>()) = 0;

  // Get this domain ID.
  int domainID() const;

  // Total number of domains.
  int numDomains() const;

  // Global number of nodes.
  int numGlobalNodes(const DataBase<Dimension>& dataBase) const;

  // Calculate the current domain decomposition, and return it as a set of 
  // DomainNode identifiers.
  std::vector<DomainNode<Dimension>> currentDomainDecomposition(const DataBase<Dimension>& dataBase,
                                                                const FieldList<Dimension, size_t>& globalNodeIDs) const;

  // Same as above, but fills in work field in the DomainNodes.
  std::vector<DomainNode<Dimension>> currentDomainDecomposition(const DataBase<Dimension>& dataBase,
                                                                const FieldList<Dimension, size_t>& globalNodeIDs,
                                                                const FieldList<Dimension, Scalar>& workPerNode) const;

  // Given a desired domain decomposition (as a vector<DomainNode>), reassign
  // nodes appropriately.
  void enforceDomainDecomposition(const std::vector<DomainNode<Dimension> >& nodeDistribution,
                                  DataBase<Dimension>& dataBase) const;

  // Test that the given domain decomposition is valid (all nodes accounted 
  // for, once and only once, etc.).
  bool validDomainDecomposition(const std::vector<DomainNode<Dimension> >& nodeDistribution,
                                const DataBase<Dimension>& dataBase) const;

  // Compute the work per node.
  FieldList<Dimension, Scalar> workPerNode(const DataBase<Dimension>& dataBase,
                                           const double Hextent) const;

  // Gather up and print the statistics of the current domain distribution based on
  // a computed work field.
  std::string gatherDomainDistributionStatistics(const FieldList<Dimension, Scalar>& work) const;

  // Flag controlling how we set the work field.
  bool computeWork() const;
  void computeWork(bool x);

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mComputeWork;

  // Pack/unpack a vector<DomainNode> as a vector<double>, for use with MPI.
  std::vector<char> packDomainNodes(const std::vector<DomainNode<Dimension> >& distribution) const;
  std::vector<DomainNode<Dimension> > unpackDomainNodes(const std::vector<char>& buf) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // No copy or assignment operations.
  RedistributeNodes(const RedistributeNodes& nodes);
  RedistributeNodes& operator=(const RedistributeNodes& rhs);
};

}

#include "RedistributeNodesInline.hh"

#endif
