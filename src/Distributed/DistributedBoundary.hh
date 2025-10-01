//---------------------------------Spheral++----------------------------------//
// DistributedBoundary -- Base class for distributed parallel boundary
// conditions, connecting NodeLists across parallel domains.
//
// Created by JMO, Thu Aug 23 21:34:32 PDT 2001
//----------------------------------------------------------------------------//

#ifndef DistributedBoundary_HH
#define DistributedBoundary_HH


#ifdef SPHERAL_ENABLE_MPI
#include <mpi.h>
#endif

#include "Boundary/Boundary.hh"

#include <vector>
#include <map>
#include <list>

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class NodeList;
  template<typename Dimension> class FieldBase;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension> class DataBase;
}

namespace Spheral {

template<typename Dimension>
class DistributedBoundary: public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  struct DomainBoundaryNodes {
    std::vector<size_t> sendNodes;
    std::vector<size_t> receiveNodes;
  };

  typedef std::map<int, DomainBoundaryNodes> DomainBoundaryNodeMap;
  typedef std::map<NodeList<Dimension>*, DomainBoundaryNodeMap> NodeListDomainBoundaryNodeMap;

  // Constructors and destructors.
  DistributedBoundary();
  virtual ~DistributedBoundary();

  // Get this domain ID.
  int domainID() const;

  // Total number of domains.
  int numDomains() const;

  // Test if the given NodeList is communicated on this domain or not.
  bool communicatedNodeList(const NodeList<Dimension>& nodeList) const;

  // Test if the given NodeList is communicated with the given domain.
  bool nodeListSharedWithDomain(const NodeList<Dimension>& nodeList,
                                int neighborDomainID) const;

  // Allow read access to the communication information.
  const NodeListDomainBoundaryNodeMap& nodeListDomainBoundaryNodeMap() const;
  const DomainBoundaryNodeMap& domainBoundaryNodeMap(const NodeList<Dimension>& nodeList) const;
  const DomainBoundaryNodes& domainBoundaryNodes(const NodeList<Dimension>&,
                                                 int neighborDomainID) const;

  // Extract the current set of processors we're communicating with.
  void communicatedProcs(std::vector<int>& sendProcs,
			 std::vector<int>& recvProcs) const;

  //****************************************************************************
  // Required Boundary interface
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;       // This one should not be used with DistributedBoundary
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override;

  // Distributed boundaries don't have "violate" nodes, so override these methods to no-ops.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the given Field.
  virtual void applyGhostBoundary(FieldBase<Dimension>& field) const override;

  //****************************************************************************
  // Require descendent Distributed Neighbors to provide the setGhostNodes method for DataBases.
  virtual void setAllGhostNodes(DataBase<Dimension>& dataBase) override = 0;

  // Override the Boundary method for culling ghost nodes.
  virtual void cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
                              FieldList<Dimension, size_t>& old2newIndexMap,
                              std::vector<size_t>& numNodesRemoved) override;

  // Override the base method to finalize ghost boundaries.
  virtual void finalizeGhostBoundary() const override;

  // We do not want to use the parallel ghost nodes as generators.
  virtual bool meshGhostNodes() const override;

  //****************************************************************************
  // Non-blocking exchanges.
  void beginExchangeFieldFixedSize(FieldBase<Dimension>& field) const;
  void beginExchangeFieldVariableSize(FieldBase<Dimension>& field) const;

  // Force the exchanges which have been registered to execute.
  void finalizeExchanges();

  // Unpack a packed set of Field values back into the Field.
  void unpackField(FieldBase<Dimension>& field,
                   const std::list< std::vector<char> >& packedValues) const;

  // Update the control and ghost nodes of the base class.
  void setControlAndGhostNodes();

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Descendent classes get read/write access to the communication maps.
  NodeListDomainBoundaryNodeMap& accessNodeListDomainBoundaryNodeMap();
  DomainBoundaryNodeMap& accessDomainBoundaryNodeMap(const NodeList<Dimension>& nodeList);
  DomainBoundaryNodes& accessDomainBoundaryNodes(const NodeList<Dimension>&,
                                                 int neighborDomainID);

  // Convenience method to return an iterator to the DomainBoundaryNodes for the
  // given NodeList and domain pair.  If there isn't an entry for this pair already,
  // this method creates one.
  DomainBoundaryNodes& openDomainBoundaryNodes(NodeList<Dimension>* nodeListPtr,
                                               const int domainID);

  // Inverse of above -- remove the DomainBoundaryNodes for a NodeList/procID pair.
  void removeDomainBoundaryNodes(NodeList<Dimension>* nodeListPtr,
                                 const int domainID);

  // Override the Boundary method for clearing the maps.
  virtual void reset(const DataBase<Dimension>& dataBase) override;

  // This handy helper method will build receive and ghost nodes on all each
  // domain based on send nodes that have already been filled in.
  void buildReceiveAndGhostNodes(const DataBase<Dimension>& dataBase);

private:
  //--------------------------- Private Interface ---------------------------//
  int mDomainID;
  NodeListDomainBoundaryNodeMap mNodeListDomainBoundaryNodeMap;
  
  // List of the fields that are currently backlogged for exchanging.
  mutable std::vector<FieldBase<Dimension>*> mExchangeFields;

  // Internal tag for MPI communiators.
  mutable int mMPIFieldTag;

#ifdef SPHERAL_ENABLE_MPI
  // Send/receive requests.
  mutable std::vector<MPI_Request> mSendRequests;
  mutable std::vector<MPI_Request> mRecvRequests;
#endif

  // TODO: Determine if this is still a thing?
#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
  mutable std::vector<int> mSendProcIDs;
  mutable std::vector<int> mRecvProcIDs;
#endif

  // Buffers for holding send/receive data.
  typedef std::list<std::list<std::vector<char>>> CommBufferSet;
  mutable CommBufferSet mSendBuffers;
  mutable CommBufferSet mRecvBuffers;

  // Maps linking fields to their communication buffers.
  typedef std::map<const FieldBase<Dimension>*, std::list<std::vector<char>>*> Field2BufferType;
  mutable Field2BufferType mField2SendBuffer;
  mutable Field2BufferType mField2RecvBuffer;

};

}

#include "DistributedBoundaryInline.hh"

#endif
