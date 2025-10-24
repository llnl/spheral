//---------------------------------Spheral++----------------------------------//
// ParmetisRedistributeNodes -- Redistribute nodes by calling the Parmetis
// library to determine the new optimal domain decomposition.
// http://www-users.cs.umn.edu/~karypis/metis/parmetis/
//
// Created by JMO, Wed Feb 12 17:29:09 PST 2003
//----------------------------------------------------------------------------//
#include "RedistributeNodes.hh"
#include "ParmetisRedistributeNodes.hh"
#include "Utilities/DomainNode.hh"
#include "BoundingVolumeDistributedBoundary.hh"
#include "DistributedBoundary.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Communicator.hh"

#include "Utilities/DBC.hh"

#include <algorithm>

#include <sstream>
#include <fstream>
#include <cstdlib>
#include <mpi.h>
namespace Spheral {


//------------------------------------------------------------------------------
// Two ways of comparing the pair<globalID, distance> for sorting.
//------------------------------------------------------------------------------
struct SortNeighborsByID {
  bool operator()(const std::pair<int, double>& lhs,
                  const std::pair<int, double>& rhs) {
    return lhs.first < rhs.first;
  }
};

struct SortNeighborsByDistance {
  bool operator()(const std::pair<int, double>& lhs,
                  const std::pair<int, double>& rhs) {
    return lhs.second < rhs.second;
  }
};

//------------------------------------------------------------------------------
// Construct with the given node extent.
//------------------------------------------------------------------------------
template<typename Dimension>
ParmetisRedistributeNodes<Dimension>::
ParmetisRedistributeNodes(double extent):
  RedistributeNodes<Dimension>(),
  mNormalizedNodeExtent(0.0) {
  setNormalizedNodeExtent(extent);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ParmetisRedistributeNodes<Dimension>::
~ParmetisRedistributeNodes() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Call on the Parmetis library to describe an
// optimal partioning of the nodes, and then apply that partitioning.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ParmetisRedistributeNodes<Dimension>::
redistributeNodes(DataBase<Dimension>& dataBase,
                  vector<Boundary<Dimension>*> boundaries) {

  const int numProcs = this->numDomains();

  // Go over each NodeList, and clear out any ghost nodes.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // If the user did not specify any boundaries, then create a Distributed
  // boundary for local use.
  // Build a BoundingVolumeDistributedBoundary, and construct the ghost nodes.
  BoundingVolumeDistributedBoundary<Dimension>& bound = BoundingVolumeDistributedBoundary<Dimension>::instance();
  if (boundaries.size() == 0) boundaries.push_back(&bound);

  // Use the specified boundaries to create the ghost nodes.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) {
    (*boundItr)->setAllGhostNodes(dataBase);
    (*boundItr)->finalizeGhostBoundary();
    for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) (*nodeListItr)->neighbor().updateNodes();
  }

  // Assign global IDs to all nodes.
  FieldList<Dimension, int> globalIDs = globalNodeIDs(dataBase);

  // Use the boundary conditions to communicate global IDs.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(globalIDs);
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Get the local description of the domain distribution.
  vector<DomainNode<Dimension> > nodeDistribution = currentDomainDecomposition(dataBase, globalIDs);
  const int numLocalNodes = nodeDistribution.size();

  // Build the Parmetis data structures.
  vector<idx_t> xadj;
  vector<idx_t> adjacency;
  vector<idx_t> vtxdist;
  vector<idx_t> vweight;
  vector<float> xyz;
  buildCSRGraph(dataBase, nodeDistribution, globalIDs,
                xadj, adjacency, vtxdist, vweight, xyz);
//   const bool ok = validCSRGraph(nodeDistribution, dataBase, xadj, adjacency, vtxdist);
//   if (!ok) cerr << "Failed valid test for CSR graph." << endl;
//   {
//     stringstream filename;
//     filename << "CSR" << this->domainID();
//     ofstream file(filename.str().c_str());
//     file << "xadj:";
//     for (vector<idx_t>::const_iterator itr = xadj.begin();
//          itr < xadj.end();
//          ++itr) 
//       file  << " " << *itr;
//     file << endl 
//          << "adj:";
//     for (vector<idx_t>::const_iterator itr = adjacency.begin();
//          itr < adjacency.end();
//          ++itr)
//       file  << " " << *itr;
//     file << endl
//          << "vtxdist:";
//     for (vector<idx_t>::const_iterator itr = vtxdist.begin();
//          itr < vtxdist.end();
//          ++itr)
//       file  << " " << *itr;
//     file << endl;
//     file.close();
//   }
//   CHECK(validCSRGraph(nodeDistribution, dataBase, xadj, adjacency, vtxdist));

  // We're done with the ghost nodes now, so eliminate them before we try 
  // rearranging nodes.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // Build the Parmetis tpwgts array, which is the distribution of vertex weights
  // between domains.  Since we want equal weight per domain, this is just 
  // an array of length numDomains, each element set to 1/numDomains.
  CHECK(numProcs > 0);
  vector<float> weightPerDomain(this->numDomains(), 1.0/numProcs);

  // Call Parmetis and get the new partitioning.
  int wgtflag = 2;              // only providing vertex weights
  int numflag = 0;              // C style numbering
  int ndims = Dimension::nDim;  // number of dimensions
  int ncon = 1;                 // number of weights per vertex
  int nparts = this->numDomains();    // number of domains we want
  float ubvec[1] = {1.05};        // imbalance tolerance per vertex weight
  int options[3] = {0, 0, 0};  // optional parameters that can be passed
  int edgecut;                 // Output -- number of edges cut in the parmetis decomp.
  vector<idx_t> part((vector<idx_t>::size_type) numLocalNodes); // Output, the parmetis decoposition.
  MPI_Barrier(Communicator::communicator());
  ParMETIS_V3_PartGeomKway(&(*vtxdist.begin()), 
                           &(*xadj.begin()), 
                           &(*adjacency.begin()),
                           &(*vweight.begin()),
                           NULL,
                           &wgtflag,
                           &numflag,
                           &ndims,
                           &(*xyz.begin()),
                           &ncon,
                           &nparts,
                           &(*weightPerDomain.begin()),
                           ubvec,
                           options,
                           &edgecut,
                           &(*part.begin()),
                           &Communicator::communicator());
  MPI_Barrier(Communicator::communicator());

  // Loop over the domain decomposition, and fill in the new domain assignments.
  CHECK(nodeDistribution.size() == part.size());
  for (int i = 0; i < numLocalNodes; ++i) {
    CHECK(part[i] >= 0 && part[i] < numProcs);
    nodeDistribution[i].domainID = part[i];
  }

  // OK, the localDistribution now holds the desired redistribution of the nodes.
  // Go ahead and redistribute them.
  CHECK(validDomainDecomposition(nodeDistribution, dataBase));
  enforceDomainDecomposition(nodeDistribution, dataBase);
}

//------------------------------------------------------------------------------
// Refine the distribution (assumes the input distribution is pretty good).
//------------------------------------------------------------------------------
template<typename Dimension>
void
ParmetisRedistributeNodes<Dimension>::
refineAndRedistributeNodes(DataBase<Dimension>& dataBase,
                           vector<Boundary<Dimension>*> boundaries) {

  const int numProcs = this->numDomains();

  // Go over each NodeList, and clear out any ghost nodes.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // If the user did not specify any boundaries, then create a Distributed
  // boundary for local use.
  // Build a BoundingVolumeDistributedBoundary, and construct the ghost nodes.
  BoundingVolumeDistributedBoundary<Dimension>& bound = BoundingVolumeDistributedBoundary<Dimension>::instance();
  if (boundaries.size() == 0) boundaries.push_back(&bound);

  // Use the specified boundaries to create the ghost nodes.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) {
    (*boundItr)->setAllGhostNodes(dataBase);
    (*boundItr)->finalizeGhostBoundary();
    for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) (*nodeListItr)->neighbor().updateNodes();
  }

  // Assign global IDs to all nodes.
  FieldList<Dimension, int> globalIDs = globalNodeIDs(dataBase);

  // Use the boundary conditions to communicate global IDs.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(globalIDs);
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Get the local description of the domain distribution.
  vector<DomainNode<Dimension> > nodeDistribution = currentDomainDecomposition(dataBase, globalIDs);
  const int numLocalNodes = nodeDistribution.size();

  // Build the Parmetis data structures.
  vector<idx_t> xadj;
  vector<idx_t> adjacency;
  vector<idx_t> vtxdist;
  vector<idx_t> vweight;
  vector<float> xyz;
  buildCSRGraph(dataBase, nodeDistribution, globalIDs,
                xadj, adjacency, vtxdist, vweight, xyz);
//   const bool ok = validCSRGraph(nodeDistribution, dataBase, xadj, adjacency, vtxdist);
//   if (!ok) cerr << "Failed valid test for CSR graph." << endl;
//   {
//     stringstream filename;
//     filename << "CSR" << this->domainID();
//     ofstream file(filename.str().c_str());
//     file << "xadj:";
//     for (vector<idx_t>::const_iterator itr = xadj.begin();
//          itr < xadj.end();
//          ++itr) 
//       file  << " " << *itr;
//     file << endl 
//          << "adj:";
//     for (vector<idx_t>::const_iterator itr = adjacency.begin();
//          itr < adjacency.end();
//          ++itr)
//       file  << " " << *itr;
//     file << endl
//          << "vtxdist:";
//     for (vector<idx_t>::const_iterator itr = vtxdist.begin();
//          itr < vtxdist.end();
//          ++itr)
//       file  << " " << *itr;
//     file << endl;
//     file.close();
//   }
//   CHECK(validCSRGraph(nodeDistribution, dataBase, xadj, adjacency, vtxdist));

  // We're done with the ghost nodes now, so eliminate them before we try 
  // rearranging nodes.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // Build the Parmetis tpwgts array, which is the distribution of vertex weights
  // between domains.  Since we want equal weight per domain, this is just 
  // an array of length numDomains, each element set to 1/numDomains.
  CHECK(numProcs > 0);
  vector<float> weightPerDomain(this->numDomains(), 1.0/numProcs);

  // Call Parmetis and get the new partitioning.
  int wgtflag = 2;              // only providing vertex weights
  int numflag = 0;              // C style numbering
  int ndims = Dimension::nDim;  // number of dimensions
  int ncon = 1;                 // number of weights per vertex
  int nparts = this->numDomains();    // number of domains we want
  float ubvec[1] = {1.05};        // imbalance tolerance per vertex weight
  int options[3] = {0, 0, 0};  // optional parameters that can be passed
  int edgecut;                 // Output -- number of edges cut in the parmetis decomp.
  vector<idx_t> part((vector<idx_t>::size_type) numLocalNodes); // Output, the parmetis decoposition.
  MPI_Barrier(Communicator::communicator());
  ParMETIS_V3_RefineKway(&(*vtxdist.begin()), 
                         &(*xadj.begin()), 
                         &(*adjacency.begin()),
                         &(*vweight.begin()),
                         NULL,
                         &wgtflag,
                         &numflag,
                         &ncon,
                         &nparts,
                         &(*weightPerDomain.begin()),
                         ubvec,
                         options,
                         &edgecut,
                         &(*part.begin()),
                         &Communicator::communicator());
  MPI_Barrier(Communicator::communicator());

  // Loop over the domain decomposition, and fill in the new domain assignments.
  CHECK(nodeDistribution.size() == part.size());
  for (int i = 0; i < numLocalNodes; ++i) {
    CHECK(part[i] >= 0 && part[i] < numProcs);
    nodeDistribution[i].domainID = part[i];
  }

  // OK, the localDistribution now holds the desired redistribution of the nodes.
  // Go ahead and redistribute them.
  CHECK(validDomainDecomposition(nodeDistribution, dataBase));
  enforceDomainDecomposition(nodeDistribution, dataBase);
}

//------------------------------------------------------------------------------
// Build the CSR adjacency information used by Parmetis.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ParmetisRedistributeNodes<Dimension>::
buildCSRGraph(const DataBase<Dimension>& dataBase,
              const vector<DomainNode<Dimension> >& nodeDistribution,
              const FieldList<Dimension, int>& globalNodeIDs,
              vector<idx_t>& xadj,
              vector<idx_t>& adjacency,
              vector<idx_t>& vtxdist,
              vector<idx_t>& vweight,
              vector<float>& xyz) const {


  // Grab that parallel info.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();

  // Create a list of flags which we will use to keep track of which nodes have
  // been calculated.
  FieldList<Dimension, int> flagNodeDone = dataBase.newGlobalFieldList(int());
  flagNodeDone = 0;

  // Get the state fields from the DataBase.
  const FieldList<Dimension, Vector> positions = dataBase.globalPosition();
  const FieldList<Dimension, SymTensor> Hfield = dataBase.globalHfield();

  // The square of the cutoff distance, which we'll use to determine which
  // nodes are neighbors to whom.
  const Scalar cutoff2 = normalizedNodeExtent()*normalizedNodeExtent();
  CHECK(cutoff2 > 0.0);

  // The map of globalNodeID connectivities.
  ConnectivityType neighbors;

  // Loop over the internal nodes in the DataBase.
  int totalNumNeighbors = 0;
  int minGlobalID = INT_MAX;
  int maxGlobalID = -INT_MAX;
  for (InternalNodeIterator<Dimension> nodeItr = dataBase.internalNodeBegin();
       nodeItr != dataBase.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (flagNodeDone(nodeItr) == 0) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      dataBase.setMasterNodeLists(positions(nodeItr), Hfield(nodeItr));

      // Now loop over all the master nodes.
      for (MasterNodeIterator<Dimension> masterItr = dataBase.masterNodeBegin();
           masterItr != dataBase.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone(masterItr) == 0);

        // State for node I.
        const Vector& ri = positions(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const int globalNodeI = globalNodeIDs(masterItr);
        minGlobalID = std::min(minGlobalID, globalNodeI);
        maxGlobalID = std::max(maxGlobalID, globalNodeI);

        // Insert an entry into the connectivity map for this node.
        CHECK(neighbors.find(globalNodeI) == neighbors.end());
        neighbors[globalNodeI] = vector<pair<int, double> >();

        // Set the refined neighbor information for this master node.
        dataBase.setRefineNodeLists(positions(masterItr), Hfield(masterItr));

        // Loop over the refined neighbors, and determine what nodes
        // are attached to this one.
        for (RefineNodeIterator<Dimension> refineItr = dataBase.refineNodeBegin();
             refineItr != dataBase.refineNodeEnd();
             ++refineItr) {

          // State for node J.
          const Vector& rj = positions(refineItr);
          const SymTensor& Hj = Hfield(refineItr);
          const int globalNodeJ = globalNodeIDs(refineItr);

          // Compute the normalized distance between these nodes.
          const Vector rij = ri - rj;
          const Scalar etai2 = (Hi*rij).magnitude2();
          const Scalar etaj2 = (Hj*rij).magnitude2();

          // If this node counts, add it to our connectivity list.
          // We force the symmetric Gather/Scatter formalism here, since ParMETIS
          // requires a symmetric graph.
          if (masterItr != refineItr &&
              (etai2 < cutoff2 || etaj2 < cutoff2)) {
            neighbors[globalNodeI].push_back(pair<int, double>(globalNodeJ, rij.magnitude2()));
            ++totalNumNeighbors;
          }
        }

        // Flag this master node as done.
        flagNodeDone(masterItr) = 1;
      }
    }
  }
  
  BEGIN_CONTRACT_SCOPE
  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (InternalNodeIterator<Dimension> itr = flagNodeDone.internalNodeBegin();
       itr != flagNodeDone.internalNodeEnd();
       ++itr) ENSURE(flagNodeDone(itr) == 1);
  END_CONTRACT_SCOPE

  // We need to ensure all nodes have at least one neighor, so if there are any
  // isolated ones assign them a neighbor.
  {
    // First find the first node on this domain that has neighbors.
    int gid = 0;
    {
      typename ConnectivityType::iterator itr = neighbors.begin();
      while (itr != neighbors.end() && itr->second.size() == 0) ++itr;
      CHECK(itr != neighbors.end());
      gid = itr->first;
    }

    // Now assign all lone nodes to have this as a neighbor.
    for (typename ConnectivityType::iterator itr = neighbors.begin();
         itr != neighbors.end();
         ++itr) {
      if (itr->second.size() == 0) {
        itr->second.push_back(pair<int, double>(gid, 1.0));
        neighbors[gid].push_back(pair<int, double>(itr->first, 1.0));
      }
    }
  }

//   // Force the connectivity to be symmetric.
//   enforceSymmetricConnectivity(neighbors);

  BEGIN_CONTRACT_SCOPE
  // Make sure there are no nodes without neighbors.
  for (typename ConnectivityType::iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) CHECK(itr->second.size() > 0);
  END_CONTRACT_SCOPE

  // Print the initial statistics about the connectivity of the nodes.
  printConnectivityStatistics(neighbors);

  // Use the full connectivity to construct the nodal weights.
  const int numLocalNodes = nodeDistribution.size();
  CHECK(numLocalNodes == neighbors.size());
  vweight = vector<idx_t>();
  vweight.reserve(numLocalNodes);
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) vweight.push_back(neighbors[itr->globalNodeID].size());

//   // Now try to cull the graph back to something small, since ParMETIS seems to have trouble
//   // with a typical full SPH graph.
//   cullConnectivity(neighbors, nodeDistribution);

  // Print the final statistics about the connectivity of the nodes.
  printConnectivityStatistics(neighbors);

  // Check that the internal connectivity map is sensible.
  BEGIN_CONTRACT_SCOPE
  {
    const bool checkConnect = validConnectivity(nodeDistribution, dataBase, neighbors, globalNodeIDs);
    CHECK(checkConnect);
  }
  END_CONTRACT_SCOPE

  // Size the output vectors.
  xadj = vector<idx_t>();
  adjacency = vector<idx_t>();
  vtxdist = vector<idx_t>();
  xyz = vector<float>();
  xadj.reserve(numLocalNodes + 1);
  adjacency.reserve(totalNumNeighbors);
  vtxdist.reserve(numProcs + 1);
  xyz.reserve(Dimension::nDim*numLocalNodes);

  // Double check that the DomainNodes are arranged properly.
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr < nodeDistribution.end() - 1;
       ++itr) CHECK(itr->globalNodeID == (itr + 1)->globalNodeID - 1);

  // Now construct the CSR format for the graph.
  xadj.push_back(0);
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) {
    const int globalNodeID = itr->globalNodeID;
    CHECK(neighbors.find(globalNodeID) != neighbors.end());
    vector<pair<int, double> >& neighborsi = neighbors[globalNodeID];

    // With boundary conditions, we may have redundant entries in the node
    // connectivity.  Remove any of those redundancies.
    sort(neighborsi.begin(), neighborsi.end(), SortNeighborsByID());
    vector<pair<int, double> >::iterator uniqueEnd = unique(neighborsi.begin(), neighborsi.end());
    for (vector<pair<int, double> >::const_iterator neighborItr = neighborsi.begin();
         neighborItr != uniqueEnd;
         ++neighborItr) adjacency.push_back(neighborItr->first);
    xadj.push_back(adjacency.size());
    CHECK(adjacency.size() == xadj.back());
  }
  CHECK(xadj.size() == numLocalNodes + 1);
  CHECK(adjacency.size() <= totalNumNeighbors);
  CHECK(vweight.size() == numLocalNodes);

  // Build the Parmetis form for the node positions.
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) {
    const Vector& ri = itr->position;
    for (typename Vector::const_iterator xItr = ri.begin();
         xItr != ri.end();
         ++xItr) xyz.push_back(float(*xItr));
  }
  CHECK(xyz.size() == Dimension::nDim*numLocalNodes);

  // Build the vtxdist array, which is just the info about the range of global
  // node IDs on each processor.
  MPI_Barrier(Communicator::communicator());
  vtxdist.push_back(0);
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    int endNodeID = nodeDistribution.back().globalNodeID + 1;
    MPI_Bcast(&endNodeID, 1, MPI_INT, sendProc, Communicator::communicator());
    vtxdist.push_back(endNodeID);
  }
  CHECK(vtxdist.size() == numProcs + 1);

  // Make sure that the CSR form of the graph is OK.
  ENSURE(validCSRGraph(nodeDistribution, dataBase, xadj, adjacency, vtxdist));
}


//------------------------------------------------------------------------------
// Cull the neighbor set down to a small number of nearest neighbors for each
// node.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ParmetisRedistributeNodes<Dimension>::
cullConnectivity(map<int, vector<pair<int, double> > >& neighbors,
                 const vector<DomainNode<Dimension> >& nodeDistribution) const {

  // The target number of neighbors for each node.
  const int ntarget = int(std::pow(4.0, Dimension::nDim) + 0.5);

  // Go over the starting connectivity, and cull down to the target number of 
  // neighbors for each point.
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) {
    const Vector& ri = itr->position;
    vector<pair<int, double> >& neighborsi = neighbors[itr->globalNodeID];

    // Can we afford to pitch any neighbors for this node?
    if (neighborsi.size() > ntarget) {

      // Sort the neighbor set in order of distance.
      sort(neighborsi.begin(), neighborsi.end(), SortNeighborsByDistance());

      // Now pick only the target number of closest neighbors as the ones we
      // want to keep.
      CHECK(neighborsi.begin() + ntarget < neighborsi.end());
      neighborsi.erase(neighborsi.begin() + ntarget, neighborsi.end());
    }
    CHECK(neighborsi.size() >= 1);
  }

  // Now we have to make sure that the connectivity is symmetric.  If any 
  // symmetric pairs were rejected (which almost certainly happened), restore them.
  enforceSymmetricConnectivity(neighbors);
}

//------------------------------------------------------------------------------
// Print the neighbor statistics for the given neighbor set.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ParmetisRedistributeNodes<Dimension>::
printConnectivityStatistics(const map<int, vector<pair<int, double> > >& neighbors) const {

  int minNeighbor = 100000;
  int maxNeighbor = 0;
  double avgNeighbor = 0.0;
  int navgNeighbor = neighbors.size();
  for (ConnectivityType::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) {
    const int n = itr->second.size();
    minNeighbor = std::min(minNeighbor, n);
    maxNeighbor = std::max(maxNeighbor, n);
    avgNeighbor += n;
  }
  const int procID = this->domainID();
  const int numProcs = this->numDomains();
  if (procID == 0) {
    for (int recvProc = 1; recvProc < numProcs; ++recvProc) {
      int minRecv, maxRecv, navgRecv;
      double avgRecv;
      MPI_Status status1, status2, status3, status4;
      MPI_Recv(&minRecv, 1, MPI_INT, recvProc, 10, Communicator::communicator(), &status1);
      MPI_Recv(&maxRecv, 1, MPI_INT, recvProc, 11, Communicator::communicator(), &status2);
      MPI_Recv(&avgRecv, 1, MPI_DOUBLE, recvProc, 12, Communicator::communicator(), &status3);
      MPI_Recv(&navgRecv, 1, MPI_INT, recvProc, 13, Communicator::communicator(), &status4);
      minNeighbor = std::min(minNeighbor, minRecv);
      maxNeighbor = std::max(maxNeighbor, maxRecv);
      avgNeighbor += avgRecv;
      navgNeighbor += navgRecv;
    }

    CHECK(navgNeighbor > 0);
    avgNeighbor /= navgNeighbor;
    cout << "ParmetisRedistributeNodes:: min connections = "
         << minNeighbor << endl
         << "                            max connections = "
         << maxNeighbor << endl
         << "                            avg connections = "
         << avgNeighbor << endl;
   
  } else {
    MPI_Send(&minNeighbor, 1, MPI_INT, 0, 10, Communicator::communicator());
    MPI_Send(&maxNeighbor, 1, MPI_INT, 0, 11, Communicator::communicator());
    MPI_Send(&avgNeighbor, 1, MPI_DOUBLE, 0, 12, Communicator::communicator());
    MPI_Send(&navgNeighbor, 1, MPI_INT, 0, 13, Communicator::communicator());
  }

}

//------------------------------------------------------------------------------
// Check the connectivity map to see if it's minimally sensible.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ParmetisRedistributeNodes<Dimension>::
validConnectivity(const vector<DomainNode<Dimension> >& nodeDistribution,
                  const DataBase<Dimension>& dataBase,
                  const map<int, vector<pair<int, double> > >& neighbors,
		  const FieldList<Dimension, int>& globalNodeIDs) const {

  // The result.
  bool valid = true;

  // First, do we have the right number of entries?
  const int numLocalNodes = nodeDistribution.size();
  const int numProcs = this->numDomains();
  const int procID = this->domainID();
  valid = valid && neighbors.size() == numLocalNodes;

  // Make sure that all the global node IDs are in a kosher range.
  const int numGlobal = numGlobalNodes(dataBase);
  int minGlobalID = INT_MAX;
  int maxGlobalID = -INT_MAX;
  for (typename ConnectivityType::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) {
    valid = valid && itr->first >= 0 && itr->first < numGlobal;
    minGlobalID = min(minGlobalID, itr->first);
    maxGlobalID = max(maxGlobalID, itr->first);
    for (vector<pair<int, double> >::const_iterator nodeItr = itr->second.begin();
         nodeItr != itr->second.end();
         ++nodeItr) {
      valid = valid && nodeItr->first >= 0 && nodeItr->first < numGlobal;
    }
  }
  valid = valid && (maxGlobalID - minGlobalID + 1) == numLocalNodes;
  if (!valid) {
    cerr << "validConnectivity: failing global id range" << endl;
    return valid;
  }

  // Make sure that the graph is symmetric.
  valid = valid && symmetricConnectivity(neighbors);
  if (!valid) cerr << "validConnectivity: failing symmetry test." << endl;

  return valid;
}

//------------------------------------------------------------------------------
// Check the given CSR adjacency information to see if it's minimally sensible.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ParmetisRedistributeNodes<Dimension>::
validCSRGraph(const vector<DomainNode<Dimension> >& nodeDistribution,
              const DataBase<Dimension>& dataBase,
              const vector<idx_t>& xadj,
              const vector<idx_t>& adjacency,
              const vector<idx_t>& vtxdist) const {

  // The result.
  bool valid = true;

  // First, do we have the right number of entries?
  const int numLocalNodes = nodeDistribution.size();
  const int numProcs = this->numDomains();
  valid = valid && xadj.size() == numLocalNodes + 1;
  valid = valid && (xadj.front() == 0 &&
                    xadj.back() == adjacency.size());
  if (!valid) {
    cerr << "validCSRGraph: failing xadj size test:  "
         << xadj.size() << " "
         << numLocalNodes + 1 << " "
         << xadj.front() << " "
         << xadj.back() << " "
         << adjacency.size() << endl;
    return valid;
  }

  // The xadj array should be monotonically increasing.
  {
    vector<idx_t>::const_iterator xadjItr = xadj.begin();
    while (xadjItr < xadj.end() - 1 && valid) {
      valid = valid && *xadjItr <= *(xadjItr + 1);
      ++xadjItr;
    }
  }
  if (!valid) {
    cerr << "validCSRGraph: failing xadj monotonicity test" << endl;
    return valid;
  }

  // Make sure that all the global node IDs are in a kosher range.
  const int numGlobal = numGlobalNodes(dataBase);
  for (vector<idx_t>::const_iterator adjItr = adjacency.begin();
       adjItr < adjacency.end();
       ++adjItr) {
    valid = valid && *adjItr >= 0 && *adjItr < numGlobal;
  }
  if (!valid) {
    cerr << "validCSRGraph: failing globalID range test" << endl;
    return valid;
  }

  // Make sure the vtxdist array is sized properly.
  valid = valid && (vtxdist.size() == numProcs + 1);
  if (!valid) {
    cerr << "validCSRGraph: failing vtxdist size:"
         << vtxdist.size() << " "
         << numProcs + 1 << endl;
    return valid;
  }

  // Make sure that the vtxdist array contains the proper global IDs for each
  // processor.
  const int procID = this->domainID();
  CHECK(vtxdist.size() > procID + 1);
  const int lowID = vtxdist[procID];
  const int highID = vtxdist[procID + 1];
  typename vector<DomainNode<Dimension> >::const_iterator
    nodeItr = nodeDistribution.begin();
  while (nodeItr < nodeDistribution.end() && valid) {
    valid = valid && (nodeItr->globalNodeID >= lowID &&
                      nodeItr->globalNodeID < highID);
    ++nodeItr;
  }
  if (!valid) {
    --nodeItr;
    cerr << "validCSRGraph: failing vtxdist range "
         << lowID << " " << highID << " " << nodeItr->globalNodeID << endl;
    return valid;
  }

  // Make sure that the graph is symmetric.
  // First we build up our standard neighbor map used for connectivity.
  ConnectivityType neighbors;
  int k = 0;
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) {
    const int i = itr->globalNodeID;
    neighbors[i] = vector<pair<int, double> >();
    CHECK(k >= 0 && k < xadj.size() - 1);
    const int firstIndex = xadj[k];
    const int lastIndex = xadj[k + 1];
    for (int j = firstIndex; j != lastIndex; ++j) {
      CHECK(j < adjacency.size());
      neighbors[i].push_back(pair<int, double>(adjacency[j], 0.0));
    }
    ++k;
  }

  valid = valid && symmetricConnectivity(neighbors);

  return valid;
}

//------------------------------------------------------------------------------
// Helper method to check the given CSR graph lists a given node as connected
// to another given node.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ParmetisRedistributeNodes<Dimension>::
verifyNeighborPresent(const int globalNeighborID,
                      const int globalID,
                      const vector<DomainNode<Dimension> >& nodeDistribution,
                      const vector<int>& xadj,
                      const vector<int>& adjacency,
                      const vector<int>& vtxdist) const {

  // Make sure that this global node is indeed on this domain.
  {
    vector<int>::const_iterator itr = upper_bound(vtxdist.begin(),
                                                  vtxdist.end(),
                                                  globalNeighborID);
    CHECK(itr != vtxdist.end());
    const int dom = distance(vtxdist.begin(), itr) - 1;
    CHECK(dom == this->domainID());
  }

  // Find the local node ID for this global ID.
  typename vector<DomainNode<Dimension> >::const_iterator 
    itr = nodeDistribution.begin();
  while (itr < nodeDistribution.end() &&
         itr->globalNodeID != globalNeighborID) ++itr;
  CHECK(itr < nodeDistribution.end());
  const int localID = itr->uniqueLocalNodeID;

  // Now look for the given global ID in the list of adjacency neighbors for 
  // this local ID.
  CHECK(localID >= 0 && localID < xadj.size());
  const int lowID = xadj[localID];
  const int highID = xadj[localID + 1];
  int i = lowID;
  while (i < highID && adjacency[i] != globalID) ++i;
  if (i < highID) {
    return true;
  } else {
    return false;
  }
}

//------------------------------------------------------------------------------
// Compute the inverse of a connectivity graph.
//------------------------------------------------------------------------------
template<typename Dimension>
map<int, vector<pair<int, double> > >
ParmetisRedistributeNodes<Dimension>::
inverseConnectivity(const map<int, vector<pair<int, double> > >& neighbors) const {

  const int numProcs = this->numDomains();
  const int procID = this->domainID();

  // Find the min and max global IDs on this domain.
  int minGlobalID = INT_MAX;
  int maxGlobalID = -INT_MAX;
  for (typename ConnectivityType::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) {
    minGlobalID = min(minGlobalID, itr->first);
    maxGlobalID = max(maxGlobalID, itr->first);
  }

  // Make sure the input connectivity has entries for all our nodes.
  BEGIN_CONTRACT_SCOPE
  for (int i = minGlobalID; i != maxGlobalID + 1; ++i) REQUIRE(neighbors.find(i) != neighbors.end());
  END_CONTRACT_SCOPE

  // Build as much of the inverse as we can on the local domain.
  ConnectivityType result;
  for (typename ConnectivityType::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) {
    const int i = itr->first;
    const vector<pair<int, double> >& js = itr->second;
    for (vector<pair<int, double> >::const_iterator jItr = js.begin();
         jItr != js.end();
         ++jItr) {
      const int j = jItr->first;
      if (result.find(j) == result.end()) result[j] = vector<pair<int, double> >();
      result[j].push_back(pair<int, double>(i, jItr->second));
    }
  }

  // Send nodes which are not local to this processor to where they belong.
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    vector<int> sendIDs;
    if (procID == sendProc) {
      for (typename ConnectivityType::const_iterator itr = result.begin();
           itr != result.end();
           ++itr) {
        const int i = itr->first;
        if (i < minGlobalID || i > maxGlobalID) sendIDs.push_back(i);
      }
    }
    int numSend = sendIDs.size();
    MPI_Bcast(&numSend, 1, MPI_INT, sendProc, Communicator::communicator());
    for (int k = 0; k != numSend; ++k) {
      int globalID = -1;
      int numNeighbors = -1;
      vector<int> neighborIDs;
      vector<double> distances;
      if (sendProc == procID) {
        globalID = sendIDs[k];
        typename ConnectivityType::iterator itr = result.find(globalID);
        CHECK(itr != result.end());
        numNeighbors = itr->second.size();
        neighborIDs.reserve(numNeighbors);
        distances.reserve(numNeighbors);
        for (typename vector<pair<int, double> >::const_iterator jItr = itr->second.begin();
             jItr != itr->second.end();
             ++jItr) {
          neighborIDs.push_back(jItr->first);
          distances.push_back(jItr->second);
        }
        CHECK(neighborIDs.size() == numNeighbors);
        CHECK(distances.size() == numNeighbors);
        result.erase(itr);
      }
      MPI_Bcast(&globalID, 1, MPI_INT, sendProc, Communicator::communicator());
      MPI_Bcast(&numNeighbors, 1, MPI_INT, sendProc, Communicator::communicator());
      CHECK(numNeighbors > 0);
      neighborIDs.resize(numNeighbors);
      distances.resize(numNeighbors);
      MPI_Bcast(&(*neighborIDs.begin()), numNeighbors, MPI_INT, sendProc, Communicator::communicator());
      MPI_Bcast(&(*distances.begin()), numNeighbors, MPI_DOUBLE, sendProc, Communicator::communicator());
      if (globalID >= minGlobalID && globalID <= maxGlobalID) {
        if (result.find(globalID) == result.end()) result[globalID] = vector<pair<int, double> >();
        vector<pair<int, double> >& currentNeighbors = result[globalID];
        currentNeighbors.reserve(currentNeighbors.size() + numNeighbors);
        for (int j = 0; j != numNeighbors; ++j) 
          currentNeighbors.push_back(pair<int, double>(neighborIDs[j], distances[j]));
      }
    }
  }

  // If any of our local IDs was not found in the inverse mapping, add them.
  for (int i = minGlobalID; i != maxGlobalID + 1; ++i) {
    if (result.find(i) == result.end()) result[i] = vector<pair<int, double> >();
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    // We should only have the number of global IDs.
    const int numIDs = maxGlobalID - minGlobalID + 1;
    ENSURE(result.size() == numIDs);

    // We should have all of our global IDs.
    for (typename ConnectivityType::const_iterator itr = result.begin();
         itr != result.end();
         ++itr) {
      ENSURE(itr->first >= minGlobalID && itr->first <= maxGlobalID);
    }
  }
  END_CONTRACT_SCOPE

  return result;
}

//------------------------------------------------------------------------------
// Force the connectivity to be symmetric.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ParmetisRedistributeNodes<Dimension>::
enforceSymmetricConnectivity(map<int, vector<pair<int, double> > >& neighbors) const {

  // Invert the connectivity map.
  ConnectivityType neighborInverse = inverseConnectivity(neighbors);

  // If we're missing any neighbors from the reverse mapping, add them!
  for (ConnectivityType::iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) {
    const int i = itr->first;
    vector<pair<int, double> >& neighborsi = itr->second;
    if (neighborInverse.find(i) != neighborInverse.end()) {
      const vector<pair<int, double> >& inversei = neighborInverse[i];
      for (vector<pair<int, double> >::const_iterator jItr = inversei.begin();
           jItr != inversei.end();
           ++jItr) {
        if (find(neighborsi.begin(), neighborsi.end(), *jItr) == neighborsi.end())
          neighborsi.push_back(*jItr);
      }
    }
  }

  // Post-conditions.
  ENSURE(symmetricConnectivity(neighbors));
}

//------------------------------------------------------------------------------
// Make sure that the given connectivity graph is symmetric.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ParmetisRedistributeNodes<Dimension>::
symmetricConnectivity(const map<int, vector<pair<int, double> > >& neighbors) const {

  // Build up a duplicate neighbor array by inverting the current
  // list.  This *should* wind up identical to our current 
  // neighbor set if it truly is symmetric.
  ConnectivityType neighborCheck = inverseConnectivity(neighbors);

  // Do both the neighbor map and it's "reverse" have the right number of entries?
  if (neighbors.size() != neighborCheck.size()) {
    cerr << "Wrong size for inverse neighbor map in symmetry check "
         << neighbors.size() << " "
         << neighborCheck.size() << " "
         << endl;
    return false;
  }

  // The sorted values for each neighbor map should be the same.
  for (typename ConnectivityType::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) {
    const int i = itr->first;
    vector<pair<int, double> >  n0 = itr->second;
    vector<pair<int, double> >& n1 = neighborCheck[i];
    sort(n0.begin(), n0.end(), SortNeighborsByID());
    sort(n1.begin(), n1.end(), SortNeighborsByID());
    if (n0 != n1) {
      return false;
    }
  }

  return true;
}

}

