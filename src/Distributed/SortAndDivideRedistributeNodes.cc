//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes -- (Re)domain decompose the nodes by using 
// a sort and divide algorithm.
//
// This is a template base class -- the actual dimension dependent objects
// are: 
//   SortAndDivideRedistributeNodes1d
//   SortAndDivideRedistributeNodes2d
//   SortAndDivideRedistributeNodes3d
//
// Created by JMO, Thu Dec 2 10:44:07 2004
//----------------------------------------------------------------------------//
#include "SortAndDivideRedistributeNodes.hh"
#include "CompareDomainNodesByPosition.hh"
#include "Field/FieldList.hh"
#include "Geometry/EigenStruct.hh"
#include "Communicator.hh"

#include <limits>
using std::vector;
using std::list;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Sort the eigen values/vectors such that the eigen values are arranged from
// max to min.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
swapIndicies(EigenStruct<nDim>& result,
             const size_t i,
             const size_t j) {
  REQUIRE(i < nDim);
  REQUIRE(j < nDim);
  REQUIRE(i != nDim);
  std::swap(result.eigenValues(i), result.eigenValues(j));
  const GeomVector<nDim> tmp = result.eigenVectors.getColumn(i);
  result.eigenVectors.setColumn(i, result.eigenVectors.getColumn(j));
  result.eigenVectors.setColumn(j, tmp);
}

inline
void
sortEigenValues(EigenStruct<1>& /*result*/) {
}
  
inline
void
sortEigenValues(EigenStruct<2>& result) {
  if (result.eigenValues(0) < result.eigenValues(1)) swapIndicies<2>(result, 0, 1);
}

inline
void
sortEigenValues(EigenStruct<3>& result) {
  if (result.eigenValues(0) < result.eigenValues(1)) swapIndicies<3>(result, 0, 1);
  if (result.eigenValues(1) < result.eigenValues(2)) swapIndicies<3>(result, 1, 2);
  if (result.eigenValues(0) < result.eigenValues(1)) swapIndicies<3>(result, 0, 1);
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SortAndDivideRedistributeNodes<Dimension>::
SortAndDivideRedistributeNodes(const double Hextent):
  mHextent(Hextent) {
    ENSURE(distinctlyGreaterThan(mHextent, 0.0));
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SortAndDivideRedistributeNodes<Dimension>::
~SortAndDivideRedistributeNodes() {
}

//------------------------------------------------------------------------------
// Sort the set of domain nodes by positions, returning the result as a list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SortAndDivideRedistributeNodes<Dimension>::
sortByPositions(list<DomainNode<Dimension> >& domainNodes,
                const int positionIndex) const {

  // Create an appropriate sorting criterion.
  const CompareDomainNodesByPosition<Dimension> cmp(positionIndex);
  domainNodes.sort(cmp);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  for (typename list<DomainNode<Dimension> >::const_iterator itr1 = domainNodes.begin();
       itr1 != domainNodes.end();
       ++itr1) {
    typename list<DomainNode<Dimension> >::const_iterator itr2 = itr1;
    ++itr2;
    if (itr2 != domainNodes.end()) {
      ENSURE(itr1->position(positionIndex) <= itr2->position(positionIndex));
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Given a sorted list of domain nodes, pop nodes off of the front building
// a return vector until the requested work is met.
//------------------------------------------------------------------------------
template<typename Dimension>
list<DomainNode<Dimension> >
SortAndDivideRedistributeNodes<Dimension>::
popFrontNodes(list<DomainNode<Dimension> >& sortedCandidateNodes,
              const double targetDomainWork,
              const int positionIndex) const {

  // Prepare the result.
  list<DomainNode<Dimension> > result;

  // Build a comparator.
  const CompareDomainNodesByPosition<Dimension> cmp(positionIndex);

  // The usual parallel IDs.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();

  // The amount of work we'll examine from each domain per iteration.
  const double chunkWork = targetDomainWork; //  / numProcs;

  // Iterate over the candidate nodes until we either fill the target
  // domains work quota or run out of nodes.
  double totalWork = 0.0;
  bool domainFinished = false;
  int numAvailableNodes = sortedCandidateNodes.size();
  int globalNumAvailableNodes;
  MPI_Allreduce(&numAvailableNodes, &globalNumAvailableNodes, 1, MPI_INT, MPI_SUM, Communicator::communicator());
  while (!domainFinished && (globalNumAvailableNodes > 0)) {

    // Have each processor select its available candidate nodes up to the chunk work size.
    vector<DomainNode<Dimension> > currentCandidates;
    double currentWork = 0.0;
    {
      typename list<DomainNode<Dimension> >::iterator itr = sortedCandidateNodes.begin();
      while (itr != sortedCandidateNodes.end() && 
             (abs(currentWork + itr->work - chunkWork) < abs(currentWork - chunkWork))) {
        currentCandidates.push_back(*itr);
        currentWork += itr->work;
        ++itr;
      }
    }

    // Reduce the current set of candidates to process 0, sort them, and then
    // communicate the result back to everyone.
    currentCandidates = reduceDomainNodes(currentCandidates, 0);
    if (procID == 0) sort(currentCandidates.begin(), currentCandidates.end(), cmp);
    MPI_Barrier(Communicator::communicator());
    currentCandidates = broadcastDomainNodes(currentCandidates, 0);

    // Step through the candidate nodes, assigning them and incrementing the 
    // work until we either fill our work quota or run out of nodes.
    {
      typename vector<DomainNode<Dimension> >::const_iterator itr = currentCandidates.begin();
      while (itr != currentCandidates.end() && 
             (!domainFinished && (globalNumAvailableNodes > 0))) {
        const double currentDeltaWork = abs(totalWork - targetDomainWork);
        const double proposedDeltaWork = abs(totalWork + itr->work - targetDomainWork);
        domainFinished = proposedDeltaWork > currentDeltaWork;
        if (!domainFinished) {
          totalWork += itr->work;
          --globalNumAvailableNodes;
          if (sortedCandidateNodes.size() > 0 &&
              itr->globalNodeID == sortedCandidateNodes.front().globalNodeID) {
            result.push_back(*itr);
            sortedCandidateNodes.pop_front();
          }
        }
        ++itr;
      }
    }

    BEGIN_CONTRACT_SCOPE
    {
      double sumTotalWork;
      int sumNumAvailableNodes;
      MPI_Allreduce(&totalWork, &sumTotalWork, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
      MPI_Allreduce(&globalNumAvailableNodes, &sumNumAvailableNodes, 1, MPI_INT, MPI_SUM, Communicator::communicator());
      CONTRACT_VAR(numProcs);
      ENSURE(fuzzyEqual(sumTotalWork, numProcs * totalWork, 1.0e-12));
      ENSURE(sumNumAvailableNodes == numProcs * globalNumAvailableNodes);
    }
    END_CONTRACT_SCOPE

  }

  return result;

}

//------------------------------------------------------------------------------
// Compute the appropriate shape tensor for a set of domain nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor::EigenStructType
SortAndDivideRedistributeNodes<Dimension>::
shapeTensor(const vector<DomainNode<Dimension> >& domainNodes) const {

  // Iterate over all the nodes to find the center of "weight" (assuming unit mass
  // per node).
  typedef vector<DomainNode<Dimension> > ContainerType;
  Vector com;
  for (typename ContainerType::const_iterator itr = domainNodes.begin();
       itr != domainNodes.end();
       ++itr) com += itr->position;
  Vector globalCOM;
  int localNumNodes = domainNodes.size();
  int globalNumNodes;
  MPI_Allreduce(&(*com.begin()), &(*globalCOM.begin()), Dimension::nDim, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
  MPI_Allreduce(&localNumNodes, &globalNumNodes, 1, MPI_INT, MPI_SUM, Communicator::communicator());
  CHECK(globalNumNodes > 0);
  globalCOM /= globalNumNodes;

  // Now find the second moment of the nodes about this center of mass.
  SymTensor J;
  for (typename ContainerType::const_iterator itr = domainNodes.begin();
       itr != domainNodes.end();
       ++itr) {
    const Vector r = itr->position - globalCOM;
    J += r.selfdyad();
  }
  SymTensor globalJ;
  const int numTensorElements = std::distance(J.begin(), J.end());
  MPI_Allreduce(&(*J.begin()), &(*globalJ.begin()), numTensorElements, MPI_DOUBLE, MPI_SUM, Communicator::communicator());

  // The shape tensor we want is the square root of the second moment.  We also
  // normalize the eigenvalues for a unit volume.
  typename SymTensor::EigenStructType result = globalJ.eigenVectors();
  for (typename Vector::iterator itr = result.eigenValues.begin();
       itr != result.eigenValues.end();
       ++itr) *itr = sqrt(*itr);
  result.eigenValues = result.eigenValues.unitVector();

  // We want the eigen values sorted from max to min.
  sortEigenValues(result);

  // That's it.
  ENSURE(fuzzyEqual(result.eigenValues.magnitude(), 1.0));
  ENSURE(fuzzyEqual(abs(result.eigenVectors.Determinant()), 1.0));
  BEGIN_CONTRACT_SCOPE
  {
    double x = result.eigenValues(0);
    for (int i = 1; i != Dimension::nDim; ++i) {
      CONTRACT_VAR(x);
      ENSURE(x >= result.eigenValues(i));
      x = result.eigenValues(i);
    }
  }
  END_CONTRACT_SCOPE

  return result;
}

//------------------------------------------------------------------------------
// Apply the necessary rotation to the positions of the domain nodes to transform into
// the primary frame of the given shape tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SortAndDivideRedistributeNodes<Dimension>::
rotateIntoShapeTensorFrame(const typename Dimension::SymTensor::EigenStructType& shapeTensor,
                           vector<DomainNode<Dimension> >& domainNodes) const {

  // Get the rotational transformation we'll be using.
  const Tensor R = shapeTensor.eigenVectors.Transpose();

  // Iterate over the nodes.
  typedef vector<DomainNode<Dimension> > ContainerType;
  for (typename ContainerType::iterator itr = domainNodes.begin();
       itr != domainNodes.end();
       ++itr) {

    // Apply the transformation to the positions.
    itr->position = R*(itr->position);

  }

}

//------------------------------------------------------------------------------
// Reduce a vector<DomainNode> to the given processor.
// Upon completion this routine returns the full reduced set of nodes on the
// target processor, and just the input per process set of nodes on all others.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<DomainNode<Dimension> >
SortAndDivideRedistributeNodes<Dimension>::
reduceDomainNodes(const std::vector<DomainNode<Dimension> >& nodes,
                  const int targetProc) const {

  // The usual parallel IDs.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();
  REQUIRE(targetProc >= 0 && targetProc < numProcs);

  // Prepare the result.
  vector<DomainNode<Dimension> > result(nodes);

  // Are we sending or receiving?
  if (procID == targetProc) {

    // We're receiving, so iterate over the other processors.
    for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
      if (sendProc != procID) {
        
        // Get the packed data.
        MPI_Status status1, status2;
        int bufferSize;
        MPI_Recv(&bufferSize, 1, MPI_INT, sendProc, 200, Communicator::communicator(), &status1);
        CHECK(bufferSize % DomainNode<Dimension>::packSize() == 0);
        if (bufferSize > 0) {
          vector<double> buffer(bufferSize);
          MPI_Recv(&(*buffer.begin()), bufferSize, MPI_DOUBLE, sendProc, 201, Communicator::communicator(), &status2);

          // Unpack the data and append it to the result.
          const int oldNumNodes = result.size();
          const int numRecvNodes = bufferSize / DomainNode<Dimension>::packSize();
          result.reserve(oldNumNodes + numRecvNodes);
          typename vector<double>::const_iterator itr = buffer.begin();
          while (itr != buffer.end()) {
            DomainNode<Dimension> node;
            node.unpack(itr);
            CHECK(itr <= buffer.end());
            result.push_back(node);
            CHECK((int)result.size() <= oldNumNodes + numRecvNodes);
          }
          CHECK((int)result.size() == oldNumNodes + numRecvNodes);

        }

      }
    }

  } else {

    // We're sending.
    // Pack up the domain nodes into a communicable buffer.
    vector<double> buffer;
    buffer.reserve(nodes.size() * DomainNode<Dimension>::packSize());
    for (typename vector<DomainNode<Dimension> >::const_iterator nodeItr = nodes.begin();
         nodeItr != nodes.end();
         ++nodeItr) {
      const vector<double> nodeBuffer = nodeItr->pack();
      for (vector<double>::const_iterator itr = nodeBuffer.begin();
           itr != nodeBuffer.end();
           ++itr) buffer.push_back(*itr);
    }

    // Now send that sucker.
    int bufferSize = buffer.size();
    CHECK(bufferSize == (int)nodes.size() * (int)DomainNode<Dimension>::packSize());
    MPI_Send(&bufferSize, 1, MPI_INT, targetProc, 200, Communicator::communicator());
    if (bufferSize > 0) MPI_Send(&(*buffer.begin()), bufferSize, MPI_DOUBLE, targetProc, 201, Communicator::communicator());

  }

  // That's it.
  MPI_Barrier(Communicator::communicator());
  return result;
}

//------------------------------------------------------------------------------
// Broadcast a vector<DomainNode> from the given processor.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<DomainNode<Dimension> >
SortAndDivideRedistributeNodes<Dimension>::
broadcastDomainNodes(const std::vector<DomainNode<Dimension> >& nodes,
                     const int targetProc) const {

  // The usual parallel IDs.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();
  CONTRACT_VAR(numProcs);
  REQUIRE(targetProc >= 0 && targetProc < numProcs);

  // Are we sending or receiving?
  vector<double> buffer;
  if (procID == targetProc) {

    // We're the sending process.  Pack the domain nodes into a send buffer.
    buffer.reserve(nodes.size() * DomainNode<Dimension>::packSize());
    for (typename vector<DomainNode<Dimension> >::const_iterator nodeItr = nodes.begin();
         nodeItr != nodes.end();
         ++nodeItr) {
      const vector<double> nodeBuffer = nodeItr->pack();
      for (vector<double>::const_iterator itr = nodeBuffer.begin();
           itr != nodeBuffer.end();
           ++itr) buffer.push_back(*itr);
    }

  }

  // Now broadcast the packed info.
  int bufferSize = buffer.size();
  MPI_Bcast(&bufferSize, 1, MPI_INT, targetProc, Communicator::communicator());
  CHECK(bufferSize >= 0);
  CHECK(bufferSize % DomainNode<Dimension>::packSize() == 0);
  buffer.resize(bufferSize);
  MPI_Bcast(&(*buffer.begin()), bufferSize, MPI_DOUBLE, targetProc, Communicator::communicator());

  // Unpack the buffer.
  const int numNodes = bufferSize / DomainNode<Dimension>::packSize();
  vector<DomainNode<Dimension> > result;
  result.reserve(numNodes);
  typename vector<double>::const_iterator itr = buffer.begin();
  while (itr != buffer.end()) {
    DomainNode<Dimension> node;
    node.unpack(itr);
    CHECK(itr <= buffer.end());
    result.push_back(node);
  }

  // That's it.
  ENSURE((int)result.size() == numNodes);
  return result;
}

}

