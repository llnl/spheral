//---------------------------------Spheral++----------------------------------//
// SpaceFillingCurveRedistributeNodes
//
// This is an abstract base for the space filling curve family of 
// repartitioners.  The assumption is that the descendent classes will provide
// the computeHashedIndices method to assign unique keys to each point in the
// order that that algorithm wants the points distributed.
//
// Created by JMO, Wed Apr  9 13:13:46 PDT 2008
//----------------------------------------------------------------------------//
#include "SpaceFillingCurveRedistributeNodes.hh"
#include "Utilities/DomainNode.hh"
#include "DistributedBoundary.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/RedistributionRegistrar.hh"
#include "allReduce.hh"
#include "Communicator.hh"

#include "Utilities/DBC.hh"

#include <algorithm>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <bitset>
using std::vector;
using std::pair;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {


//------------------------------------------------------------------------------
// Compare pairs of <Key, *> by the first element (Key).
//------------------------------------------------------------------------------
template<typename DataType>
struct SortNodesByHashedIndex {
  bool operator()(const pair<unsigned long long, DataType>& lhs,
                  const pair<unsigned long long, DataType>& rhs) {
    return lhs.first < rhs.first;
  }
};

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SpaceFillingCurveRedistributeNodes<Dimension>::
SpaceFillingCurveRedistributeNodes(double,
                                   const double minNodesPerDomainFraction,
                                   const double maxNodesPerDomainFraction,
                                   const bool workBalance,
                                   const bool localReorderOnly):
  RedistributeNodes<Dimension>(),
  mMinNodesPerDomainFraction(minNodesPerDomainFraction),
  mMaxNodesPerDomainFraction(maxNodesPerDomainFraction),
  mWorkBalance(workBalance),
  mLocalReorderOnly(localReorderOnly) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SpaceFillingCurveRedistributeNodes<Dimension>::
~SpaceFillingCurveRedistributeNodes() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Call on the SpaceFillingCurve library to describe an
// optimal partioning of the nodes, and then apply that partitioning.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
redistributeNodes(DataBase<Dimension>& dataBase,
                  vector<Boundary<Dimension>*> boundaries) {

  // The usual parallel info.
  const int numProcs = this->numDomains();
  const int procID = this->domainID();

  // The min and max range of nodes we want per domain.
  const int totalNumNodes = dataBase.globalNumInternalNodes();
  const int minNodes = int(mMinNodesPerDomainFraction * double(totalNumNodes)/numProcs);
  const int maxNodes = int(mMaxNodesPerDomainFraction * double(totalNumNodes)/numProcs);

  // Below a certain level of work per process we just punt.
  if (double(totalNumNodes)/double(Process::getTotalNumberOfProcesses()) >= 1.0) {

    // Get the global IDs.
    const FieldList<Dimension, size_t> globalIDs = globalNodeIDs(dataBase);

    // Compute the work per node.
    FieldList<Dimension, Scalar> workField(FieldStorageType::CopyFields);
    if (this->workBalance()) {

      // Enforce boundary conditions for the work computation.
      for (auto* nodeListPtr: dataBase.nodeListPtrs()) {
        nodeListPtr->numGhostNodes(0);
        nodeListPtr->neighbor().updateNodes();
      }
      for (auto* boundaryPtr: boundaries) {
        boundaryPtr->setAllGhostNodes(dataBase);
        boundaryPtr->finalizeGhostBoundary();
        for (auto* nodeListPtr: dataBase.fluidNodeListPtrs()) nodeListPtr->neighbor().updateNodes();
      }

      // Get the local description of the domain distribution, with the work per node filled in.
      workField = this->workPerNode(dataBase, 1.0);

    } else {

      // We assign a constant work of 1 to each node, effectively balancing by
      // node count.
      workField = dataBase.newGlobalFieldList(1.0, "work");

    }

    // Print the beginning statistics.
    std::string stats0 = this->gatherDomainDistributionStatistics(workField);
    if (Process::getRank() == 0) cout << "SpaceFillingCurveRedistributeNodes: INITIAL node distribution statistics:" << endl
                                      << stats0 << endl;

    // Now we can get the node distribution description.
    vector<DomainNode<Dimension> > nodeDistribution = this->currentDomainDecomposition(dataBase, globalIDs, workField);

    // Clear out any ghost nodes.
    // We won't bother to update the neighbor info at this point -- we don't need 
    // it for this algorithm, so we just update it when we're done.
    for (auto* nodeListPtr: dataBase.nodeListPtrs()) nodeListPtr->numGhostNodes(0);

    // Compute the target work per domain.
    const Scalar targetWork = workField.sumElements()/numProcs;
    if (procID == 0) cout << "SpaceFillingCurveRedistributeNodes: Target work per process " << targetWork << endl;

    // Compute the Key indices for each point on this processor.
    if (procID == 0) cout << "SpaceFillingCurveRedistributeNodes: Hashing indices" << endl;
    FieldList<Dimension, Key> indices = computeHashedIndices(dataBase);

    // Find the range of hashed indices.
    const Key indexMin = indices.min();
    const Key indexMax = indices.max();
    CHECK(indexMax < indexMax + indexMax);
    if (procID == 0) cout << "SpaceFillingCurveRedistributeNodes: Index min/max : " << indexMin << " " << indexMax << endl;

    // Build the array of (hashed index, DomainNode) pairs.
    // Note this comes back locally sorted.
    if (procID == 0) cout << "SpaceFillingCurveRedistributeNodes: sorting indices" << endl;
    vector<pair<Key, DomainNode<Dimension> > > sortedIndices = buildIndex2IDPairs(indices,
                                                                                  nodeDistribution);
    const int numLocalNodes = nodeDistribution.size();

    // Build our set of unique indices and their count.
    if (procID == 0) cout << "SpaceFillingCurveRedistributeNodes: Counting uniques and such" << endl;
    vector<Key> uniqueIndices;
    vector<int> count;
    vector<Scalar> work;
    int maxCount = 1;
    if (sortedIndices.size() > 0) {
      uniqueIndices.reserve(sortedIndices.size());
      count.reserve(sortedIndices.size());
      work.reserve(sortedIndices.size());
      uniqueIndices.push_back(sortedIndices.front().first);
      count.push_back(1);
      work.push_back(sortedIndices.front().second.work);
      {
        int j = 0;
        for (auto i = 1u; i < sortedIndices.size(); ++i) {
          CHECK(sortedIndices[i].first >= sortedIndices[i-1].first);
          if (sortedIndices[i].first == uniqueIndices[j]) {
            ++count[j];
            work[j] += sortedIndices[i].second.work;
          } else {
            uniqueIndices.push_back(sortedIndices[i].first);
            count.push_back(1);
            work.push_back(sortedIndices[i].second.work);
            ++j;
          }
          maxCount = max(maxCount, count[j]);
        }
      }
      CHECK(count.size() == uniqueIndices.size());
      CHECK(work.size() == uniqueIndices.size());
    }
    maxCount = allReduce(maxCount, SPHERAL_OP_MAX);
    if (procID == 0) cout << "SpaceFillingCurveRedistributeNodes: max redundancy is " << maxCount << endl;

    //   // DEBUG
    //   {
    //     for (int ii = 0; ii != numProcs; ++ii) {
    //       if (procID == ii) {
    //         for (size_t i = 0; i != uniqueIndices.size(); ++i) {
    //           cerr << uniqueIndices[i] << endl;
    //         }
    //       }
    //       MPI_Barrier(Communicator::communicator());
    //     }
    //   }
    //   // DEBUG

    // Figure out the range of hashed indices we want for each process.
    // Note this will not be optimal when there are degnerate indices!
    Key lowerBound = indexMin;
    vector<pair<Key, Key> > indexRanges;
    for (int iProc = 0; iProc != numProcs; ++iProc) {
      Key upperBound;
      int numNodes;
      findUpperKey(uniqueIndices,
                   count,
                   work,
                   lowerBound,
                   indexMax,
                   targetWork,
                   minNodes,
                   maxNodes,
                   upperBound,
                   numNodes);
      // cerr << "Chose indices in range: "
      //      << lowerBound << " " << upperBound << endl;
      // if (procID == 0) cerr << "  range [" 
      //                       << lowerBound << " " 
      //                       << upperBound << "], work in range "
      //                       << workInRange(uniqueIndices, work, lowerBound, upperBound)
      //                       << ", num in range "
      //                       << numIndicesInRange(uniqueIndices, count, lowerBound, upperBound)
      //                       << endl;
      indexRanges.push_back(pair<Key, Key>(lowerBound, upperBound));
      lowerBound = findNextIndex(uniqueIndices, upperBound, indexMax);
    }
    CHECK(lowerBound == indexMax);
    CHECK(indexRanges[0].first == indexMin);
    CHECK(indexRanges.back().second == indexMax);
    CHECK((int)indexRanges.size() == numProcs);

    // We now know the target index range for each domain.
    // Go through our local DomainNode set and assign them appropriately.
    nodeDistribution = vector<DomainNode<Dimension>>();
    for (auto [key, dnode]: sortedIndices) {
      nodeDistribution.push_back(dnode);
      nodeDistribution.back().domainID = domainForIndex(key, indexRanges);
      CHECK(nodeDistribution.back().domainID >= 0 and
            nodeDistribution.back().domainID < numProcs);
    }
    CONTRACT_VAR(numLocalNodes);
    CHECK((int)nodeDistribution.size() == numLocalNodes);

    // Redistribute nodes between domains.
    CHECK(this->validDomainDecomposition(nodeDistribution, dataBase));
    if (not this->localReorderOnly()) this->enforceDomainDecomposition(nodeDistribution, dataBase);

    // Sort the local indices in each NodeList to find the local ordering.
    for (auto nodeListPtr: dataBase.nodeListPtrs()) {
      const auto n = nodeListPtr->numInternalNodes();
      const auto& keys = **indices.fieldForNodeList(*nodeListPtr);
      vector<pair<Key, int>> orderedkeys(n);
      for (auto i = 0u; i < n; ++i) orderedkeys[i] = make_pair(keys(i), i);
      sort(orderedkeys.begin(), orderedkeys.end(), SortNodesByHashedIndex<int>());  // [](const pair<Key, int>& lhs, const pair<Key, int>& rhs) { return lhs.first < rhs.first; });
      vector<size_t> ordering(n);
      for (auto i = 0u; i < n; ++i) ordering[orderedkeys[i].second] = i;
      nodeListPtr->reorderNodes(ordering);
      nodeListPtr->neighbor().updateNodes();
    }

    // Notify everyone that the nodes have just been shuffled around.
    RedistributionRegistrar::instance().broadcastRedistributionNotifications();

    // Print the final statistics.
    std::string stats1 = this->gatherDomainDistributionStatistics(workField);
    if (Process::getRank() == 0) cout << "SpaceFillingCurveRedistributeNodes: FINAL node distribution statistics:" << endl
                                      << stats1 << endl;

    // Post-conditions.
    // Make sure the nodes are now sorted by the index keys.
    BEGIN_CONTRACT_SCOPE
      {
        for (auto* nodeListPtr: dataBase.nodeListPtrs()) {
          const Field<Dimension, Key>& keyField = **indices.fieldForNodeList(*nodeListPtr);
          CONTRACT_VAR(keyField);
          for (auto i = 1u; i < nodeListPtr->numInternalNodes(); ++i) {
            ENSURE2(keyField(i) >= keyField(i - 1), nodeListPtr->name()
                    << " (" << (i - 1) << " " << i << ") (" 
                    << keyField(i-1) << " " << keyField(i) << ")");
          }
        }
      }
    END_CONTRACT_SCOPE
  }
}


//------------------------------------------------------------------------------
// Based on the bounding box and our hardwired number of cells, compute the
// step size in each dimension.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
SpaceFillingCurveRedistributeNodes<Dimension>::
computeStepSize(const pair<typename Dimension::Vector, typename Dimension::Vector>& box) const {
  // We choose the number of cells per dimension to be the maximum power of 2 such that
  // mNumGridCellsPerSide**3 < 2**64, so
  // that the indices will fit into a 64 bit Key.  This gives us 2**21.
  const Vector result = (box.second - box.first)/1048576;
//   const Vector result = (box.second - box.first)/2097152;
  ENSURE(result > 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Stitch together the set of indices and DomainNodes.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<pair<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key, DomainNode<Dimension> > >
SpaceFillingCurveRedistributeNodes<Dimension>::
buildIndex2IDPairs(const FieldList<Dimension, typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
                   const vector<DomainNode<Dimension> >& domainNodes) const {

  vector<pair<Key, DomainNode<Dimension> > > result;
  result.reserve(domainNodes.size());
  size_t i = 0;
  for (InternalNodeIterator<Dimension> nodeItr = indices.internalNodeBegin();
       nodeItr != indices.internalNodeEnd();
       ++nodeItr, ++i) {
    CHECK(i < domainNodes.size());
    result.push_back(pair<Key, DomainNode<Dimension> >(indices(nodeItr),
                                                       domainNodes[i]));
  }
  ENSURE(i == domainNodes.size());
  ENSURE(result.size() == domainNodes.size());

  // Sort locally by the hashed index.
  sort(result.begin(), result.end(), SortNodesByHashedIndex<DomainNode<Dimension> >());
  return result;
}

//------------------------------------------------------------------------------
// Find the (global) hashed index the given amount of work above the specified
// lower bound.
// We assume here that the vector<indices> is unique, and the associated 
// work array contains the work per node for each index.
//------------------------------------------------------------------------------
template<typename Dimension>
typename SpaceFillingCurveRedistributeNodes<Dimension>::Key
SpaceFillingCurveRedistributeNodes<Dimension>::
findUpperKey(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
             const vector<int>& count,
             const vector<typename Dimension::Scalar>& work,
             const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
             const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key maxUpperBound,
             const typename Dimension::Scalar workTarget,
             const int minNodes,
             const int maxNodes,
             Key& upperBound1,
             int& numNodes) const {
  REQUIRE(count.size() == indices.size());
  REQUIRE(work.size() == indices.size());
  REQUIRE(lowerBound <= maxUpperBound);
  REQUIRE(workTarget > 0.0);
  
  // Now bisect in on the upper bound.
  Key upperBound0 = lowerBound;
  upperBound1 = maxUpperBound;
  while (upperBound1 - upperBound0 > KeyTraits::one) {
    const Key trialUpperBound = (upperBound0 + upperBound1)/KeyTraits::two;
    CHECK(trialUpperBound >= upperBound0 and
          trialUpperBound <= upperBound1);
    Scalar workTrial;
    workAndNodesInRange(indices, count, work, lowerBound, trialUpperBound, numNodes, workTrial);
    // if (Process::getRank() == 0) cerr << "  --> [" 
    //                                   << upperBound0 << " " << trialUpperBound 
    //                                   << "]   "
    //                                   << workTrial << " " << workTarget 
    //                                   << endl;
    if (workTrial < workTarget) {
      if (numNodes <= maxNodes) {
        upperBound0 = trialUpperBound;
      } else {
        upperBound1 = trialUpperBound;
      }
    } else {
      if (numNodes >= minNodes) {
        upperBound1 = trialUpperBound;
      } else {
        upperBound0 = trialUpperBound;
      }
    }
  }

  // That's it.
  ENSURE(upperBound1 <= maxUpperBound);
  return upperBound1;
}

//------------------------------------------------------------------------------
// Compute the (global) number of nodes in the given range of indices.
//------------------------------------------------------------------------------
template<typename Dimension>
int
SpaceFillingCurveRedistributeNodes<Dimension>::
numIndicesInRange(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
                   const vector<int>& count,
                   const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
                   const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key upperBound) const {
  REQUIRE(lowerBound <= upperBound);

  int result = 0;
  if (indices.size() > 0) {

    // Find the positions of the requested range in our local set of indices.
    const int ilower = max(0, bisectSearch(indices, lowerBound));
    const int iupper = max(0, std::min(int(indices.size()) - 1, bisectSearch(indices, upperBound)));

    // Count up how many of the suckers we have locally.
    for (int i = ilower; i != iupper + 1; ++i) {
      if (indices[i] >= lowerBound and 
          indices[i] <= upperBound) result += count[i];
    }
  }

  // Globally reduce that sucker.
  result = allReduce(result, SPHERAL_OP_SUM);

  return result;
}

//------------------------------------------------------------------------------
// Compute the (global) work for the nodes in the given range of indices.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SpaceFillingCurveRedistributeNodes<Dimension>::
workInRange(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
            const vector<typename Dimension::Scalar>& work,
            const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
            const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key upperBound) const {
  REQUIRE(lowerBound <= upperBound);

  Scalar result = 0;

  if (indices.size() > 0) {
    // Find the positions of the requested range in our local set of indices.
    const int ilower = max(0, bisectSearch(indices, lowerBound));
    const int iupper = max(0, min(int(int(indices.size()) - 1), bisectSearch(indices, upperBound)));

    // Sum our local work.
    for (int i = ilower; i != iupper + 1; ++i) {
      if (indices[i] >= lowerBound and 
          indices[i] <= upperBound) result += work[i];
    }
  }

  // Globally reduce that sucker.
  result = allReduce(result, SPHERAL_OP_SUM);

  return result;
}

//------------------------------------------------------------------------------
// Compute the (global) work and number of nodes for the nodes in the given
//  range of indices.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
workAndNodesInRange(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
                    const vector<int>& count,
                    const vector<typename Dimension::Scalar>& work,
                    const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
                    const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key upperBound,
                    int& countInRange,
                    Scalar& workInRange) const {

  REQUIRE(lowerBound <= upperBound);

  // Find the positions of the requested range in our local set of indices.
  workInRange = 0.0;
  countInRange = 0;
  if (indices.size() > 0) {
    const int ilower = max(0, bisectSearch(indices, lowerBound));
    const int iupper = max(0, min(int(int(indices.size()) - 1), bisectSearch(indices, upperBound)));

    // Sum our local work.
    for (int i = ilower; i < iupper + 1; ++i) {
      if (indices[i] >= lowerBound and 
          indices[i] <= upperBound) {
        workInRange += work[i];
        countInRange += count[i];
      }
    }
  }

  // Globally reduce that sucker.
  workInRange = allReduce(workInRange, SPHERAL_OP_SUM);
  countInRange = allReduce(countInRange, SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Compute the number of nodes we want on the given domain.
//------------------------------------------------------------------------------
template<typename Dimension>
int
SpaceFillingCurveRedistributeNodes<Dimension>::
targetNumNodes(const int numGlobal,
               const int numProcs,
               const int targetProc) const {
  const int nlocalp = numGlobal/numProcs;
  const int nremainder = numGlobal % numProcs;
  const int result = targetProc > nremainder ? nlocalp : (nlocalp + 1);
  ENSURE(result > 0 and result <= numGlobal);
  return result;
}

//------------------------------------------------------------------------------
// Find the next unique index globally following the given one.
// Note that if there is no next index, we return the maximum value.
//------------------------------------------------------------------------------
template<typename Dimension>
typename SpaceFillingCurveRedistributeNodes<Dimension>::Key
SpaceFillingCurveRedistributeNodes<Dimension>::
findNextIndex(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
              const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key index,
              const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key maxIndex) const {

  // Find the position of the given index in our local array.
  const int inext = bisectSearch(indices, index) + 1;
  CHECK(inext >= 0 and inext <= (int)indices.size());

  // The next value on this processor.
  Key result = maxIndex;
  if (indices.size() > 0) {
    if (inext < (int)indices.size() and indices[inext] > index) result = indices[inext];
    //   cerr << "  Local result: " << result << endl;

    //   cerr << "Global result: " << result << endl;
  }

  // Get the global answer.
  result = allReduce(result, SPHERAL_OP_MIN);

  // That's it.
  ENSURE(result <= maxIndex);
  return result;
}

//------------------------------------------------------------------------------
// Find the domain for the given index given the set of index ranges for 
// each processors.
//------------------------------------------------------------------------------
template<typename Dimension>
int
SpaceFillingCurveRedistributeNodes<Dimension>::
domainForIndex(const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key index,
               const vector<pair<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key, 
                                 typename SpaceFillingCurveRedistributeNodes<Dimension>::Key> >& indexRanges) const {

  // Find that sucker by bisection.
  const int numDomains = indexRanges.size();
  int ilower = 0;
  int iupper = numDomains;
  while (iupper - ilower > 1) {
    const int imid = (iupper + ilower)/2;
    if (indexRanges[imid].first > index) {
      iupper = imid;
    } else {
      ilower = imid;
    }
  }

  // Thats it.
  ENSURE(indexRanges[ilower].first <= index and
         indexRanges[ilower].second >= index);
  return ilower;
}

//------------------------------------------------------------------------------
// Range of allowed nodes per domain as a fraction of the average number per 
// domain.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SpaceFillingCurveRedistributeNodes<Dimension>::
minNodesPerDomainFraction() const {
  return mMinNodesPerDomainFraction;
}

template<typename Dimension>
double
SpaceFillingCurveRedistributeNodes<Dimension>::
maxNodesPerDomainFraction() const {
  return mMaxNodesPerDomainFraction;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
minNodesPerDomainFraction(double x) {
  mMinNodesPerDomainFraction = x;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
maxNodesPerDomainFraction(double x) {
  mMaxNodesPerDomainFraction = x;
}

//------------------------------------------------------------------------------
// Flag for whether we should compute the work per node or strictly balance by
// node count.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpaceFillingCurveRedistributeNodes<Dimension>::
workBalance() const {
  return mWorkBalance;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
workBalance(bool val) {
  mWorkBalance = val;
}

//------------------------------------------------------------------------------
// Flag that will cause us not to repartition between domains, but only sort locally on 
// each domain.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpaceFillingCurveRedistributeNodes<Dimension>::
localReorderOnly() const {
  return mLocalReorderOnly;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
localReorderOnly(bool val) {
  mLocalReorderOnly = val;
}

}
