//---------------------------------Spheral++----------------------------------//
// nodeOrdering
//
// Compute the order that a given set of nodes should be stepped through
// given a FieldList of things to sort them by.
// The FieldList returned is the one to N indexing corresponding to sorting the 
// input in increasing order.
// 
// Created by JMO, Fri Dec 19 16:13:39 PST 2008
//----------------------------------------------------------------------------//
#include "nodeOrdering.hh"
#include "KeyTraits.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Distributed/allReduce.hh"

#include <algorithm>
#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename DataType>
struct CompareTuples {
  bool operator()(const std::tuple<int, int, DataType>& lhs,
                  const std::tuple<int, int, DataType>& rhs) {
    return std::get<2>(lhs) < std::get<2>(rhs);
  }
};

template<typename Dimension, typename DataType>
FieldList<Dimension, int>
nodeOrdering(const FieldList<Dimension, DataType>& criteria) {
  typedef KeyTraits::Key Key;

  // Prepare the result.
  FieldList<Dimension, int> result(FieldStorageType::CopyFields);

  // Parallel info.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // Make a set of tuples containing the node info and indicies.
  int iNodeList = 0;
  vector<std::tuple<int, int, DataType> > sortedList;
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = criteria.begin();
       fieldItr != criteria.end();
       ++fieldItr, ++iNodeList) {
    const NodeList<Dimension>& nodeList = (**fieldItr).nodeList();
    result.appendField(Field<Dimension, int>("node indicies", nodeList, -1));
    for (auto i = 0u; i != nodeList.numInternalNodes(); ++i) {
      sortedList.push_back(std::make_tuple(iNodeList, i, (**fieldItr)(i)));
    }
  }

  // Locally sort the list.
  std::sort(sortedList.begin(), sortedList.end(), CompareTuples<DataType>());

  // Find the total number of nodes.
  const int numLocalNodes = sortedList.size();
  int numGlobalNodes = allReduce(numLocalNodes, SPHERAL_OP_SUM);

  // Iterate over the local nodes in order.
  int iLocal = 0;
  int iGlobal = 0;
  while (iGlobal < numGlobalNodes) {
    Key localKey = KeyTraits::maxKey;
    if (iLocal < numLocalNodes) localKey = std::get<2>(sortedList[iLocal]);

    // Find the next key globally.
    Key globalKey = allReduce(localKey, SPHERAL_OP_MIN);

    // If we have the next index, check for duplicates on other domains.
    int minProcID = numProcs + 1;
    if (localKey == globalKey) minProcID = procID;
    minProcID = allReduce(minProcID, SPHERAL_OP_MIN);

    // Are we the next global key?
    if (localKey == globalKey and procID == minProcID) {
      result(std::get<0>(sortedList[iLocal]), std::get<1>(sortedList[iLocal])) = iGlobal;
      ++iLocal;
    }

    // Advance to the next global ID.
    ++iGlobal;
  }
  CHECK(iGlobal == numGlobalNodes);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    int iNodeList = 0;
    for (typename FieldList<Dimension, int>::const_iterator itr = result.begin();
         itr != result.end();
         ++itr, ++iNodeList) {
      const NodeList<Dimension>& nodeList = (**itr).nodeList();
      for (int i = 0; i != (int)nodeList.numInternalNodes(); ++i) {
        ENSURE((**itr)(i) >= 0 and (**itr)(i) < numGlobalNodes);
      }
    }
    for (int iGlobal = 0; iGlobal != numGlobalNodes; ++iGlobal) {
      int countGlobal = 0;
      for (typename FieldList<Dimension, int>::const_iterator itr = result.begin();
           itr != result.end();
           ++itr, ++iNodeList) countGlobal += std::count((**itr).internalBegin(), (**itr).internalEnd(), iGlobal);
      countGlobal = allReduce(countGlobal, SPHERAL_OP_SUM);
      ENSURE(countGlobal == 1);
    }
  }
  END_CONTRACT_SCOPE

  return result;
}

}

