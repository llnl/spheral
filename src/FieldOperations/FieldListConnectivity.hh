#ifndef __Spheral_FieldListConnectivity__
#define __Spheral_FieldListConnectivity__

#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

#include <limits>
#include <utility>
#include <vector>

namespace Spheral {

// Host-side connectivity state shared by the legacy FieldOperations.  It
// gives those routines an explicit mapping between FieldList and DataBase
// indices without restoring the old per-NodeList Neighbor interface.
template<typename Dimension, typename DataType>
class FieldListConnectivity {
public:
  FieldListConnectivity(const FieldList<Dimension, DataType>& fields,
                        const typename Dimension::Scalar kernelExtent,
                        const bool ghostConnectivity = false):
    mDataBase(kernelExtent),
    mDataBaseToField(fields.numFields(), invalidIndex),
    mFieldToDataBase(fields.numFields(), invalidIndex) {
    for (auto* nodeListPtr: fields.nodeListPtrs()) mDataBase.appendNodeList(*nodeListPtr);
    mDataBase.reinitializeNeighbors();
    mDataBase.updateConnectivityMap(ghostConnectivity, false, false);
    for (auto fieldi = 0u; fieldi < fields.numFields(); ++fieldi) {
      const auto dataBasei = mDataBase.nodeListIndex(fields[fieldi]->nodeListPtr());
      REQUIRE(dataBasei < fields.numFields());
      mFieldToDataBase[fieldi] = dataBasei;
      mDataBaseToField[dataBasei] = fieldi;
    }
  }

  const ConnectivityMap<Dimension>& connectivityMap() const {
    return mDataBase.connectivityMap();
  }

  NeighborSearchType searchType() const {
    return visitNeighbor(mDataBase.neighbor(),
                         [](const auto& neighbor) { return neighbor.neighborSearchType(); });
  }

  size_t fieldIndex(const size_t dataBaseIndex) const {
    REQUIRE(dataBaseIndex < mDataBaseToField.size());
    REQUIRE(mDataBaseToField[dataBaseIndex] != invalidIndex);
    return mDataBaseToField[dataBaseIndex];
  }

  size_t dataBaseIndex(const size_t fieldIndex) const {
    REQUIRE(fieldIndex < mFieldToDataBase.size());
    REQUIRE(mFieldToDataBase[fieldIndex] != invalidIndex);
    return mFieldToDataBase[fieldIndex];
  }

  template<typename Function>
  void forEachNeighbor(const size_t fieldIndex,
                       const size_t nodeIndex,
                       Function&& function) const {
    const auto& neighbors = connectivityMap().connectivityForNode(dataBaseIndex(fieldIndex), nodeIndex);
    REQUIRE(neighbors.size() == mDataBaseToField.size());
    for (auto dataBasej = 0u; dataBasej < neighbors.size(); ++dataBasej) {
      const auto fieldj = this->fieldIndex(dataBasej);
      for (const auto j: neighbors[dataBasej]) function(fieldj, j);
    }
  }

private:
  static constexpr size_t invalidIndex = std::numeric_limits<size_t>::max();
  DataBase<Dimension> mDataBase;
  std::vector<size_t> mDataBaseToField;
  std::vector<size_t> mFieldToDataBase;
};

// Connectivity over the union of two FieldLists.  Each NodeList can belong to
// either or both roles, which is needed by sampling and splatting operations.
template<typename Dimension, typename FirstType, typename SecondType>
class FieldListPairConnectivity {
public:
  FieldListPairConnectivity(const FieldList<Dimension, FirstType>& first,
                            const FieldList<Dimension, SecondType>& second,
                            const typename Dimension::Scalar kernelExtent,
                            const bool ghostConnectivity = false):
    mDataBase(kernelExtent),
    mFirstToDataBase(first.numFields(), invalidIndex),
    mSecondToDataBase(second.numFields(), invalidIndex) {
    for (auto* nodeListPtr: first.nodeListPtrs()) mDataBase.appendNodeList(*nodeListPtr);
    for (auto* nodeListPtr: second.nodeListPtrs()) {
      if (mDataBase.nodeListIndex(nodeListPtr) == mDataBase.numNodeLists()) mDataBase.appendNodeList(*nodeListPtr);
    }
    mFirstForDataBase.resize(mDataBase.numNodeLists(), invalidIndex);
    mSecondForDataBase.resize(mDataBase.numNodeLists(), invalidIndex);
    for (auto fieldi = 0u; fieldi < first.numFields(); ++fieldi) {
      const auto dbi = mDataBase.nodeListIndex(first[fieldi]->nodeListPtr());
      mFirstToDataBase[fieldi] = dbi;
      mFirstForDataBase[dbi] = fieldi;
    }
    for (auto fieldi = 0u; fieldi < second.numFields(); ++fieldi) {
      const auto dbi = mDataBase.nodeListIndex(second[fieldi]->nodeListPtr());
      mSecondToDataBase[fieldi] = dbi;
      mSecondForDataBase[dbi] = fieldi;
    }
    mDataBase.reinitializeNeighbors();
    mDataBase.updateConnectivityMap(ghostConnectivity, false, false);
  }

  NeighborSearchType searchType() const {
    return visitNeighbor(mDataBase.neighbor(),
                         [](const auto& neighbor) { return neighbor.neighborSearchType(); });
  }

  template<typename Function>
  void forEachSecondNeighbor(const size_t firstField,
                             const size_t nodeIndex,
                             Function&& function) const {
    const auto dbi = mFirstToDataBase[firstField];
    if (mSecondForDataBase[dbi] != invalidIndex) function(mSecondForDataBase[dbi], nodeIndex);
    this->forEachRoleNeighbor(dbi, nodeIndex, mSecondForDataBase,
                              std::forward<Function>(function));
  }

  template<typename Function>
  void forEachFirstNeighbor(const size_t secondField,
                            const size_t nodeIndex,
                            Function&& function) const {
    const auto dbi = mSecondToDataBase[secondField];
    if (mFirstForDataBase[dbi] != invalidIndex) function(mFirstForDataBase[dbi], nodeIndex);
    this->forEachRoleNeighbor(dbi, nodeIndex, mFirstForDataBase,
                              std::forward<Function>(function));
  }

private:
  template<typename Function>
  void forEachRoleNeighbor(const size_t dataBaseIndex,
                           const size_t nodeIndex,
                           const std::vector<size_t>& roleForDataBase,
                           Function&& function) const {
    const auto& neighbors = mDataBase.connectivityMap().connectivityForNode(dataBaseIndex, nodeIndex);
    for (auto dbj = 0u; dbj < neighbors.size(); ++dbj) {
      if (roleForDataBase[dbj] != invalidIndex) {
        for (const auto j: neighbors[dbj]) function(roleForDataBase[dbj], j);
      }
    }
  }

  static constexpr size_t invalidIndex = std::numeric_limits<size_t>::max();
  DataBase<Dimension> mDataBase;
  std::vector<size_t> mFirstToDataBase, mSecondToDataBase;
  std::vector<size_t> mFirstForDataBase, mSecondForDataBase;
};

}

#endif
