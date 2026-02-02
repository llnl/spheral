//---------------------------------Spheral++----------------------------------//
// FieldListThreadWrapper
// 
// A lightweight thread helper for the FieldList.  Only implements the
// operator()(nodeList, i) indexing operation.
//
// Created by JMO, Mon Sep 30 16:03:30 PDT 2019
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldSpace__FieldListThreadWrapper_hh__
#define __Spheral__FieldSpace__FieldListThreadWrapper_hh__

#include "Field/FieldList.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/Hashes.hh"

#include <unordered_map>

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldListThreadWrapper {
public:
  //--------------------------- Public Interface ---------------------------//
  FieldListThreadWrapper(FieldList<Dimension, DataType>& fieldList,
                         const bool useFieldList):
    mFieldListPtr(&fieldList),
    mUseFieldList(useFieldList) {}

  FieldListThreadWrapper() = delete;

  // Non-const indexing
  DataType& operator()(const unsigned fieldIndex,
                       const unsigned nodeIndex) {
    if (mUseFieldList) {
      return (*mFieldListPtr)(fieldIndex, nodeIndex);
    } else {
      const auto key = std::make_pair(fieldIndex, nodeIndex);
      const auto itr = mValues.find(key);
      if (itr == mValues.end()) mValues[key] = DataTypeTraits<DataType>::zero();
      return mValues[key];
    }
  }

  // const indexing
  const DataType& operator()(const unsigned fieldIndex,
                             const unsigned nodeIndex) const {
    return const_cast<FieldListThreadWrapper<Dimension, DataType>&>(*this).operator()(fieldIndex, nodeIndex);
  }

private:
  //--------------------------- Private Interface ---------------------------//
  FieldList<Dimension, DataType>* mFieldListPtr;
  bool mUseFieldList;
  typedef std::unordered_map<std::pair<size_t, size_t>, DataType> StorageType;
  StorageType mValues;
};

}

#endif

