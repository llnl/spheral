#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Get the current connectivity map.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const ConnectivityMap<Dimension>&
DataBase<Dimension>::
connectivityMap() const {
  VERIFY2(mConnectivityMapPtr != 0,
          "DataBase::connectivityMap ERROR -- need to ensure ConnectivityMap is constructed before calling this method.");
  return *mConnectivityMapPtr;
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConnectivityMapPtr
DataBase<Dimension>::
connectivityMapPtr() const {
  return mConnectivityMapPtr;
}

template<typename Dimension>
inline
const ConnectivityMap<Dimension>&
DataBase<Dimension>::
connectivityMap(const bool computeGhostConnectivity,
                const bool computeOverlapConnectivity,
                const bool computeIntersectionConnectivity) const {
  if (not mConnectivityMapPtr) this->updateConnectivityMap(computeGhostConnectivity, computeOverlapConnectivity, computeIntersectionConnectivity);
  return *mConnectivityMapPtr;
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConnectivityMapPtr
DataBase<Dimension>::
connectivityMapPtr(const bool computeGhostConnectivity,
                   const bool computeOverlapConnectivity,
                   const bool computeIntersectionConnectivity) const {
  if (not mConnectivityMapPtr) this->updateConnectivityMap(computeGhostConnectivity, computeOverlapConnectivity, computeIntersectionConnectivity);
  return mConnectivityMapPtr;
}

//------------------------------------------------------------------------------
// Convenience method to construct a new FieldList with a Field for every
// NodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
FieldList<Dimension, DataType>
DataBase<Dimension>::
newGlobalFieldList(const DataType value,
                   const typename Field<Dimension, DataType>::FieldName name) const {
  FieldList<Dimension, DataType> result(FieldStorageType::CopyFields);
  for (auto* nptr: mNodeListPtrs) result.appendNewField(name, *nptr, value);
  ENSURE(result.numFields() == numNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Convenience method to construct a new FieldList with a Field for every
// FluidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
FieldList<Dimension, DataType>
DataBase<Dimension>::
newFluidFieldList(const DataType value,
                  const typename Field<Dimension, DataType>::FieldName name) const {
  FieldList<Dimension, DataType> result(FieldStorageType::CopyFields);
  for (auto* nptr: mFluidNodeListPtrs) result.appendNewField(name, *nptr, value);
  ENSURE(result.numFields() == numFluidNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Convenience method to construct a new FieldList with a Field for every
// SolidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
FieldList<Dimension, DataType>
DataBase<Dimension>::
newSolidFieldList(const DataType value,
                  const typename Field<Dimension, DataType>::FieldName name) const {
  FieldList<Dimension, DataType> result(FieldStorageType::CopyFields);
  for (auto* nptr: mSolidNodeListPtrs) result.appendNewField(name, *nptr, value);
  ENSURE(result.numFields() == numSolidNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Convenience method to construct a new FieldList with a Field for every
// DEMNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
FieldList<Dimension, DataType>
DataBase<Dimension>::
newDEMFieldList(const DataType value,
                  const typename Field<Dimension, DataType>::FieldName name) const {
  FieldList<Dimension, DataType> result(FieldStorageType::CopyFields);
  for (auto* nptr: mDEMNodeListPtrs) result.appendNewField(name, *nptr, value);
  ENSURE(result.numFields() == numDEMNodeLists());
  return result;
}


//------------------------------------------------------------------------------
// Convenience method to resize a FieldList such that it has a Field for every
// NodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeGlobalFieldList(FieldList<Dimension, DataType>& fieldList,
                      const DataType value,
                      const typename Field<Dimension, DataType>::FieldName name,
                      const bool resetValues) const {
  VERIFY((fieldList.storageType() == FieldStorageType::CopyFields));

  // First check if it's necessary to resize the FieldList.
  bool reinitialize = fieldList.numFields() != numNodeLists();
  ConstNodeListIterator nodeListItr = nodeListBegin();
  typename FieldList<Dimension, DataType>::const_iterator itr = fieldList.begin();
  while (!reinitialize && 
         nodeListItr != nodeListEnd() &&
         itr != fieldList.end()) {
    reinitialize = (*itr)->nodeListPtr() != *nodeListItr;
    ++nodeListItr;
    ++itr;
  }

  if (reinitialize) {
    fieldList = FieldList<Dimension, DataType>(fieldList.storageType());
    for (ConstNodeListIterator nodeListItr = nodeListBegin();
         nodeListItr != nodeListEnd();
         ++nodeListItr) fieldList.appendNewField(name, **nodeListItr, value);
  } else if (resetValues) {
    fieldList = value;
  }

  ENSURE(fieldList.numFields() == numNodeLists());
}

//------------------------------------------------------------------------------
// Convenience method to resize a FieldList such that it has a Field for every
// FluidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeFluidFieldList(FieldList<Dimension, DataType>& fieldList,
                     const DataType value,
                     const typename Field<Dimension, DataType>::FieldName name,
                     const bool resetValues) const {
  VERIFY((fieldList.storageType() == FieldStorageType::CopyFields));

  // First check if it's necessary to resize the FieldList.
  bool reinitialize = fieldList.numFields() != numFluidNodeLists();
  ConstFluidNodeListIterator nodeListItr = fluidNodeListBegin();
  typename FieldList<Dimension, DataType>::const_iterator itr = fieldList.begin();
  while (!reinitialize && 
         nodeListItr != fluidNodeListEnd() &&
         itr != fieldList.end()) {
    reinitialize = (*itr)->nodeListPtr() != *nodeListItr;
    ++nodeListItr;
    ++itr;
  }

  if (reinitialize) {
    fieldList = FieldList<Dimension, DataType>(fieldList.storageType());
    for (ConstFluidNodeListIterator nodeListItr = fluidNodeListBegin();
         nodeListItr != fluidNodeListEnd();
         ++nodeListItr) fieldList.appendNewField(name, **nodeListItr, value);
  } else if (resetValues) {
    fieldList = value;
  }

  ENSURE(fieldList.numFields() == numFluidNodeLists());
}


//------------------------------------------------------------------------------
// Convenience method to resize a FieldList such that it has a Field for every
// DEMNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeDEMFieldList(FieldList<Dimension, DataType>& fieldList,
                     const DataType value,
                     const typename Field<Dimension, DataType>::FieldName name,
                     const bool resetValues) const {
  VERIFY((fieldList.storageType() == FieldStorageType::CopyFields));

  // First check if it's necessary to resize the FieldList.
  bool reinitialize = fieldList.numFields() != numDEMNodeLists();
  ConstDEMNodeListIterator nodeListItr = DEMNodeListBegin();
  typename FieldList<Dimension, DataType>::const_iterator itr = fieldList.begin();
  while (!reinitialize && 
         nodeListItr != DEMNodeListEnd() &&
         itr != fieldList.end()) {
    reinitialize = (*itr)->nodeListPtr() != *nodeListItr;
    ++nodeListItr;
    ++itr;
  }

  if (reinitialize) {
    fieldList = FieldList<Dimension, DataType>(fieldList.storageType());
    for (ConstDEMNodeListIterator nodeListItr = DEMNodeListBegin();
         nodeListItr != DEMNodeListEnd();
         ++nodeListItr) fieldList.appendNewField(name, **nodeListItr, value);
  } else if (resetValues) {
    fieldList = value;
  }

  ENSURE(fieldList.numFields() == numDEMNodeLists());
}


//------------------------------------------------------------------------------
// Convenience method to resize a FieldList such that it has a Field for every
// SolidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeSolidFieldList(FieldList<Dimension, DataType>& fieldList,
                     const DataType value,
                     const typename Field<Dimension, DataType>::FieldName name,
                     const bool resetValues) const {
  VERIFY((fieldList.storageType() == FieldStorageType::CopyFields));

  // First check if it's necessary to resize the FieldList.
  bool reinitialize = fieldList.numFields() != numSolidNodeLists();
  ConstSolidNodeListIterator nodeListItr = solidNodeListBegin();
  typename FieldList<Dimension, DataType>::const_iterator itr = fieldList.begin();
  while (!reinitialize && 
         nodeListItr != solidNodeListEnd() &&
         itr != fieldList.end()) {
    reinitialize = (*itr)->nodeListPtr() != *nodeListItr;
    ++nodeListItr;
    ++itr;
  }

  if (reinitialize) {
    fieldList = FieldList<Dimension, DataType>(fieldList.storageType());
    for (ConstSolidNodeListIterator nodeListItr = solidNodeListBegin();
         nodeListItr != solidNodeListEnd();
         ++nodeListItr) fieldList.appendNewField(name, **nodeListItr, value);
  } else if (resetValues) {
    fieldList = value;
  }

  ENSURE(fieldList.numFields() == numSolidNodeLists());
}

//------------------------------------------------------------------------------
// Construct a new array<array> with an array for every NodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
std::vector<std::vector<DataType>>
DataBase<Dimension>::
newGlobalArray(const DataType value) const {
  std::vector<std::vector<DataType>> result;
  for (auto nodeListItr = nodeListBegin();
       nodeListItr != nodeListEnd();
       ++nodeListItr) {
    result.push_back(std::vector<DataType>((*nodeListItr)->numNodes(), value));
  }

  ENSURE(result.size() == numNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Construct a new array<array> with an array for every FluidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
std::vector<std::vector<DataType>>
DataBase<Dimension>::
newFluidArray(const DataType value) const {
  std::vector<std::vector<DataType>> result;
  for (auto nodeListItr = fluidNodeListBegin();
       nodeListItr != fluidNodeListEnd();
       ++nodeListItr) {
    result.push_back(std::vector<DataType>((*nodeListItr)->numNodes(), value));
  }

  ENSURE(result.size() == numNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Construct a new array<array> with an array for every SolidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
std::vector<std::vector<DataType>>
DataBase<Dimension>::
newSolidArray(const DataType value) const {
  std::vector<std::vector<DataType>> result;
  for (auto nodeListItr = solidNodeListBegin();
       nodeListItr != solidNodeListEnd();
       ++nodeListItr) {
    result.push_back(std::vector<DataType>((*nodeListItr)->numNodes(), value));
  }

  ENSURE(result.size() == numNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Construct a new array<array> with an array for every DEMNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
std::vector<std::vector<DataType>>
DataBase<Dimension>::
newDEMArray(const DataType value) const {
  std::vector<std::vector<DataType>> result;
  for (auto nodeListItr = DEMNodeListBegin();
       nodeListItr != DEMNodeListEnd();
       ++nodeListItr) {
    result.push_back(std::vector<DataType>((*nodeListItr)->numNodes(), value));
  }

  ENSURE(result.size() == numNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Resize an array<array> for every NodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeGlobalArray(std::vector<std::vector<DataType>>& array,
                  const DataType value,
                  const bool resetValues) const {
  array.resize(this->numNodeLists());
  auto k = 0;
  for (auto nodeListItr = nodeListBegin();
       nodeListItr != nodeListEnd();
       ++nodeListItr, ++k) {
    if (resetValues) {
      array[k] = std::vector<DataType>((*nodeListItr)->numNodes(), value);
    } else {
      array[k].resize((*nodeListItr)->numNodes(), value);
    }
  }

  ENSURE(array.size() == numNodeLists());
}

//------------------------------------------------------------------------------
// Resize an array<array> for every FluidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeFluidArray(std::vector<std::vector<DataType>>& array,
                 const DataType value,
                 const bool resetValues) const {
  array.resize(this->numFluidNodeLists());
  auto k = 0;
  for (auto nodeListItr = fluidNodeListBegin();
       nodeListItr != fluidNodeListEnd();
       ++nodeListItr, ++k) {
    if (resetValues) {
      array[k] = std::vector<DataType>((*nodeListItr)->numNodes(), value);
    } else {
      array[k].resize((*nodeListItr)->numNodes(), value);
    }
  }

  ENSURE(array.size() == numNodeLists());
}

//------------------------------------------------------------------------------
// Resize an array<array> for every SolidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeSolidArray(std::vector<std::vector<DataType>>& array,
                 const DataType value,
                 const bool resetValues) const {
  array.resize(this->numSolidNodeLists());
  auto k = 0;
  for (auto nodeListItr = solidNodeListBegin();
       nodeListItr != solidNodeListEnd();
       ++nodeListItr, ++k) {
    if (resetValues) {
      array[k] = std::vector<DataType>((*nodeListItr)->numNodes(), value);
    } else {
      array[k].resize((*nodeListItr)->numNodes(), value);
    }
  }

  ENSURE(array.size() == numNodeLists());
}

//------------------------------------------------------------------------------
// Resize an array<array> for every DEMNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeDEMArray(std::vector<std::vector<DataType>>& array,
                 const DataType value,
                 const bool resetValues) const {
  array.resize(this->numDEMNodeLists());
  auto k = 0;
  for (auto nodeListItr = DEMNodeListBegin();
       nodeListItr != DEMNodeListEnd();
       ++nodeListItr, ++k) {
    if (resetValues) {
      array[k] = std::vector<DataType>((*nodeListItr)->numNodes(), value);
    } else {
      array[k].resize((*nodeListItr)->numNodes(), value);
    }
  }

  ENSURE(array.size() == numNodeLists());
}


}
