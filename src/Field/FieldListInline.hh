// Includes.
#include "Geometry/MathTraits.hh"
#include "NodeIterators.hh"

#include "NodeList/FluidNodeList.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Neighbor/Neighbor.hh"
#include "Field/Field.hh"
#include "Kernel/TableKernel.hh"
#include "Distributed/allReduce.hh"

#ifdef SPHERAL_ENABLE_MPI
#include <mpi.h>
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/packElement.hh"
#include "Distributed/Communicator.hh"
#endif

#include <algorithm>
#include <limits.h>
#include <float.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::FieldList():
  FieldListBase<Dimension>(),
  mFieldPtrs(0),
  mFieldBasePtrs(0),
  mFieldCache(0),
  mStorageType(FieldStorageType::ReferenceFields),
  mNodeListPtrs(0),
  mNodeListIndexMap(),
  reductionType(ThreadReduction::SUM),
  threadMasterPtr(NULL) {
}

//------------------------------------------------------------------------------
// Construct an empty Field List, with an explict storage type.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::FieldList(FieldStorageType aStorageType):
  FieldListBase<Dimension>(),
  mFieldPtrs(0),
  mFieldBasePtrs(0),
  mFieldCache(0),
  mStorageType(aStorageType),
  mNodeListPtrs(0),
  mNodeListIndexMap(),
  reductionType(ThreadReduction::SUM),
  threadMasterPtr(NULL) {
}

//------------------------------------------------------------------------------
// Copy constructor.
// Note that the copy constructor copies the data by reference by default.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::
FieldList(const FieldList<Dimension, DataType>& rhs):
  FieldListBase<Dimension>(rhs),
  mFieldPtrs(rhs.mFieldPtrs),
  mFieldBasePtrs(rhs.mFieldBasePtrs),
  mFieldCache(),
  mStorageType(rhs.storageType()),
  mNodeListPtrs(rhs.mNodeListPtrs),
  mNodeListIndexMap(rhs.mNodeListIndexMap),
  reductionType(rhs.reductionType),
  threadMasterPtr(rhs.threadMasterPtr) {

  // If we're maintaining Fields by copy, then copy the Field cache.
#pragma omp critical (FieldList_copy)
  {
    if (storageType() == FieldStorageType::CopyFields) {
      CHECK(mFieldCache.size() == 0);
      for (typename FieldCacheType::const_iterator itr = rhs.mFieldCache.begin();
           itr != rhs.mFieldCache.end();
           ++itr) {
        auto newField = std::make_shared<Field<Dimension, DataType>>(**itr);
        mFieldCache.push_back(newField);
      }

      CHECK(this->size() == mFieldCache.size());
      auto fieldPtrItr = this->begin();
      auto fieldBasePtrItr = this->begin_base();
      auto fieldCacheItr = mFieldCache.begin();
      for(; fieldPtrItr != this->end(); ++fieldPtrItr, ++fieldBasePtrItr, ++fieldCacheItr) {
        CHECK(fieldCacheItr != mFieldCache.end());
        (*fieldPtrItr) = fieldCacheItr->get();
        (*fieldBasePtrItr) = fieldCacheItr->get();
      }

      CHECK(fieldPtrItr == this->end() &&
            fieldBasePtrItr == this->end_base() &&
            fieldCacheItr == mFieldCache.end());
    }

    NodeListRegistrar<Dimension>::sortInNodeListOrder(mFieldPtrs.begin(), mFieldPtrs.end());
    mFieldBasePtrs.clear();
    mNodeListPtrs.clear();
    for (auto fptr: mFieldPtrs) {
      mFieldBasePtrs.push_back(fptr);
      mNodeListPtrs.push_back(const_cast<NodeList<Dimension>*>(fptr->nodeListPtr()));
    }
  } // OMP critical

//   // Register this FieldList with each Field we point to.
//   for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//     registerWithField(**fieldPtrItr);

}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::~FieldList() {
  mFieldViews.free();
}

//------------------------------------------------------------------------------
// Assignment with another FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator=(const FieldList<Dimension, DataType>& rhs) {
#pragma omp critical (FieldList_assign)
  {
    if (this != &rhs) {
      mStorageType = rhs.storageType();
      mNodeListPtrs = rhs.mNodeListPtrs;
      mFieldCache = FieldCacheType();
      mNodeListIndexMap = rhs.mNodeListIndexMap;
      mFieldPtrs = std::vector<ElementType>();
      mFieldBasePtrs = std::vector<BaseElementType>();
      mFieldPtrs.reserve(rhs.size());
      mFieldBasePtrs.reserve(rhs.size());
      reductionType = rhs.reductionType;
      threadMasterPtr = rhs.threadMasterPtr;

      //     // Unregister from our current set of Fields.
      //     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
      //       unregisterFromField(**fieldPtrItr);

      switch(storageType()) {

      case FieldStorageType::ReferenceFields:
        for (auto fieldPtrItr = rhs.begin(); 
             fieldPtrItr != rhs.end();
             ++fieldPtrItr) {
          mFieldPtrs.push_back(*fieldPtrItr);
          mFieldBasePtrs.push_back(*fieldPtrItr);
        }
        break;

      case FieldStorageType::CopyFields:
        for (auto itr = rhs.mFieldCache.begin();
             itr != rhs.mFieldCache.end();
             ++itr) {
          auto newField = std::make_shared<Field<Dimension, DataType>>(**itr);
          mFieldCache.push_back(newField);
          mFieldPtrs.push_back(newField.get());
        }
        NodeListRegistrar<Dimension>::sortInNodeListOrder(mFieldPtrs.begin(), mFieldPtrs.end());
        for (auto fptr: mFieldPtrs) mFieldBasePtrs.push_back(fptr);
        CHECK(this->size() == mFieldCache.size());
        CHECK(mFieldBasePtrs.size() == mFieldCache.size());
        break;
      }

      //     // Register this FieldList with each Field we point to.
      //     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
      //       registerWithField(**fieldPtrItr);

    }
  } // OMP critical
  ENSURE(this->size() == rhs.size());
  return *this;
}

//------------------------------------------------------------------------------
// Assignment with a constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto fieldPtrItr = begin(); fieldPtrItr < end(); ++fieldPtrItr) {
    (*fieldPtrItr)->operator=(rhs);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Return the storage type.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldStorageType
FieldList<Dimension, DataType>::storageType() const {
  return mStorageType;
}

//------------------------------------------------------------------------------
// Force the FieldList to store the fields it points to.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::copyFields() {
  if (storageType() != FieldStorageType::CopyFields) {
    mStorageType = FieldStorageType::CopyFields;

//     // Unregister this FieldList from each Field we currently point to.
//     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//       unregisterFromField(**fieldPtrItr);

    // Store local copies of the fields we're pointing at.
    mFieldCache = FieldCacheType();
    auto itr = begin();
    auto baseItr = begin_base();
    for (; itr != end(); ++itr, ++baseItr) {
      auto newField = std::make_shared<Field<Dimension, DataType>>(**itr);
      mFieldCache.push_back(newField);
      *itr = mFieldCache.back().get();
      *baseItr = mFieldCache.back().get();
    }

    // Make sure the FieldPtrs are in the correct order.
    NodeListRegistrar<Dimension>::sortInNodeListOrder(mFieldPtrs.begin(), mFieldPtrs.end());
    mFieldBasePtrs.clear();
    mNodeListPtrs.clear();
    for (auto fptr: mFieldPtrs) {
      mFieldBasePtrs.push_back(fptr);
      mNodeListPtrs.push_back(const_cast<NodeList<Dimension>*>(fptr->nodeListPtr()));
    }

//     for (int i = 0; i < this->size(); ++i) {
//       mFieldCache.push_back(*((*this)[i]));
//       (*this)[i] = &(mFieldCache.back());
//     }
  }
}

//------------------------------------------------------------------------------
// Make this FieldList store copies of Fields from another.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::copyFields(const FieldList<Dimension, DataType>& fieldList) {
  mFieldPtrs.clear();
  mFieldBasePtrs.clear();
  mFieldCache.clear();
  mStorageType = FieldStorageType::CopyFields;
  mNodeListPtrs = fieldList.mNodeListPtrs;
  mNodeListIndexMap = fieldList.mNodeListIndexMap;
  reductionType = fieldList.reductionType;
  threadMasterPtr = fieldList.threadMasterPtr;

  // Store new copies of the Fields from the other FieldList
  for (const auto& fieldPtr: fieldList.mFieldPtrs) {
    auto newFieldPtr = std::make_shared<Field<Dimension, DataType>>(*fieldPtr);
    mFieldCache.push_back(newFieldPtr);
    mFieldPtrs.push_back(newFieldPtr.get());
    mFieldBasePtrs.push_back(newFieldPtr.get());
  }
}

//------------------------------------------------------------------------------
// Test if a given field is part of this field list.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
haveField(const Field<Dimension, DataType>& field) const {
  auto fieldListItr = std::find(this->begin(), this->end(), &field);
  return fieldListItr != this->end();
}

//------------------------------------------------------------------------------
// Test if there is a Field associated with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
haveNodeList(const NodeList<Dimension>& nodeList) const {
  return mNodeListIndexMap.find(&nodeList) != mNodeListIndexMap.end();
}

//------------------------------------------------------------------------------
// Set the fields of this field list equal to those of another.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
assignFields(const FieldList<Dimension, DataType>& fieldList) {
#pragma omp critical (FieldList_assignFields)
  {
    auto otherFieldListItr = fieldList.begin();
    for (auto fieldListItr = this->begin(); fieldListItr < this->end();
         ++fieldListItr, ++otherFieldListItr) {
      CHECK(otherFieldListItr < fieldList.end());
      CHECK((*fieldListItr)->nodeListPtr() == (*otherFieldListItr)->nodeListPtr());
      **fieldListItr = **otherFieldListItr;
    }
  } // OMP critical
}

//------------------------------------------------------------------------------
// Make this FieldList reference the fields of another.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
referenceFields(const FieldList<Dimension, DataType>& fieldList) {
  mFieldPtrs = fieldList.mFieldPtrs;
  mFieldBasePtrs = fieldList.mFieldBasePtrs;
  mFieldCache.clear();
  mStorageType = FieldStorageType::ReferenceFields;
  mNodeListPtrs = fieldList.mNodeListPtrs;
  mNodeListIndexMap = fieldList.mNodeListIndexMap;
}

//------------------------------------------------------------------------------
// Append the given field to the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::appendField(const Field<Dimension, DataType>& field) {
  if (haveField(field)) {
    std::cerr << "FieldList::appendField Warning: attempt to append field " << &field
              << " to FieldList " << this
              << " which already has it." << std::endl;
    return;
  }

  // Determine the order this Field should be in.
  const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
  auto orderItr = nlr.findInsertionPoint(&field,
                                         begin(),
                                         end());
  const auto delta = std::distance(begin(), orderItr);

  // Insert the field.
  switch(storageType()) {
  case FieldStorageType::ReferenceFields:
    mFieldPtrs.insert(orderItr, const_cast<Field<Dimension, DataType>*>(&field));
    mFieldBasePtrs.insert(mFieldBasePtrs.begin() + delta, const_cast<FieldBase<Dimension>*>(dynamic_cast<const FieldBase<Dimension>*>(&field)));
    break;

  case FieldStorageType::CopyFields:
    auto newField = std::make_shared<Field<Dimension, DataType>>(field);
    mFieldCache.push_back(newField);
    mFieldPtrs.insert(orderItr, newField.get());
    mFieldBasePtrs.insert(mFieldBasePtrs.begin() + delta, newField.get());
  }

//   registerWithField(*mFieldPtrs.back());

  // We also update the set of NodeListPtrs in proper order.
  NodeListRegistrar<Dimension>::sortInNodeListOrder(mFieldPtrs.begin(), mFieldPtrs.end());
  mFieldBasePtrs.clear();
  mNodeListPtrs.clear();
  for (auto fptr: mFieldPtrs) {
    mFieldBasePtrs.push_back(fptr);
    mNodeListPtrs.push_back(const_cast<NodeList<Dimension>*>(fptr->nodeListPtr()));
  }
  CHECK(mNodeListIndexMap.find(field.nodeListPtr()) == mNodeListIndexMap.end());
  buildNodeListIndexMap();

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  ENSURE(mFieldPtrs[mNodeListIndexMap[field.nodeListPtr()]] == *(fieldForNodeList(*field.nodeListPtr())));
  ENSURE(this->size() == mNodeListPtrs.size());
  END_CONTRACT_SCOPE

}

//------------------------------------------------------------------------------
// Delete the given field from the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::deleteField(const Field<Dimension, DataType>& field) {
  if (!haveField(field)) {
    std::cerr << "FieldList::deleteField Warning: attempt to delete field " << &field
              << " from FieldList " << this
              << " which does not recognize it." << std::endl;
    return;
  }

  const auto fieldPtrItr = std::find(this->begin(), this->end(), &field);
  CHECK(fieldPtrItr != this->end());
  const size_t delta = std::distance(this->begin(), fieldPtrItr);
  auto fieldItr = mFieldCache.begin();
  switch(storageType()) {
  case FieldStorageType::CopyFields:
    while (fieldItr != mFieldCache.end() && fieldItr->get() != &field) {
      ++fieldItr;
    }
    CHECK(fieldItr != mFieldCache.end());
    mFieldCache.erase(fieldItr);
    [[fallthrough]]; // C++17 for deliberate case fallthrough
  case FieldStorageType::ReferenceFields:
    mFieldPtrs.erase(fieldPtrItr);
    mFieldBasePtrs.erase(mFieldBasePtrs.begin() + delta);
    break;
  }

  // Remove the NodeList pointer.
  auto nodeListItr = std::find(mNodeListPtrs.begin(), mNodeListPtrs.end(), field.nodeListPtr());
  CHECK(nodeListItr != mNodeListPtrs.end());
  mNodeListPtrs.erase(nodeListItr);
  buildNodeListIndexMap();
}

//------------------------------------------------------------------------------
// Construct a new field and add it to the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
appendNewField(const typename Field<Dimension, DataType>::FieldName name,
               const NodeList<Dimension>& nodeList,
               const DataType value) {
  VERIFY(mStorageType == FieldStorageType::CopyFields);

  // Create the field in our cache.
  mFieldCache.push_back(std::make_shared<Field<Dimension, DataType>>(name, nodeList, value));
  Field<Dimension, DataType>* fieldPtr = mFieldCache.back().get();

  mFieldPtrs.push_back(fieldPtr);
  NodeListRegistrar<Dimension>::sortInNodeListOrder(mFieldPtrs.begin(), mFieldPtrs.end());
  mFieldBasePtrs.clear();
  mNodeListPtrs.clear();
  for (auto fptr: mFieldPtrs) {
    mFieldBasePtrs.push_back(fptr);
    mNodeListPtrs.push_back(const_cast<NodeList<Dimension>*>(fptr->nodeListPtr()));
  }

  // // Determine the order this Field should be in.
  // const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
  // auto orderItr = nlr.findInsertionPoint(fieldPtr,
  //                                            begin(),
  //                                            end());
  // const int delta = std::distance(begin(), orderItr);

  // // Insert the field.
  // mFieldPtrs.insert(orderItr, fieldPtr);
  // mFieldBasePtrs.insert(mFieldBasePtrs.begin() + delta, fieldPtr);

  // // We also update the set of NodeListPtrs in proper order.
  // mNodeListPtrs.insert(mNodeListPtrs.begin() + delta, const_cast<NodeList<Dimension>*>(&nodeList));
  CHECK(mNodeListIndexMap.find(fieldPtr->nodeListPtr()) == mNodeListIndexMap.end());
  buildNodeListIndexMap();

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  ENSURE(mFieldPtrs[mNodeListIndexMap[fieldPtr->nodeListPtr()]] == *(fieldForNodeList(*(fieldPtr->nodeListPtr()))));
  ENSURE(this->size() == mNodeListPtrs.size());
  ENSURE(mFieldBasePtrs.size() == mFieldPtrs.size());
  END_CONTRACT_SCOPE

}

//------------------------------------------------------------------------------
// Standard iterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::iterator
FieldList<Dimension, DataType>::
begin() {
  return mFieldPtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::iterator
FieldList<Dimension, DataType>::
end() {
  return mFieldPtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::reverse_iterator
FieldList<Dimension, DataType>::
rbegin() {
  return mFieldPtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::reverse_iterator
FieldList<Dimension, DataType>::
rend() {
  return mFieldPtrs.rend();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_iterator
FieldList<Dimension, DataType>::
begin() const {
  return mFieldPtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_iterator
FieldList<Dimension, DataType>::
end() const {
  return mFieldPtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_reverse_iterator
FieldList<Dimension, DataType>::
rbegin() const {
  return mFieldPtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_reverse_iterator
FieldList<Dimension, DataType>::
rend() const {
  return mFieldPtrs.rend();
}

//------------------------------------------------------------------------------
// Standard iterators for FieldListBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::iterator
FieldList<Dimension, DataType>::
begin_base() {
  return mFieldBasePtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::iterator
FieldList<Dimension, DataType>::
end_base() {
  return mFieldBasePtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::reverse_iterator
FieldList<Dimension, DataType>::
rbegin_base() {
  return mFieldBasePtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::reverse_iterator
FieldList<Dimension, DataType>::
rend_base() {
  return mFieldBasePtrs.rend();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_iterator
FieldList<Dimension, DataType>::
begin_base() const {
  return mFieldBasePtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_iterator
FieldList<Dimension, DataType>::
end_base() const {
  return mFieldBasePtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_reverse_iterator
FieldList<Dimension, DataType>::
rbegin_base() const {
  return mFieldBasePtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_reverse_iterator
FieldList<Dimension, DataType>::
rend_base() const {
  return mFieldBasePtrs.rend();
}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
operator[](const size_t index) const {
  REQUIRE2(index < mFieldPtrs.size(), "FieldList index ERROR: out of bounds " << index << " !< " << mFieldPtrs.size());
  return mFieldPtrs[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
operator[](const size_t index) {
  REQUIRE2(index < mFieldPtrs.size(), "FieldList index ERROR: out of bounds " << index << " !< " << mFieldPtrs.size());
  return mFieldPtrs[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
at(const size_t index) const {
  return (*this)[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
at(const size_t index) {
  return (*this)[index];
}

//------------------------------------------------------------------------------
// Return an iterator pointing to the Field member associated with the given
// NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::iterator
FieldList<Dimension, DataType>::
fieldForNodeList(const NodeList<Dimension>& nodeList) {
  if (haveNodeList(nodeList)) {
    return begin() + mNodeListIndexMap.find(&nodeList)->second;
  } else {
    return end();
  }
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_iterator
FieldList<Dimension, DataType>::
fieldForNodeList(const NodeList<Dimension>& nodeList) const {
  if (haveNodeList(nodeList)) {
    return begin() + mNodeListIndexMap.find(&nodeList)->second;
  } else {
    return end();
  }
}

//------------------------------------------------------------------------------
// Element access by NodeIteratorBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldList<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) {
  return this->operator()(itr.fieldID(), itr.nodeID());
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldList<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) const {
  return this->operator()(itr.fieldID(), itr.nodeID());
}

//------------------------------------------------------------------------------
// Provide a more primitive access to Field elements, based on the index of the Field
// and the node index within that field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldList<Dimension, DataType>::
operator()(const size_t fieldIndex,
           const size_t nodeIndex) {
  REQUIRE2(fieldIndex < mFieldPtrs.size(), "FieldList index ERROR: out of bounds " << fieldIndex << " !< " << mFieldPtrs.size());
  REQUIRE2(nodeIndex < mFieldPtrs[fieldIndex]->size(), "FieldList node index ERROR: out of bounds " << nodeIndex << " !< " << mFieldPtrs[fieldIndex]->size());
  return mFieldPtrs[fieldIndex]->operator()(nodeIndex);
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldList<Dimension, DataType>::
operator()(const size_t fieldIndex,
           const size_t nodeIndex) const {
  REQUIRE2(fieldIndex < mFieldPtrs.size(), "FieldList index ERROR: out of bounds " << fieldIndex << " !< " << mFieldPtrs.size());
  REQUIRE2(nodeIndex < mFieldPtrs[fieldIndex]->size(), "FieldList node index ERROR: out of bounds " << nodeIndex << " !< " << mFieldPtrs[fieldIndex]->size());
  return mFieldPtrs[fieldIndex]->operator()(nodeIndex);
}

//------------------------------------------------------------------------------
// Return the interpolated value of the FieldList for the given position.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
operator()(const typename Dimension::Vector& position,
           const TableKernel<Dimension>& W) const {

  DataType result(0.0);

  // Set the neighbor information for all NodeLists.
  std::vector<std::vector<int>> masterLists, coarseNeighbors, refineNeighbors;
  setMasterNodeLists(position, masterLists, coarseNeighbors);
  setRefineNodeLists(position, coarseNeighbors, refineNeighbors);

  // Loop over the neighbors.
  for (auto neighborItr = refineNodeBegin(refineNeighbors);
       neighborItr < refineNodeEnd();
       ++neighborItr) {

    // Ignore any ghost nodes.
    if (neighborItr.internalNode()) {

      // This is horrible, but I need the node weights to do this job.  This routine
      // needs to be rewritten!
      const Scalar& wj = (neighborItr.fluidNodeListPtr()->weight())(neighborItr);

      // Add this nodes contribution.
      const Vector rij = position - neighborItr.fluidNodeListPtr()->positions()(neighborItr);
      const SymTensor& Hj = neighborItr.fluidNodeListPtr()->Hfield()(neighborItr);
      double Wj = W(Hj*rij, Hj);
      result += (*this)(neighborItr)*wj*Wj;
    }
  }

#ifdef SPHERAL_ENABLE_MPI
  // In parallel, we need to sum up the result across all processors.
  {
    int procID, numProcs;
    MPI_Comm_rank(Communicator::communicator(), &procID);
    MPI_Comm_size(Communicator::communicator(), &numProcs);
    VERIFY(DataTypeTraits<DataType>::fixedSize());
    const size_t sizeOfElement = DataTypeTraits<DataType>::numElements()*sizeof(typename DataTypeTraits<DataType>::ElementType);
    std::vector<char> sendBuffer;
    std::vector<char> recvBuffer(sizeOfElement);
    packElement(result, sendBuffer);
    CHECK(sendBuffer.size() == sizeOfElement);
    for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
      recvBuffer = sendBuffer;
      MPI_Bcast(&(*recvBuffer.begin()), sizeOfElement, MPI_CHAR, sendProc, Communicator::communicator());
      if (procID != sendProc) {
        DataType otherValue;
        auto itr = recvBuffer.begin();
        unpackElement(otherValue, itr, recvBuffer.end());
        CHECK(itr == recvBuffer.end());
        result += otherValue;
      }
    }
  }
#endif

  return result;
}

//------------------------------------------------------------------------------
// Node ID iterators for all node IDs.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
AllNodeIterator<Dimension>
FieldList<Dimension, DataType>::nodeBegin() const {
  typename std::vector<NodeList<Dimension>*>::const_iterator 
    nodeListItr = mNodeListPtrs.begin();
  while (nodeListItr < mNodeListPtrs.end() &&
         (*nodeListItr)->numNodes() == 0) {
    ++nodeListItr;
  }
  return AllNodeIterator<Dimension>(nodeListItr,
                                    mNodeListPtrs.begin(),
                                    mNodeListPtrs.end());
}

template<typename Dimension, typename DataType>
inline
AllNodeIterator<Dimension>
FieldList<Dimension, DataType>::nodeEnd() const {
  return AllNodeIterator<Dimension>(mNodeListPtrs.end(),
                                    mNodeListPtrs.begin(),
                                    mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Node ID iterators for internal nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
InternalNodeIterator<Dimension>
FieldList<Dimension, DataType>::internalNodeBegin() const {
  typename std::vector<NodeList<Dimension>*>::const_iterator
    nodeListItr = mNodeListPtrs.begin();
  while (nodeListItr < mNodeListPtrs.end() &&
         (*nodeListItr)->numInternalNodes() == 0) {
    ++nodeListItr;
  }
  return InternalNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end());
}

template<typename Dimension, typename DataType>
inline
InternalNodeIterator<Dimension>
FieldList<Dimension, DataType>::internalNodeEnd() const {
  return InternalNodeIterator<Dimension>(mNodeListPtrs.end(),
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Node ID iterators for ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
GhostNodeIterator<Dimension>
FieldList<Dimension, DataType>::ghostNodeBegin() const {
  auto nodeListItr = mNodeListPtrs.begin();
  while (nodeListItr < mNodeListPtrs.end() &&
         (*nodeListItr)->numGhostNodes() == 0) {
    ++nodeListItr;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    CHECK((*nodeListItr)->firstGhostNode() < (*nodeListItr)->numNodes());
    return GhostNodeIterator<Dimension>(nodeListItr,
                                        mNodeListPtrs.begin(),
                                        mNodeListPtrs.end(),
                                        (*nodeListItr)->firstGhostNode());
  } else {
    return ghostNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
GhostNodeIterator<Dimension>
FieldList<Dimension, DataType>::ghostNodeEnd() const {
  return GhostNodeIterator<Dimension>(mNodeListPtrs.end(),
                                      mNodeListPtrs.begin(),
                                      mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Node ID iterators for master neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
MasterNodeIterator<Dimension>
FieldList<Dimension, DataType>::masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const {
  auto nodeListItr = mNodeListPtrs.begin();
  size_t iNodeList = 0;
  while (nodeListItr < mNodeListPtrs.end() && masterLists[iNodeList].empty()) {
    ++nodeListItr;
    ++iNodeList;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    return MasterNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end(),
                                         masterLists[iNodeList].begin(),
                                         masterLists);
  } else {
    return this->masterNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
MasterNodeIterator<Dimension>
FieldList<Dimension, DataType>::masterNodeEnd() const {
  return MasterNodeIterator<Dimension>(mNodeListPtrs.end(),
                                       mNodeListPtrs.begin(),
                                       mNodeListPtrs.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node ID iterators for coarse neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
CoarseNodeIterator<Dimension>
FieldList<Dimension, DataType>::coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const {
  auto nodeListItr = mNodeListPtrs.begin();
  size_t iNodeList = 0;
  while (nodeListItr < mNodeListPtrs.end() && coarseNeighbors[iNodeList].empty()) {
    ++nodeListItr;
    ++iNodeList;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    return CoarseNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end(),
                                         coarseNeighbors[iNodeList].begin(),
                                         coarseNeighbors);
  } else {
    return this->coarseNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
CoarseNodeIterator<Dimension>
FieldList<Dimension, DataType>::coarseNodeEnd() const {
  return CoarseNodeIterator<Dimension>(mNodeListPtrs.end(),
                                       mNodeListPtrs.begin(),
                                       mNodeListPtrs.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node ID iterators for refine neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
RefineNodeIterator<Dimension>
FieldList<Dimension, DataType>::refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const {
  auto nodeListItr = mNodeListPtrs.begin();
  size_t iNodeList = 0;
  while (nodeListItr < mNodeListPtrs.end() && refineNeighbors[iNodeList].empty()) {
    ++nodeListItr;
    ++iNodeList;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    return RefineNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end(),
                                         refineNeighbors[iNodeList].begin(),
                                         refineNeighbors);
  } else {
    return this->refineNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
RefineNodeIterator<Dimension>
FieldList<Dimension, DataType>::refineNodeEnd() const {
  return RefineNodeIterator<Dimension>(mNodeListPtrs.end(),
                                       mNodeListPtrs.begin(),
                                       mNodeListPtrs.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Set the master node lists of all the NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setMasterNodeLists(const typename Dimension::Vector& r,
                   const typename Dimension::SymTensor& H,
                   std::vector<std::vector<int>>& masterLists,
                   std::vector<std::vector<int>>& coarseNeighbors) const {
  auto etaMax = 0.0;
  for (auto nodeListItr = mNodeListPtrs.begin();
       nodeListItr != mNodeListPtrs.end();
       ++nodeListItr) etaMax = std::max(etaMax, (**nodeListItr).neighbor().kernelExtent());
  Neighbor<Dimension>::setMasterNeighborGroup(r, H,
                                              mNodeListPtrs.begin(),
                                              mNodeListPtrs.end(),
                                              etaMax,
                                              masterLists,
                                              coarseNeighbors);
}

template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setMasterNodeLists(const typename Dimension::Vector& r,
                   std::vector<std::vector<int>>& masterLists,
                   std::vector<std::vector<int>>& coarseNeighbors) const {
  this->setMasterNodeLists(r, 1e-30*SymTensor::one, masterLists, coarseNeighbors);
}

//------------------------------------------------------------------------------
// Set the refine node lists of all the NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setRefineNodeLists(const typename Dimension::Vector& r,
                   const typename Dimension::SymTensor& H,
                   const std::vector<std::vector<int>>& coarseNeighbors,
                   std::vector<std::vector<int>>& refineNeighbors) const {
  const auto numNodeLists = mNodeListPtrs.size();
  REQUIRE(coarseNeighbors.size() == numNodeLists);
  refineNeighbors = std::vector<std::vector<int>>(numNodeLists);
  auto iNodeList = 0;
  for (auto nodeListItr = mNodeListPtrs.begin();
       nodeListItr < mNodeListPtrs.end();
       ++nodeListItr, ++iNodeList) {
    (*nodeListItr)->neighbor().setRefineNeighborList(r, H, 
                                                     coarseNeighbors[iNodeList],
                                                     refineNeighbors[iNodeList]);
  }
}

template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setRefineNodeLists(const typename Dimension::Vector& r,
                   const std::vector<std::vector<int>>& coarseNeighbors,
                   std::vector<std::vector<int>>& refineNeighbors) const {
  this->setRefineNodeLists(r, 1e-30*SymTensor::one, coarseNeighbors, refineNeighbors);
}

//------------------------------------------------------------------------------
// Zero out the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::Zero() {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->Zero();
  }
}

//------------------------------------------------------------------------------
// Apply a minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyMin(dataMin);
  }
}

//------------------------------------------------------------------------------
// Apply a maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyMax(dataMax);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyScalarMin(const double dataMin) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyScalarMin(dataMin);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyScalarMax(const double dataMax) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyScalarMax(dataMax);
  }
}

//------------------------------------------------------------------------------
// Add two FieldLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator+(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != (int)numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a FieldList from another.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator-(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != (int)numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Add two FieldLists in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator+=(const FieldList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (auto i = 0u; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  for (auto i = 0u; i < numFields(); ++i) {
    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
    *((*this)[i]) += *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a FieldList from another in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator-=(const FieldList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != (int)numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  for (auto i = 0u; i < numFields(); ++i) {
    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
    *((*this)[i]) -= *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator+(const DataType& rhs) const {
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a single value from a Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator-(const DataType& rhs) const {
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Add a single value to the Field in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator+=(const DataType& rhs) {
  for (auto i = 0u; i < numFields(); ++i) {
    *((*this)[i]) += rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value from the FieldList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator-=(const DataType& rhs) {
  for (auto i = 0u; i < numFields(); ++i) {
    *((*this)[i]) -= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldList by a Scalar FieldList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator*=(const FieldList<Dimension, typename Dimension::Scalar>& rhs) {
  REQUIRE(this->numFields() == rhs.numFields());
  for (auto i = 0u; i < numFields(); ++i) {
    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
    *((*this)[i]) *= *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator*=(const Scalar& rhs) {
  for (auto i = 0u; i < numFields(); ++i) {
    *((*this)[i]) *= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::
operator/(const FieldList<Dimension, typename Dimension::Scalar>& rhs) const {
  REQUIRE(this->numFields() == rhs.numFields());
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator/(const Scalar& rhs) const {
  REQUIRE(rhs != 0.0);
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar FieldList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator/=(const FieldList<Dimension, typename Dimension::Scalar>& rhs) {
  REQUIRE(this->numFields() == rhs.numFields());
  for (auto i = 0u; i < numFields(); ++i) {
    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
    *((*this)[i]) /= *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator/=(const typename Dimension::Scalar& rhs) {
  REQUIRE(rhs != 0.0);
  for (auto i = 0u; i < numFields(); ++i) {
    *((*this)[i]) /= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Sum the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
sumElements() const {
  return allReduce(this->localSumElements(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Find the minimum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
min() const {
  return allReduce(this->localMin(), SPHERAL_OP_MIN);
}

//------------------------------------------------------------------------------
// Find the maximum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
max() const {
  return allReduce(this->localMax(), SPHERAL_OP_MAX);
}

//------------------------------------------------------------------------------
// Sum the field elements.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
localSumElements() const {
  auto result = DataTypeTraits<DataType>::zero();
  for (auto itr = begin(); itr < end(); ++itr) result += (*itr)->localSumElements();
  return result;
}

//------------------------------------------------------------------------------
// Find the minimum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
localMin() const {
  auto result = std::numeric_limits<DataType>::max();
  for (auto itr = begin(); itr != end(); ++itr) result = std::min(result, (*itr)->localMin());
  return result;
}

//------------------------------------------------------------------------------
// Find the maximum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
localMax() const {
  auto result = std::numeric_limits<DataType>::lowest();
  for (auto itr = begin(); itr < end(); ++itr) result = std::max(result, (*itr)->localMax());
  return result;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator==(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != (int)numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != (int)numFields()) {
    result = result && (*(mFieldPtrs[i]) == *(rhs[i]));
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator!=(const FieldList<Dimension, DataType>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) > *(rhs[i]));
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) < *(rhs[i]));
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>=(const FieldList<Dimension, DataType>& rhs) const {
  return operator>(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<=(const FieldList<Dimension, DataType>& rhs) const {
  return operator<(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator==(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != (int)numFields()) {
    result = result && (*(mFieldPtrs[i]) == rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != (int)numFields()) {
    result = result && (*(mFieldPtrs[i]) > rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != (int)numFields()) {
    result = result && (*(mFieldPtrs[i]) < rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  return operator>(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  return operator<(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// Return the number of Fields stored in this Field List.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t 
FieldList<Dimension, DataType>::numFields() const {
  return mFieldPtrs.size();
}

template<typename Dimension, typename DataType>
inline
size_t 
FieldList<Dimension, DataType>::size() const {
  return numFields();
}

template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::empty() const {
  return mFieldPtrs.empty();
}

template<typename Dimension, typename DataType>
inline
size_t 
FieldList<Dimension, DataType>::numElements() const {
  size_t numberOfNodes = 0;
  for (auto iter = begin(); iter != end(); ++iter) {
    numberOfNodes += (*iter)->nodeList().numNodes();
  }
  return numberOfNodes;
}

template<typename Dimension, typename DataType>
inline
size_t 
FieldList<Dimension, DataType>::numInternalElements() const {
  size_t numberOfNodes = 0;
  for (auto iter = begin(); iter != end(); ++iter) {
    numberOfNodes += (*iter)->nodeList().numInternalNodes();
  }
  return numberOfNodes;
}

template<typename Dimension, typename DataType>
inline
size_t 
FieldList<Dimension, DataType>::numGhostElements() const {
  size_t numberOfNodes = 0;
  for (auto iter = begin(); iter != end(); ++iter) {
    numberOfNodes += (*iter)->nodeList().numGhostNodes();
  }
  return numberOfNodes;
}

//------------------------------------------------------------------------------
// Return the set of NodeList pointers.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
const std::vector<NodeList<Dimension>*>&
FieldList<Dimension, DataType>::
nodeListPtrs() const {
  return mNodeListPtrs;
}

//------------------------------------------------------------------------------
// All internal values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<DataType>
FieldList<Dimension, DataType>::
internalValues() const {
  const size_t ntot = this->numInternalElements();
  std::vector<DataType> result(ntot);
  size_t offset = 0;
  for (auto itr = this->begin(); itr != this->end(); ++itr) {
    const auto n = (*itr)->numInternalElements();
    std::copy((*itr)->begin(), (*itr)->begin() + n, result.begin() + offset);
    offset += n;
  }
  CHECK(offset == ntot);
  return result;
}

//------------------------------------------------------------------------------
// All ghost values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<DataType>
FieldList<Dimension, DataType>::
ghostValues() const {
  const size_t ntot = this->numGhostElements();
  std::vector<DataType> result(ntot);
  size_t offset = 0;
  for (auto itr = this->begin(); itr != this->end(); ++itr) {
    const auto n = (*itr)->numGhostElements();
    std::copy((*itr)->begin(), (*itr)->begin() + n, result.begin() + offset);
    offset += n;
  }
  CHECK(offset == ntot);
  return result;
}

//------------------------------------------------------------------------------
// All values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<DataType>
FieldList<Dimension, DataType>::
allValues() const {
  const size_t ntot = this->numElements();
  std::vector<DataType> result(ntot);
  size_t offset = 0;
  for (auto itr = this->begin(); itr != this->end(); ++itr) {
    const auto n = (*itr)->numElements();
    std::copy((*itr)->begin(), (*itr)->begin() + n, result.begin() + offset);
    offset += n;
  }
  CHECK(offset == ntot);
  return result;
}

//------------------------------------------------------------------------------
// Internal method to build the NodeListIndexMap from scratch.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
buildNodeListIndexMap() {
  mNodeListIndexMap = HashMapType();
  int i = 0;
  for (auto itr = begin();
       itr != end();
       ++itr, ++i) mNodeListIndexMap[(*itr)->nodeListPtr()] = i;
  ENSURE(mNodeListIndexMap.size() == numFields());
}

//------------------------------------------------------------------------------
// Make a thread local copy of the FieldList -- assumes we're already in a
// threaded region.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::
threadCopy(const ThreadReduction reductionType,
           const bool copy) {
  FieldList<Dimension, DataType> result;
#pragma omp critical (FieldList_threadCopy)
  {
    if (omp_get_num_threads() == 1) {

      // In serial we can skip all the work copying
      result.referenceFields(*this);

    } else if (copy or
               reductionType == ThreadReduction::MIN or
               reductionType == ThreadReduction::MAX) {

      // For min/max operations, we need to copy the original data
      result.copyFields(*this);

    } else {    

      // Otherwise make standalone Fields of zeros
      result = FieldList<Dimension, DataType>(FieldStorageType::CopyFields);
      for (auto fitr = this->begin(); fitr < this->end(); ++fitr) result.appendNewField((*fitr)->name(),
                                                                                        (*fitr)->nodeList(),
                                                                                        DataTypeTraits<DataType>::zero());
    }
    result.reductionType = reductionType;
    result.threadMasterPtr = this;
  }
  return result;
}

template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::
threadCopy(typename SpheralThreads<Dimension>::FieldListStack& stack,
           const ThreadReduction reductionType,
           const bool copy) {
  FieldList<Dimension, DataType> result = this->threadCopy(reductionType, copy);
  stack.push_back(&result);
  return result;
}

//------------------------------------------------------------------------------
// Reduce the values in the FieldList with the passed thread-local values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
threadReduce() const {
  REQUIRE(threadMasterPtr != NULL);
  REQUIRE(threadMasterPtr->size() == this->size());
  if (omp_get_num_threads() > 1) {
    const auto numNL = this->size();
    for (auto k = 0u; k < numNL; ++k) {
      const auto n = mFieldPtrs[k]->numInternalElements();
      for (auto i = 0u; i < n; ++i) {

        switch (reductionType) {

        case ThreadReduction::SUM:
          (*threadMasterPtr)(k,i) += (*this)(k,i);
          break;

        case ThreadReduction::MIN:
          (*threadMasterPtr)(k,i) = std::min((*threadMasterPtr)(k,i), (*this)(k,i));
          break;

        case ThreadReduction::MAX:
          (*threadMasterPtr)(k,i) = std::max((*threadMasterPtr)(k,i), (*this)(k,i));
          break;

        }
      }
    }
  }
}

//****************************** Global Functions ******************************

//------------------------------------------------------------------------------
// Multiply a FieldList by another FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
inline
FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const FieldList<Dimension, DataType>& lhs,
          const FieldList<Dimension, OtherDataType>& rhs) {
  REQUIRE(lhs.numFields() == rhs.numFields());
  FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
  result.copyFields();
  for (auto i = 0u; i < lhs.numFields(); ++i) {
    CHECK2(lhs[i]->nodeListPtr() == rhs[i]->nodeListPtr(), lhs[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
    result.appendField((*(lhs[i])) * (*(rhs[i])));
  }
  return result;
}

//------------------------------------------------------------------------------
// Multiply a FieldList by a single value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
inline
FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const FieldList<Dimension, DataType>& lhs,
          const OtherDataType& rhs) {
  FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
  result.copyFields();
  for (int i = 0; i < lhs.numFields(); ++i) {
    result.appendField((*(lhs[i])) * rhs);
  }
  return result;
}

template<typename Dimension, typename DataType, typename OtherDataType>
inline
FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const DataType& lhs,
          const FieldList<Dimension, OtherDataType>& rhs) {
  FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
  result.copyFields();
  for (auto i = 0u; i < rhs.numFields(); ++i) {
    result.appendField(lhs * (*(rhs[i])));
  }
  return result;
}

//------------------------------------------------------------------------------
// View functions.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>
FieldList<Dimension, DataType>::toView() {
  auto func = [](
      const chai::PointerRecord *,
      chai::Action,
      chai::ExecutionSpace) {};

  return this->toView(func, func);
}

template<typename Dimension, typename DataType>
template<typename FL>
inline
FieldListView<Dimension, DataType>
FieldList<Dimension, DataType>::toView(FL&& extension) {
  auto func = [](
      const chai::PointerRecord *,
      chai::Action,
      chai::ExecutionSpace) {};

  return this->toView(std::forward<FL>(extension), func);
}

template<typename Dimension, typename DataType>
template<typename FL, typename F>
inline
FieldListView<Dimension, DataType>
FieldList<Dimension, DataType>::toView(FL&& extension, F&& field_extension) {
  auto callback = getFieldListCallback(std::forward<FL>(extension));

  if (mFieldViews.size() == 0 && !mFieldViews.data(chai::CPU, false)) {
    mFieldViews.allocate(size(), chai::CPU, callback);
  } else {
    mFieldViews.setUserCallback(callback);
    mFieldViews.reallocate(size());
  }

  for (size_t i = 0; i < size(); ++i) {
    mFieldViews[i] = mFieldPtrs[i]->toView(std::forward<F>(field_extension));
  }

  mFieldViews.registerTouch(chai::CPU);

  return ViewType(mFieldViews);
}

template<typename Dimension, typename DataType>
template<typename F>
inline
auto FieldList<Dimension, DataType>::getFieldListCallback(F callback) {
  return [callback](
    const chai::PointerRecord * record,
    chai::Action action,
    chai::ExecutionSpace space) {
      if (action == chai::ACTION_MOVE) {
        if (space == chai::CPU)
          DEBUG_LOG << "FieldList : MOVED to the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "FieldList : MOVED to the GPU";
      } else if (action == chai::ACTION_ALLOC) {
        if (space == chai::CPU)
          DEBUG_LOG << "FieldList : ALLOC on the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "FieldList : ALLOC on the GPU";
      } else if (action == chai::ACTION_FREE) {
        if (space == chai::CPU)
          DEBUG_LOG << "FieldList : FREE on the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "FieldList : FREE on the GPU";
      }
      callback(record, action, space);
    };
}

// //------------------------------------------------------------------------------
// // Absolute value.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// abs(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(abs(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse cosine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// acos(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(acos(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse sine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// asin(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(asin(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse tangent.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// atan(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(atan(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse tangent2.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// atan2(const FieldList<Dimension, typename Dimension::Scalar>& fieldList1,
//       const FieldList<Dimension, typename Dimension::Scalar>& fieldList2) {
//   REQUIRE(fieldList1.numFields() == fieldList2.numFields());
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList1.numFields(); ++i) {
//     result.appendField(atan2(*(fieldList1[i]), *(fieldList2[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Ceiling.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// ceil(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(ceil(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Cosine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// cos(FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(cos(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Hyperbolic cosine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// cosh(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(cosh(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Exponential.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// exp(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(exp(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Absolute value.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// fabs(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(fabs(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Floor.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// floor(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(floor(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Natural logarithm.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// log(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(log(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Log base 10.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// log10(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(log10(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // powN.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow2(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow2(*(fieldList[i])));
//   }
//   return result;
// }


// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow3(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow3(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow4(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow4(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow5(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow5(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow6(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow6(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow7(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow7(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow8(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow8(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Sine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sin(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sin(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Hyperbolic sine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sinh(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sinh(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Square.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sqr(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sqr(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Square root.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sqrt(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sqrt(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Minimum.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// min(const FieldList<Dimension, DataType>& fieldList1,
//     const FieldList<Dimension, DataType>& fieldList2) {
//   REQUIRE(fieldList1.numFields() == fieldList2.numFields());
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList1.numFields(); ++i) {
//     result.appendField(min(*(fieldList1[i]), *(fieldList2[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// min(const DataType& value,
//     const FieldList<Dimension, DataType>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(min(value, *(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// min(const FieldList<Dimension, DataType>& fieldList,
//     const DataType& value) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(min(*(fieldList[i]), value));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Maximum.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// max(const FieldList<Dimension, DataType>& fieldList1,
//     const FieldList<Dimension, DataType>& fieldList2) {
//   REQUIRE(fieldList1.numFields() == fieldList2.numFields());
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList1.numFields(); ++i) {
//     result.appendField(max(*(fieldList1[i]), *(fieldList2[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// max(const DataType& value,
//     const FieldList<Dimension, DataType>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(max(value, *(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// max(const FieldList<Dimension, DataType>& fieldList,
//     const DataType& value) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(max(*(fieldList[i]), value));
//   }
//   return result;
// }

}
