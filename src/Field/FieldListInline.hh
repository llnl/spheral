// Includes.
#include "Geometry/MathTraits.hh"
#include "NodeIterators.hh"

#include "NodeList/FluidNodeList.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Neighbor/Neighbor.hh"
#include "Field/Field.hh"
#include "Field/FieldView.hh"
#include "Field/FieldListView.hh"
#include "Kernel/TableKernel.hh"
#include "Distributed/allReduce.hh"

#ifdef USE_MPI
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
  FieldListView<Dimension, DataType>(),
  mFieldPtrs(),
  mFieldBasePtrs(),
  mFieldCache(),
  mStorageType(FieldStorageType::ReferenceFields),
  mNodeListPtrs(),
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
  FieldListView<Dimension, DataType>(),
  mFieldPtrs(),
  mFieldBasePtrs(),
  mFieldCache(),
  mStorageType(aStorageType),
  mNodeListPtrs(),
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
  FieldListView<Dimension, DataType>(),
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
      CHECK(mFieldCache.size() == 0u);
      for (auto& fptr: rhs.mFieldCache) mFieldCache.push_back(std::make_shared<Field<Dimension, DataType>>(*fptr));
      CHECK(mFieldPtrs.size() == mFieldCache.size());
      mFieldPtrs.clear();
      for (auto& fcache: mFieldCache) mFieldPtrs.push_back(fcache.get());
    }
    buildDependentArrays();
  } // OMP critical


//   // Register this FieldList with each Field we point to.
//   for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//     registerWithField(**fieldPtrItr);

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
FieldList<Dimension, DataType>::
~FieldList() {
#ifndef SPHERAL_UNIFIED_MEMORY
  std::cerr << "FIELDLIST::~FieldList()" << std::endl;
  mSpanFieldViews.free();
  std::cerr << " --> SUCCESS" << std::endl;
#endif
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
      mFieldPtrs = std::vector<ElementType>();
      mFieldCache = FieldCacheType();
      mStorageType = rhs.storageType();
      mNodeListIndexMap = rhs.mNodeListIndexMap;
      mFieldPtrs.reserve(rhs.size());
      mFieldBasePtrs.reserve(rhs.size());
      reductionType = rhs.reductionType;
      threadMasterPtr = rhs.threadMasterPtr;

      //     // Unregister from our current set of Fields.
      //     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
      //       unregisterFromField(**fieldPtrItr);

      switch(storageType()) {

      case FieldStorageType::ReferenceFields:
        mFieldPtrs = rhs.mFieldPtrs;
        break;

      case FieldStorageType::CopyFields:
        for (auto& fptr: rhs.mFieldCache) {
          auto newField = std::make_shared<Field<Dimension, DataType>>(*fptr);
          mFieldCache.push_back(newField);
          mFieldPtrs.push_back(newField.get());
        }
        break;
      }

      //     // Register this FieldList with each Field we point to.
      //     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
      //       registerWithField(**fieldPtrItr);

      buildDependentArrays();
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
  for (auto* fptr: mFieldPtrs) *fptr = rhs;
  buildDependentArrays();  // Probably not necessary
  return *this;
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
    mFieldCache.clear();
    for (auto* fptr: mFieldPtrs) mFieldCache.push_back(std::make_shared<Field<Dimension, DataType>>(*fptr));
    mFieldPtrs.clear();
    for (auto& fptr: mFieldCache) mFieldPtrs.push_back(fptr.get());
    buildDependentArrays();

//     for (int i = 0; i < this->size(); ++i) {
//       mFieldCache.push_back(*((*this)[i]));
//       (*this)[i] = &(mFieldCache.back());
//     }

    buildDependentArrays();
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
  mFieldCache.clear();
  mStorageType = FieldStorageType::CopyFields;
  reductionType = fieldList.reductionType;
  threadMasterPtr = fieldList.threadMasterPtr;

  // Store new copies of the Fields from the other FieldList
  for (const auto* fptr: fieldList.mFieldPtrs) {
    auto newFieldPtr = std::make_shared<Field<Dimension, DataType>>(*fptr);
    mFieldCache.push_back(newFieldPtr);
    mFieldPtrs.push_back(newFieldPtr.get());
  }
  buildDependentArrays();
}

//------------------------------------------------------------------------------
// Test if a given field is part of this field list.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
haveField(const Field<Dimension, DataType>& field) const {
  return std::find(this->begin(), this->end(), &field) != this->end();
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
    const auto n = this->size();
    ENSURE(n == fieldList.size());
    ENSURE(mNodeListPtrs == fieldList.mNodeListPtrs);
    for (auto i = 0u; i < n; ++i) {
      *mFieldPtrs[i] = *fieldList.mFieldPtrs[i];
    }
  } // OMP critical
  buildDependentArrays();
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
  mFieldCache.clear();
  mStorageType = FieldStorageType::ReferenceFields;
  buildDependentArrays();
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

  // Insert the field.
  switch(storageType()) {
  case FieldStorageType::ReferenceFields:
    mFieldPtrs.push_back(const_cast<Field<Dimension, DataType>*>(&field));
    break;

  case FieldStorageType::CopyFields:
    mFieldCache.push_back(std::make_shared<Field<Dimension, DataType>>(field));
    mFieldPtrs.push_back(mFieldCache.back().get());
  }

//   registerWithField(*mFieldPtrs.back());

  buildDependentArrays();

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
    break;
  }
  buildDependentArrays();
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
  auto* fieldPtr = mFieldCache.back().get();
  mFieldPtrs.push_back(fieldPtr);

  buildDependentArrays();

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  ENSURE(mFieldPtrs[mNodeListIndexMap[fieldPtr->nodeListPtr()]] == *(fieldForNodeList(*(fieldPtr->nodeListPtr()))));
  ENSURE(this->size() == mNodeListPtrs.size());
  ENSURE(mFieldBasePtrs.size() == mFieldPtrs.size());
  END_CONTRACT_SCOPE

}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::value_type
FieldList<Dimension, DataType>::
operator[](const size_t index) const {
  REQUIRE2(index < this->size(), "FieldList index ERROR: out of bounds " << index << " !< " << this->size());
  return mFieldPtrs[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::value_type
FieldList<Dimension, DataType>::
at(const size_t index) const {
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

#ifdef USE_MPI
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
  size_t iNodeList = 0u;
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
  size_t iNodeList = 0u;
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
  size_t iNodeList = 0u;
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
// Add two FieldLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator+(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    const auto n = this->numFields();
    REQUIRE(n == rhs.numFields());
    for (auto i = 0u; i < n; ++i) {
      REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
    }
  }
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
  {
    const auto n = this->numFields();
    REQUIRE(n == rhs.numFields());
    for (auto i = 0u; i < n; ++i) {
      REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
    }
  }
  END_CONTRACT_SCOPE

  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result -= rhs;
  return result;
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
// All internal values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<DataType>
FieldList<Dimension, DataType>::
internalValues() const {
  const auto ntot = this->numInternalElements();
  std::vector<DataType> result(ntot);
  size_t offset = 0u;
  for (auto* fptr: mFieldPtrs) {
    const auto n = fptr->numInternalElements();
    std::copy(fptr->begin(), fptr->begin() + n, result.begin() + offset);
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
  for (auto* fptr: mFieldPtrs) {
    const auto n = fptr->numGhostElements();
    std::copy(fptr->begin(), fptr->begin() + n, result.begin() + offset);
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
  for (auto* fptr: mFieldPtrs) {
    const auto n = fptr->numElements();
    std::copy(fptr->begin(), fptr->begin() + n, result.begin() + offset);
    offset += n;
  }
  CHECK(offset == ntot);
  return result;
}

//------------------------------------------------------------------------------
// Build all the dependent FieldList internal state assuming mFieldPtrs and
// mFieldCache are provided.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
buildDependentArrays() {
  NodeListRegistrar<Dimension>::sortInNodeListOrder(mFieldPtrs.begin(), mFieldPtrs.end());
  mFieldBasePtrs.clear();
  mFieldViewPtrs.clear();
  mNodeListPtrs.clear();
  mNodeListIndexMap.clear();
  int i = 0;
  for (auto* fptr: mFieldPtrs) {
    mFieldBasePtrs.push_back(fptr);
    mFieldViewPtrs.push_back(fptr);
    auto* nptr = const_cast<NodeList<Dimension>*>(fptr->nodeListPtr());
    mNodeListPtrs.push_back(nptr);
    mNodeListIndexMap[nptr] = i++;
  }
  CHECK(i == int(mFieldPtrs.size()));
#ifdef SPHERAL_UNIFIED_MEMORY
  mSpanFieldViews = SPHERAL_SPAN_TYPE<typename FieldListView<Dimension, DataType>::value_type>(&mFieldViewPtrs[0], mFieldViewPtrs.size());
#else
  const auto n = mFieldPtrs.size();
  if (mSpanFieldViews.size() == 0u and !mSpanFieldViews.getPointer(chai::CPU, false)) {
    std::cerr << "FIELDLIST::BDA : allocate" << std::endl;
    mSpanFieldViews.free();
    mSpanFieldViews.allocate(n, chai::CPU);
    std::cerr << " --> SUCCESS" << std::endl;
  } else {
    std::cerr << "FIELDLIST::BDA : reallocate" << std::endl;
    mSpanFieldViews.free();
    mSpanFieldViews.reallocate(n);
    std::cerr << " --> SUCCESS" << std::endl;
  }
  for (size_t i = 0; i < n; ++i) {
    std::cerr << "FIELDLIST::BDA : Getting FieldView " << i << " of " << n << std::endl;
    mSpanFieldViews[i] = &mFieldPtrs[i]->view(); // static_cast<FieldView<Dimension, DataType>*>(mFieldPtrs[i]);
    std::cerr << " --> SUCCESS" << std::endl;
  }
  std::cerr << "FIELDLIST:BDA : registerTouch" << std::endl;
  mSpanFieldViews.registerTouch(chai::CPU);
  std::cerr << " --> SUCCESS" << std::endl;
#endif
  ENSURE(mFieldBasePtrs.size() == mFieldPtrs.size());
  ENSURE(mFieldViewPtrs.size() == mFieldPtrs.size());
  ENSURE(mNodeListPtrs.size() == mFieldPtrs.size());
  ENSURE(mNodeListIndexMap.size() == mFieldPtrs.size());
  ENSURE(mSpanFieldViews.size() == mFieldPtrs.size());
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

//------------------------------------------------------------------------------
// Extract a view
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>&
FieldList<Dimension, DataType>::view() {
#ifdef SPHERAL_UNIFIED_MEMORY
  return static_cast<FieldListView<Dimension, DataType>&>(*this);
#else
  auto func = [](
      const chai::PointerRecord *,
      chai::Action,
      chai::ExecutionSpace) {};

  return this->view(func, func);
#endif
}

template<typename Dimension, typename DataType>
template<typename FL>
inline
FieldListView<Dimension, DataType>&
FieldList<Dimension, DataType>::view(FL&& extension) {
  auto func = [](
      const chai::PointerRecord *,
      chai::Action,
      chai::ExecutionSpace) {};

  return this->view(std::forward<FL>(extension), func);
}

template<typename Dimension, typename DataType>
template<typename FL, typename F>
inline
FieldListView<Dimension, DataType>&
FieldList<Dimension, DataType>::view(FL&& extension, F&& field_extension) {
#ifndef SPHERAL_UNIFIED_MEMORY
  auto callback = getFieldListCallback(std::forward<FL>(extension));

  const auto n = this->size();
  if (mSpanFieldViews.size() == 0 && !mSpanFieldViews.getPointer(chai::CPU, false)) {
    std::cerr << "FIELDLIST::view(e, e) : allocate" << std::endl;
    mSpanFieldViews.allocate(n, chai::CPU, callback);
    std::cerr << " --> SUCCESS" << std::endl;
  } else {
    std::cerr << "FIELDLIST::view(e, e) : setcallback" << std::endl;
    mSpanFieldViews.setUserCallback(callback);
    mSpanFieldViews.reallocate(n);
    std::cerr << " --> SUCCESS" << std::endl;
  }

  std::cerr << "FIELDLIST::view(e, e) : set field views" << std::endl;
  for (size_t i = 0; i < n; ++i) {
    mSpanFieldViews[i] = &(mFieldPtrs[i]->view(std::forward<F>(field_extension)));
  }
  std::cerr << " --> SUCCESS" << std::endl;

  std::cerr << "FIELDLIST::view(e, e) : register touch" << std::endl;
  mSpanFieldViews.registerTouch(chai::CPU);
  std::cerr << " --> SUCCESS" << std::endl;
#endif

  return static_cast<FieldListView<Dimension, DataType>&>(*this);
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
