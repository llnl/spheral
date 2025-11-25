//---------------------------------Spheral++----------------------------------//
// FieldList -- A list container for Fields.  This is how Spheral++ defines
//              "global" fields, or fields that extend over more than one
//              NodeList.
// A FieldList can either just hold pointers to externally stored Fields, or
// copy the Fields to internal storage.
//
// Created by JMO, Sat Feb  5 12:57:58 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldList__
#define __Spheral__FieldList__

#include "Field/FieldListBase.hh"
#include "Field/FieldListView.hh"
#include "Utilities/span.hh"
#include "Utilities/OpenMP_wrapper.hh"
#include "Utilities/Logger.hh"

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"

#include <vector>
#include <list>
#include <map>
#include <memory>
#include <functional>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class NodeIteratorBase;
template<typename Dimension> class AllNodeIterator;
template<typename Dimension> class InternalNodeIterator;
template<typename Dimension> class GhostNodeIterator;
template<typename Dimension> class MasterNodeIterator;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class NodeList;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldView;
template<typename Dimension, typename DataType> class FieldListView;

// An enum for selecting how Fields are stored in FieldLists.
enum class FieldStorageType {
  ReferenceFields = 0,
  CopyFields = 1
};

template<typename Dimension, typename DataType>
class FieldList:
    public FieldListBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  
  using FieldDimension = Dimension;
  using FieldDataType = DataType;

  using BaseElementType = FieldBase<Dimension>*;
  using ElementType = Field<Dimension, DataType>*;
  using value_type = Field<Dimension, DataType>*;    // STL compatibility
  using StorageType = std::vector<ElementType>;

  using iterator = typename StorageType::iterator;
  using const_iterator = typename StorageType::const_iterator;
  using reverse_iterator = typename StorageType::reverse_iterator;
  using const_reverse_iterator = typename StorageType::const_reverse_iterator;

  using CacheElementsType = std::vector<DataType>;
  using cache_iterator = typename CacheElementsType::iterator;
  using const_cache_iterator = typename CacheElementsType::const_iterator;

  using ViewType = FieldListView<Dimension, DataType>;
  using FieldPtrSpan = SPHERAL_SPAN_TYPE<Field<Dimension, DataType>*>;
  using FieldBasePtrSpan = SPHERAL_SPAN_TYPE<FieldBase<Dimension>*>;
  
#ifdef SPHERAL_UNIFIED_MEMORY
  using FieldViewSpan = SPHERAL_SPAN_TYPE<FieldView<Dimension, DataType>>;
#else
  using FieldViewSpan = chai::ManagedArray<FieldView<Dimension, DataType>>;
#endif

  // constructors.
  FieldList();
  explicit FieldList(FieldStorageType aStorageType);
  FieldList(const FieldList& rhs);

  // Destructor.
  virtual ~FieldList();

  // Assignment.
  FieldList& operator=(const FieldList& rhs);
  FieldList& operator=(const DataType& rhs);

  // Access the storage type of the field list.
  FieldStorageType storageType() const { return mStorageType; }

  // Force the Field storage to be Copy.
  void copyFields();

  // Store copies of Fields from another FieldList
  void copyFields(const FieldList<Dimension, DataType>& fieldList);

  // Test if the given field (or NodeList) is part of a FieldList.
  bool haveField(const Field<Dimension, DataType>& field) const;
  bool haveNodeList(const NodeList<Dimension>& nodeList) const;

  // Force the Field members of this FieldList to be equal to those of
  // another FieldList.
  void assignFields(const FieldList& fieldList);

  // Make this FieldList reference the Fields of another.
  void referenceFields(const FieldList& fieldList);

  // Convenience methods to add and delete Fields.
  void appendField(const Field<Dimension, DataType>& field);
  void deleteField(const Field<Dimension, DataType>& field);

  // Construct a new field and add it to the FieldList.
  // Note this only makes sense when we're storing fields as copies!
  void appendNewField(const typename Field<Dimension, DataType>::FieldName name,
                      const NodeList<Dimension>& nodeList,
                      const DataType value);

  // Span views
  FieldPtrSpan fieldPtrs() const                                                        { return FieldPtrSpan(mFieldPtrs); }
  FieldBasePtrSpan fieldBasePtrs() const                                                { return FieldBasePtrSpan(mFieldBasePtrs); }
  FieldViewSpan fieldViews()                                                            { return FieldViewSpan(mFieldViews); }

  // Provide the standard iterators over the Fields.
  iterator begin()                                                                      { return mFieldPtrs.begin(); } 
  iterator end()                                                                        { return mFieldPtrs.end(); }   
  reverse_iterator rbegin()                                                             { return mFieldPtrs.rbegin(); }
  reverse_iterator rend()                                                               { return mFieldPtrs.rend(); }  

  const_iterator begin()                                                          const { return mFieldPtrs.begin(); } 
  const_iterator end()                                                            const { return mFieldPtrs.end(); }   
  const_reverse_iterator rbegin()                                                 const { return mFieldPtrs.rbegin(); }
  const_reverse_iterator rend()                                                   const { return mFieldPtrs.rend(); }  

  // Iterators over FieldBase* required by base class.
  virtual typename FieldListBase<Dimension>::iterator begin_base()                      override { return mFieldBasePtrs.begin(); } 
  virtual typename FieldListBase<Dimension>::iterator end_base()                        override { return mFieldBasePtrs.end(); }   
  virtual typename FieldListBase<Dimension>::reverse_iterator rbegin_base()             override { return mFieldBasePtrs.rbegin(); }
  virtual typename FieldListBase<Dimension>::reverse_iterator rend_base()               override { return mFieldBasePtrs.rend(); }  

  virtual typename FieldListBase<Dimension>::const_iterator begin_base()          const override { return mFieldBasePtrs.begin(); } 
  virtual typename FieldListBase<Dimension>::const_iterator end_base()            const override { return mFieldBasePtrs.end(); }   
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rbegin_base() const override { return mFieldBasePtrs.rbegin(); }
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rend_base()   const override { return mFieldBasePtrs.rend(); }  

  // Index operator.
  value_type operator[](const size_t index) const;
  value_type at(const size_t index) const;

  // Provide direct access to Field elements
  DataType& operator()(const size_t fieldIndex,
                       const size_t nodeIndex) const;

  // Return an iterator to the Field associated with the given NodeList.
  iterator fieldForNodeList(const NodeList<Dimension>& nodeList);
  const_iterator fieldForNodeList(const NodeList<Dimension>& nodeList) const;

  // Provide access to the Field elements via NodeIterators.
  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  // Return the interpolated value of the FieldList at a position.
  DataType operator()(const Vector& position,
                      const TableKernel<Dimension>& W) const;

  // Provide NodeIterators on the elements of the FieldList.
  AllNodeIterator<Dimension> nodeBegin() const;
  AllNodeIterator<Dimension> nodeEnd() const;
  
  InternalNodeIterator<Dimension> internalNodeBegin() const;
  InternalNodeIterator<Dimension> internalNodeEnd() const;
  
  GhostNodeIterator<Dimension> ghostNodeBegin() const;
  GhostNodeIterator<Dimension> ghostNodeEnd() const;
  
  MasterNodeIterator<Dimension> masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const;
  MasterNodeIterator<Dimension> masterNodeEnd() const;
  
  CoarseNodeIterator<Dimension> coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const;
  CoarseNodeIterator<Dimension> coarseNodeEnd() const;
  
  RefineNodeIterator<Dimension> refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const;
  RefineNodeIterator<Dimension> refineNodeEnd() const;

  // Provide a convenience function for setting the neighbor node information
  // for all the NodeList in this FieldList.
  void setMasterNodeLists(const Vector& r, const SymTensor& H,
                          std::vector<std::vector<int>>& masterLists,
                          std::vector<std::vector<int>>& coarseNeighbors) const;
  void setMasterNodeLists(const Vector& r,
                          std::vector<std::vector<int>>& masterLists,
                          std::vector<std::vector<int>>& coarseNeighbors) const;

  void setRefineNodeLists(const Vector& r, const SymTensor& H,
                          const std::vector<std::vector<int>>& coarseNeighbors,
                          std::vector<std::vector<int>>& refineNeighbors) const;
  void setRefineNodeLists(const Vector& r,
                          const std::vector<std::vector<int>>& coarseNeighbors,
                          std::vector<std::vector<int>>& refineNeighbors) const;

  // Zero a FieldList
  void Zero();

  // Reproduce the standard Field operators for FieldLists.
  FieldList& operator+=(const FieldList& rhs);
  FieldList& operator-=(const FieldList& rhs);

  FieldList& operator+=(const DataType& rhs);
  FieldList& operator-=(const DataType& rhs);

  FieldList& operator*=(const FieldList<Dimension, Scalar>& rhs);
  FieldList& operator*=(const Scalar& rhs);

  FieldList& operator/=(const FieldList<Dimension, Scalar>& rhs);
  FieldList& operator/=(const Scalar& rhs);

  FieldList operator+(const FieldList& rhs) const;
  FieldList operator-(const FieldList& rhs) const;

  FieldList operator+(const DataType& rhs) const;
  FieldList operator-(const DataType& rhs) const;

  FieldList operator/(const FieldList<Dimension, Scalar>& rhs) const;
  FieldList operator/(const Scalar& rhs) const;

  // Some useful reduction operations.
  DataType localSumElements(const bool includeGhosts = false) const;
  DataType localMin(const bool includeGhosts = false) const;
  DataType localMax(const bool includeGhosts = false) const;

  DataType sumElements(const bool includeGhosts = false) const;
  DataType min(const bool includeGhosts = false) const;
  DataType max(const bool includeGhosts = false) const;

  // Apply limiting
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const Scalar dataMin);
  void applyScalarMax(const Scalar dataMax);

  // Comparison operators (Field-Field element wise).
  bool operator==(const FieldList& rhs) const;
  bool operator!=(const FieldList& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator> (const DataType& rhs) const;
  bool operator< (const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  // The number of fields in the FieldList
  size_t numFields() const                                                              { return mFieldPtrs.size(); } 
  size_t size() const                                                                   { return mFieldPtrs.size(); } 
  bool empty() const                                                                    { return mFieldPtrs.empty(); }

  // The number of nodes in the FieldList
  size_t numElements() const;
  size_t numInternalElements() const;
  size_t numGhostElements() const;

  // Get the NodeLists this FieldList is defined on.
  const std::vector<NodeList<Dimension>*>& nodeListPtrs()                         const { return mNodeListPtrs; }

  // Helpers to flatten the values across all Fields to a single array.
  std::vector<DataType> internalValues() const;
  std::vector<DataType> ghostValues() const;
  std::vector<DataType> allValues() const;

  //----------------------------------------------------------------------------
  // Methods to facilitate threaded computing
  // Make a local thread copy of all the Fields
  FieldList<Dimension, DataType> threadCopy(const ThreadReduction reductionType = ThreadReduction::SUM,
                                            const bool copy = false);

  // Same thing, with a "stack" object to simplify final reduction
  FieldList<Dimension, DataType> threadCopy(typename SpheralThreads<Dimension>::FieldListStack& stack,
                                            const ThreadReduction reductionType = ThreadReduction::SUM,
                                            const bool copy = false);

  // Reduce the values in the FieldList with the passed thread-local values.
  // Not all types support less than/greater than, so we have to distinguish here
  template<typename U = DataType>
  std::enable_if_t< TypeTraits::has_less_than<U, U>::value, void>
  threadReduce() const;

  template<typename U = DataType>
  std::enable_if_t<!TypeTraits::has_less_than<U, U>::value, void>
  threadReduce() const;

  //----------------------------------------------------------------------------
  // Return a view of the FieldList (appropriate for on accelerator devices)
  // The 2nd and 3rd versions are for debugging/diagnostics.
  ViewType view();
  template<typename FLCB>               ViewType view(FLCB&& fieldlist_callback);
  template<typename FLCB, typename FCB> ViewType view(FLCB&& fieldlist_callback, FCB&& field_callback);
  void setCallback(std::function<void(const chai::PointerRecord*, chai::Action, chai::ExecutionSpace)> f);

private:
  //--------------------------- Private Interface ---------------------------//
  using FieldCacheType = std::list<std::shared_ptr<Field<Dimension, DataType>>>;
  using HashMapType = std::map<const NodeList<Dimension>*, int>;

  std::vector<ElementType> mFieldPtrs;
  std::vector<BaseElementType> mFieldBasePtrs;
  FieldCacheType mFieldCache;
  FieldStorageType mStorageType;

  // Maintain a vector of the NodeLists this FieldList is defined in order to
  // construct NodeIterators.
  std::vector<NodeList<Dimension>*> mNodeListPtrs;
  HashMapType mNodeListIndexMap;

  // FieldViews for use in FieldListView
#ifdef SPHERAL_UNIFIED_MEMORY  
  std::vector<FieldView<Dimension, DataType>> mFieldViews;
#else
  FieldViewSpan mFieldViews;
#endif

  // CHAI callback functions for debugging
  std::function<void(const chai::PointerRecord*, chai::Action, chai::ExecutionSpace)> mChaiCallback;
  auto getCallback();

  // Set the internal dependent arrays based on the Field pointers.
  virtual void buildDependentArrays() override;

public:
  // A data attribute to indicate how to reduce this field across threads.
  ThreadReduction reductionType;

  // The master FieldList if this is a thread copy.
  FieldList<Dimension, DataType>* threadMasterPtr;

};

}

#include "FieldListInline.hh"

#endif
