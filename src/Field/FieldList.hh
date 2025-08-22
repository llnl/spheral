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
#include "Utilities/OpenMP_wrapper.hh"
#include "Utilities/Logger.hh"

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"

#include <vector>
#include <list>
#include <map>
#include <memory>

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
    public FieldListBase<Dimension>,
    public FieldListView<Dimension, DataType> {

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

  // Bring in various methods hidden in FieldListView
  using FieldListView<Dimension, DataType>::operator();
  // using FieldListView<Dimension, DataType>::operator[];

  // Constructors.
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

  // Index operator.
  value_type operator[](const size_t index) const;
  value_type at(const size_t index) const;

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
  virtual typename FieldListBase<Dimension>::iterator begin_base()                      { return mFieldBasePtrs.begin(); } 
  virtual typename FieldListBase<Dimension>::iterator end_base()                        { return mFieldBasePtrs.end(); }   
  virtual typename FieldListBase<Dimension>::reverse_iterator rbegin_base()             { return mFieldBasePtrs.rbegin(); }
  virtual typename FieldListBase<Dimension>::reverse_iterator rend_base()               { return mFieldBasePtrs.rend(); }  

  virtual typename FieldListBase<Dimension>::const_iterator begin_base()          const { return mFieldBasePtrs.begin(); } 
  virtual typename FieldListBase<Dimension>::const_iterator end_base()            const { return mFieldBasePtrs.end(); }   
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rbegin_base() const { return mFieldBasePtrs.rbegin(); }
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rend_base()   const { return mFieldBasePtrs.rend(); }  

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
  FieldList operator+(const FieldList& rhs) const;
  FieldList operator-(const FieldList& rhs) const;

  FieldList operator+(const DataType& rhs) const;
  FieldList operator-(const DataType& rhs) const;

  FieldList operator/(const FieldList<Dimension, Scalar>& rhs) const;
  FieldList operator/(const Scalar& rhs) const;

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Comparison operators (Field-Field element wise).
  bool operator==(const FieldList& rhs) const { return FieldListView<Dimension, DataType>::operator==(rhs); }
  bool operator!=(const FieldList& rhs) const { return FieldListView<Dimension, DataType>::operator!=(rhs); }

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const { return FieldListView<Dimension, DataType>::operator==(rhs); }
  bool operator!=(const DataType& rhs) const { return FieldListView<Dimension, DataType>::operator!=(rhs); }
  bool operator> (const DataType& rhs) const { return FieldListView<Dimension, DataType>::operator> (rhs); }
  bool operator< (const DataType& rhs) const { return FieldListView<Dimension, DataType>::operator< (rhs); }
  bool operator>=(const DataType& rhs) const { return FieldListView<Dimension, DataType>::operator>=(rhs); }
  bool operator<=(const DataType& rhs) const { return FieldListView<Dimension, DataType>::operator<=(rhs); }

  // Get the NodeLists this FieldList is defined on.
  const std::vector<NodeList<Dimension>*>& nodeListPtrs() const { return mNodeListPtrs; }

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
  void threadReduce() const;

  //----------------------------------------------------------------------------
  // Return a view of the Field (appropriate for on accelerator devices)
  ViewType& view();

  template<typename FL>
  ViewType& view(FL&& extension);

  template<typename FL, typename F>
  ViewType& view(FL&& extension, F&& field_extension);

protected:
  //--------------------------- Protected Interface ---------------------------//
  using FieldListView<Dimension, DataType>::mSpanFieldViews;

  template<typename F>
  auto getFieldListCallback(F callback);

private:
  //--------------------------- Private Interface ---------------------------//
  using FieldCacheType = std::list<std::shared_ptr<Field<Dimension, DataType>>>;
  using HashMapType = std::map<const NodeList<Dimension>*, int>;

  std::vector<ElementType> mFieldPtrs;
  std::vector<BaseElementType> mFieldBasePtrs;
  FieldCacheType mFieldCache;
  FieldStorageType mStorageType;

  // For use when building a span view of the FieldList
  std::vector<FieldView<Dimension, DataType>*> mFieldViewPtrs;

  // Maintain a vector of the NodeLists this FieldList is defined in order to
  // construct NodeIterators.
  std::vector<NodeList<Dimension>*> mNodeListPtrs;
  HashMapType mNodeListIndexMap;

  // Internal methods
  void buildDependentArrays();

public:
  // A data attribute to indicate how to reduce this field across threads.
  ThreadReduction reductionType;

  // The master FieldList if this is a thread copy.
  FieldList<Dimension, DataType>* threadMasterPtr;

};

}

#include "FieldListInline.hh"

#endif
