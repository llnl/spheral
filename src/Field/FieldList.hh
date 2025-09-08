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

#include "FieldListBase.hh"
#include "Utilities/OpenMP_wrapper.hh"
#include "Utilities/Logger.hh"
#include "FieldView.hh"
#include "FieldListView.hh"
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

// An enum for selecting how Fields are stored in FieldLists.
enum class FieldStorageType {
  ReferenceFields = 0,
  CopyFields = 1
};

template<typename Dimension, typename DataType>
class FieldList: public FieldListBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  
  typedef Dimension FieldDimension;
  typedef DataType FieldDataType;

  typedef FieldBase<Dimension>* BaseElementType;
  typedef Field<Dimension, DataType>* ElementType;
  typedef Field<Dimension, DataType>* value_type;    // STL compatibility
  typedef std::vector<ElementType> StorageType;

  typedef typename StorageType::iterator iterator;
  typedef typename StorageType::const_iterator const_iterator;
  typedef typename StorageType::reverse_iterator reverse_iterator;
  typedef typename StorageType::const_reverse_iterator const_reverse_iterator;

  typedef std::vector<DataType> CacheElementsType;
  typedef typename CacheElementsType::iterator cache_iterator;
  typedef typename CacheElementsType::const_iterator const_cache_iterator;

  // Constructors.
  FieldList();
  explicit FieldList(FieldStorageType aStorageType);
  FieldList(const FieldList& rhs);

  // Destructor.
  ~FieldList();

  // Assignment.
  FieldList& operator=(const FieldList& rhs);
  FieldList& operator=(const DataType& rhs);

  // Access the storage type of the field list.
  FieldStorageType storageType() const;

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

  // Provide the standard iterators over the Fields.
  iterator begin();
  iterator end();
  reverse_iterator rbegin();
  reverse_iterator rend();

  const_iterator begin() const;
  const_iterator end() const;
  const_reverse_iterator rbegin() const;
  const_reverse_iterator rend() const;

  // Iterators over FieldBase* required by base class.
  virtual typename FieldListBase<Dimension>::iterator begin_base();
  virtual typename FieldListBase<Dimension>::iterator end_base();
  virtual typename FieldListBase<Dimension>::reverse_iterator rbegin_base();
  virtual typename FieldListBase<Dimension>::reverse_iterator rend_base();

  virtual typename FieldListBase<Dimension>::const_iterator begin_base() const;
  virtual typename FieldListBase<Dimension>::const_iterator end_base() const;
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rbegin_base() const;
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rend_base() const;

  // Index operator.
  ElementType operator[](const size_t index);
  ElementType operator[](const size_t index) const;

  ElementType at(const size_t index);
  ElementType at(const size_t index) const;

  // Return an iterator to the Field associated with the given NodeList.
  iterator fieldForNodeList(const NodeList<Dimension>& nodeList);
  const_iterator fieldForNodeList(const NodeList<Dimension>& nodeList) const;

  // Provide access to the Field elements via NodeIterators.
  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  // Provide a more primitive access to Field elements, based on the index of the Field
  // and the node index within that field.
  DataType& operator()(const size_t fieldIndex,
                       const size_t nodeIndex);
  const DataType& operator()(const size_t fieldIndex,
                             const size_t nodeIndex) const;

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

  // Reproduce the standard Field operators for FieldLists.
  void Zero();
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const double dataMin);
  void applyScalarMax(const double dataMax);

  FieldList operator+(const FieldList& rhs) const;
  FieldList operator-(const FieldList& rhs) const;

  FieldList& operator+=(const FieldList& rhs);
  FieldList& operator-=(const FieldList& rhs);

  FieldList operator+(const DataType& rhs) const;
  FieldList operator-(const DataType& rhs) const;

  FieldList& operator+=(const DataType& rhs);
  FieldList& operator-=(const DataType& rhs);

  FieldList<Dimension, DataType>& operator*=(const FieldList<Dimension, Scalar>& rhs);
  FieldList<Dimension, DataType>& operator*=(const Scalar& rhs);

  FieldList<Dimension, DataType> operator/(const FieldList<Dimension, Scalar>& rhs) const;
  FieldList<Dimension, DataType> operator/(const Scalar& rhs) const;

  FieldList<Dimension, DataType>& operator/=(const FieldList<Dimension, Scalar>& rhs);
  FieldList<Dimension, DataType>& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const;
  DataType localMin() const;
  DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  bool operator==(const FieldList& rhs) const;
  bool operator!=(const FieldList& rhs) const;
  bool operator>(const FieldList& rhs) const;
  bool operator<(const FieldList& rhs) const;
  bool operator>=(const FieldList& rhs) const;
  bool operator<=(const FieldList& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator>(const DataType& rhs) const;
  bool operator<(const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  // The number of fields in the FieldList.
  size_t numFields() const;
  size_t size() const;
  bool empty() const;

  // The number of nodes in the FieldList.
  size_t numElements() const;
  
  // The number of internal nodes in the FieldList.
  size_t numInternalElements() const;
  
  // The number of ghost nodes in the FieldList.
  size_t numGhostElements() const;

  // Get the NodeLists this FieldList is defined on.
  const std::vector<NodeList<Dimension>*>& nodeListPtrs() const;

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

  // TODO: Good for debugging but not necessary
  using ViewType = FieldListView<Dimension, DataType>;
  ViewType toView();

  template<typename FL>
  ViewType toView(FL&& extension);

  template<typename FL, typename F>
  ViewType toView(FL&& extension, F&& field_extension);

protected:
  template<typename F>
  auto getFieldListCallback(F callback);

private:
  //--------------------------- Private Interface ---------------------------//
  typedef std::list<std::shared_ptr<Field<Dimension, DataType> > > FieldCacheType;
  typedef std::map<const NodeList<Dimension>*, int> HashMapType;

  chai::ManagedArray<FieldView<Dimension, DataType>> mFieldViews;

  std::vector<ElementType> mFieldPtrs;
  std::vector<BaseElementType> mFieldBasePtrs;
  FieldCacheType mFieldCache;
  FieldStorageType mStorageType;

  // Maintain a vector of the NodeLists this FieldList is defined in order to
  // construct NodeIterators.
  std::vector<NodeList<Dimension>*> mNodeListPtrs;
  HashMapType mNodeListIndexMap;

  // Internal method to build the NodeListIndexMap from scratch.
  void buildNodeListIndexMap();
public:
  // A data attribute to indicate how to reduce this field across threads.
  ThreadReduction reductionType;

  // The master FieldList if this is a thread copy.
  FieldList<Dimension, DataType>* threadMasterPtr;

};

}

#include "FieldListInline.hh"

#endif
