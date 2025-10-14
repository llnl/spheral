//---------------------------------Spheral++----------------------------------//
// Field -- Provide a field of a type (Scalar, Vector, Tensor) over the nodes
//          in a NodeList.
//
// This version of the Field is based on standard constructs like the STL
// vector.  This will certainly be slower at run time than the Blitz Array
// class.
//
// Created by JMO, Thu Jun 10 23:26:50 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_Field__
#define __Spheral_Field__

#include "Field/FieldBase.hh"
#include "Field/FieldView.hh"
#include "Utilities/Logger.hh"

#include "axom/sidre.hpp"
#include "chai/ExecutionSpaces.hpp"
#include "chai/ManagedArray.hpp"
#include "chai/PointerRecord.hpp"
#include "chai/Types.hpp"

#ifdef USE_UVM
#include "uvm_allocator.hh"
#endif

#include <vector>
#include <functional>

namespace Spheral {

template<typename Dimension> class NodeIteratorBase;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class NodeList;
template<typename Dimension> class TableKernel;

#ifdef USE_UVM
template<typename DataType>
using DataAllocator = typename uvm_allocator::UVMAllocator<DataType>;
#else
template<typename DataType>
using DataAllocator = std::allocator<DataType>;
#endif

template<typename Dimension, typename DataType>
class Field:
    public FieldBase<Dimension>,
    public FieldView<Dimension, DataType> {
   
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  
  using FieldName = typename FieldBase<Dimension>::FieldName;
  using FieldDimension = Dimension;
  using FieldDataType = DataType;
  using value_type = DataType;      // STL compatibility.

  using iterator = typename FieldView<Dimension, DataType>::iterator;
  using const_iterator = typename std::vector<DataType,DataAllocator<DataType>>::const_iterator;

  using ViewType = FieldView<Dimension, DataType>;

  // Bring in various methods hidden in FieldView
  using FieldView<Dimension, DataType>::operator();
  using FieldView<Dimension, DataType>::operator[];

  // Constructors.
  Field() = delete;
  explicit Field(FieldName name);
  Field(FieldName name, const Field& field);
  Field(FieldName name,
        const NodeList<Dimension>& nodeList);
  Field(FieldName name,
        const NodeList<Dimension>& nodeList,
        DataType value);
  Field(FieldName name,
        const NodeList<Dimension>& nodeList, 
        const std::vector<DataType,DataAllocator<DataType>>& array);
  Field(const NodeList<Dimension>& nodeList, const Field& field);
  Field(const Field& field);
  virtual std::shared_ptr<FieldBase<Dimension> > clone() const override;

  // Destructor.
  virtual ~Field();

  // Assignment operator.
  virtual FieldBase<Dimension>& operator=(const FieldBase<Dimension>& rhs) override;
  Field& operator=(const Field& rhs);
  Field& operator=(const std::vector<DataType,DataAllocator<DataType>>& rhs);
  Field& operator=(const DataType& rhs);

  // Comparisons
  virtual bool operator==(const FieldBase<Dimension>& rhs) const override;

  // Foward the FieldView Field-Field comparison operators
  bool operator==(const Field& rhs) const { return FieldView<Dimension, DataType>::operator==(rhs); }
  bool operator!=(const Field& rhs) const { return FieldView<Dimension, DataType>::operator!=(rhs); }

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const { return FieldView<Dimension, DataType>::operator==(rhs); }
  bool operator!=(const DataType& rhs) const { return FieldView<Dimension, DataType>::operator!=(rhs); }
  bool operator> (const DataType& rhs) const { return FieldView<Dimension, DataType>::operator> (rhs); }
  bool operator< (const DataType& rhs) const { return FieldView<Dimension, DataType>::operator< (rhs); }
  bool operator>=(const DataType& rhs) const { return FieldView<Dimension, DataType>::operator>=(rhs); }
  bool operator<=(const DataType& rhs) const { return FieldView<Dimension, DataType>::operator<=(rhs); }

  // Element access by NodeIterator
  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  // The number of elements in the field.
  virtual size_t size() const override                                      { return mDataArray.size(); }

  // Zero out the field elements.
  virtual void Zero() override;

  // Standard field additive operators.
  Field operator+(const Field& rhs) const;
  Field operator-(const Field& rhs) const;

  Field operator+(const DataType& rhs) const;
  Field operator-(const DataType& rhs) const;

  // Multiplication and division by scalar(s)
  Field operator*(const Field<Dimension, Scalar>& rhs) const;
  Field operator/(const Field<Dimension, Scalar>& rhs) const;

  Field operator*(const Scalar& rhs) const;
  Field operator/(const Scalar& rhs) const;

  // Some useful reduction operations.
  DataType sumElements(const bool includeGhosts = false) const;
  DataType min(const bool includeGhosts = false) const;
  DataType max(const bool includeGhosts = false) const;

  // Provide the standard iterator methods over the field.
  const_iterator begin() const                                              { return mDataArray.begin(); }
  const_iterator end() const                                                { return mDataArray.end(); }
  const_iterator internalBegin() const                                      { return mDataArray.begin(); }
  const_iterator internalEnd() const                                        { return mDataArray.begin() + mNumInternalElements; }
  const_iterator ghostBegin() const                                         { return mDataArray.begin() + mNumInternalElements; }
  const_iterator ghostEnd() const                                           { return mDataArray.end(); }

  // We have to explicitly redefine the non-const iterators
  iterator begin()                                                          { return FieldView<Dimension, DataType>::begin(); }
  iterator end()                                                            { return FieldView<Dimension, DataType>::end(); }
  iterator internalBegin()                                                  { return FieldView<Dimension, DataType>::internalBegin(); }
  iterator internalEnd()                                                    { return FieldView<Dimension, DataType>::internalEnd(); }
  iterator ghostBegin()                                                     { return FieldView<Dimension, DataType>::ghostBegin(); }
  iterator ghostEnd()                                                       { return FieldView<Dimension, DataType>::ghostEnd(); }

  // Required functions from FieldBase
  virtual void setNodeList(const NodeList<Dimension>& nodeList) override;
  virtual std::vector<char> packValues(const std::vector<size_t>& nodeIDs) const override;
  virtual void unpackValues(const std::vector<size_t>& nodeIDs,
                            const std::vector<char>& buffer) override;
  virtual void copyElements(const std::vector<size_t>& fromIndices,
                            const std::vector<size_t>& toIndices) override;
  virtual bool fixedSizeDataType() const override;
  virtual size_t numValsInDataType() const override;
  virtual size_t sizeofDataType() const override;
  virtual size_t computeCommBufferSize(const std::vector<size_t>& packIndices,
                                       const int sendProc,
                                       const int recvProc) const override;

  // Serialization methods
  std::vector<char> serialize() const;
  void deserialize(const std::vector<char>& buf);

  // Provide std::vector copies of the data.  This is mostly useful for the
  // python interface.
  std::vector<DataType> internalValues() const;
  std::vector<DataType> ghostValues() const;
  std::vector<DataType> allValues() const;

  // Functions to help with storing the field in a Sidre datastore.
  axom::sidre::DataTypeId getAxomTypeID() const;

  // Callback method for CHAI
  void setCallback(std::function<void(const chai::PointerRecord*, chai::Action, chai::ExecutionSpace)> f);

  // Get the view (for trivially copyable types)
  ViewType view();
  template<typename CB>  ViewType view(CB&& field_callback);

protected:
  //--------------------------- Protected Interface ---------------------------//
  virtual void resizeField(size_t size) override;
  virtual void resizeFieldInternal(size_t size, size_t oldFirstGhostNode) override;
  virtual void resizeFieldGhost(size_t size) override;
  virtual void deleteElement(size_t nodeID) override;
  virtual void deleteElements(const std::vector<size_t>& nodeIDs) override;

private:
  //--------------------------- Private Interface ---------------------------//
  // Private Data
  std::vector<DataType, DataAllocator<DataType>> mDataArray;

  using FieldView<Dimension, DataType>::mDataSpan;
  using FieldView<Dimension, DataType>::mNumInternalElements;
  using FieldView<Dimension, DataType>::mNumGhostElements;

  // Callback function for debugging CHAI
  std::function<void(const chai::PointerRecord*, chai::Action, chai::ExecutionSpace)> mChaiCallback;
  auto getCallback();

  // Helper method to keep mDataSpan and mDataArray consistent
  void assignDataSpan();
};

} // namespace Spheral

#include "FieldInline.hh"

#endif
