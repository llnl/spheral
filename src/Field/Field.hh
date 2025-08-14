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

#include "FieldBase.hh"
#include "FieldView.hh"
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
    public FieldBase<Dimension> {
   
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

  using ViewType = FieldView<Dimension, DataType>;

  using iterator = typename std::vector<DataType,DataAllocator<DataType>>::iterator;
  using const_iterator = typename std::vector<DataType,DataAllocator<DataType>>::const_iterator;

  // Constructors.
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

  // Element access.
  DataType& operator()(size_t index);
  const DataType& operator()(size_t index) const;

  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  DataType& operator[](const size_t index);
  const DataType& operator[](const size_t index) const;

  DataType& at(size_t index);
  const DataType& at(size_t index) const;

  // The number of elements in the field.
  size_t numElements() const;
  size_t numInternalElements() const;
  size_t numGhostElements() const;
  virtual size_t size() const override                                      { return mDataArray.size(); }

  // Zero out the field elements.
  virtual void Zero() override;

  // Methods to apply limits to Field data members.
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const double dataMin);
  void applyScalarMax(const double dataMax);

  // Standard field additive operators.
  Field operator+(const Field& rhs) const;
  Field operator-(const Field& rhs) const;

  Field& operator+=(const Field& rhs);
  Field& operator-=(const Field& rhs);

  Field operator+(const DataType& rhs) const;
  Field operator-(const DataType& rhs) const;

  Field& operator+=(const DataType& rhs);
  Field& operator-=(const DataType& rhs);

  // Multiplication and division by scalar(s)
  Field operator*(const Field<Dimension, Scalar>& rhs) const;
  Field operator/(const Field<Dimension, Scalar>& rhs) const;

  Field& operator*=(const Field<Dimension, Scalar>& rhs);
  Field& operator/=(const Field<Dimension, Scalar>& rhs);

  Field operator*(const Scalar& rhs) const;
  Field operator/(const Scalar& rhs) const;

  Field& operator*=(const Scalar& rhs);
  Field& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const;
  DataType localMin() const;
  DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  bool operator==(const Field& rhs) const;
  bool operator!=(const Field& rhs) const;
  bool operator>(const Field& rhs) const;
  bool operator<(const Field& rhs) const;
  bool operator>=(const Field& rhs) const;
  bool operator<=(const Field& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator>(const DataType& rhs) const;
  bool operator<(const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  // Provide the standard iterator methods over the field.
  iterator begin();
  iterator end();
  iterator internalBegin();
  iterator internalEnd();
  iterator ghostBegin();
  iterator ghostEnd();

  const_iterator begin() const;
  const_iterator end() const;
  const_iterator internalBegin() const;
  const_iterator internalEnd() const;
  const_iterator ghostBegin() const;
  const_iterator ghostEnd() const;

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

  // No default constructor.
  Field() = delete;

  // ViewType controls.
  ViewType toView();

  template<typename T=DataType, typename F>
  std::enable_if_t<std::is_trivially_copyable<T>::value, ViewType>
  toView(F&& extension);

  template<typename T=DataType, typename F>
  std::enable_if_t<!std::is_trivially_copyable<T>::value, ViewType>
  toView(F&&);

protected:
  virtual void resizeField(size_t size) override;
  virtual void resizeFieldInternal(size_t size, size_t oldFirstGhostNode) override;
  virtual void resizeFieldGhost(size_t size) override;
  virtual void deleteElement(size_t nodeID) override;
  virtual void deleteElements(const std::vector<size_t>& nodeIDs) override;

  template<typename F>
  auto getFieldCallback(F callback);

private:
  //--------------------------- Private Interface ---------------------------//
  // Private Data
  std::vector<DataType, DataAllocator<DataType>> mDataArray;

  // ManagedArray owned by the field to ensure proper lifetime of the GPU data.
  chai::ManagedArray<DataType> mManagedData;

  bool mValid;
};

} // namespace Spheral

#include "FieldInline.hh"

#endif
