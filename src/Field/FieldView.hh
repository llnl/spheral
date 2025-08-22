//---------------------------------Spheral++----------------------------------//
// FieldView -- provides a reference view (span) of the elements in an existing
// Field.
//
// Created by JMO, Mon Apr 28 15:05:15 PDT 2025
// Merged with FieldView work by Mike Davis and Landon Owen
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldView__
#define __Spheral_FieldView__

#ifdef SPHERAL_UNIFIED_MEMORY
#include "Utilities/span.hh"
#else
#include "chai/ManagedArray.hpp"
#endif

namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension, typename DataType>
class FieldView {
   
public:
  //--------------------------- Public Interface ---------------------------//
#ifdef SPHERAL_UNIFIED_MEMORY
  using ContainerType = SPHERAL_SPAN_TYPE<DataType>;
  using iterator = typename ContainerType::iterator;
  // using const_iterator = typename SPHERAL_SPAN_TYPE<DataType>::const_iterator;  // Not until C++23
#else
  using ContainerType = typename chai::ManagedArray<DataType>;
  using iterator = DataType*;
#endif

  using Scalar = typename Dimension::Scalar;

  using FieldDimension = Dimension;
  using FieldDataType = DataType;
  using value_type = DataType;      // STL compatibility.

  // Constructors, destructor
  SPHERAL_HOST_DEVICE FieldView() = default;
  SPHERAL_HOST_DEVICE FieldView(const FieldView& rhs) = default;
  SPHERAL_HOST_DEVICE FieldView(FieldView&& rhs) = default;
  SPHERAL_HOST_DEVICE ~FieldView() = default;

  // Assignment
  SPHERAL_HOST_DEVICE FieldView& operator=(FieldView& rhs) = default;
  SPHERAL_HOST_DEVICE FieldView& operator=(const DataType& rhs);

  // Element access.
  SPHERAL_HOST_DEVICE DataType& operator()(size_t index);
  SPHERAL_HOST_DEVICE const DataType& operator()(size_t index) const;

  SPHERAL_HOST_DEVICE DataType& at(size_t index);
  SPHERAL_HOST_DEVICE const DataType& at(size_t index) const;

  SPHERAL_HOST_DEVICE DataType& operator[](const size_t index);
  SPHERAL_HOST_DEVICE const DataType& operator[](const size_t index) const;

  // The number of elements in the field.
  SPHERAL_HOST_DEVICE size_t numElements()         const { return mDataSpan.size(); }
  SPHERAL_HOST_DEVICE size_t numInternalElements() const { return mNumInternalElements; }
  SPHERAL_HOST_DEVICE size_t numGhostElements()    const { return mNumGhostElements; }

  // Methods to apply limits to Field data members.
  SPHERAL_HOST_DEVICE void applyMin(const DataType& dataMin);
  SPHERAL_HOST_DEVICE void applyMax(const DataType& dataMax);

  SPHERAL_HOST_DEVICE void applyScalarMin(const Scalar& dataMin);
  SPHERAL_HOST_DEVICE void applyScalarMax(const Scalar& dataMax);

  // Standard field additive operators.
  SPHERAL_HOST_DEVICE FieldView& operator+=(const FieldView& rhs);
  SPHERAL_HOST_DEVICE FieldView& operator-=(const FieldView& rhs);

  SPHERAL_HOST_DEVICE FieldView& operator+=(const DataType& rhs);
  SPHERAL_HOST_DEVICE FieldView& operator-=(const DataType& rhs);

  // Multiplication and division by scalar(s)
  SPHERAL_HOST_DEVICE FieldView& operator*=(const FieldView<Dimension, Scalar>& rhs);
  SPHERAL_HOST_DEVICE FieldView& operator/=(const FieldView<Dimension, Scalar>& rhs);

  SPHERAL_HOST_DEVICE FieldView& operator*=(const Scalar& rhs);
  SPHERAL_HOST_DEVICE FieldView& operator/=(const Scalar& rhs);

  // Some useful reduction operations (local versions -- no MPI reductions)
  SPHERAL_HOST_DEVICE DataType localSumElements() const;
  SPHERAL_HOST_DEVICE DataType localMin() const;
  SPHERAL_HOST_DEVICE DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  SPHERAL_HOST_DEVICE bool operator==(const FieldView& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const FieldView& rhs) const;
  // bool operator> (const FieldView& rhs) const;
  // bool operator< (const FieldView& rhs) const;
  // bool operator>=(const FieldView& rhs) const;
  // bool operator<=(const FieldView& rhs) const;

  // Comparison operators (Field-value element wise).
  SPHERAL_HOST_DEVICE bool operator==(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator> (const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator< (const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const DataType& rhs) const;

  // Provide the standard iterator methods over the field.
#ifdef SPHERAL_UNIFIED_MEMORY
  SPHERAL_HOST        iterator begin() const                                              { return mDataSpan.begin(); }
  SPHERAL_HOST        iterator end() const                                                { return mDataSpan.end(); }
  SPHERAL_HOST        iterator internalBegin() const                                      { return mDataSpan.begin(); }
  SPHERAL_HOST        iterator internalEnd() const                                        { return mDataSpan.begin() + mNumInternalElements; }
  SPHERAL_HOST        iterator ghostBegin() const                                         { return mDataSpan.begin() + mNumInternalElements; }
  SPHERAL_HOST        iterator ghostEnd() const                                           { return mDataSpan.end(); }
  SPHERAL_HOST        value_type& front()                                                 { return mDataSpan.front(); }
  SPHERAL_HOST        value_type& back()                                                  { return mDataSpan.back(); }
  SPHERAL_HOST        const value_type& front() const                                     { return mDataSpan.front(); }
  SPHERAL_HOST        const value_type& back()  const                                     { return mDataSpan.back(); }
  SPHERAL_HOST        bool empty() const                                                  { return mDataSpan.empty(); }
#else
  SPHERAL_HOST        iterator begin() const                                              { return &mDataSpan[0]; }
  SPHERAL_HOST        iterator end() const                                                { return &mDataSpan[0] + mDataSpan.size(); }
  SPHERAL_HOST        iterator internalBegin() const                                      { return this->begin(); }
  SPHERAL_HOST        iterator internalEnd() const                                        { return this->begin() + mNumInternalElements; }
  SPHERAL_HOST        iterator ghostBegin() const                                         { return this->begin() + mNumInternalElements; }
  SPHERAL_HOST        iterator ghostEnd() const                                           { return this->end(); }
  SPHERAL_HOST        value_type& front()                                                 { return mDataSpan[0]; }
  SPHERAL_HOST        value_type& back()                                                  { return this->empty() ? mDataSpan[0] : mDataSpan[mDataSpan.size() - 1u]; }
  SPHERAL_HOST        const value_type& front() const                                     { return mDataSpan[0]; }
  SPHERAL_HOST        const value_type& back()  const                                     { return this->empty() ? mDataSpan[0] : mDataSpan[mDataSpan.size() - 1u]; }
  SPHERAL_HOST        bool empty() const                                                  { return mDataSpan.size() == 0u; }
#endif

#ifndef SPHERAL_UNIFIED_MEMORY
  //..........................................................................
  // These methods only make sense when we're using the ManagedArray
  SPHERAL_HOST_DEVICE DataType* data() const                                             { return mDataSpan.getActivePointer(); }
  SPHERAL_HOST        DataType* data(chai::ExecutionSpace space,
                                     bool do_move = true) const                          { return mDataSpan.data(space, do_move); }
  void move(chai::ExecutionSpace space)                                                  { mDataSpan.move(space); }
  SPHERAL_HOST_DEVICE void shallowCopy(FieldView const& other) const                     { mDataSpan.shallowCopy(other.mDataSpan); }

  SPHERAL_HOST        void touch(chai::ExecutionSpace space)                             { mDataSpan.registerTouch(space); }
  //..........................................................................
#endif

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Private Data
  ContainerType mDataSpan;
  size_t mNumInternalElements, mNumGhostElements;
};

} // namespace Spheral

#include "FieldViewInline.hh"

#endif
