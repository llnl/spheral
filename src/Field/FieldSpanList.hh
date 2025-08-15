//---------------------------------Spheral++----------------------------------//
// FieldSpanList -- A list container for FieldSpans
//
// Mostly a thin reference container for FieldSpans corresponding to the Fields
// in a normal FieldList.
//
// Created by JMO, Thu May  1 15:20:11 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldSpanList__
#define __Spheral__FieldSpanList__

#include "Utilities/span.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class FieldSpan;

template<typename Dimension, typename DataType>
class FieldSpanList {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  
  using FieldDimension = Dimension;
  using FieldDataType = DataType;

  using value_type = FieldSpan<Dimension, DataType>*;    // STL compatibility
#ifdef SPHERAL_UNIFIED_MEMORY
  using ContainerType = SPHERAL_SPAN_TYPE<value_type>;
#else
  using ContainerType = typename chai::ManagedArray<value_type>;
#endif

  using iterator = typename ContainerType::iterator;
  using reverse_iterator = typename ContainerType::reverse_iterator;

  // Constructors, destructor
  SPHERAL_HOST_DEVICE FieldSpanList() = default;
  SPHERAL_HOST_DEVICE FieldSpanList(FieldSpanList& rhs) = default;
  SPHERAL_HOST_DEVICE FieldSpanList(FieldSpanList&& rhs) = default;
  SPHERAL_HOST        virtual ~FieldSpanList();

  // Assignment
  SPHERAL_HOST_DEVICE FieldSpanList& operator=(FieldSpanList& rhs) = default;
  SPHERAL_HOST_DEVICE FieldSpanList& operator=(const DataType& rhs);

  // Provide the standard iterators over the FieldSpans
  SPHERAL_HOST_DEVICE iterator begin()                                                         { return mSpanFieldSpans.begin(); } 
  SPHERAL_HOST_DEVICE iterator end()                                                           { return mSpanFieldSpans.end(); }   
  SPHERAL_HOST_DEVICE reverse_iterator rbegin()                                                { return mSpanFieldSpans.rbegin(); }
  SPHERAL_HOST_DEVICE reverse_iterator rend()                                                  { return mSpanFieldSpans.rend(); }  

  // Index operator.
  SPHERAL_HOST_DEVICE value_type operator[](const size_t index) const;
  SPHERAL_HOST_DEVICE value_type at(const size_t index) const;

  // Provide direct access to Field elements
  SPHERAL_HOST_DEVICE DataType& operator()(const size_t fieldIndex,
                                           const size_t nodeIndex) const;

  // Reproduce the standard Field operators for FieldSpanLists.
  SPHERAL_HOST_DEVICE void applyMin(const DataType& dataMin);
  SPHERAL_HOST_DEVICE void applyMax(const DataType& dataMax);

  SPHERAL_HOST_DEVICE void applyScalarMin(const Scalar dataMin);
  SPHERAL_HOST_DEVICE void applyScalarMax(const Scalar dataMax);

  SPHERAL_HOST_DEVICE FieldSpanList& operator+=(const FieldSpanList& rhs);
  SPHERAL_HOST_DEVICE FieldSpanList& operator-=(const FieldSpanList& rhs);

  SPHERAL_HOST_DEVICE FieldSpanList& operator+=(const DataType& rhs);
  SPHERAL_HOST_DEVICE FieldSpanList& operator-=(const DataType& rhs);

  SPHERAL_HOST_DEVICE FieldSpanList& operator*=(const FieldSpanList<Dimension, Scalar>& rhs);
  SPHERAL_HOST_DEVICE FieldSpanList& operator*=(const Scalar& rhs);

  SPHERAL_HOST_DEVICE FieldSpanList& operator/=(const FieldSpanList<Dimension, Scalar>& rhs);
  SPHERAL_HOST_DEVICE FieldSpanList& operator/=(const Scalar& rhs);

  // Some useful reduction operations (local versions -- no MPI reductions)
  SPHERAL_HOST_DEVICE DataType localSumElements() const;
  SPHERAL_HOST_DEVICE DataType localMin() const;
  SPHERAL_HOST_DEVICE DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  SPHERAL_HOST_DEVICE bool operator==(const FieldSpanList& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const FieldSpanList& rhs) const;

  // Comparison operators (Field-value element wise).
  SPHERAL_HOST_DEVICE bool operator==(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const DataType& rhs) const;

  // The number of fields in the FieldSpanList.
  SPHERAL_HOST_DEVICE size_t numFields() const                                                 { return mSpanFieldSpans.size(); } 
  SPHERAL_HOST_DEVICE size_t size() const                                                      { return mSpanFieldSpans.size(); } 
  SPHERAL_HOST_DEVICE bool empty() const                                                       { return mSpanFieldSpans.empty(); }

  // The number of nodes in the FieldSpanList.
  SPHERAL_HOST_DEVICE size_t numElements() const;
  
  // The number of internal nodes in the FieldSpanList.
  SPHERAL_HOST_DEVICE size_t numInternalElements() const;
  
  // The number of ghost nodes in the FieldSpanList.
  SPHERAL_HOST_DEVICE size_t numGhostElements() const;

#ifndef SPHERAL_UNIFIED_MEMORY
  //..........................................................................
  // These methods only make sense when we're using the ManagedArray
  SPHERAL_HOST        void move(chai::ExecutionSpace space, bool recursive = true):
  SPHERAL_HOST_DEVICE value_type* data() const                                                 { return mSpanFieldSpans.data(); }
  SPHERAL_HOST        value_type* data(chai::ExecutionSpace space, bool do_move = true) const  { return mSpanFieldSpans.data(space, do_move); }
  SPHERAL_HOST        void touch(chai::ExecutionSpace space, bool recursive = true);
  //..........................................................................
#endif

protected:
  //--------------------------- Protected Interface ---------------------------//
  ContainerType mSpanFieldSpans;
};

}

#include "FieldSpanListInline.hh"

#endif
