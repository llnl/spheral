//---------------------------------Spheral++----------------------------------//
// FieldListView -- A list container for FieldViews
//
// Mostly a thin reference container for FieldViews corresponding to the Fields
// in a normal FieldList.
//
// Created by JMO, Thu May  1 15:20:11 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldListView__
#define __Spheral__FieldListView__

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

#ifdef SPHERAL_UNIFIED_MEMORY
#include "Utilities/span.hh"
#endif

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class FieldView;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension, typename DataType>
class FieldListView: public chai::CHAICopyable {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  
  using FieldDimension = Dimension;
  using FieldDataType = DataType;

  using value_type = FieldView<Dimension, DataType>;    // STL compatibility

#ifdef SPHERAL_UNIFIED_MEMORY
  using ContainerType = SPHERAL_SPAN_TYPE<value_type>;
  using iterator = typename ContainerType::iterator;
#else
  using ContainerType = typename chai::ManagedArray<value_type>;
  using iterator = value_type*;
#endif

  // Constructors, destructor
  SPHERAL_HOST_DEVICE FieldListView() = default;
  SPHERAL_HOST        FieldListView(FieldList<Dimension, DataType>& rhs);
  SPHERAL_HOST_DEVICE FieldListView(const FieldListView& rhs) = default;
  SPHERAL_HOST_DEVICE FieldListView(FieldListView&& rhs) = default;
  SPHERAL_HOST_DEVICE ~FieldListView() = default;

  // Assignment
  SPHERAL_HOST_DEVICE FieldListView& operator=(FieldListView& rhs) = default;
  SPHERAL_HOST_DEVICE FieldListView& operator=(FieldListView&& rhs) = default;
  SPHERAL_HOST_DEVICE FieldListView& operator=(const DataType& rhs);

  // Provide the standard iterators over the FieldViews
#ifdef SPHERAL_UNIFIED_MEMORY
  SPHERAL_HOST_DEVICE iterator begin()                                                         { return mFieldViews.begin(); }
  SPHERAL_HOST_DEVICE iterator end()                                                           { return mFieldViews.end(); }
  SPHERAL_HOST_DEVICE bool empty()                                                             { return mFieldViews.empty(); }
#else
  SPHERAL_HOST_DEVICE iterator begin()                                                         { return &mFieldViews[0]; }
  SPHERAL_HOST_DEVICE iterator end()                                                           { return &mFieldViews[0] + mFieldViews.size(); }
  SPHERAL_HOST_DEVICE bool empty()                                                             { return mFieldViews.size() == 0u; }
#endif

  // Index operator.
  SPHERAL_HOST_DEVICE value_type& operator[](const size_t index) const;
  SPHERAL_HOST_DEVICE value_type& at(const size_t index) const;

  // Provide direct access to Field elements
  SPHERAL_HOST_DEVICE DataType& operator()(const size_t fieldIndex,
                                           const size_t nodeIndex) const;

  // Reproduce the standard Field operators for FieldListViews.
  SPHERAL_HOST_DEVICE void applyMin(const DataType& dataMin);
  SPHERAL_HOST_DEVICE void applyMax(const DataType& dataMax);

  SPHERAL_HOST_DEVICE void applyScalarMin(const Scalar dataMin);
  SPHERAL_HOST_DEVICE void applyScalarMax(const Scalar dataMax);

  SPHERAL_HOST_DEVICE FieldListView& operator+=(const FieldListView& rhs);
  SPHERAL_HOST_DEVICE FieldListView& operator-=(const FieldListView& rhs);

  SPHERAL_HOST_DEVICE FieldListView& operator+=(const DataType& rhs);
  SPHERAL_HOST_DEVICE FieldListView& operator-=(const DataType& rhs);

  SPHERAL_HOST_DEVICE FieldListView& operator*=(const FieldListView<Dimension, Scalar>& rhs);
  SPHERAL_HOST_DEVICE FieldListView& operator*=(const Scalar& rhs);

  SPHERAL_HOST_DEVICE FieldListView& operator/=(const FieldListView<Dimension, Scalar>& rhs);
  SPHERAL_HOST_DEVICE FieldListView& operator/=(const Scalar& rhs);

  // Some useful reduction operations (local versions -- no MPI reductions)
  SPHERAL_HOST_DEVICE DataType localSumElements(const bool includeGhosts = false) const;
  SPHERAL_HOST_DEVICE DataType localMin(const bool includeGhosts = false) const;
  SPHERAL_HOST_DEVICE DataType localMax(const bool includeGhosts = false) const;

  // Comparison operators (Field-Field element wise).
  SPHERAL_HOST_DEVICE bool operator==(const FieldListView& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const FieldListView& rhs) const;

  // Comparison operators (Field-value element wise).
  SPHERAL_HOST_DEVICE bool operator==(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const DataType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const DataType& rhs) const;

  // The number of fields in the FieldListView.
  SPHERAL_HOST_DEVICE size_t numFields() const                                                 { return mFieldViews.size(); } 
  SPHERAL_HOST_DEVICE size_t size() const                                                      { return mFieldViews.size(); } 

  // The number of nodes in the FieldListView.
  SPHERAL_HOST_DEVICE size_t numElements() const;
  
  // The number of internal nodes in the FieldListView.
  SPHERAL_HOST_DEVICE size_t numInternalElements() const;
  
  // The number of ghost nodes in the FieldListView.
  SPHERAL_HOST_DEVICE size_t numGhostElements() const;

  //..........................................................................
  // CHAI/ManagedArray specific operations
  SPHERAL_HOST        void move(chai::ExecutionSpace space, bool recursive = true);
  SPHERAL_HOST_DEVICE value_type* data() const;
  SPHERAL_HOST        value_type* data(chai::ExecutionSpace space, bool do_move = true) const;
  SPHERAL_HOST        void touch(chai::ExecutionSpace space, bool recursive = true);
  //..........................................................................

private:
  //--------------------------- Private Interface ---------------------------//
  ContainerType mFieldViews;
};

}

#include "FieldListViewInline.hh"

#endif
