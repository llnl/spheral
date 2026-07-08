//---------------------------------Spheral++----------------------------------//
// FieldListViewBase -- provide a handle for FieldListView without the DataType
//
// Created by JMO, Mon Jul  6 08:54:21 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldListViewBase__
#define __Spheral__FieldListViewBase__

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

namespace Spheral {

template<typename Dimension>
class FieldListViewBase:
    public chai::CHAICopyable {
  
public:
  //--------------------------- Public Interface ---------------------------//
  using FieldDimension = Dimension;

  // Constructors, destructor
  SPHERAL_HOST_DEVICE FieldListViewBase() = default;
  SPHERAL_HOST_DEVICE FieldListViewBase(const FieldListViewBase& rhs) = default;
  SPHERAL_HOST_DEVICE FieldListViewBase(FieldListViewBase&& rhs) = default;
  SPHERAL_HOST_DEVICE virtual ~FieldListViewBase() = default;

  // Assignment
  SPHERAL_HOST_DEVICE FieldListViewBase& operator=(FieldListViewBase& rhs) = default;
  SPHERAL_HOST_DEVICE FieldListViewBase& operator=(FieldListViewBase&& rhs) = default;

  //..........................................................................
  // CHAI/ManagedArray specific operations
  SPHERAL_HOST virtual void move(chai::ExecutionSpace space, bool recursive = true) = 0;
  SPHERAL_HOST virtual void touch(chai::ExecutionSpace space, bool recursive = true) = 0;
};

}

#endif
