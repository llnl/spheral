//---------------------------------Spheral++----------------------------------//
// FieldViewBase -- much like FieldBase this class allows us to have a handle
// on FieldViews without knowing the DataType template parameter.
//
// Created by JMO, Mon Jul  6 08:54:21 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldViewBase__
#define __Spheral_FieldViewBase__

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

namespace Spheral {

template<typename Dimension>
class FieldViewBase:
    public chai::CHAICopyable {

public:
  //--------------------------- Public Interface ---------------------------//
  using FieldDimension = Dimension;

  // Constructors, destructor
  SPHERAL_HOST_DEVICE FieldViewBase() = default;
  SPHERAL_HOST_DEVICE FieldViewBase(const FieldViewBase& rhs) = default;
  SPHERAL_HOST_DEVICE FieldViewBase(FieldViewBase&& rhs) = default;
  SPHERAL_HOST_DEVICE virtual ~FieldViewBase() = default;

  // Assignment
  SPHERAL_HOST_DEVICE FieldViewBase& operator=(FieldViewBase& rhs) = default;
  SPHERAL_HOST_DEVICE FieldViewBase& operator=(FieldViewBase&& rhs) = default;

  // These methods only make sense when we're using the ManagedArray, but
  // are why we have this class so we can call them virtually
  SPHERAL_HOST        virtual void move(chai::ExecutionSpace space) = 0;
  SPHERAL_HOST        virtual void touch(chai::ExecutionSpace space) = 0;
};

} // namespace Spheral

#endif
