//---------------------------------Spheral++----------------------------------//
// ManagedViewBase -- provide a handle with virtual methods View classes can
// use to be managed by the ViewManager.
//
// Created by JMO, Mon Jul  6 08:54:21 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_ManagedViewBase__
#define __Spheral_ManagedViewBase__

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

namespace Spheral {

class ManagedViewBase:
    public chai::CHAICopyable {

public:
  //--------------------------- Public Interface ---------------------------//

  // Constructors, destructor
  ManagedViewBase() = default;
  ManagedViewBase(const ManagedViewBase& rhs) = default;
  ManagedViewBase(ManagedViewBase&& rhs) = default;
  virtual ~ManagedViewBase() = default;

  // Assignment
  ManagedViewBase& operator=(ManagedViewBase& rhs) = default;
  ManagedViewBase& operator=(ManagedViewBase&& rhs) = default;

  // These methods only make sense when we're using the ManagedArray, but
  // are why we have this class so we can call them virtually
  virtual void move(chai::ExecutionSpace space) = 0;
  virtual void touch(chai::ExecutionSpace space) = 0;
};

} // namespace Spheral

#endif
