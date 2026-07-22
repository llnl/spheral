//---------------------------------Spheral++----------------------------------//
// ManagedView -- wrap a View type and provide forwards for move and touch
//
// Created by JMO, Mon Jul 20 13:45:32 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_ManagedView__
#define __Spheral_ManagedView__

#include "Threading/ManagedViewBase.hh"

namespace Spheral {

template<typename ViewType>
class ManagedView:
    public ManagedViewBase {

public:
  //--------------------------- Public Interface ---------------------------//

  // Constructors, destructor
  ManagedView(ViewType& obj): ManagedViewBase(), mView(obj) {}
  ManagedView(const ManagedView& rhs) = default;
  ManagedView(ManagedView&& rhs) = default;
  virtual ~ManagedView() = default;

  // Assignment
  ManagedView& operator=(ManagedView& rhs) = default;
  ManagedView& operator=(ManagedView&& rhs) = default;

  // These methods only make sense when we're using the ManagedArray, but
  // are why we have this class so we can call them virtually
  virtual void move(chai::ExecutionSpace space)  { mView.move(space); }
  virtual void touch(chai::ExecutionSpace space) { mView.touch(space); }

private:
  //--------------------------- Private Interface ---------------------------//
  ViewType mView;
};

} // namespace Spheral

#endif
