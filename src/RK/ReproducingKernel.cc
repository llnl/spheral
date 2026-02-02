#include "RK/ReproducingKernel.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernel<Dimension>::
ReproducingKernel(const TableKernel<Dimension>& W,
                  const RKOrder order):
  ReproducingKernelMethods<Dimension>(order),
  mWptr(&W) {
}

//------------------------------------------------------------------------------
// Default constructor (everything NULL, invalid)
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernel<Dimension>::
ReproducingKernel():
  ReproducingKernelMethods<Dimension>(),
  mWptr(nullptr) {
}

//------------------------------------------------------------------------------
// Copy
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernel<Dimension>::
ReproducingKernel(const ReproducingKernel<Dimension>& rhs):
  ReproducingKernelMethods<Dimension>(rhs),
  mWptr(rhs.mWptr) {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernel<Dimension>&
ReproducingKernel<Dimension>::
operator=(const ReproducingKernel<Dimension>& rhs) {
  ReproducingKernelMethods<Dimension>::operator=(rhs);
  mWptr = rhs.mWptr;
  return *this;
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ReproducingKernel<Dimension>::
operator==(const ReproducingKernel<Dimension>& rhs) const {
  return (ReproducingKernelMethods<Dimension>::operator==(rhs) and
          *mWptr == *(rhs.mWptr));
}

}
