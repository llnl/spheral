//------------------------------------------------------------------------------
// Provide a wrapper for RAJA atomics
//------------------------------------------------------------------------------
#ifndef __Spheral_Atomic_wrapper__
#define __Spheral_Atomic_wrapper__

#include "config.hh"

#include "RAJA/RAJA.hpp"

namespace Spheral {
struct AtomicAddOp {
  static SPHERAL_HOST_DEVICE inline void apply(double* dst, double value) {
    RAJA::atomicAdd<RAJA::auto_atomic>(dst, value);
  }
};

struct AtomicSubOp {
  static SPHERAL_HOST_DEVICE inline void apply(double* dst, double value) {
    RAJA::atomicSub<RAJA::auto_atomic>(dst, value);
  }
};

struct AtomicMaxOp {
  static SPHERAL_HOST_DEVICE inline void apply(double* dst, double value) {
    RAJA::atomicMax<RAJA::auto_atomic>(dst, value);
  }
};

struct AtomicMinOp {
  static SPHERAL_HOST_DEVICE inline void apply(double* dst, double value) {
    RAJA::atomicMin<RAJA::auto_atomic>(dst, value);
  }
};

}
#endif
