//------------------------------------------------------------------------------
// General functions related to CHAI, RAJA, or GPU devices
//------------------------------------------------------------------------------
#ifndef __Spheral_GPUUtils__
#define __Spheral_GPUUtils__

#include "config.hh"
#include "DBC.hh"

#include "chai/ManagedArray.hpp"
#include "chai/managed_ptr.hpp"
#include "RAJA/RAJA.hpp"

//----------------------------------------------------------------------------
//                               GPU checks
// Checks specific to GPUs
//----------------------------------------------------------------------------

#ifdef SPHERAL_ENABLE_HIP
#define GPU_CHECK(expression)                  \
{                                              \
    const hipError_t status = expression;      \
    if(status != hipSuccess){                  \
        std::cerr << "HIP error "              \
                  << status << ": "            \
                  << hipGetErrorString(status) \
                  << " at " << __FILE__ << ":" \
                  << __LINE__ << std::endl;    \
    }                                          \
}
#define GPU_ERROR_CHECK GPU_CHECK(hipGetLastError())
#else
#define GPU_CHECK(expression)
#define GPU_ERROR_CHECK
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Wrappers for essential GPU device calls
//------------------------------------------------------------------------------

void initGPUs();

//------------------------------------------------------------------------------
// Wrapper for chai::ManagedArray that protects against making
// one from empty data
//------------------------------------------------------------------------------
template<typename DataType, typename ContainerType>
void
initMAView(chai::ManagedArray<DataType>& a_ma,
           ContainerType& a_dc) {
  if (a_dc.size() == 0u) {
    a_ma.free();
  } else if ((a_dc.data() != a_ma.data(chai::CPU, false) ||
              a_dc.size() != a_ma.size())) {
    a_ma.free();
    a_ma = chai::makeManagedArray(a_dc.data(), a_dc.size(), chai::CPU, false);
  }
}

//------------------------------------------------------------------------------
// Wrappers for RAJA atomics
//------------------------------------------------------------------------------
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
