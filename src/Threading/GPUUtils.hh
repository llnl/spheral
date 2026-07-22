//------------------------------------------------------------------------------
// General functions related to CHAI, RAJA, or GPU devices
//------------------------------------------------------------------------------
#ifndef __Spheral_GPUUtils__
#define __Spheral_GPUUtils__

#include "config.hh"
#include "Utilities/DBC.hh"

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"
#include "chai/managed_ptr.hpp"
#include "chai/config.hpp"
#ifdef SPHERAL_UNIFIED_MEMORY
#include "Utilities/span.hh"
#endif
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
namespace GPUUtils {
//------------------------------------------------------------------------------
// Wrappers for essential GPU device calls
//------------------------------------------------------------------------------

int deviceCount();

void initGPUs(const int stack_mult = 8);

void deviceSync();

//------------------------------------------------------------------------------
// Wrapper for chai::ManagedArray that protects against making
// one from empty data
//------------------------------------------------------------------------------
#ifndef SPHERAL_UNIFIED_MEMORY

template<typename SpanType, typename ContainerType>
void
initMAView(SpanType& a_ma, ContainerType& a_dc) {
  if (a_dc.size() == 0u) {
    a_ma.free();
  } else if ((a_dc.data() != a_ma.data(chai::CPU, false) ||
              a_dc.size() != a_ma.size())) {
    a_ma.free();
    a_ma = chai::makeManagedArray(a_dc.data(), a_dc.size(), chai::CPU, false);
  }
}

template<typename SpanType>
void
freeMAView(SpanType& a_ma) { a_ma.free(); }

template<typename SpanType>
void move(SpanType& a_ma, chai::ExecutionSpace space) { a_ma.move(space); }

template<typename SpanType>
void touch(SpanType& a_ma, chai::ExecutionSpace space) { a_ma.registerTouch(space); }

#else // SPHERAL_UNIFIED_MEMORY enabled

template<typename SpanType, typename ContainerType>
void initMAView(SpanType& a_ma, ContainerType& a_dc) { a_ma = a_dc; }

template<typename SpanType>
void freeMAView(SpanType& /*a_ma*/) { }

template<typename SpanType>
void move(SpanType& a_ma, chai::ExecutionSpace /*space*/) { }

template<typename SpanType>
void touch(SpanType& a_ma, chai::ExecutionSpace /*space*/) { }

#endif

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
} // namespace GPUUtils
} // namespace Spheral
#endif
