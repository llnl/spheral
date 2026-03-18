//------------------------------------------------------------------------------
// General functions related to CHAI, RAJA, or GPU devices
//------------------------------------------------------------------------------

#include "GPUUtils.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Wrappers for essential GPU device calls
//------------------------------------------------------------------------------

void initGPUs() {
#ifdef SPHERAL_ENABLE_HIP
  size_t limitSize;
  GPU_ERROR_CHECK(hipDeviceGetLimit(&limitSize, hipLimitStackSize));
  size_t bytes = 8*1024;
  if (limitSize < bytes) {
    GPU_ERROR_CHECK(hipDeviceSetLimit(hipLimitStackSize, bytes));
  }
#endif
}
}
