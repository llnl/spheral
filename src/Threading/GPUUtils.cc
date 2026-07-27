//------------------------------------------------------------------------------
// General functions related to CHAI, RAJA, or GPU devices
//------------------------------------------------------------------------------

#include "GPUUtils.hh"

namespace Spheral {
namespace GPUUtils {
//------------------------------------------------------------------------------
// Wrappers for essential GPU device calls
//------------------------------------------------------------------------------

int deviceCount() {
  int device_count = 0;
#ifdef SPHERAL_ENABLE_HIP
  GPU_CHECK(hipGetDeviceCount(&device_count));
#endif
  return device_count;
}

void initGPUs(const int stack_mult) {
#ifdef SPHERAL_ENABLE_HIP
  int device_count = deviceCount();
  if (device_count > 0) {
    size_t limitSize;
    GPU_CHECK(hipDeviceGetLimit(&limitSize, hipLimitStackSize));
    size_t bytes = stack_mult*1024;
    if (limitSize != bytes) {
      GPU_CHECK(hipDeviceSetLimit(hipLimitStackSize, bytes));
    }
  }
#endif
}

void deviceSync() {
#ifdef SPHERAL_ENABLE_HIP
  GPU_CHECK(hipDeviceSynchronize());
#endif
}
} // namespace GPUUtils
} // namespace Spheral
