#ifndef SPHERAL_BASIC_EXEC_POL_HH
#define SPHERAL_BASIC_EXEC_POL_HH

#include "config.hh"
#include "test-base.hh"
#include <RAJA/policy/sequential/policy.hpp>

using SEQ_EXEC_POLICY = RAJA::seq_exec;

// clang-format off
using EXEC_TYPES = camp::list<
  RAJA::seq_exec
#ifdef SPHERAL_ENABLE_OPENMP
  ,RAJA::omp_parallel_for_exec
#endif
#ifdef SPHERAL_ENABLE_CUDA
  ,RAJA::cuda_exec<GPU_BLOCK_SIZE>
#endif
#ifdef SPHERAL_ENABLE_HIP
  ,RAJA::hip_exec<GPU_BLOCK_SIZE>
#endif
  >;

// The list of execution types we want to possibly run these tests over.
using EXEC_RESOURCE_TYPES =
  ::testing::Types<camp::list<RAJA::seq_exec, camp::resources::Host>
#ifdef SPHERAL_ENABLE_OPENMP
  ,camp::list<RAJA::omp_parallel_for_exec, camp::resources::Host>
#endif
#ifdef SPHERAL_ENABLE_CUDA
  ,camp::list<RAJA::cuda_exec<GPU_BLOCK_SIZE>, camp::resources::Cuda>
#endif
#ifdef SPHERAL_ENABLE_HIP
  ,camp::list<RAJA::hip_exec<GPU_BLOCK_SIZE>, camp::resources::Hip>
#endif
  >;

using GPU_TEST_TYPE =
#ifdef SPHERAL_ENABLE_CUDA
  RAJA::cuda_exec<GPU_BLOCK_SIZE>;
#elif defined(SPHERAL_ENABLE_HIP)
  RAJA::hip_exec<GPU_BLOCK_SIZE>;
#else
  std::nullptr_t;
#endif
// clang-format on

#endif // SPHERAL_BASIC_EXEC_POL_HH
