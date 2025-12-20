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
  ,RAJA::omp_parallel_exec
#endif
#ifdef SPHERAL_ENABLE_CUDA
  ,RAJA::cuda_exec<256>
#endif
#ifdef SPHERAL_ENABLE_HIP
  ,RAJA::hip_exec<256>
#endif
  >;

// The list of execution types we want to possibly run these tests over.
using EXEC_RESOURCE_TYPES =
  ::testing::Types<camp::list<RAJA::seq_exec, camp::resources::Host>
#ifdef SPHERAL_ENABLE_OPENMP
  ,camp::list<RAJA::omp_parallel_exec, camp::resources::Host>
#endif
#ifdef SPHERAL_ENABLE_CUDA
  ,camp::list<RAJA::cuda_exec<256>, camp::resources::Cuda>
#endif
#ifdef SPHERAL_ENABLE_HIP
  ,camp::list<RAJA::hip_exec<256>, camp::resources::Hip>
#endif
  >;
// clang-format on

#endif // SPHERAL_BASIC_EXEC_POL_HH
