#ifndef SPHERAL_BASIC_EXEC_POL_HH
#define SPHERAL_BASIC_EXEC_POL_HH

#include "config.hh"
#include "gtest/gtest.h"
#include <RAJA/policy/sequential/policy.hpp>

using SEQ_EXEC_POLICY = RAJA::seq_exec;
// The list of execution types we want to possibly run these tests over.
using EXEC_TYPES =
    ::testing::Types<camp::list<RAJA::seq_exec, camp::resources::Host>
#ifdef SPHERAL_ENABLE_CUDA
                     ,
                     camp::list<RAJA::cuda_exec<512>, camp::resources::Cuda>
#endif
#ifdef SPHERAL_ENABLE_HIP
                     ,
                     camp::list<RAJA::hip_exec<512>, camp::resources::Hip>
#endif
                     >;

#endif // SPHERAL_BASIC_EXEC_POL_HH
