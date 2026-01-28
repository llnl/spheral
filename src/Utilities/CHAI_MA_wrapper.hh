//------------------------------------------------------------------------------
// Provide a wrapper for chai::ManagedArray that protects against making
// one from empty data
//------------------------------------------------------------------------------
#ifndef __Spheral_CHAI_MA_wrapper__
#define __Spheral_CHAI_MA_wrapper__

#include "config.hh"
#include "DBC.hh"

#include "chai/ManagedArray.hpp"
#include "chai/managed_ptr.hpp"
#include "RAJA/RAJA.hpp"
#include <RAJA/policy/sequential/policy.hpp>

namespace Spheral {

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

// Macros for updating managed_ptr member data
// TODO: Modify this to work on a list of member variables
#define ASSIGN_MEMBER(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE, EXEC_SPACE) \
  {                                                                     \
  const auto local_input = INPUT_VALUE;                                 \
    RAJA::forall<EXEC_SPACE>                                            \
      (RAJA::TypedRangeSegment<unsigned>(0,1),                          \
       [=] SPHERAL_HOST_DEVICE (int) {                                  \
         MANAGED_PTR->MEMBER_NAME = local_input;                        \
       });                                                              \
    HIP_ERROR_CHECK                                                     \
  }

#define ASSIGN_MEMBER_HOST(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE) ASSIGN_MEMBER(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE, RAJA::seq_exec);
#ifdef SPHERAL_ENABLE_HIP
#define ASSIGN_MEMBER_GPU(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE) ASSIGN_MEMBER(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE, EXEC_POLICY);
#else
#define ASSIGN_MEMBER_GPU(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE)
#endif

#define ASSIGN_MEMBER_ALL(MANAGED_PTR,  MEMBER_NAME, INPUT_VALUE) \
  ASSIGN_MEMBER_HOST(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE);      \
  ASSIGN_MEMBER_GPU(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE);

}
#endif
