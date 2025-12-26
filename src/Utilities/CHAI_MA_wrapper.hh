//------------------------------------------------------------------------------
// Provide a wrapper for chai::ManagedArray that protects against making
// one from empty data
//------------------------------------------------------------------------------
#ifndef __Spheral_CHAI_MA_wrapper__
#define __Spheral_CHAI_MA_wrapper__

#include "config.hh"

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
  do {                                                                  \
    /* Get the object type from the pointer (removes 'volatile' and '&' if present) */ \
    using ObjectType = std::remove_cv_t<std::remove_reference_t<decltype(*(MANAGED_PTR))>>; \
    /* Get the member type by accessing the member */                   \
    using MemberType = decltype(std::declval<ObjectType>().MEMBER_NAME); \
    /* Create local instance of INPUT_VALUE to ensure it is captured */ \
    const MemberType local_input = INPUT_VALUE;                         \
    /* Create the pointer-to-member type */                             \
    MemberType ObjectType::* member_ptr = &ObjectType::MEMBER_NAME;     \
    chai::managed_ptr<ObjectType> local_ptr(MANAGED_PTR);               \
    RAJA::forall<EXEC_SPACE>                                            \
      (RAJA::TypedRangeSegment<unsigned>(0,1),                          \
       [=] SPHERAL_HOST_DEVICE (int) {                                  \
         local_ptr.get()->*member_ptr = local_input;                    \
       });                                                              \
  } while(0)

#define ASSIGN_MEMBER_HOST(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE) ASSIGN_MEMBER(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE, RAJA::seq_exec);
#ifdef SPHERAL_ENABLE_HIP
#define ASSIGN_MEMBER_GPU(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE) ASSIGN_MEMBER(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE, RAJA::hip_exec<512>);
#else
#define ASSIGN_MEMBER_GPU(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE)
#endif

#define ASSIGN_MEMBER_ALL(MANAGED_PTR,  MEMBER_NAME, INPUT_VALUE) \
  ASSIGN_MEMBER_HOST(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE);      \
  ASSIGN_MEMBER_GPU(MANAGED_PTR, MEMBER_NAME, INPUT_VALUE);

}
#endif
