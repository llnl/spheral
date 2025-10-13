//------------------------------------------------------------------------------
// Provide a wrapper for chai::ManagedArray that protects against making
// one from empty data
//------------------------------------------------------------------------------
#ifndef __Spheral_CHAI_MA_wrapper__
#define __Spheral_CHAI_MA_wrapper__

#include "chai/ManagedArray.hpp"
#include "chai/managed_ptr.hpp"

namespace Spheral {

template<typename DataType, typename ContainerType>
void
initMAView(chai::ManagedArray<DataType>& a_ma,
           ContainerType& a_dc) {
  if (a_dc.size() == 0u) {
    a_ma.free();
  } else if (a_dc.data() != a_ma.data(chai::CPU, false) ||
             a_dc.size() != a_ma.size()) {
    a_ma.free();
    a_ma = chai::makeManagedArray(a_dc.data(), a_dc.size(), chai::CPU, false);
  }
}

template<typename DerivedClass, typename BaseClass, typename... Args>
chai::managed_ptr<BaseClass>
initMP(chai::ExecutionSpace space,
       Args&&... args) {
  BaseClass* base_ptr;
#if defined(SPHERAL_GPU_ENABLED)
  if (space == chai::GPU) {
    base_ptr = chai::make_on_device<DerivedClass>(std::forward<Args>(args)...);
  } else
#endif
    {
      base_ptr = new DerivedClass(std::forward<Args>(args)...);
    }
  return chai::managed_ptr<BaseClass> ({space}, {base_ptr});
}

}
#endif
