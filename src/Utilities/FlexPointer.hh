//---------------------------------Spheral++----------------------------------//
// Routines for instantiating and allocating classes with virtual functions
// on the device.
//
// Created by LDO, Tues Aug 19 09:28:43 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_FlexPointer__
#define __Spheral_FlexPointer__

#include "config.hh"
#include "chai/ExecutionSpaces.hpp"
#include "umpire/Umpire.hpp"
#include <utility>

namespace Spheral {

#ifdef SPHERAL_GPU_ENABLED
template<typename Derived, typename Base, typename... Args>
SPHERAL_HOST_DEVICE void newGPU(Base** d_ptr, void* mem_loc, Args&&... args) {
  RAJA::forall<RAJA::hip_exec<1>>(RAJA::TypedRangeSegment<unsigned>(0,1),
    [=] SPHERAL_HOST_DEVICE (int) {
      *d_ptr = new(mem_loc) Derived(std::forward<Args>(args)...);
    });
}

template<typename Base>
SPHERAL_HOST_DEVICE void deleteGPU(Base** d_ptr) {
  RAJA::forall<RAJA::hip_exec<1>>(RAJA::TypedRangeSegment<unsigned>(0,1),
    [=] SPHERAL_HOST_DEVICE (int) {
      delete *d_ptr;
      delete d_ptr;
    });
}
#endif

template<typename Base>
class FlexPointer {
public:

  template<typename Derived, typename... Args>
  void initialize(chai::ExecutionSpace space, Args&&... args) {
    m_space = space;
    auto& rm = umpire::ResourceManager::getInstance();
    std::string spaceStr = (m_space == chai::GPU) ? "DEVICE" : "HOST";
    auto allocator = rm.getAllocator(spaceStr);
#ifdef SPHERAL_GPU_ENABLED
    if (m_space == chai::GPU) {
      m_deviceMem = allocator.allocate(sizeof(Derived));
      m_basePtrHost = static_cast<Base**>(allocator.allocate(sizeof(Base*)));
      newGPU<Derived, Base>(m_basePtrHost, m_deviceMem, std::forward<Args>(args)...);
    } else
#endif
      {
        m_basePtrHost = static_cast<Base**>(allocator.allocate(sizeof(Base*)));
        *m_basePtrHost = new Derived(std::forward<Args>(args)...);
      }
  }

  Base* getPointer() const {
    return *m_basePtrHost;
  }

  ~FlexPointer() {
#ifdef SPHERAL_GPU_ENABLED
    if (m_space == chai::GPU) {
      auto& rm = umpire::ResourceManager::getInstance();
      auto allocator = rm.getAllocator("DEVICE");
      allocator.deallocate(m_deviceMem);
      allocator.deallocate(m_basePtrHost);
    } else
#endif
      {
        delete *m_basePtrHost;
      }
  }
private:
  chai::ExecutionSpace m_space = chai::CPU;
  Base** m_basePtrHost = nullptr;
  void* m_deviceMem = nullptr;
};

}
#endif
