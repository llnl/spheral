GPU Porting
###########

This page attempts to document some of the lessons learned and pitfalls we have experienced during the effort of porting Spheral to GPUs.

On-Device Virtual Function Lookup
=================================

Using virtual functions on-device has been a significant source of issues along the way.
At this point, ``QPiij`` in ``evaluateDerivativesImpl`` is the only virtual function call we do on the device.
The first two issues with this function call resulted in unhelpful device memory errors and disappeared if optimizations were turned off.

1. For undetermined reasons, calling RAJA atomic operations before the ``QPiij`` function in the same kernel caused a device memory error. This was fixed by simply moving all atomics below the ``QPiij`` function. We would like to know why this is necessary.

2. Any virtual function call inside the actual ``evaluateDerivativesImpl`` function call caused a device memory error. We were unable to recreate this bug with any smaller reproducer. We hypothesize this is related to register pressure and the device stack. Register pressure seems related because attempts at recreating the issue with smaller kernel calls do not show the issue. We know the device stack is part of the issue because the solution we found (other than turning optimizations off) was to increase the device stack size by calling ``hipDeviceSetLimit(hipLimitStackSize, 8*1024)``. This particular bug was nearly impossible to pin down despite communicating with many different knowledgable sources.

3. Using virtual function calls on device require careful consideration for the vtable of the object. Specifically, the object must be constructed on the device and must not be overwritten by a host instance of the object. Initial use of the ``chai::managed_ptr`` attempted to modify member data of the device object using a kernel launch. For reasons not fully understand, modifying even member data on the device caused the vtable to be made invalid or overwritten. Ultimately, the fix was to simply delete and reconstruct a new object instance on the device whenever any member data changed. Since this only occurs during problem start up, this is not expected to have much performance impact.
