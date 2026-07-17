GPU Porting
###########

This page attempts to document some of the lessons learned and pitfalls we have experienced during the effort of porting Spheral to GPUs.

On-Device Virtual Function Lookup on AMD GPUS
=============================================

Using virtual functions on-device has been a significant source of issues along the way.
At this point, ``QPiij`` in ``evaluateDerivativesImpl`` is the only virtual function call we do on the device.
The first two issues with this function call resulted in unhelpful device memory errors and disappeared if optimizations were turned off.

1. For undetermined reasons, calling RAJA atomic operations before the ``QPiij`` function in the same kernel caused a device memory error. This was fixed by simply moving all atomics below the ``QPiij`` function. We would like to know why this is necessary.

2. Any virtual function call inside the actual ``evaluateDerivativesImpl`` function call caused a device memory error. We were unable to recreate this bug with any smaller reproducer. We hypothesize this is related to register pressure and the device stack. Register pressure seems related because attempts at recreating the issue with smaller kernel calls do not show the issue. We know the device stack is part of the issue because the solution we found (other than turning optimizations off) was to increase the device stack size by calling ``hipDeviceSetLimit(hipLimitStackSize, 8*1024)``. This particular bug was nearly impossible to pin down despite communicating with many different knowledgable sources.

3. Using virtual function calls on device require careful consideration for the vtable of the object. Specifically, the object must be constructed on the device and must not be overwritten by a host instance of the object. Initial use of the ``chai::managed_ptr`` attempted to modify member data of the device object using a kernel launch. For reasons not fully understand, modifying even member data on the device caused the vtable to be made invalid or overwritten. Ultimately, the fix was to simply delete and reconstruct a new object instance on the device whenever any member data changed. Since this only occurs during problem start up, this is not expected to have much performance impact.

HIP Stack Resource Reporting
============================

ROCm assembly files saved by ``--save-temps=obj`` include per-kernel stack metadata emitted by the compiler. This is compile-time information from generated device code, not runtime profiling data. It is useful for finding which compiled targets contain kernels that use a runtime stack, and for comparing the compiler-known fixed stack portion of those kernels. It does not measure dynamic stack usage at runtime, so it cannot prove that a ``hipDeviceSetLimit(hipLimitStackSize, ...)`` value is large enough for kernels marked as using variable stack.

Configure a HIP build with::

  -DSPHERAL_HIP_REPORT_STACK_USAGE=ON

This adds ``--save-temps=obj`` for HIP compilation. The saved AMDGPU assembly metadata provides the stack fields ``.private_segment_fixed_size`` and ``.uses_dynamic_stack``. Spheral uses ``--save-temps=obj`` so parallel compiles of sources with the same basename do not race over shared temporary filenames. It also suppresses Clang's ``-Wgnu-line-marker`` diagnostic in this mode because saved HIP preprocessor outputs can contain GNU-style line markers, which otherwise fail ``ENABLE_WARNINGS_AS_ERRORS`` builds.

Interpret the report as a static inventory:

* ``largest fixed stack`` is the largest compiler-reported fixed per-thread stack/private segment size found in the selected kernels.
* ``variable-stack kernels`` counts kernels where the compiler says additional runtime stack may be used.
* ``largest fixed stack in variable kernels`` is only the known fixed portion for those variable-stack kernels. It is not the total runtime stack requirement.

The report tool does not launch kernels, sample GPU execution, or observe stack high-water marks. HIP can report the currently configured limit with ``hipDeviceGetLimit(..., hipLimitStackSize)``, but that is only the configured limit, not the stack space a kernel actually consumed. Use this report to locate targets and kernels that need runtime stack attention, then determine an adequate ``hipDeviceSetLimit(hipLimitStackSize, ...)`` value with focused runtime validation, such as a reproducer and a stack-limit sweep.

Warning: this option can increase build directory size by roughly 10x because it preserves compiler intermediate files. For example, one ROCm build grew from 6.4 GiB to 60 GiB::

  du -sh build_rocm_timers*
  6.4G    build_rocm_timers
  60G     build_rocm_timers.report

For example::

  ./scripts/devtools/host-config-build.py --host-config toss_4_x86_64_ib_cray-llvm-amdgpu@6.4.3+rocm.cmake --build-dir build_rocm_timers.report --build --nprocs 48 -DSPHERAL_ENABLE_TIMERS=ON -DENABLE_WARNINGS_AS_ERRORS=On -DSPHERAL_HIP_REPORT_STACK_USAGE=ON

The ``--host-config`` value must be the path to an existing generated host-config. If the host-config was generated in another checkout or directory, pass that relative or absolute path.

After the build, run the standalone report tool on the CMake build directory::

  python scripts/devtools/hip_stack_report.py build_rocm_timers.report/build

The report tool accepts one or more paths. Common summary modes are::

  # Whole build tree, grouped by compiled CMake target.
  python scripts/devtools/hip_stack_report.py build_rocm_timers.report/build

  # One target directory.
  python scripts/devtools/hip_stack_report.py build_rocm_timers.report/build/src/SPH/CMakeFiles/Spheral_SPH.dir

  # One saved assembly metadata file.
  python scripts/devtools/hip_stack_report.py build_rocm_timers.report/build/src/SPH/CMakeFiles/Spheral_SPH.dir/SPH_RAJAInst-hip-amdgcn-amd-amdhsa-gfx942.s

  # Compare several focused inputs in one table.
  python scripts/devtools/hip_stack_report.py \
    build_rocm_timers.report/build/tests/cpp/Loops/CMakeFiles/sph_tests.dir \
    build_rocm_timers.report/build/src/SPH/CMakeFiles/Spheral_SPH.dir/SPH_RAJAInst-hip-amdgcn-amd-amdhsa-gfx942.s \
    build_rocm_timers.report/build/src/SPH/CMakeFiles/Spheral_SPH.dir/SolidSPH_RAJAInst-hip-amdgcn-amd-amdhsa-gfx942.s

The default output is a compact summary table. Broad CMake build directories are grouped by compiled CMake targets using ``CMakeFiles/TargetDirectories.txt``, so running on the top-level build directory can produce one row per library or executable target. This includes object-library targets used by Spheral package libraries. Compiler-identification and other non-target saved assembly files are ignored. Direct target-directory inputs are summarized as one row, and specific saved assembly file inputs include the cleaned source name in the component label. When more than one component is reported, an ``OVERALL`` row is appended::

  component                              kernels  variable-stack kernels  largest fixed stack  largest fixed stack in variable kernels
  -------------------------------------  -------  ----------------------  -------------------  ---------------------------------------
  tests/cpp/Loops/sph_tests                 127                       2               5824 B                                    624 B
  src/SPH/Spheral_SPH/SPH_RAJAInst          133                       6               5824 B                                   3616 B
  src/SPH/Spheral_SPH/SolidSPH_RAJAInst     133                       6               5824 B                                   4752 B
  OVERALL                                   393                      14               5824 B                                   4752 B

  Fixed-stack columns are compiler-reported known per-thread stack use.
  Variable-stack kernels may use additional runtime stack that is not bounded by this metadata.

Use ``--details`` when detailed per-kernel rows are needed. The detailed output is tab-separated so demangled C++ template names can contain commas without splitting columns::

  python scripts/devtools/hip_stack_report.py build_rocm_timers.report/build --details > hip-stack-details.tsv
  python scripts/devtools/hip_stack_report.py build_rocm_timers.report/build --details --top 50 > hip-stack-details.tsv

If any input has variable-stack kernels, use ``--details`` to identify the specific kernels that require runtime stack attention. The required runtime stack limit still needs to be chosen from runtime validation or another measurement source.

The ``--save-temps`` flag affects build time, compiler output volume, and disk usage, but it should not change runtime performance of the generated kernels. Use ``SPHERAL_HIP_REPORT_STACK_USAGE`` only for diagnostic builds. The runtime stack limit can affect GPU memory pressure, so keep it no larger than the measured workload requires.

The AMD ROCm advanced profiling guide has more context on kernel metadata:
https://rocm.blogs.amd.com/software-tools-optimization/profiling-guide/advanced/README.html
