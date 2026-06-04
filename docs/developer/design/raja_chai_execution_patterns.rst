RAJA/CHAI Execution Patterns
============================

Purpose
-------

This document describes how Spheral uses RAJA and CHAI around current
device-capable physics kernels, especially the RAJA SPH implementations. It is
the execution companion to :doc:`value_view_and_device_execution_model` and
:doc:`value_view_conversion_case_studies`.

The object-model details live in those pages. This page focuses on launch-time
mechanics:

* selecting the RAJA execution policy;
* creating views from host owners;
* moving or touching CHAI-backed storage;
* capturing views and managed pointers in kernels;
* using atomics in pair loops;
* returning written data to CPU consumers.

Source Map
----------

Configuration and helpers:

* ``src/config.hh.in``
* ``src/Utilities/GPUUtils.hh``

Primary RAJA hydro implementations:

* ``src/SPH/SPH_RAJA.hh`` and ``src/SPH/SPH_RAJA.cc``
* ``src/SPH/SolidSPH_RAJA.hh`` and ``src/SPH/SolidSPH_RAJA.cc``

Representative captured objects:

* ``src/Field/FieldListView.hh``
* ``src/Kernel/TableKernelView.hh``
* ``src/Neighbor/NodePairListView.hh``
* ``src/Neighbor/PairwiseFieldView.hh``
* ``src/ArtificialViscosity/ArtificialViscosityView.hh``

Build-Time Execution Model
--------------------------

``src/config.hh.in`` defines the small set of macros and aliases that most
device-capable code uses:

``SPHERAL_HOST_DEVICE``
  Expands to ``RAJA_HOST_DEVICE``.

``SPHERAL_HOST`` and ``SPHERAL_DEVICE``
  Expand to RAJA host/device annotations.

``GPU_BLOCK_SIZE``
  Currently set to ``256``.

``EXEC_POLICY``
  Selects the default ``RAJA::forall`` policy:

  * HIP builds use ``RAJA::hip_exec<GPU_BLOCK_SIZE>``;
  * CUDA builds use ``RAJA::cuda_exec<GPU_BLOCK_SIZE>``;
  * OpenMP builds use ``RAJA::omp_parallel_for_exec``;
  * otherwise Spheral uses ``RAJA::seq_exec``.

``TRS_UINT``
  Alias for ``RAJA::TypedRangeSegment<size_t>``.

Most device-capable loops are written once:

::

   RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, n),
     [=] SPHERAL_HOST_DEVICE(size_t i) {
       ...
     });

The selected backend is a build configuration decision, not a physics-package
decision.

CHAI's Role
-----------

CHAI provides two related services in Spheral:

``chai::ManagedArray``
  Wraps storage that can be moved between CPU and GPU execution spaces. In
  non-UVM builds, many Spheral views contain a ``ManagedArray``.

``chai::managed_ptr``
  Holds managed objects whose methods may be called on device. Spheral uses
  this for artificial-viscosity view objects where a kernel needs polymorphic
  ``QPiij`` behavior.

Spheral wraps common CHAI operations in ``GPUUtils.hh``:

* ``initMAView`` creates or refreshes a managed-array view over owner storage;
* ``freeMAView`` frees the managed-array view;
* ``move`` and ``touch`` abstract UVM vs non-UVM behavior;
* atomic operation wrappers call RAJA atomics with ``RAJA::auto_atomic``.

In unified-memory builds, the same wrapper APIs compile to span/no-op behavior
where movement is unnecessary. This keeps most container code independent of
the memory model.

Kernel Setup Pattern
--------------------

The current RAJA hydro code follows a repeated setup pattern:

1. Read durable host objects from the package, ``DataBase``, ``State``, and
   ``StateDerivatives``.
2. Create view objects with ``view()`` or managed view pointer accessors.
3. Move or touch captured objects when the backend requires explicit movement.
4. Launch one or more ``RAJA::forall`` loops.
5. Move written views back to CPU if later host code will consume them.

A simplified version of the pattern in ``SPH_RAJA::evaluateDerivativesImpl`` is:

::

   auto W_view = W.view();
   auto WQ_view = WQ.view();

   const auto& pairs_owner = connectivityMap.nodePairList();
   const auto pairs = pairs_owner.view();

   auto mass_v = state.fields(HydroFieldNames::mass, 0.0);
   auto DvDt_v = derivs.fields(HydroFieldNames::hydroAcceleration,
                               Vector::zero());

   auto mass = mass_v.view();
   auto DvDt = DvDt_v.view();

   RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, pairs.size()),
     [=] SPHERAL_HOST_DEVICE(size_t kk) {
       const auto pair = pairs[kk];
       const auto mi = mass(pair.i_list, pair.i_node);
       ...
       DvDt(pair.i_list, pair.i_node).atomicSub(...);
     });

   DvDt.move(chai::CPU);

The actual implementation has many more fields, but the structure is the same:
owning objects stay on the host side of the launch boundary; views and managed
view pointers cross into the kernel.

SPH_RAJA Derivative Flow
------------------------

``SPH_RAJA::evaluateDerivatives`` first dispatches on the artificial-viscosity
return type:

* scalar viscosity uses
  ``chai::managed_ptr<ArtificialViscosityView<Dimension, Scalar>>``;
* tensor viscosity uses
  ``chai::managed_ptr<ArtificialViscosityView<Dimension, Tensor>>``.

The templated ``evaluateDerivativesImpl`` then performs three broad phases.

Initial setup
  Kernel views are created from ``TableKernel`` objects. The active
  ``ConnectivityMap`` supplies a ``NodePairList`` owner, whose view is captured
  by the pair loop. State and derivative field lists are looked up by key and
  converted to ``FieldListView`` objects. The code checks that field-list sizes
  match the number of node lists.

Pair loop
  A RAJA loop over ``npairs`` computes pairwise SPH contributions. Each kernel
  iteration reads one ``NodePairIdxType`` containing ``i_node``, ``i_list``,
  ``j_node``, and ``j_list``. It reads state for both nodes, evaluates kernels,
  calls ``Q->QPiij(...)``, and atomically accumulates pair contributions into
  per-node derivative fields.

Per-node finalization
  A second set of RAJA loops walks internal nodes per node list. These loops add
  self-contributions, finish velocity gradients, compute continuity
  derivatives, finish XSPH position evolution, and apply other node-local
  finalization.

After the loops, derivative views are moved back to ``chai::CPU``. This is
necessary because the integrator and later package hooks may run host-side code
that reads the derivative fields.

SolidSPH_RAJA follows the same broad shape, with more state and derivative
fields for solid mechanics. It also shows explicit ``move(chai::GPU)`` calls for
many views under HIP builds before launching kernels.

Pair Loops and Atomics
----------------------

The SPH pair loop updates both nodes in an interacting pair. Different pairs
can contribute to the same destination node concurrently, so the kernel must
use atomic accumulation for shared per-node outputs.

Spheral uses two forms of atomic update:

* scalar wrappers such as ``GPUUtils::AtomicAddOp`` and
  ``GPUUtils::AtomicMaxOp``;
* data-type methods such as ``Vector::atomicAdd`` and ``Tensor::atomicSub``
  where geometry types provide component-wise atomics.

The pair loop therefore has a deliberate asymmetry:

* each pair is visited once through ``NodePairListView``;
* contributions for both endpoints are accumulated in the same kernel
  iteration;
* output fields that receive many pair contributions must be atomic.

This design avoids separate gather/scatter passes but makes atomic placement
and data-race analysis part of kernel development.

Managed Dispatch in the Launch Path
-----------------------------------

Artificial viscosity is the current kernel path that captures a managed pointer
to a polymorphic view. The object-model contract is described in
:doc:`value_view_and_device_execution_model`; the launch path is:

1. ``SPH_RAJA::evaluateDerivatives`` asks the host viscosity owner for its
   ``QPiTypeIndex``.
2. Host code selects the scalar or tensor templated path.
3. The owner returns a
   ``chai::managed_ptr<ArtificialViscosityView<Dimension, QPiType>>`` through
   ``getScalarView()`` or ``getTensorView()``.
4. The RAJA pair kernel captures that managed base pointer by value.
5. The kernel calls ``Q->QPiij(...)`` through the device-valid view object.

The owner remains the durable holder of viscosity parameters and restart state.
The managed view is reconstructed when host-side parameters change, rather than
being modified in place, so the device-side virtual dispatch path remains
valid.

Captured View Movement
----------------------

Movement policy is distributed by captured object type:

``FieldView`` and ``FieldListView``
  Move or touch field data. ``FieldListView`` movement may be recursive because
  the outer view contains nested field views.

``TableKernelView``
  Moves nested interpolator views.

``NodePairListView``
  Moves or touches the pair array.

``PairwiseFieldView``
  Moves or touches the pairwise-value array.

``chai::managed_ptr`` views
  Managed object construction and lifetime are controlled by the owning host
  object, such as an artificial-viscosity model.

The caller that launches a kernel is responsible for ensuring all captured
views are valid in the target execution space. In some paths the first access
through CHAI may trigger movement; in others, especially HIP-focused code,
explicit movement is used before the launch.

Device-to-host movement is a consumer-boundary decision, not an automatic
post-kernel step. Written views only need to be moved back to ``chai::CPU`` when
the next consumer is non-RAJAfied host code that will read those values. If the
next consumer is another RAJA/device-capable path using views, the data can
remain in the active execution space until a host-only boundary is reached.

CPU, OpenMP, HIP, and CUDA Behavior
-----------------------------------

Because Spheral uses ``EXEC_POLICY``, the same source code can execute under
several backends:

* sequential CPU builds use ``RAJA::seq_exec``;
* OpenMP builds use ``RAJA::omp_parallel_for_exec``;
* HIP builds use ``RAJA::hip_exec``;
* CUDA builds use ``RAJA::cuda_exec``.

That portability does not make all code backend-neutral automatically. Kernel
bodies still need to obey device restrictions:

* capture simple view objects by value;
* avoid host-only APIs;
* avoid unannotated helper functions in kernel code;
* use atomics for shared output locations;
* ensure managed data has been moved or touched correctly for the backend;
* keep virtual device calls rare and carefully managed.

Extension Guidance
------------------

When adding a new RAJA kernel:

* gather all owner objects before the launch;
* convert fields, kernels, pair lists, and scratch data to views before the
  launch;
* move or touch the views needed by the target backend;
* capture views by value in the RAJA lambda;
* use ``SPHERAL_HOST_DEVICE`` helpers only inside the lambda;
* avoid STL containers, host iterators, ``FieldBase``, and ``DataBase`` access
  inside the lambda;
* use RAJA/Spheral atomics for shared outputs;
* move written data back to CPU before host code reads it;
* reacquire views after any storage resize or connectivity rebuild.

Common Failure Modes
--------------------

Missing host/device annotation
  A helper used in a RAJA lambda is host-only. Mark small kernel helpers
  ``SPHERAL_HOST_DEVICE`` or keep them outside the kernel.

Hidden host object capture
  Capturing ``this`` or a rich owner object can pull host-only state into a
  device kernel. Capture views and primitive values instead.

Non-atomic shared writes
  Pair loops write many contributions into per-node fields. Any shared
  destination must use an atomic update.

Stale CHAI view
  The owner storage changed after the view was created. Rebuild the view by
  calling ``view()`` again after the storage/layout change.

Host reads GPU-written data
  A derivative field was written in a device kernel but not moved/touched back
  before host code read it.

Device virtual dispatch instability
  Managed polymorphic objects must be constructed in a way that preserves the
  device vtable. Follow the artificial-viscosity managed-view pattern.
