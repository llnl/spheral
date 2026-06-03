RAJA/CHAI Execution Patterns
============================

Purpose
-------

This document describes how Spheral uses RAJA and CHAI in the current C++
implementation. It focuses on the execution and memory-management patterns used
by device-capable physics kernels, especially the RAJA SPH implementations.
For applied guidance on using these patterns when porting new code, see
:doc:`../dev/gpu_porting_patterns`.

The important point for new developers is that RAJA and CHAI are not isolated
utility libraries in Spheral. They shape the object model:

* data containers expose lightweight views;
* views are marked host/device and are copied into kernels;
* CHAI manages movement of view-backed storage;
* RAJA supplies the execution policy, loop launch, and atomics.

Source Map
----------

Configuration and helpers:

* ``src/config.hh.in``
* ``src/Utilities/GPUUtils.hh``
* ``docs/developer/dev/gpu_dev.rst``
* ``docs/developer/dev/gpu_porting_patterns.rst``

Primary RAJA hydro implementations:

* ``src/SPH/SPH_RAJA.hh`` and ``src/SPH/SPH_RAJA.cc``
* ``src/SPH/SolidSPH_RAJA.hh`` and ``src/SPH/SolidSPH_RAJA.cc``

View objects used in RAJA kernels:

* ``src/Field/FieldView.hh``
* ``src/Field/FieldListView.hh``
* ``src/Kernel/TableKernelView.hh``
* ``src/Neighbor/NodePairListView.hh``
* ``src/Neighbor/PairwiseFieldView.hh``
* ``src/ArtificialViscosity/ArtificialViscosityView.hh``

Representative managed-pointer implementations:

* ``src/ArtificialViscosity/MonaghanGingoldViscosity.hh``
* ``src/ArtificialViscosity/TensorMonaghanGingoldViscosity.hh``
* ``src/ArtificialViscosity/LimitedMonaghanGingoldViscosity.hh``
* ``src/ArtificialViscosity/FiniteVolumeViscosity.hh``

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

This means most device-capable loops are written once:

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
2. Create view objects with ``view()``.
3. Move or touch views when the backend requires explicit movement.
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
owning objects stay on the host side of the launch boundary; views cross into
the kernel.

SPH_RAJA Derivative Flow
------------------------

``SPH_RAJA::evaluateDerivatives`` first dispatches on the artificial viscosity
return type:

* scalar viscosity uses
  ``chai::managed_ptr<ArtificialViscosityView<Dimension, Scalar>>``;
* tensor viscosity uses
  ``chai::managed_ptr<ArtificialViscosityView<Dimension, Tensor>>``.

The templated ``evaluateDerivativesImpl`` then performs three broad phases.

Initial setup
  Kernel views are created from ``TableKernel`` objects. The active
  ``ConnectivityMap`` supplies a ``NodePairList`` view. State and derivative
  field lists are looked up by key and converted to ``FieldListView`` objects.
  The code checks that field-list sizes match the number of node lists.

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

* each pair is visited once through ``NodePairList``;
* contributions for both endpoints are accumulated in the same kernel
  iteration;
* output fields that receive many pair contributions must be atomic.

This design avoids separate gather/scatter passes but makes atomic placement
and data-race analysis part of kernel development.

Artificial Viscosity on Device
------------------------------

Artificial viscosity is a special case because the physics package needs
runtime polymorphism in the pair loop. The base view interface is
``ArtificialViscosityView<Dimension, QPiType>``. It provides the virtual
``SPHERAL_HOST_DEVICE`` method:

::

   QPiij(QPiType& QPiij, QPiType& QPiji,
         Scalar& Qij, Scalar& Qji,
         ...,
         const FieldListView<Dimension, Scalar>& fCl,
         const FieldListView<Dimension, Scalar>& fCq,
         const FieldListView<Dimension, Tensor>& DvDx) const

Managed Pointer Dispatch for Device Virtual Calls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dispatch path is deliberately split between host-side type selection and
device-side virtual behavior:

::

   SPH_RAJA::evaluateDerivatives
     asks ArtificialViscosity owner for QPiTypeIndex()
       Scalar path -> getScalarView()
       Tensor path -> getTensorView()
     receives chai::managed_ptr<ArtificialViscosityView<Dimension, QPiType>>
     calls templated evaluateDerivativesImpl(..., Q)
       RAJA pair kernel captures Q by value
       kernel calls Q->QPiij(...)
         device virtual dispatch resolves concrete viscosity view

The host object remains the durable owner of viscosity parameters and restart
state. The managed view object is the device-callable representation used by
the RAJA kernel:

::

   MonaghanGingoldViscosity
     owns host-side parameters
     owns chai::managed_ptr<MonaghanGingoldViscosityView> m_viewPtr
       |
       | getScalarView()
       v
   chai::managed_ptr<ArtificialViscosityView<Dimension, Scalar>>
       |
       | captured in RAJA lambda
       v
   Q->QPiij(...) on device

``ArtificialViscosity`` to ``ArtificialViscosityView``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This conversion follows the same owner/view idea as ``Field`` and
``FieldView``, but with one crucial difference: the view is a managed,
polymorphic object, not just a shallow data span.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Concern
     - Owner-side role
     - Managed view-side role
   * - Durable physics object
     - ``ArtificialViscosity`` descendants own host parameters, package-facing
       API, restart identity, and configuration setters.
     - ``ArtificialViscosityView`` descendants contain only the data and
       methods needed by device kernels.
   * - Runtime type selection
     - The owner reports ``QPiTypeIndex()`` so ``SPH_RAJA`` can choose the
       scalar or tensor templated path on the host.
     - The chosen base view type fixes the ``QPiType`` seen by
       ``evaluateDerivativesImpl``.
   * - Concrete behavior
     - Concrete owners such as ``MonaghanGingoldViscosity`` decide which view
       class represents them.
     - Concrete views override virtual ``QPiij`` with ``SPHERAL_HOST_DEVICE``
       behavior.
   * - Allocation and lifetime
     - The owner stores ``chai::managed_ptr<ViewType> m_viewPtr`` and lazily
       constructs it with ``chai::make_managed``.
     - The managed view is the object whose vtable must be valid on device.
   * - Mutation
     - Host setters update owner-side parameter values.
     - If a managed view already exists, it is freed and reconstructed rather
       than modified in place.
   * - Kernel use
     - The owner is never captured by the RAJA pair loop.
     - A ``chai::managed_ptr<ArtificialViscosityView<...>>`` is captured by
       value and used as ``Q->QPiij(...)``.

The structural conversion is:

::

   host ArtificialViscosity owner
     stores durable parameters
     lazily owns chai::managed_ptr<ConcreteView>
     |
     | getScalarView() or getTensorView()
     v
   chai::managed_ptr<ArtificialViscosityView<Dimension, QPiType>>
     base pointer with device-valid concrete view object
     |
     | captured by RAJA lambda
     v
   virtual QPiij dispatch on device

This differs from ``FieldView`` in two ways. First, the managed view is
constructed as an object with behavior, not only as a handle over owner
storage. Second, its virtual table is part of correctness. That is why Spheral
reconstructs the managed view when owner-side viscosity parameters change.

Concrete viscosity owners, such as ``MonaghanGingoldViscosity`` and
``TensorMonaghanGingoldViscosity``, store a
``chai::managed_ptr<ViewType> m_viewPtr``. They lazily initialize it with
``chai::make_managed<ViewType>(...)`` and return a cast pointer to the base
view interface.

When host-side viscosity parameters change, the owner frees and reconstructs
the managed view instead of modifying the existing device object in place. This
is a defensive design choice based on problems observed with device virtual
function tables, especially on AMD GPUs.

The current GPU development notes identify three constraints that matter for
new code:

* device virtual calls require objects whose vtables are valid on device;
* changing managed object member data in place can invalidate or corrupt the
  device-side virtual dispatch path;
* in ``evaluateDerivativesImpl``, RAJA atomics are intentionally placed after
  the ``QPiij`` call because earlier ordering caused device memory errors.

The practical guidance is to avoid adding new virtual calls in device kernels
unless there is no reasonable alternative. If device polymorphism is required,
use managed view objects and follow the artificial-viscosity reconstruction
pattern.

Kernel and Interpolator Views
-----------------------------

``TableKernel`` follows the same owner/view split as fields. The owning
``TableKernel`` stores tabulated data and returns a ``TableKernelView``. The
view contains interpolator views and exposes host/device methods such as:

* ``kernelValue``;
* ``gradValue``;
* ``grad2Value``;
* ``kernelAndGradValue``;
* ``equivalentNodesPerSmoothingScale``.

``TableKernelView::move`` recursively moves its interpolator views. This is the
same layered movement problem seen in ``FieldListView``: moving the outer view
is not enough if the view contains nested managed data.

Node-Pair Views
---------------

RAJA SPH kernels consume flattened pair storage, not the nested
``ConnectivityMap`` vectors. The active path is:

1. ``ConnectivityMap`` builds a ``NodePairList``.
2. ``NodePairList`` owns ``std::vector<NodePairIdxType>``.
3. ``NodePairList::initView`` exposes it through a CHAI/spanned view.
4. ``NodePairListView`` is captured by the RAJA pair loop.

This keeps kernel traversal regular: ``kk`` ranges from ``0`` to ``npairs`` and
each ``pairs[kk]`` describes one pair interaction.

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

Movement Responsibilities
-------------------------

Movement policy is distributed by container type:

``FieldView``
  Moves or touches the field's data span/managed array.

``FieldListView``
  Moves or touches the array of field views and, recursively by default, each
  field view's data.

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

When adding a new view type:

* keep it shallow and cheap to copy;
* keep ownership out of the view;
* annotate only the methods intended for device use;
* expose ``move``/``touch`` when it owns managed data;
* provide a clear owner method that rebuilds the view span after storage
  changes.

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

Relationship to the Other Design Docs
-------------------------------------

``value_view_and_device_execution_model.rst`` describes the owner/view pattern
that makes these kernels possible.

``connectivity_data_structures.rst`` describes how ``ConnectivityMap`` creates
the flattened ``NodePairList`` traversed by RAJA pair loops.

``integrator_and_state_update_model.rst`` describes when physics packages are
called and where derivative fields enter the integrator update policy system.
