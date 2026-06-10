Value/View and Device Execution Model
=====================================

Purpose
-------

This document defines the object shapes Spheral currently uses at the boundary
between long-lived host objects and device-capable RAJA kernels. It is the
contract page for the pattern; concrete class families are described in
:doc:`value_view_conversion_case_studies`, and launch mechanics are described in
:doc:`raja_chai_execution_patterns`.

In this context, a kernel-facing object is the object that setup code copies or
passes into a RAJA kernel. It may be a small view value, or it may be a managed
view object whose host/device methods are called from the kernel. The matching
host object is included when it owns the data or configuration and produces the
kernel-facing object.

The recurring split is:

* long-lived C++ host objects own storage, identity, registration, restart
  behavior, resizing, and host-side APIs;
* kernel-facing objects are lightweight, copyable, and limited to data and
  behavior that can participate in host/device kernels;
* kernels receive views or managed view pointers rather than the full host
  object graph.

Source Map
----------

Core infrastructure:

* ``src/config.hh.in``
* ``src/Utilities/GPUUtils.hh``

Representative value/view families:

* ``src/Field/Field.hh`` and ``src/Field/FieldView.hh``
* ``src/Field/FieldList.hh`` and ``src/Field/FieldListView.hh``
* ``src/Neighbor/NodePairList.hh`` and ``src/Neighbor/NodePairListView.hh``
* ``src/Neighbor/PairwiseField.hh`` and ``src/Neighbor/PairwiseFieldView.hh``
* ``src/Kernel/TableKernel.hh`` and ``src/Kernel/TableKernelView.hh``

Representative family that uses managed pointers for behavior selected at
runtime:

* ``src/ArtificialViscosity/ArtificialViscosity.hh``
* ``src/ArtificialViscosity/ArtificialViscosityView.hh``
* concrete artificial-viscosity classes and view implementations

Kernel-Facing Object Shapes
---------------------------

Device kernels need a narrower representation than the host object usually
contains. The current examples expose:

* small values that can be captured by value in a RAJA lambda;
* methods annotated with ``SPHERAL_HOST_DEVICE`` where they run in kernels;
* spans, CHAI-managed arrays, or managed pointers valid in the selected
  execution space;
* primitive metadata needed by the kernel;
* no traversal of host-only object graphs inside the kernel body.

Two related shapes appear in the current code.

Value/View Shape
~~~~~~~~~~~~~~~~

The value/view shape appears where device code needs direct access to data
owned by an object that also carries host-only responsibilities.

::

   host object
     owns storage, identity, registration, restart, resizing, host API
     rebinds the view when storage or layout changes
     |
     | view()
     v
   view
     contains spans, managed arrays, primitive metadata, small access API
     is cheap to copy
     has SPHERAL_HOST_DEVICE methods where needed
     |
     | captured by value
     v
   RAJA kernel

The host object keeps responsibilities such as:

* storage lifetime;
* registration with ``NodeList``, ``State``, or package infrastructure;
* restart and serialization;
* resizing, deletion, reordering, and layout changes;
* Python-facing semantics;
* host-only lookup maps and helper containers.

The view keeps responsibilities such as:

* indexed data access;
* primitive metadata needed by the algorithm;
* ``move``, ``touch``, and ``data`` hooks when it owns CHAI-managed data;
* small helper methods that are valid on host and device.

Managed Pointer Shape for Runtime-Selected Behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shape appears where device code needs to call behavior selected at
runtime. The artificial-viscosity family is the current example.

::

   host object
     selects concrete behavior
     owns long-lived parameters and restart state
     lazily owns chai::managed_ptr<ConcreteView>
     |
     | getScalarView() or getTensorView()
     v
   chai::managed_ptr<BaseView>
     points at device-valid concrete view object
     |
     | captured by RAJA lambda
     v
   virtual SPHERAL_HOST_DEVICE method

Prefer host-side type selection or concrete value views when the kernel can be
specialized before launch. Use a managed pointer with a virtual method call only
when the runtime choice must remain inside the device kernel. In the current
code, this is used for ``ArtificialViscosityView::QPiij``.

The managed view is an object with behavior, not just a span over host-object
storage. The C++ virtual-method machinery has to be valid on device, so the
host object constructs and makes the managed view object current when its
concrete parameters change.

Device Memory Semantics
-----------------------

Spheral abstracts host/device storage with a small number of conventions.

``SPHERAL_HOST_DEVICE``
  Expands to RAJA's host/device annotation. It appears on small methods that run
  inside RAJA kernels.

``SPHERAL_UNIFIED_MEMORY`` / ``USE_UVM``
  Select span/UVM behavior. In this mode many explicit CHAI movement operations
  become no-ops because the pointer is already unified.

``chai::ManagedArray``
  Used in non-UVM builds to wrap host-owned storage and make it movable between
  execution spaces.

``chai::managed_ptr``
  Used for managed device-callable objects, currently artificial-viscosity view
  objects.

``move(space)``
  Requests that CHAI move data to ``chai::CPU`` or ``chai::GPU``.

``touch(space)``
  Registers that data was accessed in an execution space without necessarily
  moving immediately.

``data(space, do_move)``
  Returns a pointer in the requested execution space, optionally moving first.

The host-object/view pattern does not eliminate the need to reason about data
movement. It gives each container one local place to implement movement
behavior, and it gives kernel setup code a clear set of captured objects to move
or touch before launch.

View Lifetime and Rebinding
---------------------------

Views are shallow. A previously created view can become invalid or out of date
when the host object changes storage, layout, membership, node counts, or pair
ordering. Recreate or rebind views after changes such as:

* field assignment, resizing, deletion, or deserialization;
* node-list internal or ghost count changes;
* ghost creation, culling, or deletion;
* domain redistribution;
* field-list membership changes;
* pair-related storage reallocation after connectivity changes;
* destruction of scratch host objects.

Physics code commonly creates views immediately before the kernels that consume
them, and then reacquires views after operations that change the host object
state.

Relationship to the Other Design Docs
-------------------------------------

:doc:`value_view_conversion_case_studies` is the current family reference. It
answers which host object owns the data, what object is passed to kernels, how
that object is made current, and when old views should not be reused.

:doc:`raja_chai_execution_patterns` explains how RAJA launch code gathers host
objects, creates views, moves/touches data, launches kernels, and returns
written data to CPU consumers.

:doc:`connectivity_data_structures` explains how connectivity construction
produces the flattened ``NodePairList`` consumed through ``NodePairListView``.
