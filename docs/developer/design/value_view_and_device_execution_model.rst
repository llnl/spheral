Value/View and Device Execution Model
=====================================

Purpose
-------

This document defines the object shapes Spheral currently uses at the boundary
between durable host-owned objects and device-capable RAJA kernels. It is the
contract page for the pattern; concrete class families are described in
:doc:`value_view_conversion_case_studies`, and launch mechanics are described in
:doc:`raja_chai_execution_patterns`.

In this context, a device-facing object is an object that kernel setup can copy
or pass into a RAJA kernel as the kernel-facing representation, or a managed view
object whose host/device methods are called from that kernel. The owning side of
the same family is included when it directly produces, refreshes, and
lifetime-manages that device-facing representation.

The recurring split is:

* durable C++ objects own storage, identity, registration, restart behavior,
  resizing, and host-side APIs;
* device-facing objects are lightweight, copyable, and limited to data and
  behavior that can participate in host/device kernels;
* kernels receive views or managed view pointers rather than the full owning
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

Representative managed-dispatch family:

* ``src/ArtificialViscosity/ArtificialViscosity.hh``
* ``src/ArtificialViscosity/ArtificialViscosityView.hh``
* concrete artificial-viscosity owners and view implementations

Device-Facing Object Shapes
---------------------------

Device kernels need a narrower representation than the owning object usually
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

   owner
     owns storage, identity, registration, restart, resizing, host API
     rebuilds the device-facing binding when storage or layout changes
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

The owner keeps responsibilities such as:

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

Managed View Pointer Dispatch Shape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Managed view pointer dispatch appears where device code needs runtime
polymorphism. The artificial-viscosity family is the current example.

::

   host owner
     selects concrete behavior
     owns durable parameters and restart state
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
specialized before launch. Use managed device polymorphism only when runtime
dispatch must remain inside the device kernel. In the current code, this is used
for ``ArtificialViscosityView::QPiij``.

The managed view is an object with behavior, not just a span over owner storage.
Its virtual table is part of correctness, so the host owner constructs and
refreshes the managed view object when its concrete parameters change.

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

The owner/view pattern does not eliminate the need to reason about data
movement. It gives each container one local place to implement movement
behavior, and it gives kernel setup code a clear set of captured objects to move
or touch before launch.

View Lifetime and Refresh Points
--------------------------------

Views are shallow. A previously created view can become stale when the owner
changes storage, layout, membership, node counts, or pair ordering. Common
refresh points include:

* field assignment, resizing, deletion, or deserialization;
* node-list internal or ghost count changes;
* ghost creation, culling, or deletion;
* domain redistribution;
* field-list membership changes;
* pair-related storage reallocation after connectivity changes;
* destruction of scratch owners.

Physics code commonly creates views immediately before the kernels that consume
them, and then reacquires views after operations that change the owner-side
state.

Relationship to the Other Design Docs
-------------------------------------

:doc:`value_view_conversion_case_studies` is the current family reference. It
identifies the owner side, device-facing side, refresh point, and kernel-facing
use for each family.

:doc:`raja_chai_execution_patterns` explains how RAJA launch code gathers
owners, creates views, moves/touches data, launches kernels, and returns written
data to CPU consumers.

:doc:`connectivity_data_structures` explains how connectivity construction
produces the flattened ``NodePairList`` consumed through ``NodePairListView``.
