Spheral GPU Port Design Patterns
################################

Status and Audience
===================

This page is a developer design guide for porting Spheral components toward
GPU execution. It is intended for maintainers who need to modify C++ classes,
RAJA kernels, CHAI-managed data, PYB11 bindings, and Python construction
helpers without losing the existing host behavior.

The examples are drawn from existing Spheral code and from the RAJA
``evaluateDerivatives`` work merged in commit
``548becc88a4ed60e58873326456a9559fa1db62b``. The guidance here documents
patterns used by those examples. It does not settle project-specific design
questions for future ``Neighbor``, ``TreeNeighbor``, ``NodeList``, or
``ConnectivityMap`` work.

Evidence and Limits
===================

This document distinguishes between three kinds of statements:

Observed examples
  Descriptions of how existing code such as ``FieldView``,
  ``NodePairListView``, ``ArtificialViscosityView``, ``SPH_RAJA``, and
  ``SolidSPH_RAJA`` currently work.

Design guidance
  Recommendations for future ports based on those examples and general
  maintainability concerns.

Open design choices
  Decisions that require additional analysis of the target subsystem. This page
  should not present those choices as already resolved.

When the document discusses future ``ConnectivityMap`` or ``Neighbor`` work,
read it as guidance to investigate, not as a claim that the current examples
already prove the final design.

Design Goals
============

For maintainability, GPU-porting changes should aim to:

* Preserve the existing host-facing APIs unless there is a deliberate design
  decision to change them.
* Keep owning objects responsible for lifetime, invariants, restart, Python
  exposure, and compatibility.
* Expose compact device-facing views that can be copied into RAJA kernels.
* Make CHAI data movement visible and reviewable.
* Keep host-only discovery, allocation, and polymorphic selection out of device
  loop bodies.
* Add tests at the abstraction boundary being changed, not only at the final
  application level.

Non-goals
=========

This page is not a tutorial for RAJA, CHAI, PYB11, or CMake. It also does not
require every future port to use the exact implementation shape used by the
examples below. In particular, inheritance from a view class is an existing
Spheral pattern, not a universal requirement.

Terminology
===========

Owning object
  The full host-facing C++ object. It owns storage, invariants, compatibility
  behavior, restart state, Python-facing semantics, and any host-only
  containers or lookup structures.

View
  A lightweight object containing the subset of data and operations needed by
  device code. A view should be cheap to copy and safe to capture by value in a
  RAJA lambda.

Managed view
  A view allocated with ``chai::managed_ptr`` so that device code can call
  methods through a pointer. This is used for device-side virtual dispatch.

Device-callable method
  A method marked with ``SPHERAL_HOST_DEVICE`` and limited to operations that
  are valid in host and device execution spaces.

Pattern Summary
===============

The GPU-ready examples documented on this page use this general shape:

::

   Host-facing object
     owns data, invariants, restart, Python semantics
     |
     | view()
     v
   Device-facing view
     compact data access, SPHERAL_HOST_DEVICE methods
     |
     | captured by value
     v
   RAJA kernel
     simple range, explicit inputs, atomics where required

The documented managed-dispatch example uses this shape:

::

   Concrete host object
     owns chai::managed_ptr<ConcreteView>
     |
     | dynamic_pointer_cast
     v
   chai::managed_ptr<BaseView>
     |
     | captured in RAJA kernel
     v
   virtual device-callable method

Worked Example: ``Field`` and ``FieldView``
===========================================

Problem Addressed
-----------------

``Field`` is one of the central data containers in Spheral. Existing host code
expects a rich object that knows its name, its ``NodeList``, how to resize, how
to pack values for communication, how to serialize, and how to participate in
restart and output.

The documented kernels need much less. In those kernels, the common operation
is indexed access to field values, together with a small amount of metadata
such as the number of internal and ghost elements.

Design
------

``Field`` remains the owning object. It owns the host ``std::vector`` storage
and the full host API. It derives from ``FieldBase`` and participates in
Spheral's field, restart, communication, and output infrastructure.

``FieldView`` is the device-facing access object. It contains:

* A CHAI ``ManagedArray`` or span over the field values.
* The internal and ghost element counts.
* The nodes-per-smoothing-scale value.
* ``SPHERAL_HOST_DEVICE`` element access and simple arithmetic operations.
* ``move`` and ``touch`` methods for CHAI-managed data.

The current implementation has ``Field`` inherit from ``FieldView``. That keeps
the view operations directly available on ``Field`` and lets ``Field::view()``
return a lightweight copy of the ``FieldView`` subobject:

.. code-block:: c++

   auto positionFields = state.fields(HydroFieldNames::position,
                                      Vector::zero());
   auto position = positionFields.view();

   RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, n),
     [=] SPHERAL_HOST_DEVICE (size_t i) {
       const auto& ri = position(nodeListi, i);
       // Device-callable work using ri.
     });

Important Consequences
----------------------

The kernel captures ``position``, not the owning ``State``, ``FieldList``, or
``Field``. This is the important design boundary. Host code performs lookup and
policy decisions before the loop. Device code receives a compact access object.

``Field`` also has to keep the host vector and CHAI view synchronized. When the
host storage changes, the managed span is refreshed so that the view still
points at the correct data. That synchronization belongs in the owning object,
not in arbitrary kernels.

Review Checklist
----------------

When applying this pattern to another class, check:

* Does the owning object retain all host-only behavior?
* Does the view contain only data needed by device code?
* Are all view methods used in kernels marked ``SPHERAL_HOST_DEVICE``?
* Is CHAI movement explicit and close to the algorithm that needs it?
* Is there a clear rule for invalidating or refreshing views when host storage
  changes?

Worked Example: ``FieldList`` and ``FieldListView``
===================================================

Problem Addressed
-----------------

Many Spheral algorithms operate on one field per ``NodeList``. Host code needs
the full ``FieldList`` abstraction, but kernels need fast indexed access by
``(fieldIndex, nodeIndex)``.

Design
------

``FieldListView`` packages a collection of ``FieldView`` objects into a
device-usable container. Device code can write:

.. code-block:: c++

   auto massFields = state.fields(HydroFieldNames::mass, 0.0);
   auto mass = massFields.view();

   const auto mi = mass(nodeListi, i);

The view hides the host ``FieldList`` machinery while preserving the indexing
model used throughout the C++ code.

Important Consequences
----------------------

The RAJA hydro code first extracts all state and derivative ``FieldList``
objects, converts them to views, and only then enters the pair or node loops.
That makes kernel dependencies obvious: the loop body uses captured views rather
than reaching back into ``State`` or ``StateDerivatives``.

Worked Example: ``NodePairList`` and ``NodePairListView``
=========================================================

Problem Addressed
-----------------

Pairwise hydro loops need a compact schedule of node-node interactions. Host
code also needs richer operations: sorting, insertion, pair-to-index lookup, and
domain-decomposition ordering.

Design
------

``NodePairList`` is the host-side owner. It stores a
``std::vector<NodePairIdxType>`` and maintains host-side conveniences such as
iterators and a lazy pair-to-index lookup map.

``NodePairListView`` is intentionally much smaller. It stores only a CHAI
``ManagedArray`` or span of ``NodePairIdxType`` values and exposes:

* ``operator[]`` for indexed access.
* ``size()`` and ``data()``.
* ``move()`` and ``touch()``.

This is enough for the RAJA pair loops shown in the current SPH examples:

.. code-block:: c++

   const auto& pairsValue = connectivityMap.nodePairList();
   const auto  pairs = pairsValue.view();
   const auto  npairs = pairs.size();

   RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, npairs),
     [=] SPHERAL_HOST_DEVICE (size_t kk) {
       const auto i = pairs[kk].i_node;
       const auto j = pairs[kk].j_node;
       const auto nodeListi = pairs[kk].i_list;
       const auto nodeListj = pairs[kk].j_list;

       // Compute this pair interaction.
     });

Important Consequences
----------------------

The view does not expose the pair-to-index lookup map. That omission is a
design feature, not a missing API, for loops that operate on an already-built
pair schedule. In that case, device code can consume the flattened pair list
without carrying the host lookup structure into the kernel.

This is a useful precedent, but not by itself a complete answer for future
connectivity work. ``NodePairListView`` shows that an already-built pair
schedule can be consumed on device through a flattened representation. It does
not prove that the algorithms that build or patch connectivity can avoid richer
device-side data structures. Future connectivity work should determine which
parts can use flattened inputs and which parts require additional staging or
new device-ready structures.

Worked Example: ``ArtificialViscosity`` and Managed Views
=========================================================

Problem Addressed
-----------------

``ArtificialViscosity`` is not just a data container. It is a ``Physics``
package with state registration, derivative registration, boundary behavior,
restart state, and several concrete implementations. The hydro algorithm needs
to call the pairwise ``QPiij`` operation on device, but the concrete viscosity
type is selected through the host object model.

Design
------

The host-facing ``ArtificialViscosity`` object owns the full physics-package
behavior. It also exposes ``getScalarView()`` and ``getTensorView()`` so callers
can obtain the appropriate device-facing polymorphic view.

``ArtificialViscosityView`` is the base device-facing interface. It carries the
common viscosity parameters and defines a virtual ``SPHERAL_HOST_DEVICE``
``QPiij`` method. Concrete views, such as
``MonaghanGingoldViscosityView`` and
``LimitedMonaghanGingoldViscosityView``, implement that method.

Concrete host viscosity classes allocate the concrete view with
``chai::make_managed`` and return it through a base managed pointer:

.. code-block:: c++

   chai::managed_ptr<ArtViscView> getScalarView() override {
     initView();
     return chai::dynamic_pointer_cast<ArtViscView, ViewType>(m_viewPtr);
   }

   void initView() {
     if (!m_viewPtr) {
       m_viewPtr = chai::make_managed<ViewType>(mClinear,
                                                mCquadratic,
                                                mLinearInExpansion,
                                                mQuadraticInExpansion);
     }
   }

The RAJA hydro path chooses scalar or tensor dispatch before entering the main
pair loop:

.. code-block:: c++

   auto& Qhandle = this->artificialViscosity();
   if (Qhandle.QPiTypeIndex() == std::type_index(typeid(Scalar))) {
     chai::managed_ptr<ArtificialViscosityView<Dimension, Scalar>> Q =
       Qhandle.getScalarView();
     this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
   }

Inside the pair loop, the managed base pointer is captured and used through the
virtual interface:

.. code-block:: c++

   Q->QPiij(QPiij, QPiji, Qi, Qj,
            nodeListi, i, nodeListj, j,
            ri, Hi, etai, vi, rhoi, ci,
            rj, Hj, etaj, vj, rhoj, cj,
            fClQView, fCqQView, DvDxQView);

Important Consequences
----------------------

This pattern is best reserved for cases where the device path really needs
polymorphism. A managed virtual hierarchy has stricter lifetime rules than a
plain data view.

In the viscosity implementation, when host-side parameters change, the concrete
host class frees and reconstructs the managed view rather than modifying it in
place. This is intended to preserve the validity of the device object and its
vtable.

RAJA Algorithm Structure
========================

The RAJA ``SPH_RAJA`` and ``SolidSPH_RAJA`` implementations show how the
building blocks compose into a larger algorithm.

The RAJA SPH examples use these stages. Future ports can use this sequence as a
checklist, but should still verify which stages apply to the target algorithm:

1. Gather host-side state from ``State``, ``StateDerivatives``, kernels,
   viscosity objects, and connectivity.
2. Convert owning objects to views or managed pointers.
3. Move views to the GPU when needed.
4. Launch simple range-based pair loops over ``NodePairListView``.
5. Use atomics for shared per-node outputs.
6. Launch separate final per-node loops when reduction/finalization differs
   from pair interaction.
7. Move modified views back to CPU before host code reads them.

A representative loop looks like:

.. code-block:: c++

   RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, npairs),
     [=] SPHERAL_HOST_DEVICE (size_t kk) {
       const auto i = pairs[kk].i_node;
       const auto j = pairs[kk].j_node;
       const auto nodeListi = pairs[kk].i_list;
       const auto nodeListj = pairs[kk].j_list;

       // Read state through views.
       // Compute pair contribution.
       // Accumulate results through atomics where needed.
     });
   GPU_ERROR_CHECK

In this staged pattern, keep the loop body from allocating storage, searching
host containers, or selecting concrete implementations. Those steps are easier
to review when they happen in host setup before the kernel.

Race Management
===============

The documented pair loops update both members of a pair. Different pair
iterations can therefore update the same node concurrently. The RAJA SPH path
handles this with atomics where shared outputs are unavoidable.

Examples include:

* ``GPUUtils::AtomicAddOp::apply`` for scalar accumulation.
* ``GPUUtils::AtomicMaxOp::apply`` for maxima.
* Value-type methods such as ``atomicAdd`` and ``atomicSub`` on vector and
  tensor data.

Before porting a loop, identify each output and classify it:

Private per iteration
  No atomic is needed.

One writer per element
  No atomic is needed if the one-writer property is guaranteed by the loop
  structure.

Multiple writers per element
  Use atomics, restructure the algorithm, or split the operation into staged
  reductions.

Python and PYB11 Integration
============================

Spheral GPU ports often need Python-facing plumbing, not just C++ kernels.

The RAJA SPH work added C++ classes such as ``SPH_RAJA`` and
``SolidSPH_RAJA`` rather than replacing ``SPH`` and ``SolidSPH`` in place.
Those classes were then exposed through PYB11 and selected from the Python
``SPHHydros.py`` helper when the caller passes ``RAJA=True``.

This staged approach has several advantages:

* Existing Python scripts keep their default behavior.
* New behavior is opt-in and easy to compare against the original path.
* Tests can select the RAJA path explicitly.
* C++ and Python review can happen together, so the public construction API is
  not an afterthought.

When a future port adds a new C++ implementation class, check the full exposure
path:

* Is the class added to the package CMake instantiation lists?
* Is there a PYB11 wrapper for the new class?
* Is the wrapper included and instantiated for the enabled dimensions?
* Does the Python helper expose an explicit option for selecting the new path?
* Do functional tests exercise that Python option?

Testing Strategy
================

GPU-porting tests should start at the smallest boundary that can fail.

View tests
  Verify construction, copy behavior, indexing, movement, and touch semantics.

Managed pointer tests
  Verify ``chai::managed_ptr`` capture, dynamic cast, and virtual calls inside
  RAJA loops.

Loop tests
  Verify simple pair-list traversal and atomic updates under the enabled
  execution policies.

Functional tests
  Verify that Python-facing options select the intended implementation and that
  physics behavior remains acceptable.

Performance tests
  Isolate the ported kernel or method so timing changes can be interpreted.

The RAJA ``evaluateDerivatives`` work used all of these categories: low-level
view tests, managed pointer tests, loop tests, functional Noh tests, and an
isolated ``evaluateDerivatives`` performance driver.

Design Checklist for Future Ports
=================================

Before implementing a new GPU-ready class or algorithm, write down:

* The owning object and the behavior it must preserve.
* The minimal data needed by device code.
* The proposed view type and its lifetime.
* Which methods must be ``SPHERAL_HOST_DEVICE``.
* Whether polymorphism is required on device.
* How CHAI movement will be triggered and audited.
* Which writes can race and how races are handled.
* How the C++ implementation is exposed to Python, if applicable.
* Which unit, loop, functional, and performance tests will validate the change.

Reference Implementation Map
============================

The following files are useful starting points when applying these patterns:

* ``src/Field/Field.hh`` and ``src/Field/FieldView.hh``: owning field object
  and compact field view.
* ``src/Field/FieldList.hh`` and ``src/Field/FieldListView.hh``: host
  collection of fields and device-usable collection of views.
* ``src/Neighbor/NodePairList.hh`` and ``src/Neighbor/NodePairListView.hh``:
  host pair-list owner and device pair-list view.
* ``src/ArtificialViscosity/ArtificialViscosity.hh`` and
  ``src/ArtificialViscosity/ArtificialViscosityView.hh``: host physics package
  interface and managed device-view interface.
* ``src/SPH/SPH_RAJA.cc`` and ``src/SPH/SolidSPH_RAJA.cc``: full RAJA
  algorithm examples.
* ``src/PYB11/SPH/SPH_RAJA.py`` and ``src/SPH/SPHHydros.py``: Python exposure
  and opt-in construction path.
* ``tests/cpp/Field``, ``tests/cpp/Neighbor``, ``tests/cpp/Loops``, and
  ``tests/cpp/Utilities/managedptr_tests.cc``: tests for the building blocks
  used by GPU-ready code.
