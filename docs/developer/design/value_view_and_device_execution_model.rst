Value/View and Device Execution Model
=====================================

Purpose
-------

This document describes Spheral's owning-object/view-object pattern and how it
supports host execution, accelerator execution, restart/boundary machinery, and
per-step physics kernels. The pattern is most visible in ``Field`` and
``FieldList``, but the same design appears in kernels, interpolation tables,
node-pair lists, pairwise fields, and artificial-viscosity objects.

The central idea is:

* durable C++ objects own storage, identity, registration, restart behavior,
  resizing, and host-side polymorphic interfaces;
* view objects are lightweight, copyable, and callable from host/device code;
* kernels receive views, not the full owning object graph.

This split is important because normal Spheral simulations are assembled by
long-lived host objects, while inner loops increasingly need data structures
that can be copied into RAJA kernels and moved through CHAI execution spaces.
For a practical checklist for applying this pattern during GPU ports, see
:doc:`../dev/gpu_porting_patterns`.

Source Map
----------

Core field/value-view files:

* ``src/Field/Field.hh`` and ``src/Field/FieldInline.hh``
* ``src/Field/FieldBase.hh``
* ``src/Field/FieldView.hh`` and ``src/Field/FieldViewInline.hh``
* ``src/Field/FieldList.hh`` and ``src/Field/FieldListInline.hh``
* ``src/Field/FieldListView.hh`` and ``src/Field/FieldListViewInline.hh``
* ``src/Utilities/GPUUtils.hh``
* ``src/config.hh.in``

Related view-pattern files:

* ``src/Kernel/TableKernel.hh`` and ``src/Kernel/TableKernelView.hh``
* ``src/Utilities/*Interpolator*.hh`` and ``*InterpolatorView.hh``
* ``src/ArtificialViscosity/ArtificialViscosity.hh``
* ``src/ArtificialViscosity/ArtificialViscosityView.hh``
* ``src/Neighbor/NodePairList.hh`` and ``NodePairListView.hh``
* ``src/Neighbor/PairwiseField.hh`` and ``PairwiseFieldView.hh``

Design Summary
--------------

The most important local pattern is:

::

   host owner
     owns storage and durable identity
     participates in NodeList/DataBase/State/Boundary/Restart APIs
     rebuilds its view span whenever storage or layout changes

   view
     contains spans, managed arrays, primitive metadata, and small callable API
     is cheap to copy
     is ``SPHERAL_HOST_DEVICE`` where needed
     is captured by value in RAJA kernels

For ``Field`` this is implemented with multiple inheritance:

::

   Field<Dimension, DataType>
     derives from FieldBase<Dimension>
       host-side type-erased interface
     derives from FieldView<Dimension, DataType>
       typed element-access/device-facing interface
     owns std::vector<DataType, DataAllocator<DataType>>

This means a ``Field`` is both an owning field and a valid field view. Calling
``Field::view()`` returns a shallow ``FieldView`` copy of the view subobject. It
does not copy the underlying field values.

Why This Pattern Exists
-----------------------

Spheral has several requirements that pull in different directions.

Host-side infrastructure needs rich objects:

* fields must know their ``NodeList`` and name;
* boundaries need type-erased access to many field types;
* restarts and file I/O need packing, serialization, and cloning hooks;
* node lists must resize, reorder, and delete elements across all registered
  fields;
* Python bindings and high-level algorithms need familiar owning containers.

Device-side inner loops need simple objects:

* kernels cannot traverse the full host object graph;
* virtual/type-erased interfaces are expensive or unsafe in device kernels;
* data movement has to be explicit enough for CHAI and accelerator backends;
* kernel captures need to be small, copyable, and stable for the kernel launch.

The Value/View pattern keeps both sets of requirements in the codebase without
forcing every subsystem to choose one representation.

FieldBase: Host-Side Type Erasure
---------------------------------

``FieldBase<Dimension>`` is the common host-side interface for all fields,
regardless of the stored ``DataType``. It is used when code must handle fields
generically, such as:

* boundary-condition dispatch;
* restart registration and serialization;
* packing/unpacking values for communication;
* resizing all fields on a ``NodeList``;
* cloning fields through ``FieldBase`` pointers;
* state registration where concrete field type is recovered later.

The important design boundary is that ``FieldBase`` is not the inner-loop
interface. It is intentionally rich and host-oriented. Device kernels should
work with typed views.

Field: Owning Per-Node Storage
------------------------------

``Field<Dimension, DataType>`` owns one value of ``DataType`` per node on one
``NodeList``. Its private storage is:

::

   std::vector<DataType, DataAllocator<DataType>> mDataArray

``DataAllocator`` is selected at build time:

* ``uvm_allocator::UVMAllocator`` when ``USE_UVM`` is enabled;
* ``std::allocator`` otherwise.

The field also inherits the view data members from ``FieldView``:

::

   mDataSpan
   mNumInternalElements
   mNumGhostElements
   mNodesPerSmoothingScale

The owner keeps those view members synchronized through
``Field::assignDataSpan()``. That method is called after construction,
assignment, resizing, deletion, copying, deserialization, and other operations
that can change ``mDataArray`` or the node layout.

In unified-memory builds, ``mDataSpan`` is assigned as a span over
``mDataArray``. In non-UVM builds, ``GPUUtils::initMAView`` creates or refreshes
a ``chai::ManagedArray`` view over the vector's CPU storage. The method also
records a CPU touch and installs the CHAI callback when resource management is
enabled.

``assignDataSpan`` is the key invariant for ``Field``:

* ``mDataSpan.size()`` matches ``mDataArray.size()``;
* the span/managed-array points at the vector storage;
* internal and ghost counts match the attached ``NodeList``;
* the view metadata reflects ``nodesPerSmoothingScale``.

FieldView: Kernel-Facing Field Access
-------------------------------------

``FieldView<Dimension, DataType>`` is the typed access layer used by host/device
algorithms. It derives from ``chai::CHAICopyable`` and provides:

* ``operator()(index)``, ``operator[](index)``, and ``at(index)``;
* internal/ghost count accessors;
* simple element-wise arithmetic and min/max operations;
* local reductions;
* host iterators where appropriate;
* CHAI methods: ``data(space, do_move)``, ``move(space)``, ``touch(space)``,
  and ``shallowCopy``.

In unified-memory builds the container is a Spheral span type. Otherwise it is
``chai::ManagedArray<DataType>``. Most element accessors are marked
``SPHERAL_HOST_DEVICE`` so they can be used from ordinary CPU code or from RAJA
device kernels.

``Field::view()`` is intentionally shallow:

::

   return static_cast<FieldView<Dimension, DataType>>(*this);

The returned object contains a copy of the view metadata and the span/managed
array handle. It does not own the values. The owning ``Field`` and its
``NodeList`` must outlive the view and must not resize the field while a kernel
or saved view still depends on the old storage.

FieldList: Cross-NodeList Field Aggregation
-------------------------------------------

``FieldList<Dimension, DataType>`` is the global-field abstraction. It collects
same-typed fields across node lists, usually one field per node list.

It supports two storage modes:

``ReferenceFields``
  The field list points at fields owned elsewhere. This is the common mode for
  ``DataBase`` and ``State`` views over durable node-list or package fields.

``CopyFields``
  The field list owns copies in ``mFieldCache`` through ``shared_ptr<Field>``.
  This is used for temporary/scratch field lists and state snapshots that must
  own their values.

``FieldList::buildDependentArrays()`` rebuilds all derived indexing state after
fields are appended, deleted, copied, or referenced:

* fields are sorted in ``NodeListRegistrar`` order;
* ``mFieldPtrs`` holds typed field pointers;
* ``mFieldBasePtrs`` holds corresponding ``FieldBase`` pointers;
* ``mNodeListPtrs`` records node-list order;
* ``mNodeListIndexMap`` maps node-list pointer to field-list index;
* ``mFieldViews`` stores one ``FieldView`` per field.

This method is the field-list equivalent of ``Field::assignDataSpan``. It
preserves the invariant that the field list's node-list order, pointer arrays,
and view array all describe the same collection.

FieldListView: View of Views
----------------------------

``FieldListView<Dimension, DataType>`` is a device-copyable collection of
``FieldView`` objects. It provides:

* indexed access to field views with ``operator[](fieldIndex)``;
* direct value access with ``operator()(fieldIndex, nodeIndex)``;
* field-list-wide arithmetic and local reductions;
* ``numFields()``, ``numElements()``, ``numInternalElements()``, and
  ``numGhostElements()``;
* CHAI movement through ``move(space, recursive)`` and ``touch(space,
  recursive)``.

The ``recursive`` argument is significant. Moving the array of field views is
not enough in non-UVM builds: the field views themselves contain managed arrays
for the actual field values. With ``recursive=true``, ``FieldListView::move``
moves both layers.

This is why RAJA kernels can use simple syntax such as:

::

   auto mass_v = state.fields(HydroFieldNames::mass, 0.0);
   auto mass = mass_v.view();
   ...
   const auto mi = mass(nodeListi, i);

The kernel sees a field-list view and indexes through it; the owner hierarchy
remains outside the kernel.

Conversion Anatomy: Owner to View
---------------------------------

The current code should be read as an owner/view conversion pattern, not just
as a set of paired class names. The "before" side below is the conceptual
owner-only responsibility set: everything a rich host object would need to do.
The "after" side is the current split in Spheral: host ownership remains in
the original class, while the device-capable subset is exposed through a view.

``Field`` to ``FieldView``
~~~~~~~~~~~~~~~~~~~~~~~~~~

``Field`` is the primary example because the current owner derives from both
``FieldBase`` and ``FieldView``.

.. list-table::
   :header-rows: 1
   :widths: 32 34 34

   * - Concern
     - Owner-side role
     - View-side role
   * - Storage lifetime
     - ``Field`` owns ``mDataArray`` and is resized by ``NodeList`` operations.
     - ``FieldView`` holds ``mDataSpan`` as a span or ``chai::ManagedArray``.
   * - Type-erased host API
     - ``FieldBase`` handles boundaries, restarts, packing, cloning, and
       generic field registration.
     - Not present. Kernels use typed access and avoid ``FieldBase``.
   * - Element access
     - Host code can use ``Field`` directly because it inherits the view API.
     - Kernels use ``operator()``/``operator[]`` on the shallow view copy.
   * - Layout metadata
     - ``Field`` gets counts and smoothing-scale metadata from its
       ``NodeList``.
     - ``FieldView`` stores internal count, ghost count, and
       ``nodesPerSmoothingScale`` as primitive metadata.
   * - Data movement
     - ``Field::assignDataSpan`` rebuilds the view binding after storage or
       layout changes.
     - ``FieldView::move`` and ``touch`` operate on the span/managed array.

The important mechanism is ``Field::assignDataSpan``. It is the conversion
point from owned host storage to a valid view:

::

   Field mDataArray changes
     |
     v
   assignDataSpan()
     updates mDataSpan
     records CPU touch / CHAI callback
     refreshes internal/ghost metadata
     |
     v
   Field::view()
     returns shallow FieldView copy

After conversion, the kernel-facing object is not allowed to resize, register,
serialize, or otherwise manage field lifetime. It can only access and move the
data it views.

``FieldList`` to ``FieldListView``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``FieldList`` shows the same conversion one level higher. The owner is a
collection of fields and node-list indexing metadata; the view is a collection
of field views.

.. list-table::
   :header-rows: 1
   :widths: 32 34 34

   * - Concern
     - Owner-side role
     - View-side role
   * - Field membership
     - ``FieldList`` stores ``mFieldPtrs`` and optionally owns copies in
       ``mFieldCache``.
     - ``FieldListView`` does not own fields or decide membership.
   * - Node-list order
     - ``buildDependentArrays`` sorts fields in registrar order and rebuilds
       ``mNodeListIndexMap``.
     - The view assumes the field-view array is already in the correct order.
   * - Type-erased compatibility
     - ``mFieldBasePtrs`` lets host code interact with registered fields
       generically.
     - Not present in the view.
   * - Device indexing
     - Host code prepares ``mFieldViews`` from each field's ``view()``.
     - ``operator()(fieldIndex, nodeIndex)`` indexes directly into the right
       field view.
   * - Recursive movement
     - The owner decides when the field-list membership changes and rebuilds
       ``mFieldViews``.
     - ``move(space, recursive=true)`` moves both the array of views and the
       data inside each view.

The conversion path is:

::

   FieldList membership changes
     |
     v
   buildDependentArrays()
     sort Field pointers
     rebuild NodeList index map
     assign mFieldViews[i] = field->view()
     |
     v
   FieldList::view()
     returns FieldListView over mFieldViews

This is why a RAJA kernel can treat a multi-node-list quantity as a simple
``mass(nodeListi, i)`` access even though the host owner still preserves
node-list ordering, reference-vs-copy storage mode, and type-erased field
metadata.

``NodePairList`` to ``NodePairListView``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``NodePairList`` is a smaller and cleaner example of the same pattern. It is
useful for understanding the conversion without the extra ``NodeList`` and
``FieldBase`` machinery.

.. list-table::
   :header-rows: 1
   :widths: 32 34 34

   * - Concern
     - Owner-side role
     - View-side role
   * - Pair storage
     - ``NodePairList`` owns ``std::vector<NodePairIdxType>``.
     - ``NodePairListView`` holds a span or ``chai::ManagedArray`` over that
       vector.
   * - Pair lookup by value
     - The owner has a lazy ``unordered_map<NodePairIdxType, size_t>`` for
       host-side ``index(pair)`` queries.
     - The view exposes only indexed traversal; kernels use ``pairs[kk]``.
   * - Connectivity lifetime
     - ``ConnectivityMap`` replaces the owned ``NodePairList`` when
       connectivity is rebuilt or patched.
     - Any old view is stale after that replacement.
   * - Device traversal
     - The owner prepares the view through ``initView``.
     - RAJA kernels capture the view and loop from ``0`` to ``size()``.

The conversion path is:

::

   ConnectivityMap builds vector<NodePairIdxType>
     |
     v
   NodePairList::initView()
     wraps vector storage for CHAI/span access
     |
     v
   NodePairList::view()
     returns NodePairListView
     |
     v
   RAJA pair loop reads pairs[kk]

What moved into the view is intentionally minimal: pair-array access, size,
data pointer, and movement/touch operations. Host-only conveniences such as the
pair-to-index lookup map remain outside the kernel-facing type.

``PairwiseField`` to ``PairwiseFieldView``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``PairwiseField`` applies the same conversion to values stored per pair. The
extra design point is that the owner is tied to a particular ``NodePairList``.

.. list-table::
   :header-rows: 1
   :widths: 32 34 34

   * - Concern
     - Owner-side role
     - View-side role
   * - Pair association
     - ``PairwiseField`` stores a weak pointer to the active ``NodePairList``.
     - The view does not know how pair ids map to indices.
   * - Value storage
     - The owner sizes ``mArray`` to ``numElements * pairs.size()``.
     - The view exposes strided indexed access over the value array.
   * - Lookup by pair
     - Host code can call ``operator()(NodePairIdxType)`` through the owner's
       pair-list lookup.
     - Kernels use integer pair index access, typically aligned with
       ``NodePairListView`` traversal.
   * - Validity
     - Rebuilding connectivity can orphan or invalidate the pair association.
     - The view is valid only for the current owner storage and current pair
       ordering.

This class is the clearest warning that views are not durable design objects.
They are launch-local access objects over storage whose owner defines lifetime.

Device Memory Semantics
-----------------------

Spheral abstracts host/device storage with a small number of conventions.

``SPHERAL_HOST_DEVICE``
  Expands to RAJA's host/device annotation. Use it on small methods that are
  intended to run inside RAJA kernels.

``SPHERAL_UNIFIED_MEMORY`` / ``USE_UVM``
  Select span/UVM behavior. In this mode many explicit CHAI movement operations
  become no-ops because the pointer is already unified.

``chai::ManagedArray``
  Used in non-UVM builds to wrap host-owned storage and make it movable between
  execution spaces.

``move(space)``
  Requests that CHAI move data to ``chai::CPU`` or ``chai::GPU``.

``touch(space)``
  Registers that data was accessed in an execution space without necessarily
  moving immediately.

``data(space, do_move)``
  Returns a pointer in the requested execution space, optionally moving first.

The owner/view pattern does not eliminate the need to reason about data
movement. It gives each container one local place to implement movement rules.
For fields, that is ``FieldView`` and ``Field::assignDataSpan``. For field
lists, that is ``FieldListView::move``. For node-pair lists and pairwise
fields, the equivalent methods live in ``NodePairListView`` and
``PairwiseFieldView``.

Lifecycle of a Typical Field in a Step
--------------------------------------

A common field lifecycle looks like this:

1. A ``NodeList`` constructs durable fields such as mass, position, velocity,
   and ``H``.
2. Physics packages construct additional durable or scratch fields.
3. A ``State`` or ``StateDerivatives`` object registers fields and field lists
   by key.
4. A physics package queries field lists from ``State``/``StateDerivatives``.
5. The package calls ``view()`` on each field list before the inner loop.
6. Optional CHAI movement places the views and their data in the execution
   space required by the kernel.
7. A RAJA kernel captures views by value and reads/writes values through
   ``operator()``.
8. Views that will be consumed by subsequent host code are moved or touched
   back to ``chai::CPU``.
9. The durable owning fields now contain the updated values.

The transient ``State`` object is not the owner in this sequence. It is a
registry over durable fields. Mutating a field through a view mutates the
underlying field storage.

View Invalidation Rules
-----------------------

Views are shallow. Developers should treat them as launch-local objects, not
durable handles.

A previously created view can become stale when:

* the owning ``Field`` is assigned a new vector;
* a node list changes internal or ghost node counts;
* ghost nodes are created, culled, or deleted;
* domains are redistributed;
* ``FieldList`` membership changes;
* ``ConnectivityMap`` is rebuilt and pair-related storage is reallocated;
* a scratch field list is destroyed.

The practical rule is to reacquire views after any operation that changes
storage, membership, node counts, or connectivity. Physics kernels should
create views immediately before the kernel launches that use them.

Other View-Pattern Objects
--------------------------

``TableKernel`` and ``TableKernelView``
  ``TableKernel`` owns tabulated interpolation data. ``TableKernelView`` owns
  view versions of interpolators and exposes ``SPHERAL_HOST_DEVICE`` kernel
  evaluation methods such as ``kernelValue`` and ``kernelAndGradValue``.

Interpolator views
  Quadratic and cubic-Hermite interpolators follow the same pattern: the owner
  keeps host vectors/managed arrays, while the view provides device-callable
  lookup methods and CHAI movement.

Artificial-viscosity views
  Viscosity owners expose CHAI-managed view objects for device use. This is a
  special case because ``QPiij`` is virtual on the view. See
  ``raja_chai_execution_patterns.rst`` for details and caveats.

``NodePairList`` and ``NodePairListView``
  ``NodePairList`` owns the flattened pair vector and a lazy host-side
  pair-to-index map. ``NodePairListView`` exposes indexed pair access and CHAI
  movement for kernels.

``PairwiseField`` and ``PairwiseFieldView``
  ``PairwiseField`` stores values indexed by the active ``NodePairList``. It is
  explicitly ephemeral because connectivity can change each step. Its view
  exposes pair-index storage to kernels.

Choosing the Right Interface
----------------------------

Use ``FieldBase`` when code must be type-erased and host-oriented:

* boundary dispatch;
* restart and file I/O;
* communication packing;
* generic node-list field registration.

Use ``Field`` when code owns or mutates the host field:

* construction;
* resizing;
* node deletion;
* assignment from host vectors;
* Python-facing field operations.

Use ``FieldView`` or ``FieldListView`` inside algorithms that need direct typed
element access, especially RAJA kernels.

Use ``FieldList`` for multi-node-list algorithms at the host level. Use
``FieldListView`` once the algorithm enters a device-capable inner loop.

Extension Guidance
------------------

When adding a new performance-sensitive container, follow the existing pattern.

* Put ownership, lifetime, restart behavior, and host-only indexing helpers in
  the owning object.
* Put only the data needed by kernels into the view object.
* Make view methods ``SPHERAL_HOST_DEVICE`` only when they are intended to run
  in kernels.
* Keep view copies shallow and cheap.
* Provide ``move`` and ``touch`` methods if the view owns CHAI-managed data.
* Rebuild view spans in one local helper after every storage-changing operation.
* Avoid embedding host pointers, STL containers, or host-only virtual dispatch
  in kernel-facing views.
* Document whether views remain valid after resizing, copying, connectivity
  rebuilds, or owner destruction.

Common Failure Modes
--------------------

Stale view after resize
  A field or field list was resized or rebuilt after a view was captured. Create
  the view again after the storage change.

Moving only the outer view array
  A ``FieldListView`` contains an array of ``FieldView`` objects, and each
  ``FieldView`` contains data. Use recursive movement when both layers must be
  moved.

Using host-only APIs in device code
  ``FieldBase``, STL containers, and many owner methods are not kernel-safe.
  Use the view API in RAJA lambdas.

Assuming ``State`` owns values
  ``State`` commonly stores references to durable fields. Mutations through
  ``State`` or its views usually affect the original node-list/package fields.

Retaining pairwise data across connectivity rebuilds
  ``PairwiseField`` is tied to a particular ``NodePairList``. Rebuilding or
  patching connectivity can invalidate previous pairwise fields.

Forgetting CPU movement before host consumers
  In non-UVM builds, data written on the GPU must be moved or touched for CPU
  access before later host code reads it.

Relationship to the Other Design Docs
-------------------------------------

``integrator_and_state_update_model.rst`` explains how ``State`` and
``StateDerivatives`` register field lists during a time step. This document
explains the field/view machinery those registries expose to physics packages.

``raja_chai_execution_patterns.rst`` explains how the views are used in RAJA
kernels and how CHAI movement is applied around those kernels.

``connectivity_data_structures.rst`` explains how the same owner/view pattern
is used for flattened node-pair storage and pairwise fields.
