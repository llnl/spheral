Current Device-Facing Object Families
=====================================

Purpose
-------

This document is the concrete family reference for Spheral's current
device-facing object model. It complements
:doc:`value_view_and_device_execution_model`, which defines the value/view and
managed view pointer shapes, and :doc:`raja_chai_execution_patterns`, which
describes how those objects are used around RAJA launches.

Each section identifies:

* the durable owner-side role;
* the device-facing member that kernels capture, index, or call through;
* the owner method or refresh point that makes the device-facing object current;
* lifetime or invalidation constraints that matter to kernel setup.

Source Map
----------

Field and field-list families:

* ``src/Field/FieldBase.hh``
* ``src/Field/Field.hh`` and ``src/Field/FieldInline.hh``
* ``src/Field/FieldView.hh`` and ``src/Field/FieldViewInline.hh``
* ``src/Field/FieldList.hh`` and ``src/Field/FieldListInline.hh``
* ``src/Field/FieldListView.hh`` and ``src/Field/FieldListViewInline.hh``

Pair and pairwise-value families:

* ``src/Neighbor/NodePairList.hh`` and ``src/Neighbor/NodePairList.cc``
* ``src/Neighbor/NodePairListView.hh``
* ``src/Neighbor/PairwiseField.hh`` and ``src/Neighbor/PairwiseFieldInline.hh``
* ``src/Neighbor/PairwiseFieldView.hh``
* ``src/Neighbor/PairwiseFieldElementAccessor.hh``

Kernel and interpolation families:

* ``src/Kernel/TableKernel.hh`` and ``src/Kernel/TableKernel.cc``
* ``src/Kernel/TableKernelInline.hh``
* ``src/Kernel/TableKernelView.hh`` and ``src/Kernel/TableKernelViewInline.hh``
* ``src/Utilities/QuadraticInterpolatorView.hh``
* ``src/Utilities/CubicHermiteInterpolatorView.hh``

Managed-dispatch family:

* ``src/ArtificialViscosity/ArtificialViscosity.hh``
* ``src/ArtificialViscosity/ArtificialViscosityView.hh``
* concrete artificial-viscosity owners and view implementations

``FieldBase``, ``Field``, and ``FieldView``
-------------------------------------------

``Field<Dimension, DataType>`` is the primary value/view example. It owns one
value of ``DataType`` per node on one ``NodeList``. It inherits the typed view
interface from ``FieldView<Dimension, DataType>`` while ``FieldBase<Dimension>``
provides the host-side type-erased interface.

Owner-side role:

* ``Field`` owns ``std::vector<DataType, DataAllocator<DataType>> mDataArray``;
* ``FieldBase`` supports boundary dispatch, restart registration, serialization,
  packing, cloning, and generic state registration;
* ``Field`` responds to node-list resizing, deletion, reordering, copying, and
  deserialization.

Device-facing role:

* ``FieldView`` stores ``mDataSpan`` plus primitive metadata such as internal
  and ghost counts;
* ``operator()``, ``operator[]``, and ``at`` provide typed indexed access;
* ``move``, ``touch``, and ``data`` manage CHAI-backed storage where needed;
* element accessors used by kernels are marked ``SPHERAL_HOST_DEVICE``.

Refresh point:

::

   Field storage or node layout changes
     |
     v
   Field::assignDataSpan()
     updates span/ManagedArray binding
     refreshes internal and ghost metadata
     records CPU touch and callback where needed
     |
     v
   Field::view()
     returns shallow FieldView copy

``FieldBase`` is intentionally not the inner-loop interface. Kernels use typed
``FieldView`` objects or field-list views built from them.

``FieldList`` and ``FieldListView``
-----------------------------------

``FieldList<Dimension, DataType>`` aggregates same-typed fields across node
lists. It may reference fields owned elsewhere or own copied fields in
``mFieldCache``.

Owner-side role:

* ``FieldList`` manages field membership and reference-vs-copy storage mode;
* ``buildDependentArrays`` sorts fields in ``NodeListRegistrar`` order;
* the owner maintains typed field pointers, ``FieldBase`` pointers, node-list
  pointers, and node-list index maps.

Device-facing role:

* ``FieldListView`` stores an array of ``FieldView`` objects;
* kernels use ``operator()(fieldIndex, nodeIndex)`` for direct value access;
* arithmetic, local reductions, and size/count queries operate over the view;
* ``move(space, recursive=true)`` can move both the outer array and each nested
  field view's data.

Refresh point:

::

   FieldList membership or ownership changes
     |
     v
   buildDependentArrays()
     sort field pointers
     rebuild node-list index map
     assign mFieldViews[i] = field->view()
     |
     v
   FieldList::view()
     returns FieldListView over current field views

This lets a RAJA kernel use syntax such as ``mass(nodeListi, i)`` while the
host owner preserves node-list ordering, type-erased compatibility, and field
membership semantics.

``NodePairList`` and ``NodePairListView``
------------------------------------------

``NodePairList`` is the owner/view family for an already-built flat pair
schedule. ``NodePairListView`` is the device-facing member captured by RAJA pair
kernels.

Owner-side role:

* ``NodePairList`` owns ``std::vector<NodePairIdxType>``;
* the owner has host-side lookup support for mapping a pair value to an index;
* connectivity construction replaces or refreshes the owner when pair topology
  changes.

Device-facing role:

* ``NodePairListView`` holds a span or ``chai::ManagedArray`` over the pair
  vector;
* kernels index it by integer pair position with ``pairs[kk]``;
* the view exposes pair-array size, data pointer, movement, and touch.

Refresh path:

::

   pair schedule is built or replaced
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

``ConnectivityMap`` is important context for this family because it builds and
provides the current ``NodePairList``. The device-facing family member remains
``NodePairListView``.

``PairwiseField`` and ``PairwiseFieldView``
--------------------------------------------

``PairwiseField<Dimension, Value, numElements>`` stores values indexed by the
active pair-list position rather than by node. It is intentionally ephemeral
because connectivity can change step to step.

Owner-side role:

* ``PairwiseField`` stores a weak pointer to the active ``NodePairList``;
* the owner sizes ``mArray`` to ``numElements * pairs.size()``;
* host code can access values by ``NodePairIdxType`` through the pair-list
  lookup.

Device-facing role:

* ``PairwiseFieldView`` exposes strided indexed access over the value array;
* kernels use integer pair-index access aligned with ``NodePairListView``
  traversal;
* movement and touch operate on the pairwise-value storage.

The view does not know how pair ids map to indices. Rebuilding or patching
connectivity can orphan the owner-side pair association and invalidate existing
views.

``TableKernel`` and ``TableKernelView``
----------------------------------------

``TableKernel<Dimension>`` owns tabulated interpolation data for SPH kernel
evaluation. ``TableKernelView<Dimension>`` is the device-facing member used by
RAJA hydro kernels.

Owner-side role:

* ``TableKernel`` owns interpolation tables such as ``mInterpVal``,
  ``mGradInterpVal``, and ``mGrad2InterpVal``;
* the owner initializes lookup tables from the selected kernel;
* construction and assignment refresh the embedded view members.

Device-facing role:

* ``TableKernelView`` contains interpolator views and primitive lookup metadata;
* kernels call host/device methods such as ``kernelValue``, ``gradValue``,
  ``grad2Value``, and ``kernelAndGradValue``;
* ``TableKernelView::move`` recursively moves nested interpolator views.

The nested movement is the same general issue as ``FieldListView`` movement:
moving only the outer object is not sufficient when it contains managed data in
its child views.

``ArtificialViscosity`` and ``ArtificialViscosityView``
--------------------------------------------------------

Artificial viscosity follows the owner/device-facing split but uses managed
runtime dispatch rather than a plain value view. The device-facing member is a
``chai::managed_ptr`` to an ``ArtificialViscosityView<Dimension, QPiType>`` base
object.

Owner-side role:

* ``ArtificialViscosity`` descendants own host parameters, package-facing APIs,
  restart identity, and configuration setters;
* concrete owners choose which concrete view class represents them;
* setters update owner-side values and refresh existing managed views.

Device-facing role:

* ``ArtificialViscosityView`` defines the virtual ``SPHERAL_HOST_DEVICE``
  ``QPiij`` interface;
* concrete view descendants contain only the data and methods needed by kernels;
* RAJA pair kernels capture a ``chai::managed_ptr`` to the base view and call
  ``Q->QPiij(...)``.

Dispatch path:

::

   host ArtificialViscosity owner
     reports QPiTypeIndex()
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

Host code first selects the scalar or tensor templated path. Dynamic
polymorphism remains only for the concrete artificial-viscosity behavior inside
that path. This keeps most type selection on the host while preserving the
runtime behavior that the pair kernel needs.

Shared Observations
-------------------

Across the current families:

* owners define lifetime and invariants;
* views expose the smallest useful kernel-facing API;
* host setup chooses concrete types and gathers state before launch;
* kernels capture views or managed view pointers, not package owners;
* view freshness depends on storage, membership, and pair ordering remaining
  unchanged after the view is created;
* managed polymorphism is reserved for behavior that cannot be selected entirely
  on the host side.

The important distinction is the responsibility split, not the class-name
pairing itself.
