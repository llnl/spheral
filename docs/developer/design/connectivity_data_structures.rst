Connectivity Data Structures
============================

Purpose
-------

This document describes the data structures and construction path behind
Spheral connectivity. Connectivity is the bridge between node storage and
physics pair loops: neighbor objects find candidate nodes, ``ConnectivityMap``
turns those candidates into significant neighbor lists, and ``NodePairList``
provides a flattened pair traversal for RAJA-style kernels.

The focus here is the C++ implementation: object ownership, indexing, rebuild
semantics, deterministic traversal, optional overlap/intersection structures,
and pairwise storage patterns.
For GPU-porting guidance that uses ``NodePairListView`` and flattened pair
traversal as worked examples, see :doc:`../dev/gpu_porting_patterns`.

Source Map
----------

Neighbor search:

* ``src/Neighbor/Neighbor.hh`` and ``Neighbor.cc``
* ``src/Neighbor/NestedGridNeighbor.hh`` and ``NestedGridNeighbor.cc``
* ``src/Neighbor/TreeNeighbor.hh`` and ``TreeNeighbor.cc``

Connectivity map:

* ``src/Neighbor/ConnectivityMap.hh``
* ``src/Neighbor/ConnectivityMapInline.hh``
* ``src/Neighbor/ConnectivityMap.cc``
* ``src/DataBase/DataBase.hh``, ``DataBaseInline.hh``, and ``DataBase.cc``
* ``src/Integrator/Integrator.cc``

Pair traversal and pairwise data:

* ``src/Neighbor/NodePairIdxType.hh``
* ``src/Neighbor/NodePairList.hh`` and ``NodePairList.cc``
* ``src/Neighbor/NodePairListView.hh``
* ``src/Neighbor/PairwiseField.hh`` and ``PairwiseFieldInline.hh``
* ``src/Neighbor/PairwiseFieldView.hh``
* ``src/Neighbor/PairwiseFieldElementAccessor.hh``

High-Level Design
-----------------

Connectivity is split into three layers:

``Neighbor``
  A per-``NodeList`` search object. It owns search-structure details and can
  produce master/coarse/refined candidate lists for one node or position.

``ConnectivityMap``
  A cross-node-list structure. It combines all active ``Neighbor`` objects,
  stores significant neighbors per node, creates a flattened node-pair list,
  and optionally stores overlap or intersection connectivity.

``NodePairList`` and ``PairwiseField``
  Pair-level storage. ``NodePairList`` is the flat traversal order for pair
  kernels. ``PairwiseField`` stores ephemeral values indexed by the active pair
  list.

The normal step path is:

1. Boundaries create or update ghost nodes.
2. Each node list's ``Neighbor`` updates its search structure.
3. ``DataBase::updateConnectivityMap`` rebuilds ``ConnectivityMap``.
4. ``ConnectivityMap`` stores nested neighbor vectors and a flat
   ``NodePairList``.
5. Physics packages traverse nested connectivity or flattened pairs.
6. RAJA hydro kernels usually traverse ``NodePairListView``.

Neighbor Objects
----------------

``Neighbor<Dimension>`` is the abstract base for node-list-local neighbor
search. Each ``NodeList`` owns or references one active neighbor object.

The base stores:

* ``NeighborSearchType``: ``Gather``, ``Scatter``, or ``GatherScatter``;
* kernel extent;
* pointer to the associated ``NodeList``;
* ``Field<Dimension, Vector> mNodeExtent`` for per-node search extents.

Concrete neighbor classes implement methods such as:

* ``setMasterList``;
* ``setRefineNeighborList``;
* ``updateNodes``;
* position/H overloads for scalar and tensor smoothing scales.

The base class provides shared support such as ``precullList`` and
``setMasterNeighborGroup``. ``ConnectivityMap`` uses these APIs instead of
knowing whether the underlying search is nested-grid, tree-based, or another
implementation.

ConnectivityMap Ownership and Indexing
--------------------------------------

``ConnectivityMap<Dimension>`` does not own node lists. It stores sorted
``const NodeList<Dimension>*`` pointers in ``mNodeLists``. The sorting follows
``NodeListRegistrar`` ordering so field lists, connectivity, and deterministic
traversal agree on node-list index meaning.

Its primary storage is:

::

   mOffsets
   mConnectivity

where ``mConnectivity`` is:

::

   std::vector<std::vector<std::vector<int>>>

The indexing convention is:

::

   mConnectivity[mOffsets[nodeListID] + nodeID][neighborNodeListID][k]

That entry returns the ``k``th neighbor node id in ``neighborNodeListID`` for
the source node ``(nodeListID, nodeID)``.

This layout flattens the first two dimensions, node-list id and node id, into
one vector index. ``mOffsets`` records where each node list begins in that flat
array. The last two dimensions preserve the fact that each source node has a
separate neighbor vector for every target node list.

The map also stores:

* ``mBuildGhostConnectivity``;
* ``mBuildOverlapConnectivity``;
* ``mBuildIntersectionConnectivity``;
* ``mNodePairListPtr``;
* ``mOverlapConnectivity``;
* ``mNodeTraversalIndices``;
* ``mKeys`` for domain-decomposition-independent ordering;
* ``mCouplingPtr``;
* ``mIntersectionConnectivity`` keyed by ``NodePairIdxType``.

Rebuild Semantics
-----------------

``DataBase::updateConnectivityMap`` calls:

::

   mConnectivityMapPtr->rebuild(nodeListBegin(), nodeListEnd(),
                                computeGhostConnectivity,
                                computeOverlapConnectivity,
                                computeIntersectionConnectivity);

``ConnectivityMap::rebuild`` performs the setup:

1. Requirement flags are recorded. Intersection connectivity forces ghost
   connectivity because intersections can require ghost-node information.
2. Node lists are copied into ``mNodeLists`` in registrar order.
3. The number of stored nodes per node list is selected. Domain-independent
   mode, ghost connectivity, and overlap connectivity store all nodes;
   otherwise only internal nodes are stored.
4. ``mOffsets`` is rebuilt.
5. ``computeConnectivity`` fills the actual structures.

Connectivity is therefore a per-layout/per-step structure. Anything indexed by
connectivity must be treated as invalid after a rebuild unless it is explicitly
patched by the rebuild/culling path.

computeConnectivity
-------------------

``ConnectivityMap::computeConnectivity`` is the core construction routine.
Its main stages are:

Preconditions and setup
  It checks that all neighbor objects are valid. It determines whether ghost
  connectivity is needed, builds a temporary ``DataBase`` containing the active
  node lists, and finds the maximum kernel extent.

Storage reset
  Existing ``mConnectivity`` storage is reused when possible by clearing inner
  vectors. Otherwise the full nested vector is reallocated. Intersection
  connectivity is cleared.

Domain-independent keys
  If ``NodeListRegistrar`` requests domain-decomposition-independent behavior,
  ``mKeys`` is filled with Morton-order keys. Those keys are later used to sort
  traversal indices, neighbor lists, and node pairs.

Traversal indices
  ``mNodeTraversalIndices`` is filled with the node ids to walk for each node
  list. In normal mode this is internal node order. In domain-independent mode
  it is key order over all nodes.

Candidate search
  The method gets global position and ``H`` field lists from the temporary
  database. For unfinished nodes it calls
  ``Neighbor::setMasterNeighborGroup`` to get master lists and coarse
  candidates across all node lists.

Significance test
  For each source node and each coarse candidate, the code computes
  ``eta2i = |Hi * rij|^2`` and ``eta2j = |Hj * rij|^2``. A pair is significant
  if either normalized distance is within ``kernelExtent^2``. Self-interactions
  are excluded.

Neighbor insertion
  Significant candidate node ids are appended to
  ``mConnectivity[offset + sourceNode][targetNodeList]``.

Pair-list insertion
  ``calculatePairInteraction`` decides whether the pair should be added to the
  flattened ``NodePairList``. This prevents double-counting when both nodes
  will see each other through nested connectivity.

Neighbor ordering
  Each target-node-list neighbor vector is sorted. Normal mode sorts by node
  id. Domain-independent mode sorts by Morton key.

Pair ordering
  The accumulated ``NodePairIdxType`` vector becomes a new ``NodePairList``.
  Pairs are sorted by hash in normal mode and by key-derived ordering in
  domain-independent mode.

Optional structures
  If requested, overlap connectivity and intersection connectivity are computed
  after the base connectivity and pair list are available.

Postconditions
  The method verifies that all required nodes were completed and that the map
  is valid.

NodePairIdxType
---------------

``NodePairIdxType`` is the compact value type used to identify one interacting
pair. It contains:

* ``i_node``;
* ``i_list``;
* ``j_node``;
* ``j_list``;
* ``f_couple`` for an effective coupling fraction.

It is marked ``SPHERAL_HOST_DEVICE`` and can be used in device kernels.

The ``hash`` method normalizes pair ordering before packing node-list and node
ids into a ``size_t``. This makes ``(i, j)`` and ``(j, i)`` compare as the same
pair for equality/order purposes. The packed representation assumes:

* fewer than 32 node lists;
* node ids below the maximum representable packed node index.

These limits are enforced by contract checks.

NodePairList and NodePairListView
---------------------------------

``NodePairList`` owns:

* ``std::vector<NodePairIdxType> mNodePairList``;
* a mutable lazy lookup map from pair to pair index.

It derives from ``NodePairListView``. ``initView`` wraps the vector storage with
``GPUUtils::initMAView`` so the inherited view can be copied into kernels.

``NodePairListView`` exposes:

* ``operator[](k)``;
* ``size()``;
* ``data()``;
* ``move(space)``;
* ``touch(space)``.

RAJA SPH kernels use this view directly:

::

   const auto& pairs_owner = connectivityMap.nodePairList();
   const auto pairs = pairs_owner.view();

   RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, pairs.size()),
     [=] SPHERAL_HOST_DEVICE(size_t kk) {
       const auto i = pairs[kk].i_node;
       const auto nodeListi = pairs[kk].i_list;
       ...
     });

The flattened pair list is the accelerator-friendly representation. The nested
``mConnectivity`` vectors remain useful for algorithms that need all neighbors
of one node or need set operations over neighbor lists.

calculatePairInteraction
------------------------

``ConnectivityMap::calculatePairInteraction`` determines whether a significant
neighbor relation should be represented in the flattened pair list. In the
current active logic it returns true when:

* the neighbor node list id is greater than the source node list id;
* or both nodes are in the same node list and ``j > i``;
* or the neighbor node list id is less than the source node list id and ``j``
  is a ghost node.

The first two cases keep one canonical ordering for internal/internal pairs.
The ghost-node case allows cross-boundary interactions that otherwise would be
missed when only one side of the relation is locally owned.

Domain-Independent Ordering
---------------------------

When ``NodeListRegistrar::domainDecompositionIndependent`` is active,
``ConnectivityMap`` uses Morton-order keys from ``mortonOrderIndices``. Those
keys affect:

* which nodes are traversed first;
* ordering of neighbors inside each per-node neighbor vector;
* ordering and orientation of ``NodePairList`` entries.

The goal is to make traversal less sensitive to local storage order and domain
decomposition. Developers adding connectivity-consuming algorithms should use
``ConnectivityMap::begin/end`` or ``ithNode`` when they need the stored
deterministic node order instead of raw node ids.

Overlap Connectivity
--------------------

Overlap connectivity is optional and is more expensive than basic neighbor
connectivity. When requested, ``mOverlapConnectivity`` starts as a copy of
``mConnectivity``. The builder then adds nodes whose support domains overlap
through gather/scatter relationships.

The implementation uses ``insertUnique`` so overlap neighbor vectors remain
sorted and duplicate-free. In domain-independent mode the insertion position is
computed by key order; otherwise it uses node-id order.

Use overlap connectivity only for packages that explicitly need common-support
or overlap information. The integrator collects package requirements before
building connectivity so this extra work is shared when needed and avoided when
not needed.

Intersection Connectivity
-------------------------

Intersection connectivity is also optional. It answers a different question:
for a pair ``(i, j)`` that is already an interacting pair, which nodes are in
the intersection of their neighbor sets?

When ``mBuildIntersectionConnectivity`` is true:

1. ghost connectivity is forced;
2. the base ``NodePairList`` is built;
3. ``connectivityIntersectionForNodes`` is called for every pair;
4. the result is stored in ``mIntersectionConnectivity[pair]``.

``connectivityIntersectionForNodes`` normally intersects the sorted neighbor
vectors for each node list. If a position field list is supplied, it can also
filter intersection nodes by whether the closest point lies on the segment
between the pair endpoints. The method includes the pair endpoints themselves
in the returned vectors.

PairwiseField
-------------

``PairwiseField<Dimension, Value, numElements>`` stores values indexed by the
active ``NodePairList``. It is explicitly ephemeral because meshfree
connectivity can change from step to step.

The object owns:

* a weak pointer to the ``NodePairList``;
* ``std::vector<Value> mArray`` sized to ``numElements * pairs.size()``;
* a view span/managed array inherited from ``PairwiseFieldView``.

Indexing by ``NodePairIdxType`` goes through the pair list's lazy lookup:

::

   p->index(pair)

The ``numElements`` template parameter controls the element accessor:

* ``numElements == 1`` returns one value by reference;
* ``numElements > 1`` returns a pointer to the first value in a strided block.

This is used for pair quantities such as compatible-energy pair accelerations,
where values are naturally stored per pair rather than per node.

Because ``patchConnectivity`` deliberately reallocates the ``NodePairList`` in
some paths, existing ``PairwiseField`` objects can become invalid. Pairwise
fields should be built for the current connectivity and discarded when
connectivity changes.

Integrator and Boundary Interaction
-----------------------------------

``Integrator::setGhostNodes`` is the main caller-side workflow around
connectivity. The relevant sequence is:

1. Boundaries create/update ghost nodes.
2. Neighbor structures are updated.
3. If connectivity is required, ``DataBase::updateConnectivityMap`` rebuilds
   the map with the collected ghost/overlap/intersection requirements.
4. If ghost culling is enabled and allowed, the integrator marks active nodes,
   asks boundaries to cull ghost nodes, patches the connectivity map, deletes
   culled ghost nodes from node lists, and updates neighbor structures again.

This ordering matters. Connectivity is built after ghost nodes exist, because
ghost nodes participate in pair detection and communication-boundary physics.
If ghost nodes are culled, the connectivity map is patched before the node
lists physically delete the culled nodes.

Patch Semantics
---------------

``ConnectivityMap::patchConnectivity`` is used after ghost culling. It receives:

* ``flags``: zero for deleted nodes and one for preserved nodes;
* ``old2new``: old node id to new node id mapping.

It patches:

* traversal indices;
* nested neighbor vectors;
* pair list entries;
* intersection connectivity when maintained.

It also removes references to deleted nodes and re-sorts neighbor vectors. The
method cannot fully validate the map before node lists are resized, so the
caller validates after deletion is complete.

Common Access Patterns
----------------------

Nested connectivity for one node:

::

   const auto& neighbors = cm.connectivityForNode(nodeListID, i);
   for (auto jNodeList = 0u; jNodeList < neighbors.size(); ++jNodeList) {
     for (const auto j: neighbors[jNodeList]) {
       ...
     }
   }

Flattened pair traversal:

::

   const auto& pairs = cm.nodePairList();
   for (const auto& pair: pairs) {
     ...
   }

RAJA flattened pair traversal:

::

   const auto pairs = cm.nodePairList().view();
   RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, pairs.size()),
     [=] SPHERAL_HOST_DEVICE(size_t kk) {
       const auto pair = pairs[kk];
       ...
     });

Deterministic node traversal:

::

   for (auto itr = cm.begin(nodeListID); itr != cm.end(nodeListID); ++itr) {
     const auto i = *itr;
     ...
   }

Design Invariants
-----------------

Important invariants include:

* ``mNodeLists`` are in registrar order;
* ``mOffsets.size() == mNodeLists.size()``;
* ``mConnectivity`` is sized for either internal nodes or all nodes according
  to the current ghost/overlap/domain-independent requirements;
* each ``mConnectivity[offset + nodeID]`` entry has one vector per node list;
* neighbor vectors are sorted by node id or by domain-independent key;
* ``NodePairList`` contains each pair interaction once;
* ``NodePairList`` order is deterministic for the active ordering mode;
* ``PairwiseField`` size matches the active pair-list size;
* connectivity-consuming code must not retain pairwise views across rebuilds.

Extension Guidance
------------------

When adding a new neighbor implementation:

* implement the ``Neighbor`` virtual search methods;
* keep search-structure ownership inside the neighbor class;
* ensure ``updateNodes`` leaves the object valid before connectivity rebuilds;
* preserve the meaning of gather/scatter/gather-scatter search types.

When adding a physics package that needs connectivity:

* declare the requirement through the ``Physics`` requirement hooks;
* request ghost, overlap, or intersection connectivity only when needed;
* use nested connectivity when the algorithm is node-centered;
* use ``NodePairList`` when the algorithm is pair-centered or device-oriented;
* use ``PairwiseField`` only for data tied to the current pair list.

When adding a pair kernel:

* use ``connectivityMap.nodePairList().view()`` for regular pair traversal;
* capture field-list views and pair-list views by value;
* use atomics for per-node outputs touched by multiple pairs;
* move pairwise outputs back to CPU if host code consumes them after the
  kernel.

Common Failure Modes
--------------------

Using connectivity before neighbors are updated
  ``ConnectivityMap`` assumes each node list's ``Neighbor`` is valid.

Retaining pairwise fields across rebuilds
  ``PairwiseField`` is tied to one ``NodePairList``. Rebuild or patch
  connectivity invalidates that association.

Assuming ghost connectivity exists
  Normal connectivity can store only internal-node entries. Code that queries
  ghost nodes must ensure ghost connectivity was requested or domain-independent
  mode is active.

Double-counting pair interactions
  Nested neighbor lists can contain both directions. Pair-centered algorithms
  should use ``NodePairList`` or reproduce ``calculatePairInteraction`` logic.

Breaking deterministic order
  Sorting by raw node id in domain-independent mode defeats the stored ordering.
  Use keys or existing map traversal helpers.

Ignoring overlap/intersection cost
  These structures are substantially more expensive than base connectivity.
  Request them only through package requirements that truly need them.

Relationship to the Other Design Docs
-------------------------------------

``simulation_lifecycle_and_main_loop.rst`` explains when ghost and connectivity
updates occur during startup and each step.

``integrator_and_state_update_model.rst`` explains how physics package
requirements cause the integrator to request connectivity.

``value_view_and_device_execution_model.rst`` explains the owner/view pattern
used by ``NodePairList`` and ``PairwiseField``.

``raja_chai_execution_patterns.rst`` explains how ``NodePairListView`` is used
in RAJA pair kernels.
