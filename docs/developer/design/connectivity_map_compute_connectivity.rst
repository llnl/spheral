######################
ConnectivityMap Design
######################

ConnectivityMap Object
======================

``ConnectivityMap<Dimension>`` is the ``DataBase``-owned object that records
which nodes can interact with which other nodes. It preserves ``NodeList``
ordering according to ``NodeListRegistrar``, provides the public
``connectivityForNode`` style queries used by host algorithms and Python tests,
and exposes flat pair data through ``NodePairList`` for pairwise algorithms.

This note currently focuses on the ``computeConnectivity`` builder because that
is the part of ``ConnectivityMap`` most relevant to a future GPU port. The
broader object description can be expanded as needed to cover patching,
global-connectivity queries, coupling functors, traversal APIs, and Python
binding behavior.

``computeConnectivity`` Implementation
--------------------------------------

Purpose and Scope
^^^^^^^^^^^^^^^^^

``ConnectivityMap<Dimension>::computeConnectivity`` is the internal builder for
the per-node neighbor graph used by the ``DataBase``, hydrodynamics packages,
kernel integration utilities, Python tests, and lower-level neighbor queries. It
builds several related products:

- the full per-node connectivity grouped by target ``NodeList``;
- a flat ``NodePairList`` used by pairwise algorithms such as the RAJA SPH
  derivative kernels;
- optional overlap connectivity;
- optional precomputed intersection connectivity;
- deterministic traversal and pair ordering when domain decomposition
  independence is enabled.

The current implementation is a host-side producer. GPU-enabled derivative code
does not build connectivity on the GPU; it consumes the already-built
``NodePairList`` through ``NodePairListView`` and consumes state through
``FieldView`` and ``FieldListView``.

Primary Calling Graph
^^^^^^^^^^^^^^^^^^^^^

The normal rebuild path is:

.. code-block:: text

   DataBase<Dimension>::updateConnectivityMap(...)
     -> ConnectivityMap<Dimension>::rebuild(begin, end, flags...)
        -> sort/copy NodeLists in NodeListRegistrar order
        -> compute prefix offsets
        -> ConnectivityMap<Dimension>::computeConnectivity()

Inside ``computeConnectivity`` the high-level flow is:

.. code-block:: text

   computeConnectivity()
     -> validate Neighbor state on all NodeLists
     -> create temporary DataBase over mNodeLists
     -> gather global position and H FieldLists
     -> optionally compute Morton keys for domain-independent ordering
     -> build mNodeTraversalIndices
     -> for each seed node not already processed:
        -> Neighbor<Dimension>::setMasterNeighborGroup(...)
           -> Neighbor::setMasterList(...) on each concrete Neighbor
           -> Neighbor::precullList(...) on each concrete Neighbor
        -> for each master node:
           -> walk coarse neighbor candidates
           -> test eta_i and eta_j against kernelExtent^2
           -> append neighbor IDs into mConnectivity
           -> append unique pair work into thread-local NodePair vectors
           -> sort each per-NodeList neighbor vector
           -> mark node done
        -> merge thread-local pair vectors
     -> construct mNodePairListPtr
     -> sort NodePairList
     -> optionally build mOverlapConnectivity
     -> optionally precompute mIntersectionConnectivity
     -> validate postconditions

Important downstream consumers include:

- ``DataBase::connectivityMap`` and ``DataBase::connectivityMapPtr``;
- SPH/CRKSPH/GSPH/FSISPH derivative evaluation;
- ``SPH_RAJA`` and ``SolidSPH_RAJA``, which consume
  ``connectivityMap.nodePairList().view()``;
- ``FlatConnectivity``, which flattens the nested connectivity for kernel
  integration indexing;
- Python bindings in ``src/PYB11/Neighbor/ConnectivityMap.py``;
- unit tests under ``tests/unit/Neighbor`` and C++ view tests under
  ``tests/cpp/Neighbor``.

Public Entry Points
^^^^^^^^^^^^^^^^^^^

``DataBase<Dimension>::updateConnectivityMap`` is the common explicit update
entry point. It requires an existing ``mConnectivityMapPtr`` and forwards the
``computeGhostConnectivity``, ``computeOverlapConnectivity``, and
``computeIntersectionConnectivity`` flags to ``ConnectivityMap::rebuild``.

``DataBase<Dimension>::connectivityMap(flags...)`` lazily calls
``updateConnectivityMap`` only if the pointer is missing. In current
``DataBase`` construction the pointer is normally present, so callers that need
fresh connectivity should call ``updateConnectivityMap`` before querying.

``ConnectivityMap::rebuild`` performs setup shared by all builds:

- records whether ghost, overlap, and intersection connectivity are requested;
- forces ghost connectivity when intersection connectivity is requested;
- copies ``NodeList`` pointers into the registrar-prescribed ordering;
- computes ``mOffsets`` as a prefix sum over the number of represented nodes in
  each ``NodeList``;
- calls ``computeConnectivity``.

Core Data Structures
^^^^^^^^^^^^^^^^^^^^

``mNodeLists``
   ``std::vector<const NodeList<Dimension>*>`` in ``NodeListRegistrar`` order.
   The ordering matters because public APIs address connectivity by
   ``NodeList`` index, and deterministic ordering relies on registrar ordering.

``mBuildGhostConnectivity``
   True when ghost-node connectivity is requested or when intersection
   connectivity forces it. ``computeConnectivity`` also treats domain
   decomposition independent mode as requiring ghost connectivity.

``mBuildOverlapConnectivity``
   True when overlap connectivity should be built. This also expands the
   represented node count to include ghost nodes.

``mBuildIntersectionConnectivity``
   True when intersection connectivity should be precomputed for each node pair.
   This forces ghost connectivity in ``rebuild``.

``mOffsets``
   Prefix offsets mapping ``(nodeListID, nodeID)`` to the first index of that
   node's connectivity:

   .. code-block:: text

      storageIndex = mOffsets[nodeListID] + nodeID

``mConnectivity``
   The primary nested storage:

   .. code-block:: text

      mConnectivity[mOffsets[iList] + iNode][jList][k]

   The final ``k`` index walks neighbor node IDs in ``jList``. Neighbor vectors
   are sorted by local node ID in normal mode and by Morton key in
   domain-decomposition-independent mode.

``mNodePairListPtr``
   ``std::shared_ptr<NodePairList>`` built from the connectivity walk. This is
   the key bridge to accelerator kernels. ``NodePairList`` owns a host
   ``std::vector<NodePairIdxType>`` and initializes a CHAI-backed
   ``NodePairListView`` for host/device capture.

``mOverlapConnectivity``
   Same storage shape as ``mConnectivity``. It starts as a copy of
   ``mConnectivity`` and is expanded with nodes whose kernels overlap through
   gather/scatter relationships.

``mNodeTraversalIndices``
   Per-``NodeList`` node traversal order. In normal mode this is local node ID
   order over internal nodes. In domain-decomposition-independent mode it
   includes all nodes and is sorted by Morton key.

``mKeys``
   ``FieldList<Dimension, KeyTraits::Key>`` storing Morton-order keys when
   domain-decomposition-independent mode is enabled. The keys drive deterministic
   neighbor ordering and pair ordering.

``mCouplingPtr``
   ``std::shared_ptr<NodeCoupling>``. The current ``computeConnectivity`` path
   does not apply this coupling while deciding geometric connectivity. It is
   part of the broader ``ConnectivityMap`` API.

``mIntersectionConnectivity``
   ``std::unordered_map<NodePairIdxType, std::vector<std::vector<int>>>`` keyed
   by node pair. Values are grouped by ``NodeList`` and contain common
   neighbors for that pair, optionally filtered by closest point on the segment
   between the pair positions.

Neighbor Search Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``computeConnectivity`` delegates candidate generation to the active
``Neighbor`` implementation for each ``NodeList``.

``Neighbor::setMasterNeighborGroup`` performs two stages:

1. It calls virtual ``Neighbor::setMasterList`` on each concrete neighbor
   implementation. For example, ``NestedGridNeighbor`` uses nested grid cells
   while ``TreeNeighbor`` uses tree data.
2. It computes a bounding region covering all master nodes and the seed node,
   then calls ``Neighbor::precullList`` to reduce coarse candidate lists.

The connectivity builder then performs the exact pair acceptance test itself:

.. code-block:: text

   rij = ri - rj
   eta2i = (Hi * rij).magnitude2()
   eta2j = (Hj * rij).magnitude2()
   accept if eta2i <= kernelExtent^2 or eta2j <= kernelExtent^2

This gather-or-scatter test is what turns coarse neighbor candidates into
significant connectivity entries.

Build Algorithm Details
^^^^^^^^^^^^^^^^^^^^^^^

Preconditions
"""""""""""""

Before building, ``computeConnectivity`` verifies:

- every stored ``NodeList`` has a valid ``Neighbor``;
- ``mOffsets`` has one entry per ``NodeList``.

It then determines whether ghost nodes are represented for this build:

.. code-block:: text

   ghostConnectivity =
       mBuildGhostConnectivity or
       mBuildOverlapConnectivity or
       domainDecompIndependent

The represented connectivity size is the offset of the last ``NodeList`` plus
either its total node count or internal-node count depending on that flag.

Initialization
""""""""""""""

The builder constructs a temporary ``DataBase`` over ``mNodeLists`` and reads
``globalPosition`` and ``globalHfield`` from it. It also finds the maximum kernel
extent across all ``Neighbor`` objects and squares it for the acceptance test.

Existing connectivity storage is reused if it has the correct shape; otherwise
the nested vectors are reallocated. Reused per-node neighbor vectors are cleared.
Intersection connectivity is always cleared before rebuild.

Domain-Independent Ordering
"""""""""""""""""""""""""""

If ``NodeListRegistrar::domainDecompositionIndependent()`` is true, the builder
computes ``mKeys = mortonOrderIndices(dataBase)`` and orders traversal by those
keys. Otherwise traversal is local ID order.

Neighbor lists are sorted using the same policy:

- Morton key order in domain-independent mode;
- local node ID order otherwise.

Pair lists are sorted after construction:

- ``sortPairs`` normalizes each pair orientation by comparing keys and then
  sorts by a hash of the two keys in domain-independent mode;
- ``std::sort`` over ``NodePairIdxType`` otherwise.

Main Connectivity Walk
""""""""""""""""""""""

The outer loops choose seed nodes and skip any already marked complete in
``flagNodeDone``. For each seed, ``Neighbor::setMasterNeighborGroup`` returns a
set of master nodes and coarse candidate neighbors for every ``NodeList``.

The master-node loop is OpenMP-parallelized. Each thread owns a private
``std::vector<NodePairIdxType>`` for pair output. For each master node:

1. The code gets references to that node's position, smoothing tensor, work
   field, and per-``NodeList`` neighbor storage.
2. It loops over every candidate ``j`` from every target ``NodeList``.
3. It computes the ``eta`` distance using both ``Hi`` and ``Hj``.
4. If either side is within ``kernelExtent``, it appends ``j`` to the target
   neighbor vector unless the pair is a self-interaction.
5. It appends one ``NodePairIdxType`` when ``calculatePairInteraction`` says
   this pair should be evaluated exactly once by pairwise algorithms.
6. It sorts each target neighbor vector.
7. It marks the master node done.

Thread-local pair vectors are merged into the outer ``nodePairs`` vector in an
OpenMP critical section.

Pair Uniqueness Policy
""""""""""""""""""""""

``calculatePairInteraction`` determines which accepted pairs enter
``NodePairList``:

.. code-block:: text

   nodeListj > nodeListi
   or nodeListj == nodeListi and j > i
   or nodeListj < nodeListi and j is a ghost node

This produces a single pair entry for ordinary internal/internal interactions
while preserving interactions involving ghost nodes when the target list appears
earlier in ``NodeList`` order.

Overlap Connectivity
""""""""""""""""""""

When requested, overlap connectivity starts as a copy of primary connectivity.
The builder then expands it using gather/scatter overlap logic:

1. For every node ``i``, inspect every current neighbor ``j1``.
2. If ``j1`` is a gather neighbor of ``i``, check whether ``i`` and ``j1``
   directly overlap and insert both directions.
3. Walk gather neighbors ``j2`` of ``j1``. If ``j2`` is a scatter neighbor of
   ``j1``, insert overlap between ``i`` and ``j2`` in both directions.

``insertUnique`` preserves sorted unique insertion and uses Morton key ordering
when domain-independent mode is active.

Intersection Connectivity
"""""""""""""""""""""""""

When requested, the builder loops over the completed ``NodePairList`` in
parallel. Each thread creates a private unordered map:

.. code-block:: text

   NodePairIdxType -> vector<vector<int>>

Each value is produced by ``connectivityIntersectionForNodes``. The helper
intersects the two nodes' primary connectivity lists per ``NodeList`` and adds
the pair endpoints themselves. If a position ``FieldList`` is provided, common
neighbors are filtered to those closest to the segment between the pair
positions.

Thread-local intersection maps are merged into ``mIntersectionConnectivity`` in
an OpenMP critical section.

Postconditions and Validation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``computeConnectivity`` checks that every represented node was marked done, then
calls ``valid()``.

``valid()`` verifies:

- offsets and storage size match represented node counts;
- ``NodeList`` ordering matches the registrar;
- every represented node has one per-``NodeList`` neighbor vector;
- neighbor IDs are in range;
- neighbor vectors are sorted and unique;
- connectivity is symmetric where required;
- Morton keys exist in domain-independent mode;
- traversal indices contain every represented node exactly once.

Relationship to GPU/RAJA Derivative Work
========================================

The RAJA derivative implementation consumes connectivity at a different level
than ``computeConnectivity`` produces it. In ``SPH_RAJA``, the derivative code:

1. asks the ``DataBase`` for the current ``ConnectivityMap``;
2. extracts ``connectivityMap.nodePairList()``;
3. converts it to ``NodePairListView``;
4. captures that view, ``FieldListView`` state, kernel views, and viscosity
   views in ``RAJA::forall<EXEC_POLICY>`` pair loops.

This path works because all data captured by the kernel is flat or view-based.
The full ``ConnectivityMap`` primary storage is still nested STL storage and is
not itself a device-side data structure.

Testing Notes
=============

The existing tests cover important invariants:

- ``tests/unit/Neighbor/NeighborTestBase.py`` checks
  ``connectivityForNode``, ``nodePairList``, overlap connectivity, and
  intersection behavior against direct search answers.
- ``tests/unit/Neighbor/testDistributedConnectivity.py`` compares serial and
  distributed connectivity behavior.
- ``tests/cpp/Neighbor/nodepairlistview_tests.cc`` verifies CHAI/RAJA movement
  and device capture for ``NodePairListView``.
- ``tests/cpp/Neighbor/pairwisefieldview_tests.cc`` verifies device capture for
  pairwise fields backed by a ``NodePairList``.

Any implementation that changes ``computeConnectivity`` should keep the Python
neighbor tests as correctness gates and add C++ tests for any new flat or
device-facing connectivity representation.

NodePairList-First GPU Connectivity
===================================

Author's quote:

.. epigraph::

   I was thinking for porting the current computeConnectivity to GPU it would
   make sense to build the NodePairList alone first (which should be simpler in
   a GPU kernel), and after the fact transfer that information to the
   alternative per node query data. Perhaps making that last bit only on demand
   since many of our models don't need the older per node query.

Plan:

1. Build a ``NodePairList``-only connectivity path that produces accepted,
   canonical ``NodePairIdxType`` pairs without constructing the legacy per-node
   query data in the same step.
2. Verify the ``NodePairList``-only path against the current full
   ``computeConnectivity`` result by comparing sorted pair lists across the
   existing neighbor correctness tests.
3. Derive the alternative per-node query data from the completed
   ``NodePairList`` after the fact, preferably on demand for
   ``connectivityForNode``-style consumers.
