Neighbor Family and Usage
=========================

Purpose
-------

This document analyzes ``Neighbor`` and its concrete family in the current
Spheral implementation. It focuses on ownership, base-class responsibilities,
concrete search structures, rebuild timing, and how neighbor search feeds
``ConnectivityMap`` and pair traversal.

For the broader connectivity data model, see
:doc:`connectivity_data_structures`. For the device-facing value/view and
managed view pointer model, see :doc:`value_view_and_device_execution_model`.

Source Map
----------

Neighbor family:

* ``src/Neighbor/Neighbor.hh`` and ``Neighbor.cc``
* ``src/Neighbor/NeighborInline.hh``
* ``src/Neighbor/TreeNeighbor.hh`` and ``TreeNeighbor.cc``
* ``src/Neighbor/NestedGridNeighbor.hh`` and ``NestedGridNeighbor.cc``
* ``src/Neighbor/GridCellIndex.hh``
* ``src/Neighbor/GridCellPlane.hh``

Neighbor users:

* ``src/NodeList/NodeList.hh`` and ``NodeList.cc``
* ``src/DataBase/DataBase.cc``
* ``src/Neighbor/ConnectivityMap.cc``
* ``src/PYB11/Neighbor/*.py``
* ``src/NodeList/*NodeLists.py``

Role in the Object Model
------------------------

``Neighbor<Dimension>`` is the node-list-local search abstraction. It does not
store final cross-node-list connectivity. Instead, it maintains enough
per-node-list search structure to answer questions such as:

* Which nodes are masters for this position or node?
* Which nodes are coarse candidates?
* Which coarse candidates refine to actual local neighbors?

``ConnectivityMap`` combines the answers from all active ``Neighbor`` objects,
tests pair significance across node lists, stores nested neighbor vectors, and
builds the flattened ``NodePairList`` used by pair kernels.

The normal dependency direction is:

::

   NodeList
     has one active Neighbor
     |
     | updateNodes()
     v
   Neighbor search structure
     |
     | setMasterList / setRefineNeighborList
     v
   ConnectivityMap
     stores final cross-node-list connectivity and NodePairList

Ownership and Registration
--------------------------

At the C++ level, ``NodeList`` stores a raw pointer to its active neighbor:

::

   Neighbor<Dimension>* mNeighborPtr;

``Neighbor`` registers itself with the ``NodeList`` during construction and
unregisters itself during destruction. ``NodeList::neighbor()`` checks that a
neighbor is registered and returns a reference.

This means the neighbor object is associated with exactly one active
``NodeList`` at a time, but the ``NodeList`` pointer is not an owning smart
pointer. Python construction helpers usually keep the concrete neighbor object
alive by storing it on the Python-facing node list helper, for example as
``result._neighbor``.

The base ``Neighbor`` also owns:

* ``NeighborSearchType mSearchType``;
* ``double mKernelExtent``;
* ``NodeList<Dimension>* mNodeListPtr``;
* ``Field<Dimension, Vector> mNodeExtent``.

``mNodeExtent`` is a normal ``Field`` attached to the node list. It records the
axis-aligned extent implied by each node's smoothing scale and the kernel
extent.

Base-Class Contract
-------------------

The abstract base defines the common search API:

* ``setMasterList`` for node ids, positions, scalar ``H``, tensor ``H``, and
  plane mappings;
* ``setRefineNeighborList`` for node ids and positions;
* ``updateNodes`` for all nodes or selected node ids;
* ``reinitialize`` hooks for concrete structures that need bounding boxes or a
  target smoothing length;
* ``valid`` checks.

The base also provides shared helpers:

``HExtent``
  Converts scalar or tensor smoothing-scale information into a Cartesian extent
  for search culling.

``setNodeExtents``
  Refreshes the ``mNodeExtent`` field from the node list's current ``H`` field.

``precullList``
  Applies gather, scatter, or gather-scatter bounding-box tests to a coarse
  candidate list.

``setMasterNeighborGroup``
  Static helper for multi-node-list candidate search. It calls each node
  list's concrete ``setMasterList``, accumulates the bounding box of all master
  nodes, and then preculls each coarse list with the combined master extent.

The static helper is the key cross-node-list bridge. It lets
``ConnectivityMap`` ask all node-list-local neighbor objects for candidates
without knowing the concrete search implementation behind each one.

Search Types
------------

``NeighborSearchType`` controls how ``precullList`` interprets extents:

``Gather``
  Keep candidates whose node position is inside the master support extent.

``Scatter``
  Keep candidates whose support extent overlaps the master position/extent box.

``GatherScatter``
  Keep candidates that satisfy either gather or scatter style tests.

The search type affects candidate culling only. ``ConnectivityMap`` still
performs the final significance test with normalized distances
``|Hi * rij|^2`` and ``|Hj * rij|^2``.

``TreeNeighbor``
----------------

``TreeNeighbor`` is the current tree-based neighbor implementation. It stores
an oct-tree-like hierarchy keyed by quantized coordinates:

* ``LevelKey`` identifies the tree level;
* ``CellKey`` identifies a cell within a level;
* each cell stores daughter keys, daughter pointers, and member node ids;
* ``mTree`` is a vector of hash maps, one per occupied level.

``TreeNeighbor::updateNodes`` rebuilds the tree from the current node
positions and smoothing tensors. It clears the old tree, inserts all nodes,
sorts local cell contents, merges thread-local trees, constructs daughter
pointers, and refreshes node extents.

Important properties:

* selected-node updates currently fall back to a full rebuild;
* ``reinitialize`` can reset the search box and target scale;
* tree search methods produce master/coarse/refined lists through the base
  ``Neighbor`` interface;
* serialization helpers exist for the tree cells.

``TreeNeighbor`` is usually the default neighbor type selected by Python node
list helpers.

``NestedGridNeighbor``
----------------------

``NestedGridNeighbor`` is the older nested-grid implementation. It stores a set
of grid levels and linked lists of nodes per occupied cell:

* ``mGridCellHead`` maps grid cells to the head node id;
* ``mNodeInCell`` records each node's cell on each grid level;
* ``mNextNodeInCell`` stores linked-list next pointers;
* ``mNodeOnGridLevel`` records each node's home grid level;
* ``mOccupiedGridCells`` records occupied cells per level.

``NestedGridNeighbor::updateNodes`` resizes and clears those structures,
assigns each node to a grid level and cell, links the node into the cell list,
rebuilds occupied-cell lists, and refreshes node extents.

It supports selected-node updates more directly than ``TreeNeighbor`` because
the linked-list representation can unlink and relink individual nodes. Current
Python helpers warn that ``NestedGridNeighbor`` is deprecated and suggest
``TreeNeighbor``.

ConnectivityMap Usage
---------------------

``ConnectivityMap::computeConnectivity`` assumes each node list has a valid
neighbor object. The precondition checks:

::

   (**itr).neighbor().valid()

The construction flow is:

1. Create a temporary ``DataBase`` containing the active node lists.
2. Get global position and ``H`` field lists.
3. Iterate source node lists and source nodes.
4. Call ``Neighbor::setMasterNeighborGroup`` for unfinished source nodes.
5. For each master node, walk coarse candidates from every target node list.
6. Apply the final significance test.
7. Insert nested connectivity and, when appropriate, a ``NodePairIdxType`` for
   flattened pair traversal.

This split is intentional:

* concrete ``Neighbor`` objects own local search acceleration structures;
* ``ConnectivityMap`` owns cross-node-list pair significance, ordering,
  pair-list construction, and optional overlap/intersection structures.

DataBase and Integrator Usage
-----------------------------

``DataBase::reinitializeNeighbors`` computes a global bounding box and an
average smoothing scale, calls ``neighbor.reinitialize(xmin, xmax, havg)`` on
each node list, and then calls ``neighbor.updateNodes()``.

``Integrator::setGhostNodes`` drives the per-step workflow around boundaries
and connectivity:

1. boundaries create or update ghost nodes;
2. neighbor structures are updated;
3. ``DataBase::updateConnectivityMap`` rebuilds connectivity;
4. optional ghost culling patches connectivity;
5. node lists delete culled ghosts;
6. neighbor structures are updated again.

The ordering matters because neighbor structures depend on current node
positions, smoothing tensors, and ghost-node layout. Connectivity is invalid
once the neighbor inputs change.

Python Exposure and Construction
--------------------------------

The PYB11 layer exposes:

* ``NeighborSearchType``;
* ``Neighbor``;
* ``NestedGridNeighbor``;
* ``TreeNeighbor``;
* ``ConnectivityMap``;
* pair-list and pairwise-field types.

The Python node-list helper modules choose and attach a concrete neighbor. The
common path constructs a ``TreeNeighbor`` and stores it on the returned node
list wrapper. ``NestedGridNeighbor`` remains available but emits deprecation
warnings in those helpers.

Design Invariants
-----------------

Important invariants include:

* every active ``NodeList`` must have a valid neighbor before connectivity is
  built;
* ``Neighbor`` extents must reflect the current ``H`` field and kernel extent;
* concrete search structures must be refreshed after node positions,
  smoothing tensors, or node counts change;
* ``ConnectivityMap`` owns final pair significance and ordering, not the
  concrete neighbor implementation;
* ghost-node layout changes require neighbor and connectivity refreshes;
* neighbor-derived connectivity and pair views become stale across a rebuild.

Future Device-Facing Questions
------------------------------

The current ``Neighbor`` family is host-oriented. Making
``NodePairListView`` device-facing does not by itself make neighbor search
device-facing.

Questions to answer before porting neighbor search include:

* Which phase should run on device: search-structure rebuild, candidate search,
  final significance testing, pair-list construction, or only pair traversal?
* Can the target algorithm consume a flattened schedule, or does it need
  dynamic per-node candidate lists?
* What replaces host containers such as nested ``std::vector`` and tree hash
  maps if construction moves to device?
* Is polymorphism required on device, or can the host choose a concrete
  neighbor path before launch?
* What view or data refresh points follow ghost culling, redistribution, or
  smoothing-scale updates?

Those questions map to the object shapes in
:doc:`value_view_and_device_execution_model`: a future neighbor port may need a
value view for compact traversal data, a staged device structure for
construction, host-side dispatch over concrete neighbor implementations, or a
combination of those choices.
