Simulation Lifecycle and Main Loop
==================================

Purpose
-------

This document explains the execution path that advances a Spheral simulation.
It focuses on the C++ implementation and the abstractions used at that level,
while also showing how the Python controller drives the compiled step engine.

The most important point for new developers is that Spheral does not have one
monolithic C++ ``main`` routine as the primary user-facing entry point. A
typical run is driven from Python. Python owns run control, restart/viz
scheduling, and user-facing workflow. The C++ implementation owns the time-step
mechanics: state assembly, boundary handling, connectivity, derivative
evaluation, state update, and package finalization.

Source Map
----------

The main files involved in the lifecycle are:

* ``src/SimulationControl/SpheralController.py``: Python run controller. It
  performs startup orchestration, calls ``integrator.step(goalTime)`` in a loop,
  and runs periodic work.
* ``src/Integrator/Integrator.hh`` and ``src/Integrator/Integrator.cc``: Base
  C++ step engine. It owns the common package, boundary, connectivity, state,
  derivative, time-step selection, retry, and restart logic.
* ``src/Integrator/SynchronousRK2.cc``, ``src/Integrator/PredictorCorrector.cc``,
  and other concrete integrators: time-centering implementations. They decide
  when derivatives are evaluated and how ``State::update`` is applied.
* ``src/Physics/Physics.hh`` and ``src/Physics/Physics.cc``: Base package
  contract. Physics packages register state, register derivative fields,
  compute derivatives, vote on ``dt``, and participate in boundary and
  finalization hooks.
* ``src/DataBase/DataBase.hh``: Runtime collection of node lists and shared
  connectivity.
* ``src/DataBase/State.hh``, ``src/DataBase/State.cc``,
  ``src/DataBase/StateDerivatives.hh``, and
  ``src/DataBase/StateDerivatives.cc``: Per-step state and derivative
  registries assembled from the active physics packages.
* ``src/DataBase/UpdatePolicyBase.hh`` and related policy classes: Rules used
  by ``State::update`` to evolve, replace, or recompute registered fields.
* ``src/Boundary/Boundary.hh`` and ``src/Distributed/DistributedBoundary.hh``:
  Ghost-node and violation-node boundary interfaces. Distributed boundaries
  implement MPI ghost communication through the same boundary abstraction.
* ``src/Neighbor/Neighbor.hh`` and ``src/Neighbor/ConnectivityMap.hh``:
  Neighbor-search and pair-connectivity abstractions used by physics packages.

Runtime Layers
--------------

The lifecycle is layered as follows:

::

   User problem script
     |
     v
   Python imports and helper factories
     |
     v
   SpheralController
     |
     v
   C++ Integrator base class
     |
     v
   Concrete C++ integrator
     |
     v
   Physics packages
     |
     v
   DataBase, State, StateDerivatives, Boundary, ConnectivityMap, Fields

Each layer has a different responsibility.

Problem Script
  Constructs the physical problem: materials, node lists, initial fields,
  kernels, viscosity, physics packages, boundaries, database, integrator, and
  controller.

``SpheralController``
  Controls runs. It initializes package dependencies, inserts support packages
  and boundaries, handles restarts and visualization, calls the integrator, and
  schedules periodic work.

``Integrator``
  Owns one simulation step at the C++ level. It assembles transient
  ``State``/``StateDerivatives`` objects, updates ghost nodes, builds
  connectivity, applies boundaries, picks a time step, calls package hooks, and
  delegates time-centering to a concrete integrator.

Concrete Integrator
  Implements a numerical time-integration scheme such as RK2, predictor
  corrector, forward Euler, Verlet, or implicit methods. It decides how many
  derivative evaluations happen and how state updates are staged.

``Physics``
  Owns equations and package-specific durable fields. A package registers the
  state it needs, computes derivative fields, votes on a stable ``dt``, and
  maintains its own boundary/finalization behavior.

``DataBase`` and Lower Objects
  Own and expose the simulation data layout: node lists, fields, field lists,
  neighbor objects, and connectivity.

What Is Durable State?
----------------------

The durable simulation state is not primarily stored in ``State``. Durable
state lives in:

* ``NodeList`` and derived node-list fields, such as mass, position, velocity,
  ``H``, density, thermal energy, stress, damage, and DEM fields.
* physics-package member fields, such as hydro scratch fields that must persist
  across calls or be restartable;
* integrator bookkeeping, such as current time, current cycle, last ``dt``, and
  the retry multiplier;
* controller bookkeeping, such as total controller steps and restart/viz
  schedules.

``State`` and ``StateDerivatives`` are per-step registries over durable objects.
They provide a uniform interface to integrators and policies. Most fields inside
them are references to existing fields. Copying a ``State`` object copies
references unless ``copyState`` is explicitly called.

Startup Lifecycle
-----------------

The setup path starts in the user script and then moves through
``SpheralController.__init__`` and ``SpheralController.reinitializeProblem``.
The exact script varies, but the normal sequence is:

1. Import a Spheral dimension module or ``Spheral``.
2. Create equations of state, strength models, and other material objects.
3. Create ``NodeList`` objects and initialize node fields.
4. Append node lists to a ``DataBase``.
5. Create kernels and artificial viscosity/conduction objects.
6. Create physics packages.
7. Attach user-specified boundaries to physics packages.
8. Create an ``Integrator`` around the database and package list.
9. Create a ``SpheralController`` around the integrator.

During controller construction, the controller performs several normalization
steps before normal advancement begins.

Dimension and Kernel Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The controller determines the problem dimension from ``integrator.dataBase.nDim``.
It also finds the base interpolation kernel. The kernel may be passed directly,
extracted from a physics package that exposes ``kernel``, or synthesized for DEM
cases that only need an extent.

This matters because later startup operations, such as H iteration, Voronoi
support, reproducing kernels, and redistribution, need a kernel extent and a
dimension-specific type name.

Boundary Insertion
~~~~~~~~~~~~~~~~~~

The controller can insert boundaries that the user did not explicitly attach:

* spherical-origin boundaries for spherical coordinates;
* RZ axis boundaries for cylindrical coordinates;
* a distributed boundary when running with more than one MPI rank.

The distributed boundary is inserted into every physics package's boundary list.
The controller also enforces boundary ordering. Some boundaries, such as
constant and inflow/outflow boundaries, must precede the distributed boundary
because their ghost values must be communicated to other ranks.

Physics Package Organization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SpheralController.organizePhysicsPackages`` can insert support packages ahead
of packages that request shared derived structures:

* ``VoronoiCells`` is inserted before the first package that requires Voronoi
  cells, or when Voronoi computation is forced.
* ``RKCorrections`` is inserted before the first package that requests
  reproducing kernels.

This keeps physics packages declarative. A hydro package says what shared
structures it requires through ``Physics`` requirement hooks; the controller
arranges the package list so those structures are computed.

Initial Ghosts, Connectivity, and Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``reinitializeProblem`` then prepares the initial C++ state:

1. Call ``setGlobalFlags``.
2. Reset the controller step count and timers.
3. Set the integrator current time.
4. Call ``initializeProblemStartup(False)`` on all unique boundaries.
5. Reinitialize node-list neighbor structures.
6. Call ``integrator.setGhostNodes()``.
7. Build an initial connectivity map with ``db.updateConnectivityMap(False)``.
8. Call ``initializeProblemStartup`` on each physics package.
9. Build initial ``State`` and ``StateDerivatives`` objects.
10. If packages require connectivity, rebuild neighbors/connectivity with the
    requested ghost, overlap, and intersection options, and enroll the map in
    ``State``.
11. Call ``initializeProblemStartupDependencies`` on each physics package.
12. Optionally check and iterate initial ``H``.
13. Rebuild ghost nodes and connectivity again after ``H`` changes.
14. Optionally initialize derivatives for stateful boundaries or user requests.
15. Give stateful boundaries a final startup pass.
16. Register periodic work, restart, visualization, redistribution, and
    conservation bookkeeping.

The repeated neighbor/ghost/connectivity work is intentional. Startup can
change ``H``, add support packages, and initialize dependent fields. Those
changes alter neighbor extents, ghost-node positions, and connectivity.

Advance Loop
------------

After setup, normal run control is in ``SpheralController.advance``:

::

   while current_time < goal_time and step budget not exceeded:
       start step timer
       integrator.step(goal_time)
       optionally iterate H between cycles
       stop step timer
       increment controller step counters
       run periodic work

   force final periodic work, except redistribution
   print global node statistics
   print timing diagnostics

``SpheralController.step(steps)`` is a convenience wrapper that calls
``advance(1e40, steps)``.

Periodic work is not part of the C++ step engine. It is controller-level work
run after a completed step or at forced startup/end points. Default periodic
work includes:

* cycle status printing;
* garbage collection;
* conservation-history updates;
* restart dumps;
* neighbor reinitialization;
* domain redistribution;
* visualization dumps when configured;
* user-provided periodic callbacks.

C++ Integrator Entry Point
--------------------------

The common C++ entry point is ``Integrator<Dimension>::step(maxTime)``. This is
the method Python usually calls through bindings.

Its job is not to perform a specific RK or predictor-corrector update. Its job
is to prepare the common state needed by any concrete integrator, then call the
derived implementation.

The base entry point does the following:

1. Reset package-wide connectivity requirement flags.
2. Query every physics package for:

   * ``requireConnectivity``;
   * ``requireGhostConnectivity``;
   * ``requireOverlapConnectivity``;
   * ``requireIntersectionConnectivity``.

3. If the current cycle matches ``updateBoundaryFrequency``, call
   ``setGhostNodes``.
4. Build a fresh ``State`` object from the active physics packages.
5. Build a fresh ``StateDerivatives`` object from the active physics packages.
6. Apply ghost boundaries to the registered state and derivatives.
7. Try the concrete integrator step up to ten times.
8. If a concrete step reports failure, cut the internal ``dt`` multiplier in
   half and retry.
9. Disable interim ``dt`` checking on the final retry.
10. Reset the ``dt`` multiplier and return success/failure.

In pseudo-code:

::

   bool Integrator::step(maxTime):
       gather connectivity requirements from packages

       if currentCycle % updateBoundaryFrequency == 0:
           setGhostNodes()

       State state(dataBase, packages)
       StateDerivatives derivs(dataBase, packages)

       applyGhostBoundaries(state, derivs)

       for retry in 1..10:
           success = derived.step(maxTime, state, derivs)
           if success:
               break
           dtMultiplier *= 0.5

       dtMultiplier = 1.0
       return success

This method is where package requirements first influence the shared data
structures for the step. It is also where transient state containers are built.

Objects Used During This Step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main Spheral objects used during the base C++ step are:

``Integrator<Dimension>``
  The orchestration object. It stores the active package order, time
  bookkeeping, retry multiplier, boundary cadence, and connectivity requirement
  flags.

``DataBase<Dimension>``
  The shared collection of node lists and the current ``ConnectivityMap``. It
  also supplies database-wide field lists such as the per-node work field.

``Physics<Dimension>`` packages
  The active equation and model packages. The step queries them for
  connectivity requirements and ``dt`` votes, asks them to register state and
  derivative fields, and calls their initialization, derivative, boundary, and
  finalization hooks.

``Boundary<Dimension>`` objects
  The unique boundary conditions gathered from the physics packages. They
  define ghost nodes, update ghost-node positions and ``H`` tensors, finalize
  ghost communication, and identify violation nodes for enforcement.

``NodeList`` objects
  The durable node storage, including internal and ghost nodes. Fluid and DEM
  node lists own the per-node fields that ``State`` usually references.

``Neighbor<Dimension>`` and ``ConnectivityMap<Dimension>``
  The neighbor-search and pair-connectivity structures. They are refreshed as
  ghost nodes change and are rebuilt when package requirements ask for
  connectivity.

``Field`` and ``FieldList`` objects
  Typed per-node data containers registered into ``State`` and
  ``StateDerivatives``. They are the field-level views used by packages and
  update policies.

``State<Dimension>``
  A per-step registry over mutable state and update policies. Most entries are
  references to durable ``NodeList`` or package-owned fields rather than owned
  copies.

``StateDerivatives<Dimension>``
  A per-step registry for derivative, scratch, and pairwise data. Concrete
  integrators zero it before derivative evaluations and physics packages fill
  it.

``UpdatePolicyBase<Dimension>`` and concrete policies
  The field-update rules applied by ``State::update`` inside concrete
  integrators. They decide whether fields are incremented, replaced, bounded,
  or recomputed from dependencies.

Concrete ``Integrator<Dimension>`` subclasses
  The numerical time-centering implementation, such as ``SynchronousRK2``,
  ``PredictorCorrector``, ``ForwardEuler``, ``Verlet``, or an implicit method.
  The base step prepares common objects, then delegates to one of these
  subclasses.

State Construction
------------------

Constructing ``State`` and ``StateDerivatives`` is an active operation.

``State(dataBase, packageBegin, packageEnd)``
  Calls ``registerState(dataBase, state)`` on each physics package, then enrolls
  the database's current connectivity map if one exists.

``StateDerivatives(dataBase, packageBegin, packageEnd)``
  Calls ``registerDerivatives(dataBase, derivs)`` on each physics package, then
  enrolls the database's current connectivity map if one exists.

This gives packages a chance to register:

* existing node-list fields;
* package-owned field lists;
* scratch or derived state;
* update policies for state fields;
* derivative fields that will be zeroed and filled by derivative evaluation;
* shared objects such as connectivity, mesh, or RK data.

The constructors are called every step, so registration must be deterministic
and relatively cheap. Long-lived storage should be owned by node lists or
physics packages, not allocated as durable state inside ``State`` construction.

The Hook Order Inside a Concrete Step
-------------------------------------

Concrete integrators call a common set of base-class helper hooks. The exact
order depends on the integration method, but the hook meanings are consistent.

``preStepInitialize(state, derivs)``
  Called once at the beginning of a concrete step. Delegates to
  ``Physics::preStepInitialize`` for every package.

``selectDt(dtMin, dtMax, state, derivs)``
  Asks each package for a ``dt`` vote and chooses the smallest positive value,
  subject to ``dtGrowth``, ``dtMin``, ``dtMax``, and the retry multiplier.
  Optionally performs a global minimum reduction.

``initializeDerivatives(t, dt, state, derivs)``
  Clears the database ``work`` field and calls ``Physics::initialize`` on every
  package. If any package returns ``true``, ghost boundaries are reapplied.

``derivs.Zero()``
  Clears all registered derivative storage. This uses a type visitor over the
  ``StateDerivatives`` registry and handles fields, primitive values, and
  pairwise fields.

``evaluateDerivatives(t, dt, db, state, derivs)``
  Calls ``Physics::evaluateDerivatives`` on every package. This is where most
  package-specific numerical kernels run.

``finalizeDerivatives(t, dt, db, state, derivs)``
  Calls ``Physics::finalizeDerivatives`` on every package. Packages use this to
  finish reductions, apply compatible-energy adjustments, or complete
  derivative fields after all package derivative contributions are available.

``state.update(derivs, multiplier, t, dt)``
  Advances registered state according to update policies.

``applyGhostBoundaries(state, derivs)``
  Updates ghost node positions and ``H`` values, updates neighbor data, asks
  each package to apply ghost values for package-registered fields, then
  finalizes all boundaries.

``postStateUpdate(t, dt, db, state, derivs)``
  Calls ``Physics::postStateUpdate`` on every package after a state update. If
  any package returns ``true``, ghost boundaries are reapplied.

``postStepFinalize(t, dt, state, derivs)``
  Calls ``Physics::finalize`` on every package once at the end of a completed
  step.

``enforceBoundaries(state, derivs)``
  Finds violation nodes, updates them through boundaries, then asks packages to
  enforce boundary constraints on their fields.

Representative Concrete Step: SynchronousRK2
--------------------------------------------

``SynchronousRK2`` is a useful representative because it shows the common step
concepts clearly. Its implementation does:

1. Read current time ``t`` and database reference.
2. Call ``preStepInitialize``.
3. Select ``dt`` and compute ``hdt = 0.5 * dt``.
4. Initialize derivatives at ``t`` over the half-step interval.
5. Zero derivative storage.
6. Evaluate and finalize derivatives at the beginning-of-step state.
7. Copy the beginning state into ``state0`` and call ``state0.copyState``.
8. Trial-advance ``state`` by ``hdt``.
9. Set current time to ``t + hdt``.
10. Apply and finalize ghost boundaries.
11. Call ``postStateUpdate`` for the midpoint trial state.
12. Initialize, zero, evaluate, and finalize derivatives at the midpoint.
13. Optionally reselect ``dt`` and reject the step if the new vote is too small.
14. Restore durable fields from ``state0``.
15. Advance from the beginning state by the full ``dt`` using midpoint
    derivatives.
16. Set current time to ``t + dt``.
17. Apply and finalize ghost boundaries.
18. Call ``postStateUpdate`` for the completed state.
19. Call ``postStepFinalize``.
20. Enforce boundaries.
21. Increment the integrator cycle and store ``lastDt``.

In pseudo-code:

::

   bool SynchronousRK2::step(maxTime, state, derivs):
       t = currentTime()
       db = accessDataBase()

       preStepInitialize(state, derivs)
       dt = selectDt(clamped dtMin, clamped dtMax, state, derivs)
       hdt = 0.5 * dt

       initializeDerivatives(t, hdt, state, derivs)
       derivs.Zero()
       evaluateDerivatives(t, hdt, db, state, derivs)
       finalizeDerivatives(t, hdt, db, state, derivs)

       state0 = State(state)
       state0.copyState()

       state.update(derivs, hdt, t, hdt)
       currentTime(t + hdt)
       applyGhostBoundaries(state, derivs)
       finalizeGhostBoundaries()
       postStateUpdate(t + hdt, hdt, db, state, derivs)

       initializeDerivatives(t + hdt, hdt, state, derivs)
       derivs.Zero()
       evaluateDerivatives(t + hdt, hdt, db, state, derivs)
       finalizeDerivatives(t + hdt, hdt, db, state, derivs)

       if allowDtCheck and selectDt(...) < dtCheckFrac * dt:
           currentTime(t)
           state.assign(state0)
           return false

       state.assign(state0)
       state.update(derivs, dt, t, dt)
       currentTime(t + dt)
       applyGhostBoundaries(state, derivs)
       finalizeGhostBoundaries()
       postStateUpdate(t + dt, dt, db, state, derivs)
       postStepFinalize(t + dt, dt, state, derivs)
       enforceBoundaries(state, derivs)
       currentCycle(currentCycle + 1)
       lastDt(dt)
       return true

This is not the only concrete step sequence. ``PredictorCorrector`` evaluates
beginning and end derivatives and applies a corrector using half-step
contributions from both derivative sets. Implicit integrators override more of
the base behavior to solve for independent state convergence. The important
constant is the contract: concrete integrators orchestrate the same
``State``/``StateDerivatives``/``Physics``/``Boundary`` machinery.

State Update Policy Mechanics
-----------------------------

``State::update`` is where durable fields actually change. The integrator does
not directly know how to update most fields. It passes the derivative registry,
time multiplier, current time, and ``dt`` to ``State``. ``State`` then applies
registered policies.

The update algorithm:

1. Build a map of policy-controlled state fields.
2. Track which field names remain to be updated.
3. Iterate until all fields are updated.
4. For each policy, check whether its declared dependencies have already been
   satisfied.
5. Ensure field-list-wide wildcard policies fire before node-list-specific
   policies for the same field name.
6. Call either ``policy->update`` or ``policy->updateAsIncrement`` depending on
   ``timeAdvanceOnly``.
7. Detect circular dependencies if an iteration completes without updating any
   remaining field.

This design decouples time integration from field-specific semantics. For
example, a hydro package can register density with an increment policy, pressure
with a replacement policy depending on density and thermal energy, and sound
speed with a package-specific derived-state policy. The integrator just calls
``state.update``.

Boundary Lifecycle
------------------

Boundaries participate in two related but distinct operations:

``setGhostNodes``
  Defines the ghost-node population and usually changes node-list sizes.

``applyGhostBoundaries``
  Updates ghost-node positions and field values after state changes.

``enforceBoundaries``
  Finds internal nodes that violate boundary constraints and fixes their
  internal state.

The base integrator's ``setGhostNodes`` does:

1. Collect unique boundaries across all physics packages.
2. Remove old ghost nodes from fluid and DEM node lists.
3. Temporarily enlarge neighbor kernel extents if overlap connectivity is
   required.
4. Update neighbor structures.
5. Ask every boundary to set all ghost nodes on the database.
6. Finalize each boundary and refresh neighbor structures.
7. Restore neighbor kernel extents if they were enlarged.
8. If connectivity is required, build the database connectivity map.
9. Optionally cull unused ghost nodes when domain-decomposition independence,
   ghost connectivity, and overlap connectivity do not require keeping them.
10. Patch connectivity and delete culled ghost nodes from node lists.

``applyGhostBoundaries`` is called more often. It updates ghost-node positions
and ``H`` tensors, refreshes neighbor structures, then delegates field-level
ghost value updates to each physics package. Physics packages, not the base
integrator, know which package-owned fields need boundary application.

``enforceBoundaries`` is called after state updates in concrete integrators. It
asks boundaries to identify violation nodes and then delegates field-specific
enforcement to packages.

Connectivity Lifecycle
----------------------

Connectivity is built only when some package asks for it. The request is made
through ``Physics`` requirement hooks and collected by the base integrator at
the beginning of ``Integrator::step(maxTime)``.

The main flags are:

``requireConnectivity``
  Build ordinary significant-neighbor connectivity.

``requireGhostConnectivity``
  Include ghost-node connectivity.

``requireOverlapConnectivity``
  Build overlap connectivity. During ghost-node setup the integrator
  temporarily doubles kernel extent so the ghost population is sufficient.

``requireIntersectionConnectivity``
  Store intersection connectivity for node pairs.

Once ``setGhostNodes`` calls ``DataBase::updateConnectivityMap``, the
``DataBase`` owns a current ``ConnectivityMap`` pointer. New ``State`` and
``StateDerivatives`` objects enroll that connectivity pointer during
construction.

The connectivity map then becomes available to physics packages through:

* direct database access;
* enrolled state access;
* field-list and pair-list traversal helpers.

This document intentionally stops at the lifecycle level. The internal
connectivity data structure design, including master/coarse/refine neighbor
lists, node-pair construction, overlap connectivity, intersection connectivity,
and deterministic traversal, is covered in
:doc:`connectivity_data_structures`. The ``Neighbor`` class family itself is
covered in :doc:`neighbor_family_and_usage`.

Physics Package Lifecycle
-------------------------

A ``Physics`` package participates in several phases.

Startup
  ``initializeProblemStartup`` sizes and initializes durable package fields once
  node lists and materials exist.

Startup Dependencies
  ``initializeProblemStartupDependencies`` initializes derived state that needs
  an initial ``State`` and ``StateDerivatives`` object.

State Registration
  ``registerState`` enrolls fields and update policies into ``State``.

Derivative Registration
  ``registerDerivatives`` enrolls derivative fields into ``StateDerivatives``.

Step Initialization
  ``preStepInitialize`` runs once per concrete step. ``initialize`` runs before
  each derivative evaluation and can request boundary reapplication.

Derivative Evaluation
  ``evaluateDerivatives`` fills derivative fields from the current state,
  database, and connectivity.

Derivative Finalization
  ``finalizeDerivatives`` completes derivative state after all packages have had
  a chance to contribute.

Post-State Update
  ``postStateUpdate`` reacts after an integrator applies a state update and can
  request boundary reapplication.

Step Finalization
  ``finalize`` runs once after a completed step.

Boundary Handling
  ``applyGhostBoundaries`` and ``enforceBoundaries`` apply package-specific
  boundary behavior to registered fields.

Most hooks default to no-op behavior in ``Physics.cc``. Concrete packages
override only the hooks they need.

Object Roles at the Main-Loop Level
-----------------------------------

The following table summarizes the primary objects visible in the lifecycle.

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Object
     - Represents
     - Main-loop role
   * - ``SpheralController``
     - Python run-control object
     - Calls the integrator in a loop and runs restart, viz, redistribution,
       diagnostics, and user periodic work.
   * - ``Integrator``
     - C++ time-step engine
     - Builds transient state, manages boundaries/connectivity, selects ``dt``,
       retries unstable steps, and delegates time-centering.
   * - Concrete integrator
     - Numerical time-centering scheme
     - Orders derivative evaluations and state updates.
   * - ``Physics``
     - A package of equations or forces
     - Registers state, fills derivatives, votes on ``dt``, and handles
       package-specific boundary/finalization behavior.
   * - ``DataBase``
     - Active node-list collection
     - Exposes node-list and field-list views, neighbor updates, and
       connectivity maps.
   * - ``NodeList``
     - Local nodes of one material/body type
     - Owns node counts, fundamental fields, attached fields, and a neighbor
       object.
   * - ``Field``
     - Typed values over one node list
     - Stores durable or derivative values.
   * - ``FieldList``
     - Same-typed fields across node lists
     - Gives packages global views without flattening node-list ownership.
   * - ``State``
     - Per-step registry of state fields and policies
     - Applies update policies to durable fields.
   * - ``StateDerivatives``
     - Per-step registry of derivative/change fields
     - Is zeroed and filled by physics packages.
   * - ``UpdatePolicyBase``
     - Field update rule
     - Defines how ``State::update`` evolves or replaces one registered field.
   * - ``Boundary``
     - Ghost/violation-node rule
     - Creates ghost nodes, fills ghost values, and enforces constraints.
   * - ``Neighbor``
     - Per-node-list search structure
     - Provides candidate neighbor lists for connectivity construction.
   * - ``ConnectivityMap``
     - Cross-node-list pair connectivity
     - Provides neighbor and pair traversal to physics packages.

Common Developer Questions
--------------------------

Where do I add a new time integration method?
  Add a concrete ``Integrator`` subclass. Reuse base-class hooks for package
  initialization, derivative evaluation, boundaries, and finalization. Keep the
  field-specific update semantics in ``State`` policies rather than embedding
  them directly in the integrator.

Where do I add a new equation or force?
  Add a ``Physics`` package or derive from a more specific base such as
  ``GenericHydro`` or ``GenericBodyForce``. Register fields and policies,
  compute derivatives, and declare shared-structure requirements through
  ``Physics`` hooks.

Where does a persistent field belong?
  If the field is intrinsic to nodes, put it on a node-list specialization. If
  it is owned by a physics model, store it as a package member field or field
  list. Register it with ``State`` each step. Do not rely on a transient
  ``State`` object to own durable state.

Where do ghost values get filled?
  Geometric ghost-node creation is in ``Boundary`` and ``Integrator``. Field
  values are normally filled when each package's ``applyGhostBoundaries`` calls
  boundary methods on the fields it owns or registered.

Where is connectivity built?
  The integrator asks packages what connectivity they need, then
  ``setGhostNodes`` calls ``DataBase::updateConnectivityMap`` after ghost nodes
  and neighbor structures are current.

Why is ``State`` copied and then ``copyState`` called?
  A plain ``State`` copy is reference-based. Concrete integrators need a durable
  snapshot for trial states and rejection/retry logic, so they call
  ``copyState`` when they need owned copies of registered fields.

Why do packages vote on ``dt`` instead of the integrator computing it?
  Stability limits are physics-specific. The integrator can apply global
  policies such as growth limits and global reductions, but each package knows
  its own wave speeds, accelerations, strain limits, solver convergence
  criteria, or force constraints.

Why does startup build ghosts and connectivity more than once?
  Startup changes dependent state and often changes ``H``. Those changes alter
  node extents, neighbor search, and ghost-node placement. Rebuilding keeps the
  first real step consistent with the initialized state.

Design Boundaries to Preserve
-----------------------------

When extending the main loop, preserve these boundaries:

* The controller should manage run-level workflow, not field-level numerical
  updates.
* The base integrator should manage shared step mechanics, not package-specific
  equations.
* Concrete integrators should manage time-centering, not package-owned field
  semantics.
* Physics packages should own their equations, derivative fields, and
  package-specific boundary behavior.
* ``State`` policies should own field update semantics.
* ``DataBase`` should expose layout and connectivity, not implement physics.
* Boundaries should own ghost and violation-node mappings.

Related Deep Dives
------------------

This lifecycle document is the foundation for several deeper design documents:

* integrator functionality, concrete time-centering schemes, and the
  supporting ``State``/policy object model, covered in
  :doc:`integrator_and_state_update_model`;
* device-facing value/view design, managed view pointer dispatch, and the
  field execution model, covered in
  :doc:`value_view_and_device_execution_model`;
* current device-facing object families, covered in
  :doc:`value_view_conversion_case_studies`;
* RAJA loop structure, CHAI movement, and managed device views, covered in
  :doc:`raja_chai_execution_patterns`;
* connectivity data structures and traversal patterns, covered in
  :doc:`connectivity_data_structures`.
* the ``Neighbor`` class family and its usage, covered in
  :doc:`neighbor_family_and_usage`.

Those documents reference this lifecycle rather than repeat it.
