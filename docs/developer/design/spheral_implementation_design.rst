Spheral Implementation Design
=============================

Purpose
-------

This document describes the implementation structure of Spheral as it exists in
the source tree. It is intended for developers who need to understand where the
major abstractions live, how data moves through a simulation, and how new
physics, data structures, or bindings should fit into the existing design.

Spheral is primarily a C++20 library with explicit template instantiations for
one, two, and three spatial dimensions. The Python interface is the normal
application layer: users assemble materials, node lists, kernels, physics
packages, boundaries, and an integrator from Python, while the performance
critical algorithms run in compiled C++.

The implementation is organized around a small set of concepts:

* ``Dimension`` selects the scalar, vector, tensor, and faceted-volume types
  for 1-D, 2-D, or 3-D code.
* ``NodeList`` owns the nodes and the fundamental fields defined on them.
* ``Field`` stores one value per node on one ``NodeList``.
* ``FieldList`` groups same-named fields across multiple ``NodeList`` objects.
* ``Neighbor`` is the per-``NodeList`` spatial search object used to find
  candidate interactions.
* ``DataBase`` owns the active collection of node lists and produces global
  field-list views and connectivity.
* ``State`` and ``StateDerivatives`` are transient registries built from the
  active physics packages.
* ``Physics`` packages register state, compute derivatives, vote on time steps,
  and apply package-specific boundary behavior.
* ``Integrator`` implementations orchestrate time advancement by repeatedly
  building state, updating neighbors and boundaries, evaluating physics
  derivatives, and applying update policies.

This page is the architectural overview. Detailed subsystem documents live
beside it:

* :doc:`simulation_lifecycle_and_main_loop` explains the controller and C++
  integrator execution path.
* :doc:`integrator_and_state_update_model` describes concrete integrator
  mechanics, ``State``/``StateDerivatives``, and update policies.
* :doc:`value_view_and_device_execution_model` describes the value/view and
  managed view pointer patterns used when host objects become device-facing.
* :doc:`value_view_conversion_case_studies` identifies current device-facing
  object families such as ``Field``/``FieldView`` and
  ``NodePairList``/``NodePairListView``.
* :doc:`raja_chai_execution_patterns` describes RAJA loop structure, CHAI data
  movement, managed views, and current device-kernel constraints.
* :doc:`connectivity_data_structures` describes neighbor search,
  ``ConnectivityMap``, node-pair lists, overlap/intersection connectivity, and
  pairwise fields.
* :doc:`neighbor_family_and_usage` analyzes ``Neighbor``, ``TreeNeighbor``,
  ``NestedGridNeighbor``, and their role in connectivity construction.

Repository Layout
-----------------

The source tree is package-oriented. Each package generally has a
``CMakeLists.txt``, C++ headers and sources, explicit-instantiation templates
named ``*Inst.cc.py``, and optionally Python helper files.

``src/``
  The compiled Spheral implementation. Most directories map directly to C++
  packages and Python binding submodules.

``src/PYB11/``
  PYB11Generator definitions used to produce pybind11 bindings. The generated
  extension modules are assembled into ``SpheralCompiledModules``.

``src/SimulationControl/``
  Python application layer: imports, dimension aliases, the controller, command
  line helpers, restart/viz orchestration, and high-level simulation workflow.

``cmake/``
  Build options, third-party dependency handling, BLT integration, explicit
  instantiation helpers, and Spheral-specific CMake macros.

``tests/``
  C++ unit tests built through CMake plus Python unit and functional tests.

``docs/``
  Sphinx documentation.

Major Source Packages
---------------------

The following table summarizes the important implementation packages. The list
is intentionally architectural rather than exhaustive.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Package
     - Responsibility
   * - ``Geometry``
     - Dimensional traits, vector/tensor math, planes, polygons, polyhedra,
       clipping utilities, eigenvalue helpers, and faceted-volume types.
   * - ``Kernel``
     - SPH interpolation kernels and ``TableKernel`` tabulation/view support.
   * - ``NodeList``
     - Node ownership, fluid/solid/DEM node-list specializations, and the
       global ``NodeListRegistrar`` ordering service.
   * - ``Field``
     - Per-node field storage, field views, field lists, iterators, and
       memory-management hooks for host/device execution.
   * - ``Neighbor``
     - Neighbor search abstractions and ``ConnectivityMap`` construction for
       pair interactions.
   * - ``DataBase``
     - Central registry for node lists and cached global field-list views.
   * - ``DataBase/State*``
     - Transient state and derivative containers plus update policies.
   * - ``Physics``
     - Base physics interfaces, generic hydro infrastructure, and body-force
       abstractions.
   * - ``Integrator``
     - Explicit and implicit time integration algorithms.
   * - ``Boundary``
     - Boundary-condition base classes and concrete geometric/stateful
       boundaries.
   * - ``Distributed``
     - MPI process abstraction, distributed boundaries, node redistribution,
       and domain-decomposition strategies.
   * - ``SPH``, ``CRKSPH``, ``FSISPH``, ``GSPH``, ``SVPH``, ``FVCRKH``
     - Hydrodynamics implementations and supporting policies/utilities.
   * - ``ArtificialViscosity`` and ``ArtificialConduction``
     - Dissipation models used by hydro packages.
   * - ``Material`` and ``SolidMaterial``
     - Equations of state, physical units, and strength models.
   * - ``Strength``, ``Damage``, ``Porosity``
     - Solid-material state policies and physics models for plasticity,
       fracture/damage, and porous compaction.
   * - ``DEM``
     - Discrete-element particle mechanics.
   * - ``Gravity`` and ``ExternalForce``
     - N-body, tree, polyhedral gravity, and other body-force packages.
   * - ``Mesh`` and ``VoronoiCells``
     - Tessellation support, Voronoi-cell physics, and mesh output helpers.
   * - ``RK`` and ``KernelIntegrator``
     - Reproducing-kernel corrections and kernel-integration machinery.
   * - ``FileIO`` and ``DataOutput``
     - Restart registration and storage backends such as HDF5, Silo, Sidre,
       flat files, and Python-backed file I/O.
   * - ``Utilities``
     - Contracts, timers, interpolation utilities, global-node IDs, space
       filling curves, GPU helpers, and assorted numerical support code.

Build-Time Design
-----------------

Top-level configuration starts in ``CMakeLists.txt`` and delegates most project
setup to ``cmake/SetupSpheral.cmake``. That file configures the C++ standard,
BLT, compiler and accelerator options, third-party dependencies, build
definitions, documentation, tests, and installation metadata.

Package assembly happens in ``src/CMakeLists.txt``:

1. The build chooses a list of internal packages from core packages plus
   feature-controlled optional packages such as gravity, GSPH, SVPH, LEOS, or
   SUNDIALS.
2. Each package adds a ``Spheral_<Package>`` target using
   ``spheral_add_obj_library``.
3. In normal builds those package targets are object libraries. Their objects
   are folded into the exported ``Spheral_CXX`` library.
4. In developer builds, ``ENABLE_DEV_BUILD`` makes package targets separate
   shared libraries and makes ``Spheral_CXX`` an interface target over them.
   This reduces rebuild cost while preserving the public target name.
5. If Python is enabled, ``src/PYB11`` builds generated pybind11 submodules and
   optionally the aggregate ``SpheralCompiledModules`` module.

Dimensional explicit instantiation is a central build pattern. Packages list
template bases such as ``SPHBase`` or ``DataBase`` in ``<Package>_inst``. The
``instantiate`` CMake helper runs ``src/helpers/InstantiationGenerator.py`` on
the corresponding ``*Inst.cc.py`` file and produces one combined
``*Inst.cc`` file, or one file per dimension if
``SPHERAL_COMBINE_INSTANTIATIONS`` is disabled. The enabled dimensions are
controlled by ``SPHERAL_ENABLE_1D``, ``SPHERAL_ENABLE_2D``, and
``SPHERAL_ENABLE_3D``.

Core Data Model
---------------

The core data model is deliberately particle-centric. Nodes are the independent
degrees of freedom, but a node is not represented by a standalone C++ ``Node``
object in this layer. A node is a node id, or row index, into the aligned
``Field`` objects attached to one ``NodeList``. Fields attach typed values to
those node ids. Field lists present a global view over several node lists
without forcing a single monolithic array.

The most important ownership relationships are:

::

   DataBase
     owns/references active NodeLists
     owns/caches the active ConnectivityMap
       |
       +-- NodeList
       |     owns fundamental Fields: mass, position, velocity, H, work
       |     registers all Fields defined on that NodeList
       |     has one registered Neighbor
       |       |
       |       +-- TreeNeighbor or another concrete Neighbor implementation
       |
       +-- FluidNodeList
       |     adds mass density, specific thermal energy, and an EOS
       |
       +-- SolidNodeList
       |     adds deviatoric stress, plastic strain, damage, fragments,
       |     particle types, and a strength model
       |
       +-- DEMNodeList
             adds discrete-element particle state

   Field
     stores one value per node on one NodeList

   FieldList
     stores references to, or copies of, Fields across NodeLists

   Neighbor
     stores per-NodeList search acceleration data and returns candidate node ids

   ConnectivityMap
     combines all active Neighbors into final cross-NodeList connectivity

   State / StateDerivatives
     stores registered Fields, FieldLists, connectivity, meshes, kernels,
     and other per-step objects by string key

``NodeList``
~~~~~~~~~~~~

``NodeList<Dimension>`` is the base container for simulation nodes. It stores:

* node counts split into internal and ghost nodes;
* the required state fields ``mass``, ``positions``, ``velocity``, and
  ``Hfield``;
* a mutable ``work`` field used for load-balancing and diagnostics;
* smoothing-scale bounds and neighbor-count targets;
* registered ``FieldBase`` references for all fields defined on the node list;
* a pointer to the active ``Neighbor`` implementation;
* restart registration.

Code that works with a node normally carries two pieces of identity: the
``NodeList`` and the node id within that list. For node id ``i``, the node's
state is assembled from field entries such as ``mass(i)``, ``positions(i)``,
``velocity(i)``, and ``Hfield(i)``. Node iterators package this same pair of
``NodeList`` plus node id so fields and field lists can be indexed without
constructing node objects.

The internal/ghost split is a major invariant. Internal nodes are locally owned
degrees of freedom. Ghost nodes are generated by boundaries or received from
other MPI ranks and are regenerated or culled as the simulation evolves.

``FluidNodeList`` adds thermodynamic fields and an ``EquationOfState`` pointer.
``SolidNodeList`` adds solid mechanics fields and a ``StrengthModel``. These
specializations keep material state close to the nodes while still allowing the
same ``DataBase`` and ``FieldList`` machinery to operate on all node lists.

``NodeListRegistrar`` is a singleton that keeps node lists in name-sorted order.
That ordering is used when building field lists and when deterministic traversal
is requested. The integrator exposes this through
``domainDecompositionIndependent``.

``Field`` and ``FieldList``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Field<Dimension, DataType>`` owns a typed array of values for one
``NodeList``. It derives from ``FieldBase`` and ``FieldView``:

* ``FieldBase`` supplies the type-erased interface needed by boundaries,
  restarts, packing/unpacking, and state registration.
* ``FieldView`` supplies lightweight element access and view-style operations.
* ``Field`` owns the backing storage and keeps the view span synchronized.

Fields know their node list and register with it. This lets the node list resize
or reorder all attached fields when ghost nodes are created, nodes are deleted,
or domains are redistributed.

``FieldList<Dimension, DataType>`` is the global-field abstraction. A field list
can either reference existing fields or own copies. It supports element-wise
arithmetic, reductions, min/max limiting, interpolation through a kernel, and
node iterators spanning multiple node lists. ``DataBase`` factory and query
methods are the common way to obtain field lists for standard quantities.

Memory and Device Views
~~~~~~~~~~~~~~~~~~~~~~~

The field implementation is structured to support both CPU and accelerator
execution:

* ``Field`` storage uses a standard allocator by default, and a UVM allocator
  when unified memory support is enabled.
* Views can be moved through CHAI execution spaces where supported.
* GPU-capable kernels use RAJA policies, ``SPHERAL_HOST_DEVICE`` methods, CHAI
  managed pointers, and explicit view movement.
* RAJA-specific SPH implementations, such as ``SPH_RAJA`` and
  ``SolidSPH_RAJA``, keep pair loops and per-node finalization loops in
  accelerator kernels.

The device design keeps host-side polymorphic objects mostly out of inner GPU
loops. Where virtual behavior is needed on device, such as artificial viscosity
views, the code uses CHAI-managed view objects and constructs device-valid
instances carefully.

DataBase
~~~~~~~~

``DataBase<Dimension>`` is the central runtime registry for node lists. It keeps
separate views of all node lists, fluid node lists, solid node lists, and DEM
node lists. It also owns or caches the active ``ConnectivityMap``.

Its primary responsibilities are:

* append and delete node lists;
* return counts of local and global internal/ghost nodes;
* expose node iterators over all nodes or only fluid nodes;
* reinitialize all neighbor objects;
* build and patch connectivity maps;
* construct field lists for standard global, fluid, solid, and DEM fields;
* create new field lists matching the active node-list layout.

The database does not own every physics-specific field. Physics packages keep
their own scratch and derived field lists. The database provides the common
layout and access layer that lets those packages operate over the same node
population.

State and Update Policies
-------------------------

``StateBase`` is a key-value registry used by both ``State`` and
``StateDerivatives``. Keys are strings. Field keys encode both the field name
and node-list name so same-named fields from different node lists can coexist.
The registry can hold:

* fields and field lists;
* arbitrary C++ objects through ``std::any``;
* a shared ``ConnectivityMap``;
* an optional mesh;
* reproducing-kernel and other per-step support objects.

``State`` represents the values to be advanced. ``StateDerivatives`` represents
changes, rates, or replacement values computed by physics packages.

``State`` adds an update-policy map. A policy declares dependencies and defines
how a registered field is advanced when ``State::update`` is called. Common
policies include:

* ``IncrementState``: ``x1 = x0 + multiplier * derivative``.
* ``ReplaceState``: replace state from derivative data, but support increment
  behavior when forced by an integrator.
* ``PureReplaceState``: direct replacement without derivative semantics.
* bounded variants: apply floors/ceilings while updating.
* package-specific policies: pressure, sound speed, stress, damage, porosity,
  H evolution, volume, and other derived-state updates.

This policy system is a key extension point. A physics package does not need to
teach every integrator how its state evolves. Instead, it registers fields and
policies; integrators call the common ``State::update`` API.

Runtime Simulation Lifecycle
----------------------------

Most runs are built from Python but execute the same C++ lifecycle. At a high
level:

::

   Python problem script
     creates EOS/strength/material objects
     creates NodeLists and fills Fields
     appends NodeLists to DataBase
     creates kernels, viscosity, physics packages, boundaries
     creates Integrator
     creates SpheralController

   SpheralController startup
     initializes global flags and Axom
     inserts coordinate/distributed boundaries as needed
     organizes dependent physics packages
     creates ghost nodes and initial connectivity
     calls initializeProblemStartup on each Physics package
     builds State and StateDerivatives
     calls initializeProblemStartupDependencies
     optionally iterates initial H
     restores restart data if requested

   Per-cycle Integrator::step
     computes package connectivity requirements
     sets ghost nodes and updates connectivity
     builds fresh State and StateDerivatives
     applies ghost boundaries
     selects dt from physics package votes
     evaluates derivatives through package hooks
     advances State through update policies
     enforces and reapplies boundaries
     lets packages finalize state and derivatives
     updates cycle/time bookkeeping

State and derivative containers are normally rebuilt each step. The durable
state lives in node lists and physics package fields. The transient containers
provide a consistent, per-step view for integrators and update policies.

Integrator Design
-----------------

``Integrator<Dimension>`` is the abstract base for time advancement. It stores a
reference to the ``DataBase``, an ordered vector of ``Physics`` packages,
current time/cycle, time-step limits, growth limits, boundary update controls,
and restart state.

The base class provides shared orchestration:

* append physics packages, including pre- and post-subpackages;
* choose a time step by taking the minimum positive package vote;
* reduce the selected time step globally when configured;
* call ``preStepInitialize``, ``initialize``, ``evaluateDerivatives``,
  ``finalizeDerivatives``, ``postStateUpdate``, and ``finalize`` on packages;
* collect unique boundaries from all packages;
* set ghost nodes, apply ghost values, finalize boundaries, and enforce
  violation-node corrections;
* update or patch connectivity maps;
* retry a step with a reduced ``dt`` when an integrator reports instability.

Concrete integrators implement the actual time-centering scheme in
``step(maxTime, state, derivs)``. The tree includes explicit integrators such
as ``ForwardEuler``, ``PredictorCorrector``, ``SynchronousRK2``,
``CheapSynchronousRK2``, ``SynchronousRK4``, and ``Verlet``. It also includes
implicit infrastructure through ``ImplicitIntegrator``, ``BackwardEuler``, and
``CrankNicolson``.

For example, ``SynchronousRK2`` performs:

1. beginning-of-step package initialization;
2. time-step selection;
3. derivative evaluation at the current state;
4. a trial half-step update;
5. boundary application and package post-state updates;
6. derivative evaluation at the midpoint;
7. optional stability recheck;
8. full-step update from the original state using midpoint derivatives;
9. final package hooks, boundary enforcement, and cycle/time updates.

``PredictorCorrector`` follows the same base hooks but uses beginning and end
derivatives to apply its corrector stage.

Physics Package Contract
------------------------

``Physics<Dimension>`` is the interface all physics packages implement. Required
methods are:

``evaluateDerivatives``
  Compute package contributions to ``StateDerivatives`` for the current state.

``dt``
  Vote on a stable time step and return a reason string.

``registerState``
  Enroll state fields and update policies.

``registerDerivatives``
  Enroll derivative fields that the package will compute.

``label``
  Provide a restart and diagnostic label.

The base class also defines hooks for startup, dependency initialization,
per-step initialization, derivative finalization, post-update work, boundary
application, and visualization-state registration. It exposes requirement flags
that tell the integrator whether the package needs:

* node connectivity;
* ghost connectivity;
* overlap connectivity;
* intersection connectivity;
* Voronoi cells;
* reproducing-kernel corrections;
* reproducing-kernel Hessians.

The integrator combines those requirements across all packages before it creates
ghost nodes and connectivity. This lets packages request expensive shared
structures without coupling the integrator to any specific physics type.

Hydrodynamics Architecture
--------------------------

Hydro packages share a common base through ``GenericHydro``. That class stores
the artificial viscosity object, CFL settings, time-step diagnostics, and
neighbor statistics. It implements hydro-style ``dt`` voting and implicit
residual support.

``GenericHydro`` introduces ``MassDensityType`` to select density treatment:

* summed density;
* rigorous or hybrid summed density;
* integrated continuity density;
* Voronoi-cell density;
* summed Voronoi-cell density;
* corrected summed density.

SPH
~~~

``SPHBase`` is the common base for standard SPH/ASPH packages. It stores the
interpolation kernels, hydro switches, tensile correction parameters, optional
Voronoi bounding boxes, and a large set of scratch and derivative field lists:
pressure, sound speed, omega-grad-h, entropy, density sums, velocity-gradient
fields, acceleration, density derivative, energy derivative, XSPH fields, and
related normalization/moment fields.

Concrete packages such as ``SPH``, ``PSPH``, ``SolidSPH``, spherical/RZ
variants, and RAJA variants provide the actual derivative loops. The base class
handles state registration, derivative registration, startup dependencies,
boundary application, derivative finalization, and restart persistence.

CRKSPH
~~~~~~

``CRKSPHBase`` implements common Conservative Reproducing Kernel SPH behavior.
It requests reproducing kernels, tracks the correction order, stores hydro
state/derivative field lists, and supports compatible energy evolution, total
energy evolution, XSPH, and density-update variants. Concrete CRKSPH and solid
CRKSPH packages specialize the derivative evaluation.

FSISPH
~~~~~~

``SolidFSISPH`` implements fluid-solid interaction SPH. It derives from
``GenericHydro`` rather than ``SPHBase`` and manages interface-aware fields and
slide-surface behavior. Supporting functions compute FSI-specific density
estimates, including H-weighted and interface-pressure-corrected variants.

GSPH, MFM, and MFV
~~~~~~~~~~~~~~~~~~

The GSPH package family uses ``GenericRiemannHydro`` for Godunov-style,
Riemann-solver-based methods. Concrete ``GSPH``, ``MFM``, and ``MFV`` packages
share gradient initialization, volume estimates, and Riemann hydro state while
implementing their own derivative algorithms.

SVPH and FVCRKH
~~~~~~~~~~~~~~~

SVPH packages use mesh/Voronoi face information and correction policies for
finite-volume-like hydro behavior. FVCRKH adds finite-volume conservative
reproducing-kernel hydro support.

Artificial Viscosity and Conduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Artificial viscosity is modeled as a separate package hierarchy. Hydro packages
hold an ``ArtificialViscosity`` reference and request scalar or tensor viscosity
views for inner loops. CPU implementations can use virtual C++ calls directly;
GPU implementations use device-valid view objects.

Artificial conduction is optional and is included only when
``SPHERAL_ENABLE_ARTIFICIAL_CONDUCTION`` is enabled.

Materials, Strength, Damage, Porosity, and DEM
----------------------------------------------

The material model is split between fluid thermodynamics, solid EOS support,
strength, damage, and porosity.

``Material``
  Defines ``EquationOfState`` and fluid EOS implementations such as gamma-law,
  isothermal, polytropic, stiffened gas, and Helmholtz-backed EOS support.
  EOS classes set pressure, pressure derivatives, temperature, sound speed,
  gamma, bulk modulus, entropy, and inverse pressure-energy estimates on fields.

``SolidMaterial``
  Adds solid equations of state and strength models, including Gruneisen,
  Tillotson, Murnaghan, Osborne, ANEOS, constant strength, Johnson-Cook,
  Steinberg-Guinan, Geodyn, and related material libraries/shadow Python
  helpers.

``Strength``
  Provides policies for derived solid fields such as bulk modulus, shear
  modulus, yield strength, deviatoric stress, plastic strain, and melt energy.

``Damage``
  Provides damage models, flaw distributions, damage policies, node-coupling
  modifiers, and fragment identification helpers. Damage fields are normal
  state fields carried by solid node lists and updated through policies and
  physics hooks.

``Porosity``
  Provides porosity models such as P-alpha and strain porosity, plus policies
  that connect porous solid density and alpha fields to the rest of the state.

``DEM``
  Provides discrete-element mechanics, pair-field storage policies, linear
  spring DEM, contact storage, and unique particle identity utilities.

Gravity and Body Forces
-----------------------

Body-force packages use the same ``Physics`` contract as hydro packages.
``TreeGravity`` is an instructive example. It derives from
``GenericBodyForce``, opts out of regular hydrodynamic connectivity, registers
potential and acceleration state, builds an internal cell tree during
initialization/evaluation, computes gravitational acceleration and potential,
votes on a time step, and reports extra potential energy.

Other gravity support includes direct N-body gravity, compatible gravitational
velocity policies, polyhedral gravity, approximate polyhedral models, and
Laplacian matrix support.

Neighbors and Connectivity
--------------------------

Neighbor search is separated into two layers:

``Neighbor``
  A per-node-list object that knows how to find candidate neighbors using a
  search algorithm and kernel extent.

``ConnectivityMap``
  A cross-node-list structure built from all active neighbor objects. It stores
  significant neighbors, node-pair lists, optional overlap connectivity, and
  optional intersection connectivity.

The typical path is:

1. Boundaries set or update ghost nodes.
2. Each node list's ``Neighbor`` updates its internal search structure.
3. ``DataBase::updateConnectivityMap`` asks all relevant neighbor objects for
   master/coarse/refined lists.
4. ``ConnectivityMap`` stores the resulting neighbor lists and node pairs.
5. Physics packages traverse connectivity to compute pair interactions.

The connectivity map also provides deterministic traversal helpers. When
domain-decomposition-independent mode is enabled, packages can traverse nodes in
the connectivity map's stored order instead of raw local storage order.

Boundary Conditions
-------------------

Boundary conditions define three sets of nodes per node list:

``controlNodes``
  Real nodes used as the source for ghost values.

``ghostNodes``
  Ghost nodes created or updated by the boundary.

``violationNodes``
  Internal nodes that violate the boundary and need correction.

``Boundary`` requires concrete classes to set and update ghost nodes, identify
violation nodes, and update violation nodes. It also provides type-erased and
typed ``applyGhostBoundary`` and ``enforceBoundary`` overloads for fields and
field lists.

The integrator owns the boundary workflow:

1. collect unique boundaries from all physics packages;
2. clear old ghost nodes;
3. ask each boundary to create ghost nodes;
4. update neighbor structures;
5. build connectivity if required;
6. apply ghost values for all package-registered fields;
7. finalize boundaries;
8. later, enforce violation-node corrections after state updates.

Distributed Boundaries and Redistribution
-----------------------------------------

MPI support is implemented as boundary and redistribution machinery rather than
as a separate data model. ``DistributedBoundary`` derives from ``Boundary`` and
connects node lists across ranks by maintaining, for each node list and peer
domain, the nodes to send and the ghost nodes to receive.

Important responsibilities include:

* building receive and ghost nodes from send-node maps;
* packing fixed-size and variable-size field values;
* launching nonblocking exchanges;
* unpacking received field data into ghost fields;
* culling inactive ghost nodes after connectivity has been built;
* opting distributed ghosts out of mesh generation when appropriate.

``RedistributeNodes`` and derived classes repartition internal nodes across MPI
ranks. The design works by computing a desired vector of ``DomainNode`` records
and then enforcing that decomposition by deleting, packing, sending, receiving,
and appending node field values. Implementations include position-based,
sort-and-divide, nested-grid, Morton, Peano-Hilbert, ParMETIS, and Voronoi
strategies.

Mesh, Voronoi, RK, and Kernel Integration
-----------------------------------------

Mesh and correction packages provide shared geometric infrastructure requested
by physics packages through ``Physics`` requirement hooks or state registration.

``Mesh``
  Provides line, polygonal, and polyhedral mesh types plus conversion to/from
  polytope structures, mesh generation, and Silo dump support.

``VoronoiCells``
  Provides physics packages and update policies for Voronoi-cell computation,
  sub-point pressure hourglass control, and Voronoi volume calculation.

``RK``
  Provides reproducing-kernel correction coefficients, correction methods,
  RK interpolation/gradient/hessian utilities, and volume policies.

``KernelIntegrator``
  Provides quadrature/integration support over kernels, flat connectivity,
  manufactured solutions, integration kernels, and integral coefficient
  machinery.

These systems are intentionally not hard-coded into the base integrator. A
physics package advertises what it needs, and the controller/integrator/state
assembly path provides the required shared structures.

Python Interface
----------------

Spheral's Python interface is built from two layers:

``src/PYB11``
  Binding-generation definitions. Each binding package calls
  ``spheral_add_pybind11_library``. Generated modules are named
  ``Spheral<Package>`` and are normally compiled as submodules of
  ``SpheralCompiledModules``.

``src/SimulationControl``
  User-facing Python modules. ``Spheral.py`` imports the compiled packages,
  initializes GPUs and Axom, imports material and node-list helpers, imports
  configured hydro constructors, registers pickle helpers, and exposes the
  controller and command-line parser. ``Spheral1d.py``, ``Spheral2d.py``,
  ``Spheral3d.py``, RZ, spherical, and solid convenience modules create
  dimension-specific aliases.

The generated ``SpheralCompiledPackages.py`` imports all compiled submodules
from ``SpheralCompiledModules`` and re-exports their public symbols. This is why
Python problem scripts can usually import from ``Spheral`` or a dimension module
without caring which compiled submodule originally defined a class.

``SpheralController`` is the main Python orchestration class. It is not an
integrator. Instead, it wraps an integrator with run-control services:

* restart restore and dump scheduling;
* visualization dump scheduling;
* periodic work callbacks;
* neighbor reinitialization;
* initial H iteration;
* initial derivative evaluation when needed;
* automatic insertion of axis and distributed boundaries;
* physics-package organization, including support packages requested by hydro
  methods;
* conservation and timing diagnostics.

Restart and File I/O
--------------------

Restart support is registration-based. Any restartable object provides
``dumpState``, ``restoreState``, and ``label``. The templated ``Restart`` handle
adapts that object to the global restart registrar. Core objects such as
``NodeList``, physics packages, integrators, and controllers register
themselves.

``FileIO`` defines the storage backend interface. It supports primitive values,
vectors, dimensional geometry types, fields, field lists, faceted volumes,
uniform random generators, and Python objects/bytes when Python is enabled.
Concrete backends include HDF5, Silo, Sidre, flat files, gzip-backed Python I/O,
and Python-overridable ``PyFileIO``.

The design keeps serialization responsibility with the owning object. For
example, a node list dumps its node fields, an integrator dumps time/cycle
state, and a physics package dumps its package-specific fields.

Testing Design
--------------

The test tree reflects the same layering as the source tree:

* ``tests/cpp`` contains compiled unit tests for low-level C++ packages such as
  fields, geometry, kernels, neighbors, loops, artificial viscosity, and
  utilities.
* ``tests/unit`` contains Python-driven unit tests for bindings and package
  behavior.
* ``tests/functional`` contains end-to-end physics problems, restart/reference
  checks, generators, hydro problems, gravity problems, damage problems, and
  material behavior.
* ``*.ats`` files define larger regression suites for the ATS harness.

For implementation changes, the lowest useful test level is usually the package
unit test. Shared changes in ``Field``, ``DataBase``, ``State``, ``Neighbor``,
``Boundary``, or ``Integrator`` should be validated with both C++ tests and at
least one Python functional path because those packages sit across the compiled
and user-facing layers.

Extension Patterns
------------------

Adding a Physics Package
~~~~~~~~~~~~~~~~~~~~~~~~

1. Add a new package under ``src/<Package>`` with a ``CMakeLists.txt``.
2. Derive from ``Physics`` or a more specific base such as ``GenericHydro`` or
   ``GenericBodyForce``.
3. Implement ``evaluateDerivatives``, ``dt``, ``registerState``,
   ``registerDerivatives``, and ``label``.
4. Store durable package fields as member ``FieldList`` objects.
5. Register state fields with appropriate update policies.
6. Register derivative fields with clear names, usually using field-name
   constants.
7. Override requirement hooks only for shared structures actually needed.
8. Implement boundary handling for package-owned fields.
9. Add restart support for durable fields and parameters.
10. Add explicit instantiation templates and CMake ``instantiate`` entries.
11. Add PYB11 bindings if the package is user-facing.
12. Add unit or functional tests.

Adding a State Field
~~~~~~~~~~~~~~~~~~~~

1. Define a stable field name, usually in a ``*FieldNames.hh`` file.
2. Decide ownership: node-list-owned durable state, physics-package scratch
   state, or transient derivative state.
3. Register the field in ``registerState`` with a policy if it is advanced.
4. Register the derivative or replacement field in ``registerDerivatives``.
5. Apply and enforce boundaries for the field if ghost values matter.
6. Include it in restart output if it is durable.
7. Include Python bindings or visualization registration when needed.

Adding a Dimension-Templated Class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Template the class on ``Dimension`` and use
   ``typename Dimension::Scalar``, ``Vector``, ``Tensor``, and related aliases.
2. Place nontrivial definitions in ``.cc`` files when possible.
3. Add ``<Class>Inst.cc.py`` with the explicit instantiation template.
4. Add the base name to the package ``*_inst`` list.
5. Ensure unsupported dimensions are filtered with the ``dimensions`` variable
   in the instantiation template if needed.
6. Add PYB11 bindings for each enabled dimension if the class is exposed to
   Python.

Adding Python Bindings
~~~~~~~~~~~~~~~~~~~~~~

1. Add or update the relevant ``src/PYB11/<Package>/<Package>_PYB11.py`` file.
2. Ensure the package is present in ``src/PYB11/CMakeLists.txt`` or gated by the
   right CMake option.
3. Use existing binding patterns for dimensional class names and constructor
   factories.
4. Expose Python helper/shadow modules from ``SimulationControl`` only when they
   are part of the user-facing API.

Important Invariants
--------------------

Several invariants are relied on across package boundaries:

* Node lists have a stable name and are globally ordered by
  ``NodeListRegistrar``.
* Every field belongs to zero or one node list and has a size matching that
  node list's internal-plus-ghost node count.
* Ghost nodes appear after internal nodes.
* Boundaries own the mapping between control, ghost, and violation nodes.
* Neighbor data must be refreshed after node positions, smoothing tensors,
  ghost-node counts, or domain ownership change.
* Connectivity must be rebuilt or patched after neighbor topology or ghost
  culling changes.
* ``State`` copies are reference copies unless ``copyState`` is called.
* Durable simulation state lives in node lists and physics-package members, not
  in transient ``State`` objects.
* Physics packages should request shared structures through requirement hooks
  rather than constructing them independently.
* Restartable objects must dump enough data to restore their externally visible
  state and internal durable fields.
* Python dimension aliases depend on generated class names containing ``1d``,
  ``2d``, or ``3d``.

Performance Design Considerations
---------------------------------

Spheral's performance model is built around local, typed arrays and explicit
pair traversal:

* ``Field`` data is contiguous per node list.
* ``FieldList`` avoids flattening across node lists unless explicitly needed.
* ``ConnectivityMap`` precomputes neighbor lists and node pairs so physics
  packages can focus on pair-wise math.
* The C++ template dimension model avoids runtime branching on dimension in
  inner loops.
* Explicit instantiation keeps compile/link behavior predictable for the
  enabled dimensions.
* Optional RAJA implementations put the most expensive pair and node loops on
  accelerators.
* ``work`` fields and redistribution support allow load balancing to account
  for computational cost, not just node count.

The main tradeoff is that many shared abstractions are template-heavy and
string-keyed at state-registration boundaries. The template model keeps inner
loops efficient; the state registry keeps integrators and physics packages
decoupled.

Configuration and Optional Features
-----------------------------------

The main CMake feature options are documented in ``cmake/SpheralOptions.cmake``.
Architecturally important options include:

``SPHERAL_ENABLE_PYTHON``
  Build Python bindings and install Python simulation-control modules.

``SPHERAL_BUILD_MAIN_PYTHON_MODULE``
  Build the aggregate ``SpheralCompiledModules`` pybind11 module.

``SPHERAL_ENABLE_1D``, ``SPHERAL_ENABLE_2D``, ``SPHERAL_ENABLE_3D``
  Select dimensions for explicit instantiation and bindings.

``SPHERAL_COMBINE_INSTANTIATIONS``
  Generate one instantiation file per class rather than one per class per
  dimension.

Optional physics package flags
  Include optional physics packages, such as ``SPHERAL_ENABLE_GRAVITY``,
  ``SPHERAL_ENABLE_GSPH``, ``SPHERAL_ENABLE_FSISPH``, and
  ``SPHERAL_ENABLE_SVPH``.

``ENABLE_MPI`` and ``ENABLE_OPENMP``
  Enable distributed and threaded execution through BLT dependencies.

``ENABLE_CUDA`` and ``ENABLE_HIP``
  Enable accelerator builds. These set ``SPHERAL_GPU_ENABLED`` and pull in the
  relevant compiler/runtime dependencies.

``SPHERAL_UNIFIED_MEMORY``
  Select unified-memory behavior for GPU builds.

``ENABLE_DEV_BUILD``
  Build internal packages as separate shared libraries for faster development
  rebuilds.

Conceptual Data-Flow Summary
----------------------------

The shortest complete summary of a Spheral time step is:

::

   NodeLists + package fields hold durable state
       |
       v
   DataBase builds field-list views and connectivity
       |
       v
   Physics packages register State and StateDerivatives
       |
       v
   Integrator applies boundaries and asks packages for dt
       |
       v
   Physics packages fill derivatives from fields/connectivity
       |
       v
   State update policies advance durable fields
       |
       v
   Boundaries and packages repair/finalize updated state
       |
       v
   Controller handles restart, visualization, periodic work, and loop control

This division is the main design principle to preserve. Node lists own data,
field lists expose layout-compatible views, physics packages define equations,
integrators define time-centering, and the Python controller defines run
workflow.
