.. _Sedov_initial_conditions:

=============================================
Setting initial conditions (point properties)
=============================================

The next section of ``Sedov-demo.py`` sets the initial conditions for our problem::

  #-------------------------------------------------------------------------------
  # Set the node properties.
  # The above geometry initialization set the node positions, masses, mass
  # density, and H tensors.  We still need to initialize the energy for the Sedov
  # blast energy.
  # In this case we're going to simply pick the innermost ring of points closest
  # to the origin and dump the initial spike evenly across them.
  #-------------------------------------------------------------------------------
  pos = nodes.positions()
  mass = nodes.mass()
  eps = nodes.specificThermalEnergy()

  # Find the points closest to the origin.
  dr = rmax/nRadial
  spikePoints = [i for i in range(nodes.numInternalNodes) if pos[i].magnitude() < 0.75*dr]
  print("Selected a total of {} points for energy deposition.".format(mpi.allreduce(len(spikePoints))))

  # Distribute the energy evenly in energy/mass.
  Mspike = mpi.allreduce(sum([mass[i] for i in spikePoints]), mpi.SUM)
  eps0 = Espike/Mspike
  for i in spikePoints:
      eps[i] = eps0

  # Verify we actually initialized the expected energy.
  # We are using multiplication of Spheral Fields here, and the fact that Field.sumElements
  # sums across all MPI domains by default
  Esum = (mass*eps).sumElements()
  print("Initialized a total energy of {}".format(Esum))
  assert fuzzyEqual(Esum, Espike)

In general the generator/distributor stage that precedes this block will set the node positions, masses, smoothing scales, and mass density. Any other quantities (like the specific thermal energy, velocity, etc.) needs to be explicitly set following that generation stage.

.. note::

   It is worth emphasizing that until the the distributor is called with the ``(generator, NodeList)`` pairs on line 106 of ``Sedov-demo.py`` there are no nodes in the problem, so there are no positions to check against when setting initial conditions, and all Fields of physical values are of zero size.  It is only after this distribution stage that the NodeLists are populated and the points assigned positions to check in order to set other physical state values.

The first couple of lines in this section grab the Fields for the position, mass, and specific thermal energy of the nodes in the NodeList ``nodes``.  Fields are essentially arrays of values associated with a given NodeList.  Fields are always sized correctly for the NodeList they are associated with: however many nodes are in the NodeList, that is how many elements will be in each Field associated with that NodeList.

The initial conditions for the classic Taylor-Sedov (or blastwave) problem :cite:`1959:Sedov` consist of uniform pressureless, mostionless gas into which a delta function spike of energy is introduced.
