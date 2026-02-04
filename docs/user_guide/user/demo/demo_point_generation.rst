.. _point_generation:

=======================================================
Generating and distributing our nodes across processors
=======================================================

This section of ``Sedov-demo.py`` is concerned with how we generate the point distribution that represent the geometry/initial conditions of the problem, as well as how those points should be distributed across MPI parallel domains:

.. code-block::
   :linenos:
   :lineno-start: 93

   #-------------------------------------------------------------------------------
   # Generate the initial node geometry (lay down the points).
   #-------------------------------------------------------------------------------
   generator = GenerateNodeDistribution2d(nRadial = nRadial,
                                          nTheta = 1,
                                          rho = rho0,
                                          distributionType = "constantDTheta",
                                          rmin = rmin,
                                          rmax = rmax,
                                          xmin = (0.0, 0.0),
                                          xmax = (rmax, rmax),
                                          theta = 0.5*pi,
                                          nNodePerh = nPerh)
   distributeNodes2d((nodes, generator))
   print("Point distribution across MPI ranks:")
   output("  mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
   output("  mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
   output("  mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

This is the section of our script describing how the points in our NodeList created in :ref:`nodelist_construction` will be created, and once this section is completed our initally empty NodeList ``nodes`` will be filled with points.

The concept of how points are created and initialized is broken into two pieces in Spheral: we use ``NodeGenerators`` (``GenerateNodeDistribution2d`` in this case) to describe how points should be laid down and the most basic properties (typically mass, mass density, and smoothing scale in addition to position) assigned to those points.  Think of the NodeGenerator as a mathematical prescription for how points should be generated, but it does not do the actual assignment of that information to the NodeLists. That task falls to the ``Distributor`` chosen (in this case ``distributeNodes2d`` imported from ``PeanoHilbertDistributeNodes`` at the start of the script).

There are many generators canned and available in Spheral which should fulfill the most basic sorts of point geometry descriptions. However, one of our intentions with making Spheral a scriptable Python application is that the user should be able to generate arbitrarily complex rules for how they would like to lay down points and create their problem, which can be tailored to the exact needs of the problem at hand.  Therefore for complex problems it is sensible (and expected) that the user script will include a custom NodeGenerator which has been specialized to describe the particular problem in question.  Therefore it is the NodeGenerator we would expect to vary between problems depending on the geometry being modeled and the users needs.

The Distributor however is a different beast. The Distributor takes a set of ``(NodeGenerator, NodeList)`` pairs (i.e., one generator per NodeList), and decides how those points described by the generators should be divvied up between processor domains to best load balance the problem.  There are several such distributors available in Spheral, though the Peano-Hilbert algorithm (based on the Peano-Hilbert space filling curve) used here is a good generic choice.  This sort of load balancing problem is a finicky and complicated beast, and probably does not need to vary wildly between problems.  So unlike the geometry description, it is probably less useful for the user to try and custom craft this distribution algorithm themselves (though the adventurous user is certainly allowed to!)

This reasoning is why Spheral has adopted this split in how to describe and effect the point generation.  The user can pick from a wide variety of possible NodeGenerator descriptions for how they would like to lay down points (or concoct their own from whole cloth).  This information is then fed to one of a few sorts of node distributors to actually use the generator description to fill in the NodeLists.
