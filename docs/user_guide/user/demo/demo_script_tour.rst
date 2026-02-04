.. Spheral_script_tour:

#################################
A tour of a Spheral script
#################################

Running a problem with Spheral involves assembling the necessary objects in Python: typically a set of physics packages, equations of state, materials, time integrator, etc.  These can then be handed to an object designed to coordinate and run the problem (the ``SpheralController``), generating output such as visualization files or invoking custom analysis routines the user might provide.  Each of these concepts will be discussed in depth later in this manual, but first let's take a tour of a simple example designed to run a classic hydrodynamics test: the Taylor-Sedov blastwave :cite:`1959:Sedov`.

The analytical description of this problem starts with an infinitesimal energy spike deposited in a uniform pressureless gas.  This results in a shock or blastwave that propagates radially away from the original energy spike, and can be described by a self-similar analytic solution for a simple analytical material (such as a :math:`\gamma` law gas).

.. figure:: Sedov_demo.gif
   :width: 500px

   Time sequence of the mass density evolution in the Sedov problem using default Spheral visualization via `Visit <https://visit-dav.github.io/visit-website/>`_.  The black dots show the positions of the SPH nodes, while the spatial color map of density is produced by coloring polygonal Voronoi cells containing and constructed from the node positions.

First let's look at the example Python script in its entirety, and then we'll go through it section by section in more detail:

.. literalinclude:: Sedov-demo.py
   :language: Python
   :linenos:

Now let's go through each section of this script in some detail.

.. toctree::
   :maxdepth: 1

   demo_imports.rst
   demo_command_line.rst
   demo_directories.rst
   demo_materials.rst
   demo_kernel.rst
   demo_nodes.rst
   demo_point_generation.rst
   demo_initial_conditions.rst
   demo_database.rst
   demo_hydro.rst
   demo_boundary_conditions.rst
   demo_time_integrator.rst
   demo_controller.rst
   demo_advance.rst
