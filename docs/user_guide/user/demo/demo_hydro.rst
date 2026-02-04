.. _Sedov_hydro:

==================================
Building the hydrodynamics package
==================================

Finally we come to building a Spheral physics package, in this case one of the hydrodynamics packages:

.. code-block::
   :linenos:
   :lineno-start: 151

   #-------------------------------------------------------------------------------
   # Construct the hydro physics object.
   #-------------------------------------------------------------------------------
   hydro = FSISPH(dataBase = db,
                  W = WT)
   packages = [hydro]

   output("hydro")
   output("  hydro.cfl")
   output("  hydro.compatibleEnergyEvolution")
   output("  hydro.densityUpdate")
   output("  hydro.HEvolution")

In this case we are building an FSISPH :cite:`2022:Pearl_fsisph` hydrodynamics object, which is a variation on Smoothed Particle Hydrodynamics.  Typically there are many options to these Physics objects, for most of which we are using the defaults here.  The required arguments are the DataBase of NodeLists (``db``), and the chosen interpolation kernel (``WT``).  It's generally good practice to print (``output``) some of the major options the package is using to capture in the script output, such as we do in this example on lines 158--162.

.. note::

   In many Spheral scripts you may have multiple physics packages (for instance hydrodynamics, gravity, porosity, etc.)  A standard pattern is to build up a Python list of these packages (such as is started here on line 156 in the variable ``packages`` -- then as we build any other Physics packages we simply append them to this list. In this example there is only the one Physics model running though so the situation is simplified.  The utility of this ``packages`` list will come up when we build the time integrator in :ref:`Sedov_time_integrator`.

