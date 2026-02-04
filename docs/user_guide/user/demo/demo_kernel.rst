.. _kernel_generation:

====================
Interpolation kernel
====================

The next line constructs the choice of interpolation kernel for this simulation:

.. code-block::
   :linenos:
   :lineno-start: 79

   #-------------------------------------------------------------------------------
   # Create our interpolation kernel
   #-------------------------------------------------------------------------------
   WT = TableKernel(WendlandC4Kernel())
   output("WT")

The interpolation kernel is a fundamental piece of the numerics for an SPH type simulation, and the choice of this kernel is intimately tied with the choice of how many points per smoothing scale (the parameter ``nperh`` in this example script) should be aimed for during smoothing scale evolution (see :ref:`interpolation_kernels` and :ref:`smoothing_scales` for further discussion). That said, the combination of a Wendland kernel with an ``nperh=4.01`` as used here is a good starting point for most SPH calculations in Spheral.
