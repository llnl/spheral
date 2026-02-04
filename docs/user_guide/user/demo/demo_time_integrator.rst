.. _Sedov_time_integrator:

===========================
The time integration choice
===========================

In order to advance our model in time to some final answer we're looking for, we need to tell Spheral what sort of time integration scheme to use. Spheral supports a number of time integration methods, all with different compromises in terms of speed, accuracy, preservation of invariants, implicit vs. explicit in time, etc.  There's not one correct answer for the choice of time integration, and Spheral lets the user decide what to try.  In this case we are building a variation of a predictor-corrector (or similarly a second-order Runge-Kutta time integrator where we resuse the mid-point derivatives for the next steps predictor) as follows:

.. code-block::
   :linenos:
   :lineno-start: 175

   #-------------------------------------------------------------------------------
   # Construct a time integrator, and add the one physics package.
   #-------------------------------------------------------------------------------
   integrator = CheapSynchronousRK2Integrator(db, packages)
   integrator.lastDt = dt
   integrator.allowDtCheck = True
   output("integrator")
   output("  integrator.havePhysicsPackage(hydro)")
   output("  integrator.dtGrowth")
   output("  integrator.lastDt")
   output("  integrator.dtMin")
   output("  integrator.dtMax")
   output("  integrator.allowDtCheck")

The time integrator is responsible for coordinating the time advancement for a single step across all physics packages, which is why we hand it the set of physics packages here (as well as the DataBase). On line 179 we also specify an initial guess for the timestep to be used in the first cycles advancement.  As is often the case there are many options to the time integrators, but in most cases the default values are good choices. Here we specify a minimum set of parameters to the integrator but use ``output`` to record the settings of several key parameters in the output.
