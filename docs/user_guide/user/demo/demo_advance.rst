.. _Sedov_advance:

===================================
Advancing and analysing the problem
===================================

The final lines in our example actually run the problem and plot the results.

.. code-block::
   :linenos:
   :lineno-start: 204

   #-------------------------------------------------------------------------------
   # Finally run the problem and plot the results.
   #-------------------------------------------------------------------------------
   control.advance(goalTime)

   #-------------------------------------------------------------------------------
   # Plot the final state if desired.
   #-------------------------------------------------------------------------------
   if graphics:
       rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)

Line 207 tells our SpheralController instance ``control`` to advance the simulation to the time in the variable ``goalTime``. Recall we set goalTime in the commandLine arguments block on (line 36); this means it defaults to ``goalTime=0.5``, but we can override this on the command line starting the run should we wish to.  The controller will then take as many steps as necessary to reach this goal time, periodically dropping restarts and visualization files as it progresses.

Upon completion of this advancement, the final commands in this script (line 213) plot scatter plot profiles (vs. radius) of the final state of the simulation (for mass density, radial velocity, specific thermal energy, pressure, and smoothing scale in this case).  These are plotted using the Python Matplotlib graphics package, using some Spheral utility methods imported from ``SpheralMatplotlib`` on line 10.

.. note::

   While this concludes what this scripted example does, the fact that Spheral is a Python application allows the user to continue exploring the model results however they wish. For instance running the ``spheral`` executable with the ``-i`` flag will cause the process to come back to an interactive prompt upon completion of the script rather than simply exit. At this point the user can grab any state Fields they like from the model (such as ``nodes.massDensity()`` or ``nodes.specificThermalEnergy()`` for example), and then query or numerically analyze further any interesting metrics they like. In this way Spheral is it's own data exploration utility as much as a modeling tool, depending on what the user in interested in.

