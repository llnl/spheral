.. _Sedov_controller:

======================
The Spheral controller
======================

Nearly done! The next stage is to build a SpheralController object. This is the object the user will use to advance a problem in time, specify when and how to output files such as restart checkpoint or visualization files, do any other optional periodic work, etc.

.. code-block::
   :linenos:
   :lineno-start: 189

   #-------------------------------------------------------------------------------
   # Build the controller.
   #-------------------------------------------------------------------------------
   control = SpheralController(integrator = integrator,
                               kernel = WT,
                               restartStep = restartStep,
                               restartBaseName = restartBaseName,
                               restoreCycle = restoreCycle,
                               vizBaseName = "Sedov-cylindrical-{}".format(nRadial),
                               vizDir = vizDir,
                               vizStep = vizCycle,
                               vizTime = vizTime,
                               SPH = True)
   output("control")

The SpheralController has a few required arguments, such as the time integrator (``integrator``) and interpolation kernel (``WT``). Note we have also specified when, how, and where to output restart files (``restartStep`` and ``restartBaseName``); similarly we specify where and when to drop visualization files (``vizDir``, ``vizStep``, and ``vizTime``).  There are many other options we can give the SpheralController, include the ability to have it run arbitrary user specified Python methods at regular intervals as the problem progresses.  This is the fundamental job of the SpheralController: take steps using the provided ``integrator``, and between steps do other necessary jobs including restarts, visualizations, analysis, or anything else we would like to have happen on some frequency.  See :ref:`SpheralController` for more information.
