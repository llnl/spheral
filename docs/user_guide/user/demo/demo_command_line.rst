.. _demo_command_line:

====================================
Command line definitions (optional)
====================================

The next section in this script uses a handy function to allow us to override some variables from the command-line:

.. code-block::
   :linenos:
   :lineno-start: 14

   #-------------------------------------------------------------------------------
   # Generic problem parameters
   #-------------------------------------------------------------------------------
   commandLine(
       # Geometry of the problem
       rmin = 0.0,
       rmax = 1.0,
    ...
    )

The ``commandLine()`` function is a wrapper around a standard Python method (`argparse <https://docs.python.org/3/library/argparse.html>`_) of allowing variables to be overridden on the command line, and designed to be a simpler interface to those more generic Python command line utilities.  Each entry in the ``commandLine()`` function is given a default value, and if not otherwise overridden results in a Python variable in the script that is defined to be that value.  For instance in this case following this block there will be a floating point variable named ``rmax`` with the value ``1.0``.  Because ``rmax`` is availble in the ``commandLine()`` statement we could choose to override its value to be 2 by adding ``--rmax 2.0`` on our Spheral run line.

The advantage of using this sort of ``commandLine()`` block to define variables like this is that it is then quite simple to run variations on a Spheral script by overriding these values on the command line without directly editing the script or making new versions.  This is particularly useful when we want to run a suite of calculations varying some set of parameters.  In that sort of case we can simply make sure the variables we want to vary are defined in the ``commandLine()`` statement, and then set up runs that vary this value without making whole new Spheral scripts with such minor changes.

For example, in this script the value of the energy spike that creates the blastwave is defined by the variable ``Espike``, with a default value of 0.25.  Say we wanted to run several variations of the Sedov problem with different initial energies, this can be very simply accomplished by running this script multiple times and just overriding this one value::

  spheral Sedov-demo.py --Espike 0.1
  spheral Sedov-demo.py --Espike 0.5
  spheral Sedov-demo.py --Espike 1.0
  spheral Sedov-demo.py --Espike 10.0

If you were running this on a multicore computer or cluster, these variations could all be run in parallel to one and other from our single ``Sedov-demo.py`` script.

The use of this ``commandLine()`` function is entirely optional -- one can just as easily write a Spheral script and make these variables have fixed values in this script.  This is entirely up to the script developer, as is the choice of what is useful to provide for command line options like this.
