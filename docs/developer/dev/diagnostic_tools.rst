Code Debugging and Diagnostics
##############################

Valgrind
========

We advise using Valgrind to check memory leaks when doing development on Spheral.
When using Valgrind to check Spheral, be sure to use the provided suppression file and call the Spheral python executable directly:
::

   valgrind --suppressions=./scripts/devtools/valgrind_python_suppression .venv/bin/python


Using Caliper
=============

Spheral uses Caliper to perform code diagnostics, such as timing. To enable this functionality in the code, Spheral needs to be configured with ``ENABLE_TIMER=ON``. Otherwise, the timing regions are no-ops for improved performance.
::

  ./scripts/devtools/host-config-build.py <sys_type>-<spec>.cmake -DENABLE_TIMER=ON

Caliper is configured and started through the ``cali::ConfigManager``.
The ``cali::ConfigManager`` is wrapped in a ``TimerMgr`` singleton class, which has a python interface.

.. note::
   ``TimerMgr`` is initialized in ``src/SimulationControl/SpheralTimingParser.py`` which is called during ``commandLine()`` in ``src/SimulationControl/SpheralOptionParser.py``. This is because ``commandLine()`` is almost always invoked directly near the start of a problem. However, if ``commandLine()`` is not called, the timer manager would need to be configured and started directly using the ``TimerMgr`` class. See :ref:`below <manual_caliper>` for more details.

By default, the Caliper configuration is set to ``spot`` and outputs Caliper files (``.cali``).

There are many different Caliper configurations to view various information. Here are some extra links for those who want to read or experiment with other features in Caliper that can be incorporated into Spheral:

  * `Configuration basics <https://software.llnl.gov/Caliper/CaliperBasics.html#more-on-configurations>`_
  * `Builtin Configuration <https://software.llnl.gov/Caliper/BuiltinConfigurations.html>`_
  * `Manual Configuration <https://software.llnl.gov/Caliper/configuration.html>`_
  * `Output Format <https://software.llnl.gov/Caliper/OutputFormats.html>`_

Caliper and Adiak Options
-------------------------

.. option:: --caliperFilename

            Name of Caliper timing file. Should include file extensions. Optional, default: ``name_of_file_YEAR_MONTH_DATE_TIME.cali``.

.. option:: --caliperConfig CONFIG_STR

            Specify a built-in Caliper configuration or turn off timers with ``none``. Optional, default: ``spot``.

            **Example**:
            ::

               ./spheral ex_prog.py --caliperConfig 'runtime-report(output=time.txt),calc.inclusive,region.count'

.. note::
   The configuration in the example above produces timing results similar to the previous ``Spheral::Timer`` method. This results in a file named ``time.txt`` with cumulative times for the nested regions as well as a count of how many times each region ran.

.. option:: --caliperConfigJSON JSON_FILE

            Specify a JSON file containing a non-default Caliper configuration. Optional.

.. option:: --adiakData ADIAK_DATA_STR

            Specify any Adiak data directly in the command line. Must be a string in key:value format, separated by commas. Optional.

            **Example**:
            ::

               ./spheral ex_prog.py --adiakData "test_name: the_cheat, test_num:10"

.. note::
   The ``--adiakData`` input is useful for specifying metadata but leaving the python script unchanged, such as when running tests through ATS. In most cases, adding Adiak metadata through the python interface is preferred. See :ref:`below <python_adiak>` for more details.

Adding Region Timers in C++
---------------------------

The following macros are used to create timing regions in the Spheral C++ interface:

- ``TIME_FUNCTION`` can be added to the very beginning of a function and creates a region for the entire function using the function's name. ``TIME_FUNCTION`` uses just the function name and no class or parameter information, so be careful when using this method with functions that could share names.

- ``TIME_BEGIN("timer_name")`` and ``TIME_END("timer_name")`` create a region between the two different calls and use the string (in this case ``timer_name``) as the name.


Adding Region Timers in Python
------------------------------

Region timers can be added inside the python code using the following function calls:
::

   from SpheralUtilities import TimerMgr
   TimerMgr.timer_begin("timer_name")
   some_function_call()
   TimerMgr.timer_end("timer_name")

.. note::
   All timers must have both a start and end call. Otherwise, memory issues will occur.

.. _python_adiak:

Adding Adiak Metadata in Python
-------------------------------

Adiak metadata can be added inside python code using the following function calls:

.. code-block:: python

                adiak_values("value_name", value)

Below is a list of some of the metadata that Spheral adds to Adiak by default:

======================== ==========================
Adiak Key                Description
======================== ==========================
``user``                 User
``cluster``              Hostname (ie rzgenie)
``jobsize``              Number of ranks
``numhosts``             Number of allocated nodes
``total_internal_nodes`` Number of SPH nodes
``total_steps``          Number of time steps
``dim``                  Number of dimensions
``threads_per_rank``     Number of threads per rank
``git_hash``             Short Git hash of source
``git_branch``           Git branch of source
======================== ==========================

.. _manual_caliper:

Starting Caliper Manually
-------------------------

As mentioned above, the Caliper timing manager is normally configured and started in the ``commandLine()`` routine. However, Caliper can be directly configured and started through the python interface, if desired. This can be done by putting the following into the python file:
::

   from SpheralUtilities import TimerMgr
   caliper_config = "some_configuration(output=some_filename.txt)"
   TimerMgr.add(caliper_config)
   TimerMgr.start()

Performance Regression Testing
==============================

.. note::
   The following is currently only applicable for use on LC machines at LLNL.

``tests/performance.py`` contains a set of performance regression tests. These tests allow a developer to estimate the performance implications of code under development and compare it to the current development branch of Spheral.
When a merge to the develop branch occurs, the CI runs this regression test multiple times to accumulate benchmark timing data in a centralized directory (``/usr/WS2/sduser/Spheral/benchmark``).
The general procedure to comparing performance regression tests is:

* Run the performance regression tests from an installation using 2 nodes (number of nodes used in benchmark run by CI):
   ::

      ./spheral-ats --log test_dir_name --numNodes 2 tests/performance.py

   There is also a ``--threads`` option to specify a given number of threads per rank.
   The test above will create an ATS python file, ``test_dir_name/atsr.py``, as well as

* In general, to compare the performance between two performance results, use:
   ::

      ./spheral ./scripts/performance_analysis.py --perfdata1 /path/to/perf/data --perfdata2 /path/to/other/perf/data

* To compare newly run times with reference times for regression testing purposes, use:
   ::

      ./spheral ./scripts/performance_analysis.py --perfdata test_dir_name --ref /directory/of/reference/caliper/files/

   The input to ``--ref`` can be also be an ATS directory created from running ``performance.py`` or just a directory of Caliper files.
   Removing the ``--ref`` input will default to comparing to benchmark timings in ``/usr/WS2/sduser/Spheral/benchmark``.
   Timing trees can be displayed using the ``--display`` flag.
   The script above computes the mean (:math:`\mu`) and standard deviation (:math:`\sigma`) of the inclusive average time per rank (``Avg time/rank``) timers for each test in the reference (or benchmark) data.
   It computes a timing threshold using:

   .. math::

      \delta_{\mathrm{thresh}} = 0.08 \mu + 2 \sigma

   If the :math:`t_c - \mu > \delta_{\mathrm{thresh}}` for the ``main`` region, where :math:`t_c` is the new performance time, the test is considered to have failed and the timing tree of the exclusive average time per rank (``Avg time/rank (exc)``) will be displayed.
   If the test configurations, like the number of time steps differed between the runs, or the hardware/install configurations did not match, it will consider the test skipped.
   Otherwise, the test is considered to have passed.
   If the :math:`t_c - \mu < -\delta_{\mathrm{thresh}}`, the performance improved significantly and the timing tree will be displayed.

.. note::

   If ``performance.py`` is run on a non-MPI Spheral build, it will only use 1 rank and will thread all other available cores.
