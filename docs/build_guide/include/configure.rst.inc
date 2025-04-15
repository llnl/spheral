..
   CONFIGURE SPHERAL INTRO
   ----------------------------------------

[ex_intro-start]
Configuring Spheral 
###################

After running ``tpl-manager`` you will see a file in your Spheral root directory following the format ``<sys_type>-<spec>.cmake``. 

For example if you run ``tpl-manager`` with only a single ``--spec`` e.g. ``spheral~mpi%gcc``, you will only see:

::

  <sys_type>-gcc~mpi.cmake

[ex_intro-end]

[lc_intro-start]
Configuring Spheral 
###################

After running ``tpl-manager`` you will see a file in your Spheral root directory following the format ``<sys_type>-<spec>.cmake``. 

For example if you run ``tpl-manager`` on a ``toss_4_x86_64_ib`` system with ``--spec spheral+mpi%clang``, you will see:

::

  toss_4_x86_64_ib-clang@9.0.0.cmake

If no ``--spec`` is provided to ``tpl-manager``, the TPLs for all specs in an environment will be built, but no host config file will be generated.
A ``--spec`` must be provided to generate a host config file.

[lc_intro-end]



..
   HOST-CONFIG-BUILD.PY
   ----------------------------------------

[ex_hcb-start]
host-config-build.py
====================

The ``host-config-build`` tool is located at ``scripts/devtools/host-config-build.py``. ``host-config-build`` takes a host-config file and sets up Spheral's CMake with the appropriate TPLs. ``host-config-build`` by default also sets up a basic build/install directory structure. 

::

  python3 scripts/devtools/host-config-build.py --host-config <sys_type>-gcc.cmake"

``--host-config`` is a **required** argument of the tool, by default this will create two directories in the format of:

::

  build_<host-config>/build
  build_<host-config>/install

.. note::
   ``host-config-build.py`` is simply a wrapper around CMake. It is possible to directly run CMake rather than ``host-config-build.py``; ``host-config-build.py`` reports the CMake line it is using, so this can be a good starting point if you need to run CMake manually yourself.  See :ref:`manual_config` for more details.

[ex_hcb-end]

[lc_hcb-start]
host-config-build.py
====================

The ``host-config-build`` tool is located at ``scripts/devtools/host-config-build.py``. ``host-config-build`` takes a host-config file and sets up Spheral's CMake with the appropriate TPLs. ``host-config-build`` by default also sets up a basic build/install directory structure. 

::

  ./scripts/devtools/host-config-build.py --host-config <sys_type>-gcc.cmake"

``--host-config`` is a **required** argument of the tool, by default this will create two directories in the format of:

::

  build_<host-config>/build
  build_<host-config>/install

.. note::
   ``host-config-build.py`` is simply a wrapper around CMake. It is possible to directly run CMake rather than ``host-config-build.py``; ``host-config-build.py`` reports the CMake line it is using, so this can be a good starting point if you need to run CMake manually yourself.  See :ref:`manual_config` for more details.

[lc_hcb-end]


..
   HELP
   ----------------------------------------

[help-start]
Help
----

``host-config-build`` supports ``-h`` or ``--help`` if you need to reference the available options.
[help-end]



..
   HOST-CONFIG-BUILD OPTIONS
   ----------------------------------------

[hcb_options-start]
Custom Build Directory
======================

If you wish your build directory to live somewhere else, run ``host-config-build`` from that directory and use ``--source-dir`` to point at the root Spheral dir.

Custom Install Directory
========================

You can setup a custom install location by using ``--install-dir``.

Option ``--build``
--------------------

The ``--build`` option will configure CMake as usual and then execute a build and install in parallel. 
[hcb_options-end]
[ex_hcb_options_note-start]

.. note::
   When building on non LC systems and using ``--build`` you can ignore any warning pertaining to ``lc-modules`` not set.

[ex_hcb_options_note-end]
[lc_hcb_options_note-start]

.. warning::
   If running on LC and using ``--build`` be sure to launch in a resource allocation, as ``--build`` will take advantage of parallel compilation.

--lc-modules
............

If you use build you may need some system modules in your environment during the build and install step. you can pass these to ``host-config-build`` with ``--lc-modules`` as so:

::

  ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake --build --lc-modules "gcc/8.1.0"

If ``--build`` is not passed ``--lc-modules`` will not do anything, you will need to ensure the correct modules are in your path before building manually.

[lc_hcb_options_note-end]




..
   ADDITIONAL CMAKE OPTIONS
   ----------------------------------------

[add_cmake_options-start]
Additional CMake Options
========================

With ``host-config-build`` we are still able to pass and override CMake arguments (See: `Spheral / CMake Configurations <cmake_config.html>`_). To do this add your CMake ``-D<XXXXX>`` options to your ``host-config-build`` command. This is particularly useful if you want to change the ``CMAKE_BUILD_TYPE`` or use a TPL that was not installed by ``tpl-manager``.

The example below shows how you can take the ``gcc`` host-config from above, and configure with ``Release`` and a custom ``PYB11Generator`` install.

[add_cmake_options-end]
[ex_add_cmake_example-start]
::

  python3 scripts/devtools/host-config-build.py --host-config <sys_type>-gcc.cmake" -DCMAKE_BUILD_TYPE=Release -Dpyb11generator_DIR=<PYB11generator_install_prefix>/lib/python3.9/site-packages/

[ex_add_cmake_example-end]
[lc_add_cmake_example-start]
::

  ./scripts/devtools/host-config-build.py --host-config <sys_type>-gcc.cmake" -DCMAKE_BUILD_TYPE=Release -Dpyb11generator_DIR=<PYB11generator_install_prefix>/lib/python3.9/site-packages/

[lc_add_cmake_example-end]




..
   MANUAL CMAKE CONFIGURATION
   ----------------------------------------

[manual_cmake-start]
.. _manual_config:

Manually Configure CMake
========================

``host-config-build.py`` is a tool for convenience, if you prefer to use CMake manually and set up your own build/install directory structure that is still very easy to do.

::

  mkdir -p Spheral_release/BUILD && cd Spheral_release/BUILD
  cmake -C ../../Spheral/<sys_type>-gcc.cmake \
        -DCMAKE_INSTALL_PREFIX=`cd ..; pwd` ../../Spheral

In this example we create our build directory ``Spheral_release/BUILD``, and will install Spheral in ``Spheral_release``.

The CMake option ``-C ../../Spheral/<sys_type>-gcc.cmake`` is how we tell the CMake system to use the TPLs we installed with ``tpl-manager.py`` for ``gcc``.

The somewhat obscure command ``-DCMAKE_INSTALL_PREFIX=`cd ..; pwd``` specifies the install directory as the full path to ``Spheral_release``. Alternatively, you can specify this path explicitly, such as ``-DCMAKE_INSTALL_PREFIX=/usr/local/Spheral_release``.
[manual_cmake-end]
