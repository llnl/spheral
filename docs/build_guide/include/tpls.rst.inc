
Spheral provides a tool (``tpl-manager.py``) in an attempt to simplify the Spack process for the user, which does the following:

 - Clone and setup a local Spack instance for Spheral.
 - Activate (and create, if necessary) a Spack environment for the current system architecture.
 - Generate a dependency tree of third party libraries (relative to the provided configuration).
 - Build and install all dependent libraries in the local Spack instance.
 - Generate a `CMake host-config <https://llnl-blt.readthedocs.io/en/develop/tutorial/host_configs.html>`_ file for configuring Spheral builds. This file will be named ``<sys_type>-<spec>.cmake``.

``tpl-manager.py`` is located in ``scripts/devtools`` and must be run from the Spheral root directory.
The most common usage of ``tpl-manager.py`` will look like:

::

  ./scripts/devtools/tpl-manager.py --spec spheral+mpi%gcc

The command above tells the TPL manager to build and install TPLs for Spheral with MPI enabled and using the GCC compiler.
Overall, the following options can be provided to ``tpl-manager.py``:

.. option:: --spec <SPEC>

            Spheral Spack spec to build. Optional, default behavior is to build all spec permutations in the current environment. Required if running on non-LC systems. See the Spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_ for more information.

.. option:: --show-specs

            Display specs found in current environment, if environment exists.

.. option:: --show-info

            Display Spheral Spack info, including variants and dependencies.

.. option:: --spack-dir <SPACK_DIR>

            Where to create the Spack instance to use. Optional, default ``./../spheral-spack-tpls``.

.. option:: --clean

            Cleans Spack repo instance and forces Spack to reconcretize environment. Use if running into odd behavior with the script.
