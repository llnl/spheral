[intro-start]
The Spheral build system works in two stages. 

 - Stage 1: Building and setting up Third Party Libraries (TPL)s.
 - Stage 2: Building and installing Spheral.

Spheral uses `Spack <https://github.com/llnl/spack>`_ under the hood to handle Third Party Library dependencies. Spack will track dependencies between TPLs and version constraints of TPLs as the project develops. Spack is also particularly good at handling TPL setup across various configurations, compilers and build time options.
[intro-end]


[tpl_manager-start]
tpl-manager.py
==============

Spheral provides a tool (``tpl-manager.py``) in an attempt to simplify the spack process for the user. ``tpl-manager`` 's primary responsibilities are to:

 - Set up a local ``Spack`` instance for Spheral.
 - Generate a dependency tree of third party libraries (relative to the provided configuration).
 - Build and install all dependent libraries in the local ``Spack`` instance.
 - Generate a `CMake host-config <https://llnl-blt.readthedocs.io/en/develop/tutorial/host_configs.html>`_  file for configuring Spheral builds (if a spec is given).

``tpl-manager`` is located at ``scripts/devtools/tpl-manager.py``. ``tpl-manager`` can be used in two ways:

 - Build TPL's for a single compiler configuration (ideal for users).
 - Build and setup the full range of TPLs for all of our supported compiler configurations (ideal for developers).

[tpl_manager-end]

[ex_tpl_manager_notes-start]

 .. note:
    The External User Guide does not cover building the full range of tpl configurations. Please see the LC User Guide for more information.

[ex_tpl_manager_notes-end]
[lc_tpl_manager_notes-start]

 .. note:
    The External User Guide does not cover building the full range of tpl configurations. Please see the LC User Guide for more informtaion.

[lc_tpl_manager_notes-end]



[ex_running_tpl_manager-start]
Running tpl-manager
===================

``tpl-manager`` requires ``python3`` and should be run from the root Spheral directory.

``tpl-manager`` takes a ``--spec`` argument to determine what compiler to use and what configuration we want to build Spheral in. For example,

::

  ./scripts/devtools/tpl-manager.py --spec spheral+mpi%gcc

To see what specs are available for the current environment, run

::

  ./scripts/devtools/tpl-manager.py --show-specs

This will download a new Spack instance in a directory adjacent to your Spheral root directory called ``spheral-spack-tpls/spack``. You can use ``--spack-dir`` if you would like to setup the Spack instance somewhere else. The ``%`` symbol denotes the compiler to use.

When run for the first time, ``tpl-manager`` will do a few things:

#. Create and activate a Spack environment in ``spheral/scripts/spack/environments`` based on the output of ``spack arch``. All subsequent Spack commands will modify this environment by adding specs, compilers, and external packages to it.

#. Run ``spack compiler find`` to find system compilers. Spack searches the ``$PATH`` environment variable for compilers and packages. Add any paths to ``$PATH`` to ensure Spack will find them before running the ``tpl-manager``.

#. Run ``spack external find --not-buildable`` It will run Spack commands to try to find existing system installs for things like CMake, Git, Python, and MPICH (when the spec has ``+mpi``).

#. Add the current spec to the environment with ``spack add <spec>`` and concretize using ``spack concretize -U``.

#. Install the TPLs and create the host config file for the given spec with ``spack install -u initconfig <spec>``.

.. note::
   By default we have ``python`` bindings enabled (``+python``) and docs disabled (``~docs``). Therefore the spec ``spheral+mpi+python~docs%gcc`` will build the same TPL set as just ``spheral+mpi%gcc``.For more information on ``spec`` syntax please see the spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.

.. note::
   Spheral minimally requires a C++14 compliant compiler.

[ex_running_tpl_manager-end]

[lc_running_tpl_manager-start]
Running tpl-manager
===================

``tpl-manager`` requires ``python3`` and should be run from the root Spheral directory.

``tpl-manager`` takes a ``--spec`` argument to determine what compiler to use and what configuration we want to build Spheral in.

::

  ./scripts/devtools/tpl-manager.py --spec spheral+mpi%gcc

This will install the local Spheral Spack instance into the adjacent directory of your Spheral root dir. You can use ``--spack-dir`` if you would like to setup the spack instance somewhere else. 

.. warning::
   If running on LC be sure to launch ``tpl-manager`` in a resource allocation, as ``tpl-manager`` will take advantage of parallel builds when compiling libraries.  An example appropriate for use on rzgenie would be::

    srun --exclusive ./scripts/devtools/tpl-manager.py --spec spheral+mpi%gcc

Above we are telling ``tpl-manager`` to build our TPLs with ``gcc``. The Spack compiler file for the current system (in ``./scripts/spack/configs/$SYS_TYPE/compilers.yaml``) only has one gcc compiler version, so that is the one that will be used. If the compiler file is aware of multiple versions of gcc, a specific version should be given, ie ``spheral+mpi%gcc@10.3.1``. Which variants should be provided depend on the root specs in the current environment. For example, the TOSS4 machine environment has three specs: ``spheral+mpi%gcc``, ``spheral~mpi%gcc``, and ``spheral+mpi%clang``. To turn off mpi with the gcc compiler, use the following spec:
::

  ./scripts/devtools/tpl-manager.py --spec spheral~mpi%gcc

To see what specs are available for the current environment, run
::

  ./scripts/devtools/tpl-manager.py --show-specs

.. note::
   By default we have ``python`` bindings enabled (``+python``) and docs disabled (``~docs``). Therefore the spec ``spheral+mpi+python~docs%gcc`` will build the same TPL set as just ``spheral+mpi%gcc``.For more information on ``spec`` syntax please see the spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.


CUDA Support
------------

To build Spheral TPLs and configure Spheral to build with CUDA on LC systems we need to add the spack variant ``+cuda`` to our spec. We also need to define the cuda hardware architecture we are targetting Spheral to run on. This can be done with ``cuda_arch=XX``. On most LC NVIDIA systems we will be targetting ``sm_70`` architecture. Therefore a gcc spec with cuda support will look as such.

::

  ./scripts/devtools/tpl-manager.py --spec "spheral+cuda cuda_arch=70%gcc"

.. note::
   We need to use "" here as spack requires spacing to define ``cuda_arch``. This ensures tpl-manager doesn't register the spec as separate arguments.

Setup Full TPL List
-------------------

The simplest way to have ``tpl-manager`` build all TPLs for all supported Spheral configurations is by running this command from your Spheral root directory:

::

  ./scripts/devtools/tpl-manager.py

Note, this will build TPLs for all permutations of the specs in the current SYS_TYPE environment. It will not produce any host config files.

[lc_running_tpl_manager-end]

[help-start]
Help
----

``tpl-manager`` supports ``-h`` or ``--help`` if you need to reference the available options.
[help-end]

[manual_tpl-start]
Manually installing TPLs
========================

Although it is HIGHLY recommended, you do not need to use ``tpl-manager`` to setup TPLs for Spheral. TPLs can be built manually or built through your own spack installation and passed to the Spheral CMake system. See :ref:`manual_config` for more details.
[manual_tpl-end]
