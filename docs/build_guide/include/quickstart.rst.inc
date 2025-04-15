..
   SYSTEM UPDATE
   ----------------------------------------

[ex_update_sys-section-start]
This page details commands for setting up and running Spheral on an Ubuntu 20.04 system. 

Update and install necessary package dependencies.
::

  sudo apt update
  sudo apt upgrade
  sudo apt install build-essential git gfortran mpich autotools-dev autoconf sqlite pkg-config uuid gettext cmake libncurses-dev libgdbm-dev libffi-dev libssl-dev libexpat-dev libreadline-dev liblapack-dev libbz2-dev locales python python3 unzip libtool wget curl tk-dev

[ex_update_sys-section-end]

..
   DIRECTORY STRUCTURE
   ----------------------------------------

[create_dir-section-start]
Create a directory structure and clone spheral.
::

  mkdir -p SPHERAL && cd SPHERAL
  git clone --recursive https://github.com/llnl/spheral
  cd spheral

[create_dir-section-end]

..
   TPLS
   ----------------------------------------

[ex_tpl-section-start]
Build our TPL dependencies from source with the Spheral tpl-management tool (``tpl-manager.py``).
::

  ./scripts/devtools/tpl-manager.py --spec spheral%gcc

.. note::
   This command sequence assumes ``gcc`` is installed and will use the version in your system path. If you wish to use a different compiler (such as ``clang``) ensure it is in your path and replace ``gcc`` with the compiler name of your choice (e.g. ``clang``).

.. note::
  This command will generate a ``.cmake`` file with the naming convention ``<system-type>-<compiler-spec>``. The following commands will refer to this format as ``<host-config>`` for generalization across operating systems and architectures. You will need to substitute the correct format in the following commands.

.. warning::
   For operating systems other than Ubuntu 20.04 you should get an error to the effect of:
   ``[ERROR: invalid spack config dir: /<path>/scripts/spack/configs/<OperatingSystem><Version> ]``
   You will **need** to follow the steps outlined in :ref:`ERROR: invalid spack config dir` before proceeding with this quickstart guide.


[ex_tpl-section-end]

[lc_tpl-section-start]
Build our TPL dependencies from source with the Spheral tpl-management tool (``tpl-manager.py``).
::

  ./scripts/devtools/tpl-manager.py --spec spheral+mpi%gcc

.. warning::
  This command needs to be run under an allocation as any TPLs that need to be built will be built in parallel.

.. note::
  This command will generate a ``.cmake`` file with the naming convention ``<system-type>-<compiler-spec>``. The following commands will refer to this format as ``<host-config>`` for generalization across operating systems and architectures. You will need to substitute the correct format in the following commands. 

[lc_tpl-section-end]

..
   BUILD
   ----------------------------------------

[ex_build-section-start]
Set up a build/install directory structure and configure cmake.
::

  python3 scripts/devtools/host-config-build.py --host-config <host-config>.cmake

[ex_build-section-end]

[lc_build-section-start]
Set up a build/install directory structure and configure cmake.
::

  ./scripts/devtools/host-config-build.py --host-config <host-config>.cmake

[lc_build-section-end]

..
   INSTALL
   ----------------------------------------

[ex_install-section-start]
Build and install Spheral.
::

  cd build_<host-config>/build
  make -j <N> install

where ``<N>`` should be replaced with the number of parallel make processes you want to launch (typically the number of cores on your CPU).

[ex_install-section-end]

[lc_install-section-start]
Build and install Spheral.
::

  cd build_<host-config>/build
  make -j <N> install

.. warning::
 ``make`` should be run under an allocation (i.e., ``srun`` or ``salloc`` on Slurm) as it will take considerable resources on LC to build and install Spheral.  An exammple appropriate for rzgenie would be::

  srun --exclusive make -j 36 install

[lc_install-section-end]


..
   TESTING
   ----------------------------------------

[test-section-start]
Run a basic smoke test for Spheral
::

  cd ../install
  ./spheral -c "import Spheral"

Run our full test suite.
::

  ./spheral-ats test/integration.ats

These commands are explained in further sections.
[test-section-end]
