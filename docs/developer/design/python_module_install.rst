Python Module Installation Strategy
===================================

The full list of Spherals python module dependencies:

::

  python==2.7.18
  pip==9.0.1
  setuptools
  wheel
  virtualenv==20.2.2
  pyb11generator==1.0.12
  sphinx==1.8.5
  sphinx_rtd_theme==0.5.0
  docutils==0.17.1
  numpy==1.16.6
  numpy-stl==2.11.2
  matplotlib==2.2.5
  decorator==4.4.2
  h5py==2.10.0
  twine==1.15.0
  cython==0.29.21
  sobol==0.9
  scipy==1.2.3
  pipreqs==0.4.10
  importlib_metadata==2.1.1
  mpi4py==3.3.0
  gnuplot==1.8
  ats==5.2

However many of these are not required to actually build Spheral, they are purely dependencies of the Spheral python interface.

----

Build Time Packages
-------------------

These packages are **always** required when building the Spheral python interface:

::

  python==2.7.18
  pip==9.0.1
  setuptools==?
  wheel==?
  virtualenv==20.2.2
  pyb11generator==1.0.12

Optionally if we are building Spheral with ``SPHERAL_ENABLE_DOCS=On`` we will need:

::

  sphinx==1.8.5
  sphinx_rtd_theme==0.5.0
  docutils==0.17.1

Implementation
--------------

- These packages are used during the build as they are necessary utilities to compile Spheral source code and other targets.
- Each package is designed to be install / used from a directory outside of the python ``site-packages`` directory. Therefore during the build Spheral creates a ``PYTHONENV`` list of each of the packages respective directories and calls the python executable with the ``PYTHONENV`` variable as necessary.

CMake Interface
---------------

- With the new python package management system we expose these python packages the same way we expose other TPL 	dependencies in Spheral. They can be enabled and disabled through CMake options ``<package_name>_DIR`` and ``<package_name>_BUILD`` this is useful as we move towards a Spack based TPL system where python packages are treated as separate dependencies from a given python installation.

e.g. Using a locally built / installed version of a ``pyb11generator`` package

::

  cmake ..... -Dpyb11generator_BUILD=Off -Dpyb11generator_DIR=/<local_pyb11gen_site_package_dir>

----

Runtime Packages
----------------

::

  numpy==1.16.6
  numpy-stl==2.11.2
  matplotlib==2.2.5
  decorator==4.4.2
  h5py==2.10.0
  twine==1.15.0
  cython==0.29.21
  sobol==0.9
  scipy==1.2.3
  pipreqs==0.4.10
  importlib_metadata==2.1.1

The above packages are always required by Spheral at runtime and are added to a ``requirements.txt`` file that is generated during the CMake configuration stage. If ``ENABLE_MPI=On`` then ``mpi4py`` is also added to this list of packages.

The generated ``requirements.txt`` file is used in the ``spheral-setup-env.sh`` script to install these directly to the Spheral virtual environment we create.

- Note :  To continue to support off-line installations, all package ``.tar`` files are still saved to the ``${CACHE_DIR}`` and all install commands check this directory for previously downloaded ``.tar`` files.

-----

Custom Runtime Packages
-----------------------

::

  gnuplot==1.8
  ats==5.2

These packages do not have Pypi versions of themselves and are therefore installed alongside the Runtime Packages in the Spheral virtual environment setup script.

----

Spack Integration
-----------------

As we move to a Spack based TPL system it is useful to mimic the separation of python and python-package installs. This system allows us to expose our build time python dependencies for Spack to handle. However once we are past the build stage of Spheral we do not want Spack to manage our runtime dependencies as that requires users to use Spack as an interface to Spheral: launching the spack environment, loading the spack/spheral module, dirtying the users environment etc. This approach lets us control the Spheral runtime environment by continuing to providing a ``./spheral`` command at the install directory.


