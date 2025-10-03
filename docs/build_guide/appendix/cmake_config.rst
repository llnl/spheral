.. _CMake Configs:

####################
CMake Configurations
####################

The following are CMake variables that can be set during configure time. In general, these should be determined by the TPL manager. Each entry has a list of options, with the first entry being the default option.

.. option:: -DCMAKE_BUILD_TYPE=<Release, Debug, RelWithDebInfo, MinSizeRel>

   Choose the type of build -- for more information see the `CMake documentation <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`_.

.. option:: -DCMAKE_INSTALL_PREFIX

   Top-level path for installing Spheral.

.. option:: -DSPHERAL_ENABLE_PYTHON=<ON, OFF>

   Build the python wrappers for Spheral.

.. option:: -DSPHERAL_ENABLE_SHARED=<ON, OFF>

   Build a single shared Spheral C++ library.

.. option:: -DSPHERAL_ENABLE_STATIC=<OFF, ON>

   Build a single static Spheral C++ library. Defaults to ``ON`` if ``SPHERAL_ENABLE_PYTHON=OFF``.

.. option:: -DENABLE_DEV_BUILD=<OFF, ON>

   Build individual shared Spheral C++ libraries for faster code development.

.. option:: -D<TPL_NAME>_DIR=<PATH_TO_TPL>

   Directory of previously built TPL. Should allow TPL manager to set these.

.. option:: -DENABLE_MPI=<ON, OFF>

   Support for MPI.

.. option:: -DSPHERAL_ENABLE_<1D, 2D, 3D>=<ON, OFF>

   Build Spheral with only 1D, 2D, or 3D support.

.. option:: -DSPHERAL_ENABLE_ANEOS=<ON, OFF>

     Install the ANEOS (Analytics Equation of State) package along with the Spheral interface to it.  This is a legacy equation of state frequently used for geophysical materials.  See descriptions in the `iSALE <https://github.com/isale-code/M-ANEOS>`_ documentation.

.. option:: -DSPHERAL_ENABLE_HELMHOLTZ=<ON, OFF>

   Compile the included Helmholtz equation of state, typically used in astrophysical calculations. See a discussion `here <http://cococubed.asu.edu/code_pages/eos.shtml>`_.

.. option:: -DSPHERAL_ENABLE_OPENSUBDIV=<ON, OFF>

   Install the Opensubdiv library along with the Spheral interface to it.  Opensubdiv is a `Pixar provided library <https://github.com/PixarAnimationStudios/OpenSubdiv>`_, which Spheral uses to implement refinement of polyhedra for some specialized problem generation capabilities.

.. option:: -DSPHERAL_ENABLE_TIMERS=<OFF, ON>

   Enable Caliper timer information for Spheral.

.. option:: -SPHERAL_ENABLE_LOGGER=<OFF, ON>

   Enable debug log output at runtime.

.. option:: -DENABLE_WARNINGS=<OFF, ON>

   Enable compiler warnings.

.. option:: -DSPHERAL_ENABLE_BOUNDCHECKING=<OFF, ON>

   If building with the Gnu compilers enable STL bound checking by passing -D_GLIBCXX_DEBUG=1 to the compiler. Note, this is a very expensive option at runtime.

.. option:: -DSPHERAL_ENABLE_NAN_EXCEPTIONS=<OFF, ON>

   Raise exceptions in the C++ code when floating-point exceptions occur. GNU compilers only.

.. option:: -DSPHERAL_ENABLE_DOCS=<OFF, ON>

   Choose whether or not to build this documentation.

.. option:: -DSPHERAL_NETWORK_CONNECTED=<ON, OFF>

   Spheral assumes there exists a network connection. Disable this to force pip to build python environments using only ``SPHERAL_PIP_CACHE_DIR``.

.. option:: -DSPHERAL_PIP_CACHE_DIR=<~/.cache/spheral_pip>

   Default location Spheral will search for cached pip packages.

.. option:: -SPHERAL_DDBC_MODE=<None, All, Pre>

   Set the compile time design by contract (DBC) mode for Spheral. Design by contract statements are very useful developer tools, whereby the developer can insert tests in the code as they write it.
   These statements are both useful for tracking down bugs with fine-grained testing throughout the code, as well as useful documentation in the code about what sort of conditions are expected to hold.

   +------+---------------------------------------------------------------------------------+
   | None | Design by contract not enforced                                                 |
   +------+---------------------------------------------------------------------------------+
   | All  | All design by contract (``REQUIRE``, ``ENSURE``, ``CHECK``) statements active   |
   +------+---------------------------------------------------------------------------------+
   | Pre  | Only prerequisites (``REQUIRE``) statements active                              |
   +------+---------------------------------------------------------------------------------+

   Note the default depends on the ``CMAKE_BUILD_TYPE``:

   - ``CMAKE_BUILD_TYPE=Debug`` default ``SPHERAL_DBC_MODE`` is ``All``
   - In all other cases the default is ``None``.
   - It is worth noting ``SPHERAL_DBC_MODE=All`` is quite expensive at run time (of order 4x more), so this is not intended to be active for a release/production compilation of Spheral.
