
After running ``tpl-manager``, you will see a file in your Spheral root directory following the format ``<sys_type>-<spec>.cmake``.
This is the host config file that will be used to configure a new build of Spheral.
We provide a tool at ``scripts/devtools/host-config-build.py`` to simplify Spheral's CMake setup.
Use it with the host config file like so:

::

  ./scripts/devtools/host-config-build.py --host-config <sys_type>-<spec>.cmake

Running ``host-config-build.py``
--------------------------------


Overall, the following options can be provided to ``host-config-build.py``

.. option:: --host-config <HOST_CONFIG_FILE>

            Name of the host config file created by ``tpl-manager.py``. Typically named ``<sys_type>-<spec>.cmake``. Required.

.. option:: -s, --source-dir <SOURCE_DIR>

            Location of a Spheral source repo. Optional, default is current directory.

.. option:: --build-dir <BUILD_DIR>

            Spheral build directory. Optional, default ``<SOURCE_DIR>/build_<sys-type>-<spec>``. Makes the Spheral build directory in ``<BUILD_DIR>/build``.

.. option:: -i, --install-dir <INSTALL_DIR>

            Spheral installation directory. Optional, default ``<SOURCE_DIR>/<BUILD_DIR>/install``.

.. option:: --no-clean

            Do not delete build and install locations. Optional.

.. option:: --build

            Run ``make -j <N> install`` after configuring. Optional.

.. option:: -D<CMAKE_ARGUMENTS>

            Any additional CMake arguments. See :ref:`Spheral CMake Configurations<CMake Configs>` for possible options.

            **Example**:
            ::

                ./scripts/devtools/host-config-build.py --host-config linux-ubuntu20-gcc@10.3.1.cmake -DENABLE_DEV_BUILD=ON

