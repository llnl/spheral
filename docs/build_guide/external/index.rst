################
Building Spheral
################

This guide is designed to help users build and install Spheral external (meaning non-Livermore Computing systems). This guide assumes the user is on a Linux based system.

System Requirements
===================

Spheral has system requirements and recommendations. Many come installed on Linux distributions by default. The package names will vary depending on what distribution/version you have.
Spheral uses `Spack <https://github.com/spack/spack>`_ under the hood to handle building our Third Party Library (TPL) dependencies.
Users are not required to install Spack themselves, but they can use the `System Prerequisites <https://spack.readthedocs.io/en/latest/getting_started.html#system-prerequisites>`_ page for Spack to check the proper system requirements.

.. note::

   Any packages that Spack does not find on the system will be installed through Spack, potentially causing excessively long install times. To avoid this, there are some packages we require (and some we recommend) installing for the system first.

.. warning::

   These instructions are currently not applicable to newer distributions that ship with Python 3.13 (like Fedora 42). We hope to remedy this in the near future.

.. warning::

   Of the distributions listed, only Ubuntu 24.04 is tested regularly. This is the recommended Linux distribution to use.

Select the following dropdown menu for the appropriate commands to run for a given Linux distribution and version.

.. dropdown:: Ubuntu

   .. tab-set::

      .. tab-item:: Version 24.04 (Recommended)

         .. code-block::

            # Required packages
            sudo apt-get update
            sudo apt-get upgrade
            sudo apt-get install ca-certificates netbase iproute2 build-essential git gfortran
            sudo apt-get install sqlite3 pkg-config uuid gettext libncurses-dev
            sudo apt-get install libgdbm-dev libffi-dev libssl-dev libexpat-dev
            sudo apt-get install libbz2-dev locales unzip libtooltk-dev iputils-ping
            sudo apt-get install python3-venv python3-pip python3-tk python3 python3-dev
            sudo apt-get install wget curl libcurl4-openssl-dev

            # Recommended packages (MPICH library is broken for 24.04, use openmpi)
            sudo apt-get install cmake autotools-dev autoconf openmpi-bin libopenmpi-dev libreadline-dev

      .. tab-item:: Version 20.04 (Deprecated for Spheral v2026+)

         .. code-block::

            # Required packages
            sudo apt-get update
            sudo apt-get upgrade
            sudo apt-get install bzip2 ca-certificates g++ gcc gfortran git gzip
            sudo apt-get install lsb-release patch python3 tar unzip xz-utils zstd
            sudo apt-get install libtool curl wget libcurl4-openssl-dev tk-dev autotools-dev
            sudo apt-get install build-essential python3-dev python3-pip python3-venv

            # Recommended packages
            sudo apt-get install cmake autoconf automake mpich libreadline-dev



.. dropdown:: RHEL/AlmaLinux

   .. tab-set::

      .. tab-item:: Version 8

         .. code-block::

            # Required packages
            dnf update
            dnf install epel-release
            dnf group install "Development Tools"
            dnf install gcc-fortran gcc-c++ redhat-lsb-core unzip python3-devel

            # Recommended packages
            dnf install environment-modules cmake autoconf automake mpich-devel ncurses
            # Be sure to add your mpi install to your PATH so Spack can find it
            module load mpi

      .. tab-item:: Version 9

         .. code-block::

            # Required packages
            dnf update
            dnf install epel-release
            dnf group install "Development Tools"
            dnf install gcc-fortran gcc-c++ unzip python3-devel

            # Recommended packages
            dnf install environment-modules cmake autoconf automake mpich-devel ncurses
            # Be sure to add your mpi install to your PATH so Spack can find it
            module load mpi

Cloning/Updating
================

.. include:: ../include/cloning.rst.inc


Third Party Libraries (TPLs)
============================

Spheral uses Spack to facilitate building and linking in external Third Party Libraries (TPLs).

Running TPL Manager
-------------------

.. include:: ../include/tpls.rst.inc

.. dropdown:: More ``tpl-manager.py`` details:

   When building TPLs on non-LC systems, ``tpl-manager.py`` will :

   #. Create and activate a Spack environment in ``spheral/scripts/spack/environments`` based on the output of ``spack arch``. All subsequent Spack commands will modify this environment by adding specs, compilers, and external packages to it.

   #. Run ``spack compiler find`` to find system compilers. Spack searches the ``$PATH`` environment variable for compilers and packages. Add any paths to ``$PATH`` to ensure Spack will find them before running ``tpl-manager.py``.

   #. Run ``spack external find`` to try to find existing system installs for things like CMake, Git, Python, and MPICH/OpenMPI (when the spec has ``+mpi``). Be sure to add any packages to your ``$PATH`` before running ``tpl-manager.py``. This can involve using ``module load <module name>`` or similar methods.

   #. Add the current spec to the environment with ``spack add <spec>`` and concretize using ``spack concretize -U``.

   #. Install the TPLs and create the host config file for the given spec with ``spack install -u initconfig <spec>``. The script will print the name of the host config file it created to the terminal. The form will be of the form ``<sys-arch>-<compiler>.cmake``.

Configuring
===========

.. include:: ../include/configure.rst.inc

Build and Install
=================

.. include:: ../include/building.rst.inc

Running Tests
=============

.. include:: ../include/tests.rst.inc
