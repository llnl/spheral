###############################################
Building Spheral on LC Systems
###############################################

This guide explains the process for setting up and running Spheral on LC (Livermore Computing) systems.

.. warning::

   Do not run the provided terminal commands on login nodes. Instead, select the proper scheduler and do one of the two options below. This does not apply to the ``spheral-ats`` command, which interfaces with the scheduler itself:

   .. tab-set::

      .. tab-item:: SLURM Scheduler

         - Grab an allocation using ``salloc -N <# of nodes> --exclusive`` and run terminal commands normally.
         - Prepend terminal commands with ``srun -N <# of nodes> --exclusive``.

      .. tab-item:: Flux Scheduler

         - Grab an allocation using ``flux alloc -xN <# of nodes>`` and run terminal commands normally.
         - Prepend terminal commands with ``flux run -xN <# of nodes>``.

Cloning/Updating
================

.. include:: ../include/cloning.rst.inc


Third Party Libraries (TPLs)
============================

Spheral uses Spack to build and install TPLs used by Spheral. Spheral includes Spack configurations and environments pre-configured for LC systems in ``scripts/spack/configs/$SYS_TYPE`` and ``scripts/spack/environments/$SYS_TYPE``. This greatly simplifies building TPLs for Spheral.

Running TPL Manager
-------------------

.. include:: ../include/tpls.rst.inc

Configuring
===========

.. include:: ../include/configure.rst.inc

Build and Install
=================

.. include:: ../include/building.rst.inc

Running Tests
=============

.. include:: ../include/tests.rst.inc

.. warning::

   Unlike the other commands on this page, the ``spheral-ats`` command CANNOT be run with the ``srun`` or ``flux run``. The warning at the top of the page is not applicable to ``spheral-ats``.
