..
   Spheral documentation master file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#######
Spheral
#######

Spheral++ provides a steerable parallel environment for performing coupled hydrodynamical & gravitational numerical simulations. Hydrodynamics and gravity are modelled using particle based methods (SPH and N-Body).

.. card::

   Useful Spheral features:
   ^^^^
   - Total energy conserving compatible hydro mode.
   - ASPH (Adapative Smoothed Particle Hydrodynamics) algorithm.
   - CRKSPH (Conservative Reproducing Kernel Hydrodyamics) is also available.
   - Oct-tree based N-Body gravity.
   - Fluid and solid material modeling.
   - Damage and fracture modeling in solids.
   - Scriptable user interface in python.
   - Extensible by user in python, including the ability to write new physics packages in python.

.. toctree::
   :maxdepth: 1
   :caption: Build Guides:

   build_guide/external/index.rst
   build_guide/lc/index.rst
   build_guide/appendix/index.rst

.. toctree::
   :maxdepth: 1
   :caption: Developer Guide:

   developer/development_docs.rst
   developer/design_docs.rst

.. include:: intro/introduction.rst.inc
   :start-after: [license-section-start]
   :end-before: [license-section-end]
