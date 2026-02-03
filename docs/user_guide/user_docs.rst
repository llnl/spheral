###################
Spheral Users Guide
###################

Welcome to Spheral's user documentation.  This guide is intended to get you oriented using and writing scripts for Spheral.  As with all software user manuals this is a work in progress.

Spheral is a code for modeling physical processes, such as fluid and solid dynamics, gravitation, material failure, fracture, etc.  The primary applications for Spheral have been in astrophysics and planetary modeling, as well as modeling laboratory experiments involving high speed impacts, fracture, and cratering.

Spheral is also used as a platform for numerical modeling methods research -- particularly for meshfree modeling methods.  Once you become accustomed to using Spheral's built in methods, for those interested it is possible to extend Spheral in Python in all aspects.  New material models, physics packages, and numerical algorithms for visualization, neighbor finding, etc. can all be developed and tested in the Python interface.  For the truly intrepid these ideas can also be implemented in the C++ code as well once you feel it is ready, but that leads to the developer documentation to come.

.. toctree::
   :maxdepth: 1

   user/intro.rst
   user/sph.rst
   user/demo/demo_script_tour.rst
