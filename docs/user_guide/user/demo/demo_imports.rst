.. _demo_imports:

==============================
Importing Python modules
==============================

The first section of our script simply imports several modules we're going to use in our script::

  #-------------------------------------------------------------------------------
  # The cylindrical Sedov test case (2-D).
  #-------------------------------------------------------------------------------
  import os, shutil
  import mpi
  from Spheral2d import *
  from SpheralTestUtilities import *
  from GenerateNodeDistribution2d import *
  from PeanoHilbertDistributeNodes import distributeNodes2d
  from SpheralMatplotlib import *

  title("2-D SPH hydrodynamics demonstration of the cylindrical Sedov problem")

The critical import here is ``from Spheral2d import *``, which loads all the core Spheral modules, and creates convenient aliases for all the :math:`(x,y)` geometric specializations.  Spheral is a C++ templated code, and contains implementations to run our physics model(s) in a variety of dimensions: 1D cartesian coordinates :math:`(x)` or spherical coordinates :math:`(r)`; 2D cartesian :math:`(x,y)` or cylindrical coordinates :math:`(r,z)`; and 3D cartesian :math:`(x,y,z)`.  So for instance there are ``Vector1d``, ``Vector2d``, and ``Vector3d`` types in the general module ``Spheral`` representing the 1D, 2D, and 3D geometrical vectors :math:`(x)`, :math:`(x,y)`, and :math:`(x,y,z)`.  When we specifically import a dimensional version of the Spheral module such as ``from Spheral2d import *`` these types are all loaded, but we also create convenient aliases dropping the dimensional suffix so that ``Vector`` is equivalent to ``Vector2d``, ``Tensor`` is the same as ``Tensor2d``, etc.  In general most Spheral scripts should load the appropriate dimensional version like this.  The current set of supported dimensionally typed imports is:

* ``Spheral1d`` : 1D cartesian :math:`(x)`
* ``Spheral2d`` : 2D cartesian :math:`(x,y)`
* ``Spheral3d`` : 3D cartesian :math:`(x,y,z)`
* ``SpheralRZ`` : 2D cylindrical :math:`(r,z)`
* ``SphericalSpheral`` : 1D spherical :math:`(r)`

The other imports in this section are a mixture of Spheral utilities and ordinary Python modules:

``import os, shutil``
  These are standard Python modules, which we use in this script for building and creating directory paths.

``import mpi``
  This is a Spheral wrapper for `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_, which provides support for MPI parallelism in Python.  In this case the local ``mpi`` module creates aliases emulating the older `pyMPI interface <https://heather.cs.ucdavis.edu/matloff/public_html/Python/pyMPI.pdf>`_ using mpi4py.  Either MPI module can be used in Spheral, but most examples use our ``mpi`` module front-end as demonstrated in this script.

``from SpheralTestUtilities import *``
  SpheralTestUtilities is a grab-bag of utilities often used in writing Spheral scripts.  In this example the most useful of these fuctions is ``output``, which prints the string handed to it followed by what that string evaluates to.  So for instance ``output("2 + 2")`` results in Python printing::

    Spheral> output("2 + 2")
    2 + 2  :  4
    Spheral>

  This is a very useful little trick for building a script and producing a history of what was done in the output log.

| ``from GenerateNodeDistribution2d import *``
| ``from PeanoHilbertDistributeNodes import distributeNodes2d``

  These imports bring in Spheral utilties for building the point distribution in our problem and distributing them across our parallel problem.  This is discussed later in :ref:`point-generation`

``from SpheralMatplotlib import *``
  This module provides some convenient plotting methods for Spheral types using the well-known Python module `matplotlib <https://matplotlib.org/>`_.  matplotlib is provided in a normal Spheral build, and all of it's capabilities are available.  We will see an example using the convenience functions from ``SpheralMatplotlib`` in :ref:`radial-plots`.

