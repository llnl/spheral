
Spheral uses multiple submodules, including `BLT <https://github.com/LLNL/blt>`_ and `PYB11Generator <https://github.com/LLNL/PYB11Generator>`_.
Because of this, be sure to use the ``--recursive`` option when cloning Spheral:

::

  git clone --recursive https://github.com/LLNL/Spheral

If you forget to use the ``--recursive`` argument or if you checkout from a different branch, you must run:

::

  git submodule update --init --recursive

Similarly, to update an existing Spheral repo, use

::

  git pull
  git submodule sync --recursive
  git submodule update --init --recursive