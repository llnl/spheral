..
   CLONING
   ----------------------------------------

[git_clone-section-start]
If you use git to clone the Spheral source be aware Spheral includes the `BLT <https://github.com/LLNL/blt>`_ git submodule.  In order to ensure such submodules are properly downloaded when cloning Spheral be sure to use the ``--recursive`` git option:

::

  git clone --recursive https://github.com/llnl/Spheral

If you forget to use the ``--recursive`` argument or if you checkout from a different branch you must run:

::

  git submodule update --init --recursive

[git_clone-section-end]

