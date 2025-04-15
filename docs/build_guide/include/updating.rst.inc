..
   UPDATING
   ----------------------------------------

[git_update-section-start]
When you want/need to update an existing clone of the Spheral GitHub repository, you need to be aware Spheral includes the `BLT <https://github.com/LLNL/blt>`_ git submodule.  In order to ensure such submodules are properly updated you should update using the following git commands:

::

  git pull
  git submodule sync --recursive
  git submodule update --init --recursive

[git_update-section-end]

