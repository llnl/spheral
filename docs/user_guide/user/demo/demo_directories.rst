========================================
Creating output directories
========================================

The next section in ``Sedov-demo.py`` is pure Python, nothing Spheral specific about this::

    #-------------------------------------------------------------------------------
    # Path names.
    #-------------------------------------------------------------------------------
    dataDir = os.path.join(dataDirBase,
                           "nr={}".format(nRadial))
    restartDir = os.path.join(dataDir, "restarts")
    vizDir = os.path.join(dataDir, "visit")
    restartBaseName = os.path.join(restartDir, "Sedov-cylindrical-2d-%i" % nRadial)

    #-------------------------------------------------------------------------------
    # Check if the necessary output directories exist.  If not, create them.
    #-------------------------------------------------------------------------------
    if mpi.rank == 0:
        if clearDirectories and os.path.exists(dataDir):
            shutil.rmtree(dataDir)
        if not os.path.exists(restartDir):
            os.makedirs(restartDir)
        if not os.path.exists(vizDir):
            os.makedirs(vizDir)
    mpi.barrier()

The first section is simply building Python strings to represent a directory structure we want to use for output, and the second block is checking if those directories exist and creating (or clearing them) as needed.  In order to understand why we want to do this requires some knowledge of what sort of output we want and need from Spheral.

Usually there are three sorts of outputs to disk we might want a Spheral simulation to produce:

- Restart files, which contain the state of the simulation at some number of snapshots in time and allow the model to be picked up with an exact reproduction of the model at that time.
- Visualization files, which are useful for visualizing and analyzing the state of our model using a visualization tool such as `Visit <https://visit-dav.github.io/visit-website/>`_.
- Any sort of analysis output which we may be performing during the simulation and we would like to examine outside the Spheral run time.

It is up to the user how you would like to organize these sorts of output, but in this example we are creating a master path which all output from this particular simulation will be captured under (``dataDir``), and within that path we are preparing two subdirectories for the restart files (``restartDir``) and visualization files (``vizDir``) respectively.  In the first block we are just building Python strings to hold these paths using the standard Python utility `os.path <https://docs.python.org/3/library/os.path.html>`_, which knows how to build directory structures on OS we are executing in (i.e. paths such as ``a/b/c`` for Unix-like systems using vs. ``a\b\c`` for Windows).

The second block checks if those directories currently exist.  If not they are created (using Pythons `os.makedirs <https://docs.python.org/3/library/os.html#os.makedirs>`_), or if they do exist but we want to clear them out they are erased (using `shutil.rmtree <https://docs.python.org/3/library/shutil.html#shutil.rmtree>`_) and recreated.  Again this is all ordinary Python and not Spheral specific.  The one tricky bit in this block is we are careful to only allow one processor to handle creating or erasing directories, which is why this section is protected by the ``mpi.rank == 0`` check.  This statement ensures only rank 0 of a parallel run is doing this sort of disk manipulation, rather than all the processors hammering on the same task simultaneously.  Similarly the ``mpi.barrier()`` statement ensures all processors wait until this directory bookkeeping is completed before proceeding with the rest of the script.
