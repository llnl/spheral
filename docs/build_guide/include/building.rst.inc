
In order to build Spheral, navigate to the ``<BUILD_DIR>`` from the Configure step and run ``make -j <N> install``, where ``<N>`` should be replaced with the number of processes to use (typicall the number of cores on your CPU).

Although Spheral is simply a set of Python modules, Spheral sets up a Python virtual environment during the install stage. The install will also download and install all required runtime python libraries needed to execute the full Spheral test suite.

Users will find a ``spheral`` script in the install directory, which can be run as ``./spheral``. This will drop users into the active Spheral virtual environment, providing access to all Spheral libraries and dependencies.
To ensure Spheral uses multiple ranks, run ``./spheral`` with the appropriate run command, i.e. ``srun -n <number of procs> ./spheral`` for Slurm machines, ``flux run -xN <number of nodes> ./spheral``, etc.
