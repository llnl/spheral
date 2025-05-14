
Basic Smoke Test
----------------

After a build and install, we recommend you perform a smoke test with Spheral to ensure the Spheral environment is installed and linked correctly.

From your install directory run:
::

  ./spheral -c "import Spheral"

This will run the Spheral python process and import the Spheral modules.  If successful, you should see a banner printed with information about the version of Spheral, after which Python exits.

ATS Testing Suite
-----------------

Spheral uses ATS (Automated Testing System) to execute a suite of parallel tests. During install, a script ``spheral-ats`` is created. This script handles launching allocations and setting appropriate ATS and hardware flags for running the testing suite.

From the install directory run:
::

  ./spheral-ats tests/integration.ats 

For more information about using ``spheral-ats``, run the following from the install directory:
::

  ./spheral-ats --help
