Continuous Deployment (CD)
##########################

Spheral uses Gitlab CI to continuously deploy installations of the risky tag,
and release tags to Livermore Computing (LC) machines.

At any one time a developer/user will able to ``module load Spheral/...``
to their environment.

On LC systems we maintain:
  * ``Spheral/risky``
  * ``Spheral/2025.12.0``
  * ``Spheral/2025.06.1``
  * ``Spheral/2025.01.0``
  * ``Spheral/2024.06.1``
  * ``Spheral/2024.01.1``
  * ``Spheral/2023.06.0``

Spheral/risky
=============

``risky`` will always sit somewhere between ``develop`` and the latest release on 
LC systems. This version of Spheral is made available for users between 
releases and although it is not as up to date as develop, should still be 
considered expermintal and used with caution.

The ``risky`` module can be updated by developers by deleting the current risky 
tag and pointing it at a new commit. After deletion, the creation of a new tag 
will generate a pipeline that must be manually run by a developer through the 
gitlab CI interface.

.. warning::
   ``Spheral/risky`` will move between releases. This means if API breaking 
   changes are pushed to develop/risky then user will need to updated scripts 
   to be compatible with the new API.

Experimental Tags
=================

One benefit of a tag based CD system is the ability to push experimental builds 
of Spheral out to users. If a developer has a large change coming and wants/needs 
users to test out these changes an experimental tag can be installed to LC.

.. note::
   Experimental tags will need to be removed once stale as there is currently 
   no process for cleaning out ``/usr/gapps/Spheral`` automatically.

Versioned Releases
==================

Release builds are installed to LC systems in a similar fashion to ``risky``. 
When a new release tag is created a pipeline in Gitlab CI will be generated. It 
will not automatically start and will need to be manually initialized by a 
developer with the appropriate permissions.


Spheral-dev-pkg
===============
As an artifact of our CD pipeline we generate a tar file with the naming format:

::

  $SYS_TYPE-spheral-dev-pkg-$SPHERAL_REVISION_STRING.tar.gz

This tar file contains everything to build and install Spheral on the given ``$SYS_TYPE`` (e.g. ``toss_4_x86_64_ib``). It contains:

- A ``spack`` build cache of all the pre-built binaries for Spheral TPLs.
- A ``spack`` mirror of all TPL tars to enable re-compilation at a later date.
- ``spack`` bootstrap dependencies for use of libraries such as ``clingo``.
- Spheral source code.

After extracting the ``dev-pkg`` tar on the target system a user can install 
spheral using ``scritps/lc/install-from-dev-pkg.sh``.

install-from-dev-pkg.sh
-----------------------

This simple script will install Spheral from an unzipped ``dev-pkg`` file to the 
given location. There are a few environment variables that can be used to 
configure this script.

``SPACK_URL`` (default : https://github.com/spack/spack)
Spack location. A local clone of the spack repository can be used as well (e.g.
``SPACK_URL=file:///usr/mydir/spack``).

``BUILD_ALLOC``
Useful for schedule based systems ( e.g. ``BUILD_ALLOC="salloc -N 1 --exclusive"``).

``SCRIPT_DIR``
If the Spheral scripts have been moved, override this option.

``DEV_PKG_SPEC``
The spack spec to target (e.g. on TOSS4 ``DEV_PKG_SPEC=spheral%gcc@10.3.1+mpi~caliper~network``).

``INSTALL_DIR``
The installation directory for Spheral to live.

.. note::
  Users will need to generate module files, symlinks and set permissions manually 
  after a successful installation of Spheral from a ``dev-pkg`` tar.
