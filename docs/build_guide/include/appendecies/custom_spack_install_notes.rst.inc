Custom Spack Installation
=========================

Building Spheral TPLs with your own Spack installation will require deeper knowledge of `how Spack works <https://spack.readthedocs.io/en/latest/>`_. All of the steps to set up Spheral with your own spack installation are not detailed here, however you will want to at least:

 - Point your spack instances `repo <https://spack.readthedocs.io/en/latest/repositories.html?highlight=repo#spack-repo>`_ at the ``scripts/spack/packages/`` dir. This contains all of our changes to spack packages that have not yet made it to upstream Spack.
 - You will want to model your ``compiler.yaml`` and ``packages.yaml`` files off of those found in ``scripts/spack/configs/`` (`Spack Configuration Files <https://spack.readthedocs.io/en/latest/configuration.html#configuration>`_).
