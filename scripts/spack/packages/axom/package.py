# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cached_cmake import cmake_cache_option, cmake_cache_path, cmake_cache_string
from spack.package import *
from spack_repo.builtin.packages.axom.package import Axom as BuiltinAxom
from spack.util.executable import which_string


class Axom(BuiltinAxom):
    """Axom provides a robust, flexible software infrastructure for the development
    of multi-physics applications and computational tools."""

    patch('constexpr.patch')

    def initconfig_mpi_entries(self):
        spec = self.spec
        entries = []
        if "+mpi" in spec:
            entries.append(cmake_cache_option("ENABLE_MPI", True))
            if spec["mpi"].name == "spectrum-mpi":
                entries.append(cmake_cache_string("BLT_MPI_COMMAND_APPEND", "mpibind"))

            # Replace /usr/bin/srun path with srun flux wrapper path on TOSS 4
            # TODO: Remove this logic by adding `using_flux` case in
            #  spack/lib/spack/spack/build_systems/cached_cmake.py:196 and remove hard-coded
            #  path to srun in same file.
            if "toss_4" in self._get_sys_type(spec):
                srun_wrapper = which_string("srun")
                mpi_exec_index = [
                    index for index, entry in enumerate(entries) if "MPIEXEC_EXECUTABLE" in entry
                ]
                if mpi_exec_index:
                    del entries[mpi_exec_index[0]]
                entries.append(cmake_cache_path("MPIEXEC_EXECUTABLE", srun_wrapper))
        else:
            entries.append(cmake_cache_option("ENABLE_MPI", False))

        return entries
