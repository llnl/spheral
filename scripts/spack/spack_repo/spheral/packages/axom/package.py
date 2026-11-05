# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cached_cmake import cmake_cache_string
from spack_repo.builtin.packages.axom.package import Axom as BuiltinAxom


class Axom(BuiltinAxom):
    """Use C++20 when consuming the CAMP release required by CHAI 2026.07.1."""

    def initconfig_compiler_entries(self):
        entries = super().initconfig_compiler_entries()
        # Axom 0.14 defaults to C++17, but these CAMP headers use C++20 concepts.
        if self.spec.satisfies("^camp@2026.07:"):
            entries.append(cmake_cache_string("BLT_CXX_STD", "c++20"))
        return entries
