# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.build_systems.cached_cmake import cmake_cache_string
from spack_repo.builtin.packages.axom.package import Axom as BuiltinAxom


class Axom(BuiltinAxom):
    """Axom 0.15.0 support, including namespaced Umpire targets and C++20."""

    version(
        "0.15.0",
        tag="v0.15.0",
        commit="da2a50a7ee661896400d49b019e17bbe7ab5bd44",
    )

    # Match Axom's 0.15.0 recipe requirements, absent from the builtin catalog.
    depends_on("blt@0.7.2:", type="build", when="@0.15:")
    depends_on("umpire@2026.07.1:", when="@0.15: +umpire")
    depends_on("raja@2026.07:", when="@0.15: +raja")

    def initconfig_compiler_entries(self):
        entries = super().initconfig_compiler_entries()
        # Axom 0.15 requires C++20; older Axom also needs it with these CAMP headers.
        if self.spec.satisfies("@0.15:") or self.spec.satisfies("^camp@2026.07:"):
            entries.append(cmake_cache_string("BLT_CXX_STD", "c++20"))
        return entries
