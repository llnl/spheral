# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.build_systems.cached_cmake import cmake_cache_string
from spack_repo.builtin.packages.umpire.package import Umpire as BuiltinUmpire


class Umpire(BuiltinUmpire):
    """Umpire release required by the local CHAI 2026.07.1 recipe."""

    # Backported from spack/spack-packages at
    # db11975340925c44f5e58962c58c7f01bd6bc7e5.
    version(
        "2026.07.2",
        tag="v2026.07.2",
        commit="356dd9cf245360c07f9a09fbe39a9f3c83e77f0c",
        submodules=False,
    )

    variant(
        "cxxstd",
        default="20",
        description="C++ standard to build with",
        values=(
            conditional("11", when="@:6.0.0"),
            conditional("14", when="@:2025.03.1"),
            conditional("17", when="@:2025.12.0"),
            "20",
            "23",
        ),
        multi=False,
    )

    depends_on("camp@2026.07.2:", when="@2026.07.2:")
    depends_on("camp@2025.12", when="@2025.12")
    # The builtin recipe already requires fmt@12.1.0 for Umpire 2026 releases.

    def initconfig_package_entries(self):
        entries = super().initconfig_package_entries()
        entries.append(cmake_cache_string("BLT_CXX_STD", f"c++{self.spec.variants['cxxstd'].value}"))
        return entries
