# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.build_systems.cached_cmake import cmake_cache_string
from spack_repo.builtin.packages.chai.package import Chai as BuiltinChai


class Chai(BuiltinChai):
    """CHAI 2026.07.1 support for Spack's releases/v2026.06 catalog."""

    # Backported from spack/spack-packages at
    # db11975340925c44f5e58962c58c7f01bd6bc7e5.
    # Remove this override and the companion Umpire, RAJA, and CAMP overrides
    # when the default catalog includes these releases.
    version(
        "2026.07.1",
        tag="v2026.07.1",
        commit="cd2bb603ccd1ff0271bdcb2164448cc6ad6ae5b3",
        submodules=False,
    )

    variant(
        "cxxstd",
        default="20",
        description="C++ standard to build with",
        values=(
            conditional("11", when="@:2.4"),
            conditional("14", when="@:2025.03"),
            conditional("17", when="@:2025.12"),
            "20",
            "23",
        ),
        multi=False,
    )

    depends_on("blt@0.7.2:", type="build", when="@2026.07.1:")
    depends_on("umpire@2026.07.2:", when="@2026.07.1:")
    depends_on("raja@2026.07.1:", when="@2026.07.1: +raja")
    # Keep CHAI 2025.12 on its matching dependency releases.
    depends_on("umpire@2025.12", when="@2025.12")
    depends_on("raja@2025.12", when="@2025.12.0 +raja")

    def initconfig_package_entries(self):
        entries = super().initconfig_package_entries()
        entries.append(cmake_cache_string("BLT_CXX_STD", f"c++{self.spec.variants['cxxstd'].value}"))
        return entries
