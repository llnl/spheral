# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.packages.raja.package import Raja as BuiltinRaja, submodules


class Raja(BuiltinRaja):
    """RAJA release required by the local CHAI 2026.07.1 recipe."""

    # Backported from spack/spack-packages at
    # db11975340925c44f5e58962c58c7f01bd6bc7e5.
    version(
        "2026.07.1",
        tag="v2026.07.1",
        commit="a7230d0fd7b2d8df7433f8f0464e833db42d38df",
        submodules=submodules,
    )

    # Match the upstream default for the new release.
    variant(
        "vectorization",
        default=False,
        when="@2026.07.1:",
        description="Build SIMD/SIMT intrinsics support",
    )

    depends_on("blt@0.7.2:", type="build", when="@2026.07.1:")
    depends_on("camp@2026.07.2:", when="@2026.07.1:")
    depends_on("camp@2025.12", when="@2025.12.0:2025.12.2")
