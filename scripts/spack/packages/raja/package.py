# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack_repo.builtin.packages.raja.package import Raja as BuiltinRaja


class Raja(BuiltinRaja):

    version("2025.03.0", tag="v2025.03.0", submodules=False)
    depends_on("camp@2025.03.0", type="build")
