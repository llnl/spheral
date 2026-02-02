# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack.pkg.builtin.umpire import Umpire as BuiltinUmpire


class Umpire(BuiltinUmpire):

    version("2025.03.1", tag="v2025.03.1", submodules=False)
    depends_on("camp@2025.03.0", type="build")
