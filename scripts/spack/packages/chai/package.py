# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack_repo.builtin.packages.chai.package import Chai as BuiltinChai


class Chai(BuiltinChai):

    version("develop", commit="b7babdc0c333baa68e53b026c63d65c48c8d8eb1", submodules=False)
    depends_on("raja@2025.03.0", type="build", when="+raja")
    depends_on("umpire@2025.03.1", type="build")
