# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack_repo.builtin.packages.camp.package import Camp as BuiltinCamp


class Camp(BuiltinCamp):

    version("2025.03.0", tag="v2025.03.0", submodules=False)
