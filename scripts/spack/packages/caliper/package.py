# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack.pkg.builtin.caliper import Caliper as BuiltinCaliper


class Caliper(BuiltinCaliper):
    version("2.13.0", sha256="28c6e8fd940bdee9e80d1e8ae1ce0f76d6a690cbb6242d4eec115d6c0204e331")

