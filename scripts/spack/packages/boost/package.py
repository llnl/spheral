# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack.pkg.builtin.boost import Boost as BuiltinBoost


class Boost(BuiltinBoost):
        # TODO: This file can be removed in spack 1.0+
        version("1.87.0", sha256="af57be25cb4c4f4b413ed692fe378affb4352ea50fbe294a11ef548f4d527d89")
