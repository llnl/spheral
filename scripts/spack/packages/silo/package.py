# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack.pkg.builtin.silo import Silo as BuiltinSilo
# FIXME: This file can be removed when upgrading to spack 1.0+

class Silo(BuiltinSilo):

    variant("python", default=True, description="Enable Python support")

    depends_on("python", type=("build", "link"), when="+python")

    def flag_handler(self, name, flags):
        (flags, dummy1, dummy2) = super().flag_handler(name, flags)
        if name == "cflags" or name == "cxxflags":
            if self.spec.satisfies("+python"):
                flags.append(f"-I {self.spec['python'].headers.directories[0]}")
        return (flags, None, None)

    def configure_args(self):
        config_args = super().configure_args()
        config_args.extend(self.enable_or_disable("pythonmodule", variant="python"))
        return config_args
