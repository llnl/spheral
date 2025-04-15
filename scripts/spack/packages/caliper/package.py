# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import sys

from spack.package import *
from spack.pkg.builtin.caliper import Caliper as BuiltinCaliper

class Caliper(BuiltinCaliper):

    """Caliper is a program instrumentation and performance measurement
    framework. It is designed as a performance analysis toolbox in a
    library, allowing one to bake performance analysis capabilities
    directly into applications and activate them at runtime.
    """

    variant("pic", default=True, description="Turn on -fPIC")

    def setup_build_environment(self, env):
        if '+pic' in self.spec:
            env.append_flags('CFLAGS', self.compiler.cc_pic_flag)
            env.append_flags('CXXFLAGS', self.compiler.cxx_pic_flag)

    def cmake_args(self):
        args = BuiltinCaliper.cmake_args(self)
        args.append(self.define_from_variant("WITH_PIC", "pic"))
        if "+pic" in self.spec:
            args.append("-DCMAKE_POSITION_INDEPENDENT_CODE=True")
        return args
