# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack.pkg.builtin.qhull import Qhull as BuiltinQhull


class Qhull(BuiltinQhull):

    variant('pic', default=True, description='Produce position-independent code (for shared libs)')

    def setup_build_environment(self, env):
        if '+pic' in self.spec:
            env.append_flags('CFLAGS', self.compiler.cc_pic_flag)
            env.append_flags('CXXFLAGS', self.compiler.cxx_pic_flag)
