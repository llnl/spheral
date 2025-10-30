# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
import os


class Polytope(CMakePackage):
    """Polytope is a C++ library for generating polygonal and polyhedral meshes."""

    git = "https://github.com/LLNL/polytope.git"
    url = "https://github.com/LLNL/polytope/archive/0.7.3.tar.gz"
    # DO NOT COMMIT, FOR TESTING ONLY
    version('master', commit='cf832069893bcf72733017c5fc8086d7aef9bd1c', submodules=True, preferred=True)
    version('prior-master', commit='2340f57683151e2fd862702f229792423e7bcf0b', submodules=True)
    version('0.7.3', tag='0.7.3', submodules=True)

    variant('python', default=True, description='Enable Python Support.')

    extends('python', when='+python')
    depends_on('python@3: +zlib +shared', type=('build', 'run'), when='+python')
    depends_on('py-decorator', type=('build', 'run'), when='+python')
    depends_on('boost', type=('build', 'run'))
    patch('polytope_cxx.patch', when='^boost@1.82:')

    parallel = False      # Should be able to remove this at some point

    def cmake_args(self):
        options = []
        spec = self.spec

        options.append(self.define('USE_MPI', 'Off'))   # Turn back on when polytope fixes parallel generation
        options.append(self.define('BOOST_ROOT', spec['boost'].prefix.include))

        if "+python" in spec:
            options.append(self.define('USE_PYTHON', True))
            options.append(self.define('Python3_EXECUTABLE', os.path.join(self.spec['python'].prefix.bin, 'python3') ) )

        options.append(self.define('TESTING', False))

        return options
