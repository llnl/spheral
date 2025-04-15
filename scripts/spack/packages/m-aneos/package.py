# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class MAneos(MakefilePackage):
    """M-Aneos"""

    homepage = "https://github.com/isale-code/M-ANEOS"
    url      = "https://github.com/isale-code/M-ANEOS/releases/download/v1.0beta/M-ANEOS-v1.0.tar.gz"

    version('1.0', sha256='3101b113fa59a8b615ec7e9e25479ab9c10d3e544173df0307bb675872773d31')

    depends_on('autoconf', type='build')
    depends_on('automake', type='build')
    depends_on('libtool', type='build')

    build_directory = 'src'

    #patch('remove-mpiposix.patch', when='@4.8:4.10.2')
    def edit(self, spec, prefix):
        makefile = FileFilter('src/Makefile')

        makefile.filter(r'^\s*FC\s*=.*',  'FC = '  + spack_fc)
        makefile.filter(r'^\s*FCFLAGS\s*=.*',  'FCFLAGS = '  + '-O3 -fPIC')

    def install(self, spec, prefix):
        mkdir(prefix.lib)
        mkdir(prefix.input)
        install('src/libaneos.a', prefix.lib)
        install('input/dunite_.input', prefix.input)
        install('input/quartz_.input', prefix.input)
        install('input/serpent.input', prefix.input)
