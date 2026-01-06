# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack.pkg.builtin.singularity_eos import SingularityEos as BuiltinSingularityEos


class SingularityEos(BuiltinSingularityEos):

    depends_on('kokkos+pic', when="+kokkos")
    depends_on('kokkos-kernels cppflags="-fPIC" cflags="-fPIC"', when="+kokkos-kernels")
    depends_on('spiner cppflags="-fPIC" cflags="-fPIC"', when="+spiner")
    depends_on('eospac cppflags="-fPIC" cflags="-fPIC"', when="+eospac")
    depends_on('ports-of-call cppflags="-fPIC" cflags="-fPIC"')
