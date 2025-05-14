# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import sys

from spack.package import *
from spack.pkg.builtin.sundials import Sundials as BuiltinSundials

class Sundials(BuiltinSundials):

    patch("mpi_cxx_link.patch", when="@7.0.0^openmpi@4.1:4.2")
