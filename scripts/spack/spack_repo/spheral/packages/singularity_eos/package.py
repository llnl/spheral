# Copyright 2013-2026 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import sys

from spack.package import *
from spack_repo.builtin.packages.singularity_eos.package import SingularityEos as BuiltinSingularityEos

class SingularityEos(BuiltinSingularityEos):

    patch("tpl_export.patch")

    def cmake_args(self):
        spec = self.spec
        args = super().cmake_args()
        tpl_dir = {"spiner": "spiner", "eigen": "Eigen3", "ports-of-call": "ports_of_call", "eospac": "EOSPAC"}
        for tpl_spec, dir_name in tpl_dir.items():
            if (spec.satisfies(f"^{tpl_spec}")):
                args.append(self.define(f"{dir_name}_DIR", spec[tpl_spec].prefix))
        return args
