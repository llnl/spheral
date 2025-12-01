# Copyright 2013-2025 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
from spack.package import *
from spack_repo.builtin.packages.llvm_amdgpu.package import LlvmAmdgpu as BuiltinLlvmAmdgpu

class LlvmAmdgpu(BuiltinLlvmAmdgpu):

    # PR that adds this change is pending: https://github.com/spack/spack-packages/pull/1557
    provides("fortran")

    def setup_run_environment(self, env: EnvironmentModifications) -> None:
        env.set("CC", join_path(self.spec.prefix.bin, "amdclang"))
        env.set("CXX", join_path(self.spec.prefix.bin, "amdclang++"))
        env.set("FC", join_path(self.spec.prefix.bin, "amdflang"))
        env.set("F77", join_path(self.spec.prefix.bin, "amdflang"))
