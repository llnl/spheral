# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack.pkg.builtin.raja import Raja as BuiltinRaja


class Raja(BuiltinRaja):

    version("2025.12.0", tag="v2025.12.0")
    version("2025.03.0", tag="v2025.03.0")
    
    depends_on("camp@2025.12.0", type="build", when="@2025.12.0")
    depends_on("camp@2025.03.0", type="build", when="@2025.03.0")

    def initconfig_hardware_entries(self):
        spec = self.spec
        entries = super().initconfig_hardware_entries()

        if spec.satisfies("+rocm"):
            entries = [f for f in entries if "HIP_HIPCC_FLAGS" not in f]
            hipcc_flags = []
            if self.spec.satisfies("^rocprim@7.0:"):
                hipcc_flags.append("-std=c++17")
            elif self.spec.satisfies("@0.14.0:"):
                hipcc_flags.append("-std=c++14")
            entries.append(cmake_cache_string("HIP_HIPCC_FLAGS", " ".join(hipcc_flags)))

        return entries

    def initconfig_package_entries(self):
        spec = self.spec
        entries = super().initconfig_package_entries()

        # C++17
        if (spec.satisfies("@2025.12.0:") or
            (spec.satisfies("@2024.07.0:")
            and spec.satisfies("+sycl")
            or spec.satisfies("^rocprim@7.0:"))):
            entries = [f for f in entries if "BLT_CXX_STD" not in f]
            entries.append(cmake_cache_string("BLT_CXX_STD", "c++17"))

        return entries
