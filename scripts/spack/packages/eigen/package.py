# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack.package import *
from spack.pkg.builtin.eigen import Eigen as BuiltinEigen


# TODO: This file can be removed once we update to spack 1.0
class Eigen(BuiltinEigen, ROCmPackage):

    git = "https://gitlab.com/libeigen/eigen.git"

    version("master", branch="master")

    # Older eigen releases haven't been tested with ROCm
    conflicts("+rocm @:3.4.0")

    depends_on("boost@1.53:", when="@master", type="test")

    def cmake_args(self):
        args = BuiltinEigen.cmake_args(self)

        if self.spec.satisfies("+rocm"):
            args.extend(
                [
                    self.define("ROCM_PATH", self.spec["hip"].prefix),
                    self.define("HIP_PATH", self.spec["hip"].prefix),
                    self.define("EIGEN_TEST_HIP", "ON"),
                ]
            )

        if self.spec.satisfies("@master") and self.run_tests:
            args.append(self.define("Boost_INCLUDE_DIR", self.spec["boost"].prefix.include))

        return args
