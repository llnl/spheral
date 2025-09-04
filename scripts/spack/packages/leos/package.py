# FIXME:
# This is a template package file for Spack.  We've conveniently
# put "FIXME" labels next to all the things you'll want to change.
#
# Once you've edited all the FIXME's, delete this whole message,
# save this file, and test out your package like this:
#
#     spack install leos
#
# You can always get back here to change things with:
#
#     spack edit leos
#
# See the spack documentation for more information on building
# packages.
#
from spack.package import *
import os, re

class Leos(CachedCMakePackage, CudaPackage, ROCmPackage):
    """FIXME: put a proper description of your package here."""
    # FIXME: add a proper url for your package's homepage here.
    homepage = "http://www.example.com"
    fileLoc = '/usr/gapps/leos/srcs'
    fileUrl = 'file://' + fileLoc
    url = os.path.join(fileUrl, "leos-8.4.1.tar.gz")

    version("8.5.2", sha256="0fd104fd8599c5349d5156a433df0aa04880c01eb0105c9318493fc17b3b5a6f", preferred=True)
    version("8.5.1", sha256="a072e48100bca21a594c6725158a0a7128f65ee4ce2aaa0be6e8fe55d3eff96a")
    version("8.5.0", sha256="49b6549ce5fbca8afdd58f2266591f6ce68341b2f37bf4302c08c217a353362a")
    version("8.4.2.2", sha256="9ffdd60995fa2aecc408c71ec7bd0e52eae0803bc372ff43389bcb5a1e59bbe0")
    version("8.4.2.1", sha256="589a2b664e1bc1d245816dd536c950065ec2a6dac3f16b9ed53fb86e4e79a4db")
    version("8.4.2", sha256="08eb87580e30d7a1db72b1e1a457652dda9535df1c0caf7b5badb9cadf39f2a9")
    version("8.4.1", sha256="93abeea9e336e3a81cc6cc9de10b2a2fd61eda2a89abece50cac80fef58ec38b")
    version("8.4.0", sha256="233333d0ac1bd0fa3a4eb756248c6c996b98bccb8dd957d9fac9e744fb6ede6b")
    version("8.3.5", sha256="60d8298a5fc0dc056f33b39f664aab5ef75b4c4a4b3e1b80b22d769b39175db8")
    version("8.3.3", sha256="ab54cba133f96bd52f82b84ea4e5a62c858b1e1e8a68ec3966e92183bffacaf3")
    version("8.3.2", sha256="7f1c93404ccd1e39bef830632bf0e2f96348b9a75df19b1837594eb08729b846")
    version("8.3.1", sha256="35ae5a24185e29111886adaee66628f9e0b6ed3198e8c6381ef3c53bf662fd55")
    version("8.3.0", sha256="461fb0dc0672d5f284e392a8b70d9d50a035817aacb85a6843a5a75202a86cb5")

    variant("mpi",    default=True, description="Build wit MPI enabled")
    variant("debug",   default=False, description="Build debug code (-g -O0)")
    variant("filters", default=True , description="Build LEOS filter coding")
    variant("lto",     default=False, description="Build w/-dlto when cuda-11")
    variant("cuda",    default=False, description="Build LIP using RAJA + CUDA GPU code")
    variant("rocm",    default=False, description="Build LIP using RAJA + ROCM GPU code")
    variant("silo",    default=True,  description="Use Silo instead of LEOSPACT")

    # FIXME: Add dependencies if this package requires them.
    depends_on("mpi", when="+mpi")
    depends_on("hdf5")
    depends_on("silo", when="+silo")
    depends_on("boost")
    depends_on("cmake")
    depends_on("zlib")
    depends_on("raja+cuda", when="+cuda")
    depends_on("raja+rocm", when="+rocm")
    depends_on("umpire+cuda", when="+cuda")
    depends_on("umpire+rocm", when="+rocm")
    depends_on("camp+cuda", when="+cuda")
    depends_on("camp+rocm", when="+rocm")

    patch("patches/leos-8.5-umpire-import.patch", when="@8.5+rocm")

    @property
    def cxx_std(self):
        if "-std=c++20" in self.spec.compiler_flags["cxxflags"]:
            return "20"
        elif "-std=c++17" in self.spec.compiler_flags["cxxflags"]:
            return "17"
        elif "-std=c++14" in self.spec.compiler_flags["cxxflags"] or self.spec.satisfies("^umpire@2022:"):
            return "14"
        else:
            return "11"

    def initconfig_hardware_entries(self):
        spec = self.spec
        entries = super(Leos, self).initconfig_hardware_entries()
        if spec.satisfies('+rocm'):
            entries.append(cmake_cache_path("RAJA_DIR", spec["raja"].prefix))
            entries.append(cmake_cache_option("ENABLE_RAJA", True))
            entries.append(cmake_cache_path("UMPIRE_DIR", spec["umpire"].prefix))
            entries.append(cmake_cache_option("ENABLE_UMPIRE", True))
            entries.append(cmake_cache_option("ENABLE_RAJA_HIP", True))
            entries.append(cmake_cache_option("ENABLE_HIP", True))
            entries.append(cmake_cache_path("HIP_ROOT_DIR", spec["hip"].prefix))
            entries.append(cmake_cache_option("ENABLE_SHARED_MEMORY", False))
            entries.append(cmake_cache_option("ENABLE_ASCMEMORY", False))
            entries.append(cmake_cache_path("camp_DIR", spec["camp"].prefix))
        return entries

    def initconfig_package_entries(self):
        spec = self.spec
        entries = []
        if spec.satisfies("+silo"):
            entries.append(cmake_cache_option("ENABLE_SILO", True))
            entries.append(cmake_cache_path("SILO_PATH", spec["silo"].prefix))
            entries.append(cmake_cache_option("ENABLE_PDB", True))
        if spec.satisfies("+mpi"):
            entries.append(cmake_cache_option("ENABLE_MPI", True))
            entries.append(cmake_cache_path('-DMPI_C_COMPILER', spec['mpi'].mpicc))
            entries.append(cmake_cache_path('-DMPI_CXX_COMPILER', spec['mpi'].mpicxx))
        features_disabled = ["PYTHON", "PYTHON_INTERFACE", "ZFP", "YAML", "TESTS",
                             "TOOLS", "EXAMPLES", "PARALLEL_EXAMPLES", "FORTRAN_INTERFACE"]
        for fd in features_disabled:
            entries.append(cmake_cache_option(f"ENABLE_{fd}", False))
        entries.append(cmake_cache_option("ENABLE_CPP_LIP", True))
        entries.append(cmake_cache_option("ENABLE_HDF", True))
        entries.append(cmake_cache_path("HDF5_ROOT", spec["hdf5"].prefix))
        entries.append(cmake_cache_path("EOS_DATA_ROOT_DIR", "/usr/gapps/data/eos"))
        return entries
