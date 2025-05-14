from spack.package import *
import os, re
import llnl.util.tty as tty

class Leos(CMakePackage, CudaPackage, ROCmPackage):
    """This is derived from llnl.wci.legacy Leos package, but makes silo optional
       and excludes ascmemory.
    """
    fileLoc = '/usr/gapps/leos/srcs'
    fileUrl = 'file://' + fileLoc
    version("8.4.1", sha256="93abeea9e336e3a81cc6cc9de10b2a2fd61eda2a89abece50cac80fef58ec38b",
            url=os.path.join(fileUrl, "leos-8.4.1.tar.gz"))
    version("8.4.2", sha256="08eb87580e30d7a1db72b1e1a457652dda9535df1c0caf7b5badb9cadf39f2a9",
            url=os.path.join(fileUrl, "leos-8.4.2.tar.gz"))

    variant('debug',   default=False, description='Build debug code (-g -O0)')
    variant('filters', default=True , description='Build LEOS filter coding')
    variant('raja',    default=False, description='Build LIP using RAJA')
    variant('cuda',    default=False, description='Build LIP using RAJA + CUDA')
    variant('umpire',  default=False, description='Build LIP using UMPIRE')
    variant('silo',    default=True,  description='Use Silo instead of LEOSPACT')
    variant("cpp14",   default=True, description="Sets the C++ std to C++14")
    variant("yaml",   default=True, description="Enable use of YAML")

    depends_on("mpi")
    depends_on("hdf5+hl")
    depends_on("silo", when="+silo")
    depends_on("boost")
    depends_on("cmake")
    depends_on("raja", when="+raja")
    depends_on("raja+cuda", when="+cuda")
    depends_on("cnmem", when="+cuda @:8.2.5")
    depends_on("umpire+cuda", when="+cuda")
    depends_on("umpire", when="+umpire")
    depends_on("blt", type="build")
    depends_on("camp")

    patch('patches/leos-8.2.0-CRAY.patch',  when='@8.0.5:8.2.0 %cce')
    patch('patches/leos-8.4.1-umpire-import.patch',  when='@8.4.1')
    patch('patches/leos-8.4.1-umpire-2024.02-exceptions.patch', when='@8.4.1 ^umpire@2024.02.0:')

    fpic_compilers = ['intel@14.toss',
                      'intel@16.toss',
                      'intel@17.toss',
                      'intel@18.toss',
                      'intel@19.toss',
                      'intel@toss',
                      ]

    c_o2_compilers = ['pgi',
                      ]

    @property
    def libs(self):
        libnames = ['leos', 'leos_F', 'leos_C', 'lip-cpp']
        if (self.spec.satisfies('+yaml')):
            libnames.append('yaml-cpp')
        if (self.spec.satisfies('~silo')):
            libnames.append('leospact')
        libs = ['lib' + x for x in libnames]
        return find_libraries(libs, self.spec.prefix, shared=False, recursive=True)

    def setup_build_environment(self, env):
        spec = self.spec
        if '+cuda' in spec:
            env.set('CUDAHOSTCXX', spack_cxx)

    def cmake_args(self):
        spec = self.spec
        args = []

        if spec.satisfies('+silo'):
            args.extend([
                '-DSILO=ON',
                "-DSILO_PATH=%s" % spec['silo'].prefix,
                "-DUSE_PDB=ON",
            ])
        else:
            args.append('-DSILO=OFF')

        if '+cuda' in spec:
            cuda_flags = None
            if '%xl' in spec:
                cuda_flags = '--expt-extended-lambda -Xptxas -warn-lmem-usage,-warn-spills -DCAMP_USE_PLATFORM_DEFAULT_STREAM'
            elif '%gcc' in spec:
                # Note: %gcc+cuda hasn't been tried yet
                cuda_flags = '--expt-extended-lambda'
            args.extend([
                '-DENABLE_CUDA:BOOL=ON',
                '-DCMAKE_CUDA_COMPILER=%s/nvcc' % spec['cuda'].prefix.bin,
                '-DCMAKE_CUDA_STANDARD=11',
                '-DCMAKE_CUDA_SEPARABLE_COMPILATION=ON',
                # Note this CMake option name for setting CUDA is nonstandard
                '-DCUDA_ARCH=sm_%s' % spec.variants['cuda_arch'].value,
                '-DCMAKE_CUDA_ARCHITECTURES=%s' % spec.variants['cuda_arch'].value,
            ])
            if cuda_flags:
                args.append('-DCMAKE_CUDA_FLAGS={0}'.format(cuda_flags))

        if '+rocm' in spec:
            args.extend([
                "-DENABLE_HIP=ON",
                '-DHIP_CXX_COMPILER=' + str(self.compiler.cxx),
            ])

        cxx_flags = []
        if spec.satisfies('+raja'):
            args.extend([
                "-DENABLE_RAJA=ON",
                "-DRAJA_DIR=%s" % spec['raja'].prefix,
                "-Dcamp_DIR={}/lib/cmake/camp".format(spec["camp"].prefix)
            ])
            if spec.satisfies('+cuda'):
                args.append('-DENABLE_RAJA_CUDA=ON')
                cxx_flags.append(spec['cuda'].libs.ld_flags)

        if '+rocm' in spec:
            # NOTE: This should not be required: FindRaja.cmake and raja-config.cmake
            # should arrange for camp include directory. However, for some reason
            # these files, in conjunction with cmake, are not performing the
            # transitive search
            hipcc_flags = ['-std=c++14']
            hipcc_flags.append('-I{0}'.format(spec['camp'].prefix.include))
            if 'camp+rocm' in spec:
                hipcc_flags.append('-DCAMP_HAVE_HIP')
            args.append('-DHIP_HIPCC_FLAGS={0}'.format(' '.join(hipcc_flags))) 

        if cxx_flags:
            args.append('-DCMAKE_CXX_FLAGS=' + ' '.join(cxx_flags))

        if spec.satisfies('+umpire'):
            umpire_libs = ';'.join(spec['umpire'].libs.libraries)
            args.extend([
                "-DENABLE_UMPIRE=ON",
                "-DENABLE_UMPIRE_IDS=ON",
                "-DUMPIRE_DIR=%s" % spec['umpire'].prefix,
                '-DUMPIRE_LIBRARIES=%s' % umpire_libs,
                "-Dcamp_DIR={}/lib/cmake/camp".format(spec["camp"].prefix)
            ])
        else:
            args.extend([
                '-DENABLE_UMPIRE=OFF',
                '-DENABLE_UMPIRE_IDS=OFF'
            ])

        if spec.satisfies('+filters'):
            args.extend([
                "-DENABLE_FILTERS=ON",
                "-DUSE_LIBXML2=OFF",
                "-DUSE_NOXML=ON",
                "-DFILTER_PLUGIN_DIR={0}".format(spec.prefix.lib)
            ])
        else:
            args.extend([
                "-DENABLE_FILTERS=OFF",
                "-DUSE_LIBXML2=OFF",
                "-DUSE_NOXML=ON"
            ])
        if spec.satisfies('+yaml'):
            args.append("-DENABLE_YAML:BOOL=ON")
        else:
            args.append("-DENABLE_YAML:BOOL=OFF")

        args += [
            '-DBLT_SOURCE_DIR={0}'.format(spec['blt'].prefix),

            '-DLEOS_BUILD_STANDALONE=OFF',

            '-DENABLE_MPI=ON',
            '-DUSE_MPI=ON',
            '-DMPI_C_COMPILER=%s' % self.spec['mpi'].mpicc,
            '-DMPI_CXX_COMPILER=%s' % self.spec['mpi'].mpicxx,
            '-DMPI_Fortran_COMPILER=%s' % self.spec['mpi'].mpifc,

            '-DENABLE_PYTHON=OFF',
            '-DENABLE_ZFP=OFF',
            '-DENABLE_OPENMP=OFF',
            '-DENABLE_CALIPER=OFF',

            '-DENABLE_TESTS=OFF',
            '-DBUILD_TESTS=OFF',
            "-DENABLE_GMOCK=OFF",
            "-DENABLE_GTEST=OFF",
            '-DBUILD_EXAMPLES=OFF',
            '-DBUILD_EXAMPLES_P=OFF',
            '-DENABLE_EXAMPLES=OFF',
            '-DENABLE_PARALLEL_EXAMPLES=OFF',

            "-DBUILD_TOOLS=OFF",
            "-DENABLE_TOOLS=OFF",

            "-DUSE_TCALCFILTER=ON",
            "-DENABLE_V7=ON",
            "-DSTATIC_FILTERS=ON",
            "-DENABLE_COMPRESSION=OFF",
            "-DENABLE_SM_LIP=OFF",
            "-DENABLE_C_LIP=ON",
            "-DENABLE_CPP_LIP=ON",
            "-DLIP_SHARED_MEMORY=OFF",
            "-DENABLE_ASCMEMORY=OFF",
            "-DENABLE_FPE_HANDLING=OFF",
            "-DENABLE_UNCRUSTIFY=OFF",
            "-DENABLE_VALGRIND=OFF",

            "-DUSE_HDF=ON",
            "-DHDF_ROOT=%s" % spec['hdf5'].prefix,
            "-DHDF5_ROOT=%s" % spec['hdf5'].prefix,

            "-DENABLE_FORTRAN=ON",
            "-DENABLE_FORTRAN_INTERFACE=ON",
            "-DENABLE_C=ON",
            "-DENABLE_C_INTERFACE=ON",

            '-DEOS_DATA_ROOT_DIR=/usr/gapps/data/eos',
        ]

        if "+cpp14" in spec:
            args.extend(["-DBLT_CXX_STD=c++14", "-DCMAKE_CXX_STANDARD=14"])

        # If present, pass in language-specific flags to CMake
        cppflags = " ".join(spec.compiler_flags["cppflags"])
        if cppflags:
            # avoid always ending up with ' ' with no flags defined
            cppflags += " "
        cflags = cppflags + " ".join(spec.compiler_flags["cflags"])
        if cflags:
            args.append(f'-DBLT_C_FLAGS="{cflags}"')

        cxxflags = cppflags + " ".join(spec.compiler_flags["cxxflags"])
        if cxxflags:
            args.append(f'-DBLT_CXX_FLAGS="{cxxflags}"')

        fflags = " ".join(spec.compiler_flags["fflags"])
        if spec.satisfies("%cce"):
            fflags += " -ef"
        if spec.satisfies("%rocmcc"):
            fflags += " -O3 -Mstandard -Mfreeform -Mextend -Mbackslash -fno-default-real-8 -cpp"
        if fflags:
            args.append(f'-DBLT_FORTRAN_FLAGS={fflags}')

        return args

    def build(self, spec, prefix):
        fflags_path = join_path(self.build_directory, 'interfaces/fortran/CMakeFiles/leos_F.dir/flags.make')
        tty.msg("Fortran flags path: " + fflags_path)
        filter_file('-DAND', '', fflags_path)
        super(Leos, self).install(spec, prefix)

    def patch(self):
        if self.spec.satisfies('@8.4.1:  ^raja@2024.02.2:'):
            filter_file('loop_', 'seq_', "LIP/cpp_source/raja_api/SimpleRajaDefinitions.hxx", backup=True)

        if self.spec.satisfies('@8.2.0') and self.spec.satisfies('%gcc@10:'):
            filter_file(r'functional>', 'functional>\n#include<string>', 'src/implementation/LEOS_MemoryManager.hxx')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/fortran/LEOS_Function_F.f90')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.h')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.cc')
        elif self.spec.satisfies('@8.1.1p1') and ("%gcc" in self.spec or ("%clang" in self.spec and not self.spec.target.family=='ppc64le')): # gfortran issue
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/fortran/LEOS_Function_F.f90')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.h')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.cc')
        if self.spec.satisfies('@8.0.5:'):
            filter_file(r'USE_MPI','USE_MPI AND BUILD_TESTS', 'interfaces/C/CMakeLists.txt')
            filter_file(r'USE_MPI','USE_MPI AND BUILD_TESTS', 'interfaces/fortran/CMakeLists.txt')
            filter_file(r'include\(leos-mpi\)','add_definitions(-DLEOS_USE_MPI)', 'CMakeLists.txt')
            if self.spec.satisfies('@8.2:'):
                filter_file(r'NOT \${HDF5_C_HL_LIBRARIES}', 'NOT "${HDF5_C_HL_LIBRARIES}"', 'cmake/Setup3rdParty.cmake')

            if "%pgi" in self.spec:
                filter_file(r' tolower', ' ::tolower', 'src/io/ReaderFactory.cpp')
                filter_file(r' -cpp','', 'interfaces/fortran/CMakeLists.txt')
                if self.spec.satisfies('@8.2:'):
                    filter_file(r'__x86_64__', '__FOO__', 'cereal/include/cereal/external/rapidjson/internal/biginteger.h')
                    filter_file(r'__x86_64__', '__FOO__', 'cereal/include/cereal/external/rapidjson/internal/diyfp.h') # avoid uint128 undefined

            if "%clang@bgq" in self.spec or "%gcc@bgq" in self.spec:
                filter_file(r'-WF,-D', '-D', 'interfaces/fortran/CMakeLists.txt')
        elif self.spec.satisfies('@8.0.0:'):
            if (self.spec.satisfies('platform=linux') and self.spec.target.family=='ppc64le') or ("%arm@sandia" in self.spec and self.spec.satisfies('@8.0.2:')):
                filter_file(r'typedef int LEOS_MPI_Comm', '#include "mpi.h"\n typedef MPI_Comm LEOS_MPI_Comm', 'interfaces/C/LEOS_types_C.h')

            if self.spec.satisfies("@8.0.4:"):
                datafile="src/io/ReaderFactory.cpp"
            else:
                datafile="src/io/databaseNames.h"
            filter_file(r'/usr/gapps/data/eos', self.datadir, datafile)
        if self.spec.satisfies('@8.0.5:8.1.0') and "%cce" in self.spec:
            filter_file(r'\${TMP_VERS}','"#define BOOST_LIB_VERSION \"1_54\""', 'CMAKE/leos-boost.cmake')
        if self.spec.satisfies('@8.0:8.2') and "%pgi" in self.spec:
            filter_file(r'LIP_setup.h \*/', 'LIP_setup.h */\n#ifndef LIPSETUPH\n#define LIPSETUPH', 'LIP/source/LIP_setup.h')
            with open("LIP/source/LIP_setup.h", 'a') as f:
                f.write("#endif")
            filter_file(r'LIP_proto.h \*/', 'LIP_proto.h */\n#ifndef LIPPROTOH\n#define LIPPROTOH', 'LIP/source/LIP_proto.h')
            filter_file(r'\}\n', '}\n#endif\n', 'LIP/source/LIP_proto.h')
            filter_file(r'LIP_setup.h', 'LIP_setup.h"\n#include "LIP_PBH.h', 'LIP/source/lip_setup_interp_shmem.cc')

    def old_flag_handler(self, name, flags):
        my_flags = flags
        if name in ('cflags', 'cxxflags'):
            if any(self.spec.compiler.satisfies(c) for c in self.fpic_compilers):
                my_flags.append(self.compiler.cc_pic_flag)
                my_flags.append('-g -fp-model precise -fp-model source -unroll-aggressive -finline-functions -nolib-inline -axCORE-AVX2 -xAVX -no-fma -ip -no-ftz -prec-div -prec-sqrt -diag-disable cpu-dispatch')
            else:
                my_flags.append('-g')
        if name == 'cxxflags' and any(self.spec.compiler.satisfies(c) for c in self.c_o2_compilers):
            my_flags = list(map(lambda x: "" if x == "-O2" else x, my_flags))
        return (None, my_flags, None)

    @property
    def datadir(self):
        return '/usr/gapps/data/eos'
