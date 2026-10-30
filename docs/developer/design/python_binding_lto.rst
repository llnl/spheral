Python bindings and link-time optimization
==========================================

Spheral builds its Python bindings without link-time optimization (LTO). The
rest of Spheral builds normally; LTO settings for the native libraries are
unchanged.

Some binding modules use different pybind11 conversion rules for the same C++
type. For example, CXXTypes exposes ``std::vector<std::vector<int>>`` as a
bound Python vector object, while Geometry converts the same type to and from
a Python list of lists. Each rule produces a different definition of the same
pybind11 converter, which violates C++'s One Definition Rule (ODR).

LTO can be enabled in the bindings without Spheral asking for it. In the tested
Ubuntu 26.04 environment, MPICH supplies LTO compile flags that propagate
through dependency targets into the Python bindings. With LTO, the compiler
sees all of the definitions participating in a link and can report these ODR
violations as diagnostics, which prevent the build when warnings are treated
as errors. Spheral therefore disables LTO for the binding targets explicitly
and filters LTO-related flags inherited from dependencies. This filtering does
not remove LTO flags added globally or injected directly by a compiler
wrapper; those flags must be handled separately.

Optimized builds can work in practice when inlining keeps the code for bound
Python vector objects separate from the code that converts Python lists of
lists. Spheral's CTest and ATS test suites pass in this configuration, with
some cases skipped or filtered. Passing tests do not guarantee that this
separation will persist across compiler versions or build settings. The
problem was reproduced with binding LTO disabled in a GCC 13.3.1 build with
unoptimized Debug bindings: Geometry's ``initializePolygon`` binding used
initialization code from the incompatible converter and crashed with an
invalid free while converting a Python list of lists into a C++ vector of
vectors.

Disabling LTO keeps the bindings building and preserves the existing Python
interfaces, but it is a mitigation, not a fix. The conflicting definitions
still exist; excluding the bindings from LTO avoids the associated LTO
diagnostics without resolving the underlying ODR violations.
