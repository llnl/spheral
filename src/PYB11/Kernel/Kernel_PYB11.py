"""
Spheral Kernel module.

Provide the standard SPH/CRK interpolation kernels.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/GeomPlane.hh"',
                  '"Kernel/BSplineKernel.hh"',
                  '"Kernel/W4SplineKernel.hh"',
                  '"Kernel/GaussianKernel.hh"',
                  '"Kernel/SuperGaussianKernel.hh"',
                  '"Kernel/PiGaussianKernel.hh"',
                  '"Kernel/HatKernel.hh"',
                  '"Kernel/SincKernel.hh"',
                  '"Kernel/NSincPolynomialKernel.hh"',
                  '"Kernel/NBSplineKernel.hh"',
                  '"Kernel/QuarticSplineKernel.hh"',
                  '"Kernel/QuinticSplineKernel.hh"',
                  '"Kernel/TableKernel.hh"',
                  '"Kernel/WendlandC2Kernel.hh"',
                  '"Kernel/WendlandC4Kernel.hh"',
                  '"Kernel/WendlandC6Kernel.hh"',
                  '"Kernel/ExpInvKernel.hh"',
                  '"Kernel/SphericalKernel.hh"',
                  '"Kernel/SphericalBiCubicSplineKernel.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from Kernel import *
from SphericalKernel import *
from SphericalBiCubicSplineKernel import *

# TODO: change template_parameters to match what is used for TableKernelView
for ndim in dims:
    for KT in ("BSpline","W4Spline", "Gaussian", "SuperGaussian", "PiGaussian",
               "Hat", "Sinc", "NSincPolynomial", "NBSpline", "QuarticSpline",
               "QuinticSpline", "WendlandC2", "WendlandC4", "WendlandC6", "ExpInv"):
        exec('''
_Kernel%(ndim)id_%(KT)s = PYB11TemplateClass(Kernel,
                                             template_parameters = ("Dim<%(ndim)i>", "%(KT)sKernel<Dim<%(ndim)i>>"))
%(KT)sKernel%(ndim)id = PYB11TemplateClass(%(KT)sKernel,
                                           template_parameters = ("Dim<%(ndim)i>", "%(KT)sKernel<Dim<%(ndim)i>>"),
                                            cppname = "%(KT)sKernel<Dim<%(ndim)i>>")
''' % {"ndim" : ndim,
       "KT"   : KT})
    exec(f'''
_Kernel{ndim}d_TableKernelView = PYB11TemplateClass(Kernel,
                                                    template_parameters = ("Dim<{ndim}>", "TableKernelView<Dim<{ndim}>>"))
TableKernelView{ndim}d = PYB11TemplateClass(TableKernelView,
                                            template_parameters = "Dim<{ndim}>",
                                            cppname = "TableKernelView<Dim<{ndim}>>")
''')
    exec(f'''
TableKernel{ndim}d = PYB11TemplateClass(TableKernel,
                                        template_parameters = "Dim<{ndim}>",
                                        cppname = "TableKernel<Dim<{ndim}>>")
''')
