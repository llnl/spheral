"""
Spheral SmoothignScale module.

Provides the physics packages for updating the smoothing scale
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"SmoothingScale/SmoothingScaleBase.hh"',
                  '"SmoothingScale/FixedSmoothingScale.hh"',
                  '"SmoothingScale/SPHSmoothingScale.hh"',
                  '"SmoothingScale/ASPHSmoothingScale.hh"',
                  '"SmoothingScale/ASPHClassicSmoothingScale.hh"',
                  '"SmoothingScale/ASPHSmoothingScaleUserFilter.hh"',
                  '"SmoothingScale/ASPHRadialFunctor.hh"',
                  '"SmoothingScale/polySecondMoment.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Kernel/TableKernel.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
HEvolutionType = PYB11enum(("IdealH", "IntegrateH", "FixedH"), export_values=True,
                     doc="The choices for updating the smoothing scale.")

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from SmoothingScaleBase import SmoothingScaleBase
from FixedSmoothingScale import FixedSmoothingScale
from SPHSmoothingScale import SPHSmoothingScale
from ASPHSmoothingScale import ASPHSmoothingScale
from ASPHClassicSmoothingScale import ASPHClassicSmoothingScale
from ASPHSmoothingScaleUserFilter import ASPHSmoothingScaleUserFilter
from ASPHRadialFunctor import ASPHRadialFunctor

for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    Vector = f"{Dimension}::Vector"
    exec(f'''
SmoothingScaleBase{ndim}d = PYB11TemplateClass(SmoothingScaleBase, template_parameters="{Dimension}")
FixedSmoothingScale{ndim}d = PYB11TemplateClass(FixedSmoothingScale, template_parameters="{Dimension}")
SPHSmoothingScale{ndim}d = PYB11TemplateClass(SPHSmoothingScale, template_parameters="{Dimension}")
ASPHSmoothingScale{ndim}d = PYB11TemplateClass(ASPHSmoothingScale, template_parameters="{Dimension}")
ASPHClassicSmoothingScale{ndim}d = PYB11TemplateClass(ASPHClassicSmoothingScale, template_parameters="{Dimension}")
ASPHSmoothingScaleUserFilter{ndim}d = PYB11TemplateClass(ASPHSmoothingScaleUserFilter, template_parameters="{Dimension}")
ASPHRadialFunctor{ndim}d = PYB11TemplateClass(ASPHRadialFunctor, template_parameters="{Dimension}")

@PYB11cppname("polySecondMoment")
def polySecondMoment{ndim}d(poly = "const {Dimension}::FacetedVolume&",
                            center = "const {Dimension}::Vector&"):
    "Return the second moment of a convex polytope"
    return "{Dimension}::SymTensor"
''')

        
