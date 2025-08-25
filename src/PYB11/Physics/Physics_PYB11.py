"""
Spheral Physics module.

Provides the Physics, GenericHydro, and other physics base classes
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Physics/Physics.hh"',
                  '"Physics/GenericHydro.hh"',
                  '"Physics/GenericBodyForce.hh"',
                  '"Boundary/Boundary.hh"',
                  '"ArtificialViscosity/ArtificialViscosityView.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Kernel/TableKernel.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
MassDensityType = PYB11enum(("SumDensity", 
                             "RigorousSumDensity",
                             "HybridSumDensity",
                             "IntegrateDensity",
                             "VoronoiCellDensity",
                             "SumVoronoiCellDensity",
                             "CorrectedSumDensity"), export_values=True)

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from Physics import *
from GenericHydro import *
from GenericBodyForce import *

for ndim in dims:
    exec('''
Physics%(ndim)id = PYB11TemplateClass(Physics, template_parameters="%(Dimension)s")
GenericHydro%(ndim)id = PYB11TemplateClass(GenericHydro, template_parameters="%(Dimension)s")
GenericBodyForce%(ndim)id = PYB11TemplateClass(GenericBodyForce, template_parameters="%(Dimension)s")

vector_of_Physics%(ndim)id = PYB11_bind_vector("Physics<%(Dimension)s>*", opaque=True, local=False)
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})
