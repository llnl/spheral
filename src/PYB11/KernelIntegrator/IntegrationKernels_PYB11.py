"""
Spheral IntegrationKernels module
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from BilinearIndex import *
from FlatConnectivity import *
from KernelIntegrationData import *
from IntegrationCoefficient import *
from IntegrationKernel import *
from KernelIntegrator import *
from ManufacturedSolution import *
from RKIntegrationKernel import *
from SPHIntegrationKernel import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['<memory>',
                  '<string>',
                  '<tuple>',
                  '<vector>',
                  '<unordered_map>',
                  '"Boundary/Boundary.hh"',
                  '"DataBase/DataBase.hh"',
                  '"DataBase/State.hh"',
                  '"Field/FieldList.hh"',
                  '"Geometry/GeomPlane.hh"',
                  '"Hydro/HydroFieldNames.hh"',
                  '"KernelIntegrator/BilinearIndex.hh"',
                  '"KernelIntegrator/FlatConnectivity.hh"',
                  '"KernelIntegrator/IntegrationCoefficient.hh"',
                  '"KernelIntegrator/IntegrationKernel.hh"',
                  '"KernelIntegrator/KernelIntegrationData.hh"',
                  '"KernelIntegrator/KernelIntegrator.hh"',
                  '"KernelIntegrator/ManufacturedSolution.hh"',
                  '"KernelIntegrator/RKIntegrationKernel.hh"',
                  '"KernelIntegrator/SPHIntegrationKernel.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiations
#-------------------------------------------------------------------------------
for ndim in dims:
    # Dimension-dependent
    exec('''
BilinearIndex%(ndim)id = PYB11TemplateClass(BilinearIndex, template_parameters="Dim<%(ndim)i>")
FlatConnectivity%(ndim)id = PYB11TemplateClass(FlatConnectivity, template_parameters="Dim<%(ndim)i>")
IntegrationKernel%(ndim)id = PYB11TemplateClass(IntegrationKernel, template_parameters="Dim<%(ndim)i>")
SPHIntegrationKernel%(ndim)id = PYB11TemplateClass(SPHIntegrationKernel, template_parameters="Dim<%(ndim)i>")
KernelIntegrator%(ndim)id = PYB11TemplateClass(KernelIntegrator, template_parameters="Dim<%(ndim)i>")
KernelIntegrationData%(ndim)id = PYB11TemplateClass(KernelIntegrationData, template_parameters="Dim<%(ndim)i>")
ManufacturedFunction%(ndim)id = PYB11TemplateClass(ManufacturedFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedSteadyStateFunction%(ndim)id = PYB11TemplateClass(ManufacturedSteadyStateFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedConstantFunction%(ndim)id = PYB11TemplateClass(ManufacturedConstantFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedSinusoidalFunction%(ndim)id = PYB11TemplateClass(ManufacturedSinusoidalFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedWaveFunction%(ndim)id = PYB11TemplateClass(ManufacturedWaveFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedTransportSolution%(ndim)id = PYB11TemplateClass(ManufacturedTransportSolution, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})

    # Dependent on primitives
    dim_types = (("Dim<%i>::Scalar" % ndim, "Scalar"),
                 ("Dim<%i>::Vector" % ndim, "Vector"),
                 ("Dim<%i>::SymTensor" % ndim, "SymTensor"),
                 ("Dim<%i>::Tensor" % ndim, "Tensor"),
                 ("std::vector<Dim<%i>::Scalar>" % ndim, "StdVectorScalar"),
                 ("std::vector<Dim<%i>::Vector>" % ndim, "StdVectorVector"))
    for (value, label) in dim_types:
        exec('''
%(label)sConstantIntegrationCoefficient%(ndim)id = PYB11TemplateClass(ConstantIntegrationCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sIntegrationCoefficient%(ndim)id = PYB11TemplateClass(IntegrationCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sFieldListIntegrationCoefficient%(ndim)id = PYB11TemplateClass(FieldListIntegrationCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    # Dependent on order
    for order in (0, 1, 2, 3, 4, 5, 6, 7):
        exec('''
RKIntegrationKernel%(ndim)id%(order)s = PYB11TemplateClass(RKIntegrationKernel, template_parameters=("Dim<%(ndim)i>", "%(order)s"))
''' % {"ndim" : ndim,
       "order" : order})
